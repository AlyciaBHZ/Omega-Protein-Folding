#!/usr/bin/env python3
"""
Download mmCIFs for an existing PDB ID list (best-effort) and write a manifest.

Why this exists:
- `pdb_sample_and_download.py` *samples* IDs from RCSB then downloads.
- For full1k bench runs we already have a fixed ID list (source-of-truth),
  and we only want to ensure the local cache is populated + auditable.

Outputs (tracked-friendly):
- docs/runs/<run_name>/pdb_download_manifest_{tag}.csv
- docs/runs/<run_name>/pdb_download_summary_{tag}.md

Cache:
- data/raw/pdb_cache/<PDB>.cif
"""

from __future__ import annotations

import argparse
import concurrent.futures as cf
import csv
import time
from pathlib import Path
from typing import Dict, List, Tuple

import requests

from bench_utils import ensure_dir
from pdb_sample_and_download import RCSB_DOWNLOAD_URL_TMPL, download_cif, longest_chain_ca_len


def _process_one(
    pdb_id: str,
    *,
    root: Path,
    cache_dir: Path,
    min_len: int,
    max_len: int,
    sleep_s: float,
) -> Dict[str, str]:
    pdb_id = str(pdb_id).strip().upper()
    out_path = cache_dir / f"{pdb_id}.cif"
    url = RCSB_DOWNLOAD_URL_TMPL.format(pdb_id=pdb_id)

    if out_path.exists():
        ok = True
        status = "cached"
    else:
        with requests.Session() as sess:
            ok, status = download_cif(pdb_id, out_path, session=sess, sleep_s=float(sleep_s))

    row: Dict[str, str] = {
        "pdb_id": pdb_id,
        "status": status if ok else status,
        "path": str(out_path.relative_to(root)).replace("\\", "/") if ok else "",
        "url": url,
        "chain": "",
        "N": "",
        "selected": "0",
    }

    if not ok:
        return row

    try:
        chain_id, N = longest_chain_ca_len(out_path)
        row["chain"] = str(chain_id)
        row["N"] = str(int(N))
        if int(min_len) <= int(N) <= int(max_len):
            row["selected"] = "1"
            if row["status"] == "cached":
                row["status"] = "cached_selected"
            elif row["status"] == "ok":
                row["status"] = "downloaded_selected"
        else:
            row["status"] = "len_out_of_range"
    except Exception:  # noqa: BLE001
        row["status"] = "parse_fail"

    return row


def _read_ids(path: Path) -> List[str]:
    ids = [x.strip().upper() for x in path.read_text(encoding="utf-8").splitlines() if x.strip()]
    # de-dup while preserving order
    seen = set()
    out: List[str] = []
    for x in ids:
        if x not in seen:
            out.append(x)
            seen.add(x)
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--ids", required=True, help="IDs file (one PDB ID per line).")
    ap.add_argument("--tag", default="n1000_seed0", help="Tag suffix used for manifest filenames.")
    ap.add_argument("--run-name", default="full1k_hybrid_auric", help="docs/runs/<run-name>/ output folder.")
    ap.add_argument("--workers", type=int, default=16)
    ap.add_argument("--sleep-s", type=float, default=0.05)
    ap.add_argument("--min-len", type=int, default=60)
    ap.add_argument("--max-len", type=int, default=350)
    args = ap.parse_args()

    root = Path(__file__).resolve().parents[2]
    ids_path = Path(args.ids)
    if not ids_path.is_absolute():
        ids_path = (root / ids_path).resolve()

    cache_dir = root / "data" / "raw" / "pdb_cache"
    ensure_dir(cache_dir)

    out_dir = root / "docs" / "runs" / str(args.run_name)
    ensure_dir(out_dir)

    ids = _read_ids(ids_path)
    if not ids:
        raise ValueError(f"Empty IDs file: {ids_path}")

    manifest_path = out_dir / f"pdb_download_manifest_{args.tag}.csv"
    summary_path = out_dir / f"pdb_download_summary_{args.tag}.md"

    t0 = time.perf_counter()
    rows: List[Dict[str, str]] = []
    max_workers = max(1, int(args.workers))

    with cf.ThreadPoolExecutor(max_workers=max_workers) as ex:
        futs = [
            ex.submit(
                _process_one,
                pdb_id,
                root=root,
                cache_dir=cache_dir,
                min_len=int(args.min_len),
                max_len=int(args.max_len),
                sleep_s=float(args.sleep_s),
            )
            for pdb_id in ids
        ]
        for fut in cf.as_completed(futs):
            rows.append(fut.result())

    # stable ordering like input ids
    row_by_id = {r["pdb_id"]: r for r in rows if "pdb_id" in r}
    rows_ordered = [row_by_id.get(pdb_id, {"pdb_id": pdb_id, "status": "missing", "path": "", "url": "", "chain": "", "N": "", "selected": "0"}) for pdb_id in ids]

    fieldnames = ["pdb_id", "status", "path", "url", "chain", "N", "selected"]
    with manifest_path.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        for r in rows_ordered:
            w.writerow(r)

    n_total = len(ids)
    n_cached = sum(1 for r in rows_ordered if r.get("status", "").startswith("cached"))
    n_downloaded = sum(1 for r in rows_ordered if r.get("status", "").startswith("downloaded"))
    n_selected = sum(1 for r in rows_ordered if r.get("selected", "") == "1")
    n_fail = sum(1 for r in rows_ordered if r.get("status", "") in {"404"} or r.get("status", "").startswith("error") or r.get("status", "") in {"parse_fail", "missing"})
    dt = time.perf_counter() - t0

    summary_path.write_text(
        "\n".join(
            [
                f"# PDB mmCIF cache manifest ({args.tag})",
                "",
                f"- IDs file: `{ids_path.relative_to(root).as_posix()}`",
                f"- Cache dir: `{cache_dir.relative_to(root).as_posix()}`",
                f"- Total IDs: {n_total}",
                f"- Cached: {n_cached}",
                f"- Downloaded: {n_downloaded}",
                f"- Selected (len in [{int(args.min_len)},{int(args.max_len)}]): {n_selected}",
                f"- Failures (404/error/parse/missing): {n_fail}",
                f"- Manifest: `{manifest_path.relative_to(root).as_posix()}`",
                f"- Runtime: {dt:.1f}s (workers={max_workers})",
                "",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    print(f"Wrote: {manifest_path}")
    print(f"Wrote: {summary_path}")


if __name__ == "__main__":
    main()

