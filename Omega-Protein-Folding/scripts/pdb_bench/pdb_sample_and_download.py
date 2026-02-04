#!/usr/bin/env python3
"""
Sample and download real protein structures from RCSB PDB (mmCIF).

This script is intentionally conservative:
- It queries for **Protein (only)** experimental structures (X-ray/EM) with
  resolution <= --max-resolution
- It downloads mmCIF files and then filters locally by chain length

Outputs:
- data/raw/pdb_cache/*.cif
- data/processed/pdb_bench/pdb_ids_n{N}_seed{seed}.txt
- data/processed/pdb_bench/pdb_download_manifest_n{N}_seed{seed}.csv

Example:
  python scripts/pdb_bench/pdb_sample_and_download.py --n 200 --seed 0
"""

from __future__ import annotations

import argparse
import concurrent.futures as cf
import csv
import json
import random
import time
from pathlib import Path
from typing import List, Tuple

import requests
from Bio.PDB.MMCIFParser import MMCIFParser


RCSB_SEARCH_URL = "https://search.rcsb.org/rcsbsearch/v2/query"
RCSB_DOWNLOAD_URL_TMPL = "https://files.rcsb.org/download/{pdb_id}.cif"


def ensure_dir(p: Path) -> None:
    p.mkdir(parents=True, exist_ok=True)


def rcbs_search_entry_ids(
    *,
    max_resolution: float,
    rows: int,
) -> List[str]:
    """
    Return a (potentially long) list of entry IDs that match basic quality filters.
    We intentionally do *not* filter by chain length here; we do it after download.
    """
    query = {
        "type": "group",
        "logical_operator": "and",
        "nodes": [
            {
                "type": "terminal",
                "service": "text",
                "parameters": {
                    "attribute": "rcsb_entry_info.selected_polymer_entity_types",
                    "operator": "exact_match",
                    "value": "Protein (only)",
                },
            },
            {
                "type": "terminal",
                "service": "text",
                "parameters": {
                    "attribute": "exptl.method",
                    "operator": "in",
                    "value": ["X-RAY DIFFRACTION", "ELECTRON MICROSCOPY"],
                },
            },
            {
                "type": "terminal",
                "service": "text",
                "parameters": {
                    "attribute": "rcsb_entry_info.resolution_combined",
                    "operator": "less_or_equal",
                    "value": float(max_resolution),
                },
            },
        ],
    }

    # RCSB enforces limits on paginate size; page defensively.
    page_size = min(5000, int(rows))
    ids: List[str] = []
    start = 0
    while len(ids) < rows:
        payload = {
            "query": query,
            "return_type": "entry",
            "request_options": {
                "results_verbosity": "minimal",
                "scoring_strategy": "combined",
                "paginate": {"start": int(start), "rows": int(page_size)},
            },
        }
        r = requests.post(RCSB_SEARCH_URL, json=payload, timeout=60)
        r.raise_for_status()
        data = r.json()
        result_set = data.get("result_set", [])
        batch = [x["identifier"].upper() for x in result_set if "identifier" in x]
        if not batch:
            break
        ids.extend(batch)
        start += page_size
        # Stop early if server returned fewer than requested
        if len(batch) < page_size:
            break
        # Politeness for the API
        time.sleep(0.1)
    ids = ids[: int(rows)]
    # De-dup while preserving order.
    seen = set()
    out = []
    for i in ids:
        if i not in seen:
            out.append(i)
            seen.add(i)
    return out


def download_cif(pdb_id: str, out_path: Path, *, session: requests.Session, sleep_s: float) -> Tuple[bool, str]:
    url = RCSB_DOWNLOAD_URL_TMPL.format(pdb_id=pdb_id.upper())
    try:
        resp = session.get(url, timeout=60)
        if resp.status_code == 404:
            return False, "404"
        resp.raise_for_status()
        out_path.write_bytes(resp.content)
        if sleep_s > 0:
            time.sleep(sleep_s)
        return True, "ok"
    except Exception as e:  # noqa: BLE001 (keep script dependency-light)
        return False, f"error:{type(e).__name__}"

def longest_chain_ca_len(cif_path: Path) -> Tuple[str, int]:
    """
    Return (chain_id, N) for the longest CA trace in the first model.
    """
    parser = MMCIFParser(QUIET=True)
    s = parser.get_structure(cif_path.stem, str(cif_path))
    model = next(iter(s.get_models()))

    best_chain = None
    best_n = 0
    for chain in model:
        n = 0
        for res in chain:
            if "CA" in res:
                n += 1
        if n > best_n:
            best_n = n
            best_chain = chain.id

    if best_chain is None or best_n <= 0:
        raise ValueError("No CA trace found")
    return str(best_chain), int(best_n)


def process_one(
    pdb_id: str,
    *,
    root: Path,
    cache_dir: Path,
    min_len: int,
    max_len: int,
    sleep_s: float,
) -> dict:
    """
    Download (if needed), parse longest-chain CA length, and decide selection.
    Returns a manifest row dict.
    """
    out_path = cache_dir / f"{pdb_id}.cif"
    url = RCSB_DOWNLOAD_URL_TMPL.format(pdb_id=pdb_id)

    with requests.Session() as sess:
        if out_path.exists():
            ok = True
            status = "cached"
        else:
            ok, status = download_cif(pdb_id, out_path, session=sess, sleep_s=sleep_s)

    row = {
        "pdb_id": pdb_id,
        "status": status,
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
        row["chain"] = chain_id
        row["N"] = str(N)
        if min_len <= N <= max_len:
            row["selected"] = "1"
            row["status"] = "selected" if status != "cached" else "cached_selected"
        else:
            row["status"] = "len_out_of_range"
    except Exception:
        row["status"] = "parse_fail"

    return row


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=200, help="Number of PDB entries to download (best effort).")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--max-resolution", type=float, default=2.5)
    ap.add_argument("--candidate-rows", type=int, default=5000, help="How many candidate IDs to fetch from search API.")
    ap.add_argument("--sleep-s", type=float, default=0.1, help="Politeness delay between downloads.")
    ap.add_argument("--min-len", type=int, default=60)
    ap.add_argument("--max-len", type=int, default=350)
    ap.add_argument("--workers", type=int, default=8, help="Parallel workers for download+parse.")
    args = ap.parse_args()

    root = Path(__file__).resolve().parents[2]
    cache_dir = root / "data" / "raw" / "pdb_cache"
    out_dir = root / "data" / "processed" / "pdb_bench"
    ensure_dir(cache_dir)
    ensure_dir(out_dir)

    # 1) Query candidate IDs
    candidates = rcbs_search_entry_ids(max_resolution=args.max_resolution, rows=args.candidate_rows)
    rnd = random.Random(args.seed)
    rnd.shuffle(candidates)

    # 2) Download+parse until we have n selected entries (length-filtered)
    manifest_path = out_dir / f"pdb_download_manifest_n{args.n}_seed{args.seed}.csv"
    ids_path = out_dir / f"pdb_ids_n{args.n}_seed{args.seed}.txt"
    meta_chain_path = out_dir / f"pdb_chain_meta_n{args.n}_seed{args.seed}.csv"

    rows: List[dict] = []
    selected: List[str] = []

    max_workers = max(1, int(args.workers))
    with cf.ThreadPoolExecutor(max_workers=max_workers) as ex:
        pending = set()
        it = iter(candidates)

        # Submit initial batch
        for _ in range(max_workers * 2):
            try:
                pdb_id = next(it)
            except StopIteration:
                break
            pending.add(
                ex.submit(
                    process_one,
                    pdb_id,
                    root=root,
                    cache_dir=cache_dir,
                    min_len=args.min_len,
                    max_len=args.max_len,
                    sleep_s=args.sleep_s,
                )
            )

        while pending and len(selected) < args.n:
            done, pending = cf.wait(pending, return_when=cf.FIRST_COMPLETED)
            for fut in done:
                row = fut.result()
                rows.append(row)
                if row.get("selected") == "1":
                    selected.append(str(row["pdb_id"]))
                if len(selected) >= args.n:
                    break
                try:
                    pdb_id = next(it)
                except StopIteration:
                    continue
                pending.add(
                    ex.submit(
                        process_one,
                        pdb_id,
                        root=root,
                        cache_dir=cache_dir,
                        min_len=args.min_len,
                        max_len=args.max_len,
                        sleep_s=args.sleep_s,
                    )
                )

    fieldnames = ["pdb_id", "status", "path", "url", "chain", "N", "selected"]
    with manifest_path.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        for r in rows:
            w.writerow(r)

    ids_path.write_text("\n".join(selected) + "\n", encoding="utf-8")

    with meta_chain_path.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=["pdb_id", "chain", "N"])
        w.writeheader()
        for r in rows:
            if r.get("selected") == "1":
                w.writerow({"pdb_id": r["pdb_id"], "chain": r.get("chain", ""), "N": r.get("N", "")})

    meta_path = out_dir / f"pdb_ids_n{args.n}_seed{args.seed}.meta.json"
    meta_path.write_text(
        json.dumps(
            {
                "n_requested": args.n,
                "n_selected": len(selected),
                "seed": args.seed,
                "max_resolution": args.max_resolution,
                "candidate_rows": args.candidate_rows,
                "min_len": args.min_len,
                "max_len": args.max_len,
                "workers": max_workers,
                "cache_dir": str(cache_dir.relative_to(root)).replace("\\", "/"),
            },
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )

    print(f"Selected {len(selected)}/{args.n} entries (length-filtered).")
    print(f"- IDs: {ids_path}")
    print(f"- Manifest: {manifest_path}")
    print(f"- Chain meta: {meta_chain_path}")


if __name__ == "__main__":
    main()

