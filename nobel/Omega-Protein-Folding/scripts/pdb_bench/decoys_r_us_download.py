#!/usr/bin/env python3
"""
Download and extract Decoys 'R' Us dataset(s) into the local decoy cache.

We intentionally do NOT commit the decoy binaries. This script:
- downloads a tarball (default: 4state_reduced) from Decoys 'R' Us
- extracts it under data/raw/decoy_cache/decoys_r_us/
- writes a lightweight manifest under docs/runs/decoys_4state_auric/

Dataset page (for reference): http://compbio.buffalo.edu/dd/download.shtml
"""

from __future__ import annotations

import argparse
import csv
import re
import tarfile
from pathlib import Path
from typing import Iterable, List, Tuple

import requests

from bench_utils import ensure_dir


# Direct tarball URL (preferred). The CGI page provides a link to this.
# Dataset index: http://compbio.buffalo.edu/dd/download.shtml
DEFAULT_4STATE_URL = "http://dd.compbio.org/4state_reduced.tgz"


def _iter_pdbs(root: Path) -> Iterable[Path]:
    for p in root.rglob("*.pdb"):
        if p.is_file():
            yield p


def _guess_pairs(extract_root: Path) -> List[Tuple[str, Path, Path]]:
    """
    Best-effort pairing for Decoys 'R' Us layouts.
    Returns list of (target_id, native_path, decoy_path).

    Many Decoys 'R' Us sets are organized as:
      <target>/native.pdb and <target>/decoys/*.pdb
    But the exact naming varies; we keep this permissive and let later
    scripts validate length/parse.
    """
    rows: List[Tuple[str, Path, Path]] = []

    # Dataset-specific common layout (4state_reduced, etc):
    #   doc/pdb_orig/<target>.pdb are natives
    #   <target>/*.pdb are decoys
    native_dir = extract_root / "doc" / "pdb_orig"
    native_map = {}
    if native_dir.is_dir():
        for p in native_dir.glob("*.pdb"):
            native_map[p.stem.lower()] = p

    for target_dir in sorted([d for d in extract_root.iterdir() if d.is_dir()]):
        if target_dir.name.lower() in {"doc"}:
            continue
        pdbs = list(_iter_pdbs(target_dir))
        if not pdbs:
            continue

        native = None
        if native_map:
            native = native_map.get(target_dir.name.lower())
        if native is None:
            # Fallback heuristic: if a folder contains one pdb with 'native' in name, treat it as native
            natives = [p for p in pdbs if "native" in p.name.lower()]
            if natives:
                native = natives[0]
            else:
                native = sorted(pdbs, key=lambda x: len(x.name))[0]

        for decoy in pdbs:
            if decoy.resolve() == native.resolve():
                continue
            rows.append((target_dir.name, native, decoy))

    return rows


def _maybe_follow_html_download_landing(url: str, body: bytes) -> str:
    """
    Decoys 'R' Us sometimes serves an HTML landing page with a .tgz link.
    If this looks like HTML, extract the first .tgz URL and return it.
    """
    head = body[:2048].lstrip()
    if not (head.startswith(b"<") or b"<html" in head.lower()):
        return url
    try:
        txt = body.decode("utf-8", errors="ignore")
    except Exception:
        return url
    m = re.search(r"(https?://[^\s\"']+?\.tgz)", txt)
    return m.group(1) if m else url


def download(url: str, out_path: Path) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    r = requests.get(url, timeout=180)
    r.raise_for_status()
    data = r.content
    url2 = _maybe_follow_html_download_landing(url, data)
    if url2 != url:
        r2 = requests.get(url2, timeout=180)
        r2.raise_for_status()
        data = r2.content
    out_path.write_bytes(data)


def extract_tgz(tgz_path: Path, out_dir: Path) -> Path:
    out_dir.mkdir(parents=True, exist_ok=True)
    with tarfile.open(tgz_path, "r:gz") as tf:
        tf.extractall(out_dir)

    # Some archives unpack into a single top-level folder; normalize by returning it if present
    kids = [p for p in out_dir.iterdir() if p.is_dir()]
    if len(kids) == 1:
        return kids[0]
    return out_dir


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--url", default=DEFAULT_4STATE_URL, help="Decoy tgz URL.")
    ap.add_argument("--name", default="4state_reduced", help="Dataset short name (folder naming).")
    ap.add_argument("--run-name", default="decoys_4state_auric", help="docs/runs/<run-name>/ output folder.")
    ap.add_argument("--force", action="store_true", help="Re-download and re-extract even if present.")
    args = ap.parse_args()

    root = Path(__file__).resolve().parents[2]
    cache_root = root / "data" / "raw" / "decoy_cache" / "decoys_r_us" / str(args.name)
    tgz_path = root / "data" / "raw" / "decoy_cache" / "decoys_r_us" / f"{args.name}.tgz"

    if args.force and cache_root.exists():
        # Keep it simple: user can delete cache dir manually if needed.
        pass

    if args.force or not tgz_path.exists():
        print(f"Downloading {args.url} -> {tgz_path}")
        download(str(args.url), tgz_path)
    else:
        print(f"Using cached archive: {tgz_path}")

    if args.force or not cache_root.exists():
        print(f"Extracting -> {cache_root}")
        extract_root = extract_tgz(tgz_path, cache_root)
    else:
        extract_root = cache_root

    # Many decoy sets unpack under dd/(multiple|single|loop)/<setname>/...; normalize.
    name = str(args.name)
    candidates = [
        extract_root / "dd" / "multiple" / name,
        extract_root / "dd" / "single" / name,
        extract_root / "dd" / "loop" / name,
        extract_root / "dd" / name,
        extract_root / name,
        extract_root / "dd" / "multiple",
        extract_root / "dd" / "single",
        extract_root / "dd" / "loop",
        extract_root / "dd",
    ]
    for c in candidates:
        if c.is_dir():
            extract_root = c
            break

    # Manifest (lightweight, committed)
    run_dir = root / "docs" / "runs" / str(args.run_name)
    ensure_dir(run_dir)
    manifest_path = run_dir / "decoy_manifest.csv"
    summary_path = run_dir / "decoy_download_summary.md"

    pairs = _guess_pairs(extract_root)
    with manifest_path.open("w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["dataset", "target_id", "native_path", "decoy_path"])
        for target_id, native, decoy in pairs:
            w.writerow(
                [
                    str(args.name),
                    str(target_id),
                    str(native.relative_to(root)).replace("\\", "/"),
                    str(decoy.relative_to(root)).replace("\\", "/"),
                ]
            )

    uniq_targets = sorted({t for (t, _, _) in pairs})
    n_targets = len(uniq_targets)

    summary_path.write_text(
        "\n".join(
            [
                f"# Decoys 'R' Us download summary ({args.name})",
                "",
                f"- URL: `{args.url}`",
                f"- Cache root: `{cache_root.relative_to(root).as_posix()}`",
                f"- Extract root: `{extract_root.relative_to(root).as_posix()}`",
                f"- Targets (heuristic): {n_targets}",
                f"- (native,decoy) pairs: {len(pairs)}",
                f"- Manifest: `{manifest_path.relative_to(root).as_posix()}`",
                "",
                "Note: pairing is best-effort; the downstream Auric bench validates parse/length and reports any skipped entries.",
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

