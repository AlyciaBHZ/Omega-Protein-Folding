from __future__ import annotations

import argparse
import csv
import hashlib
import sys
import time
import urllib.error
import urllib.request
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np


def ensure_dir(p: Path) -> None:
    p.mkdir(parents=True, exist_ok=True)


def afdb_pdb_url(uniprot_id: str, *, version: str = "v4") -> str:
    uid = str(uniprot_id).strip()
    if not uid:
        raise ValueError("Empty UniProt ID")
    # AlphaFold DB filename convention:
    #   AF-<UniProt>-F1-model_v4.pdb
    return f"https://alphafold.ebi.ac.uk/files/AF-{uid}-F1-model_{version}.pdb"


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def download(url: str, out_path: Path, *, force: bool = False, timeout_s: float = 120.0) -> None:
    ensure_dir(out_path.parent)
    if out_path.exists() and not force:
        return

    req = urllib.request.Request(
        url,
        headers={
            "User-Agent": "Omega-Protein-Folding/afdb_download.py (+https://github.com/AlyciaBHZ/Omega-Protein-Folding)",
            "Accept": "text/plain,*/*",
        },
        method="GET",
    )
    try:
        with urllib.request.urlopen(req, timeout=float(timeout_s)) as resp:
            tmp = out_path.with_suffix(out_path.suffix + ".tmp")
            with tmp.open("wb") as f:
                while True:
                    chunk = resp.read(1024 * 1024)
                    if not chunk:
                        break
                    f.write(chunk)
            tmp.replace(out_path)
    except urllib.error.HTTPError as e:
        raise RuntimeError(f"HTTP error {e.code} downloading {url}") from e
    except urllib.error.URLError as e:
        raise RuntimeError(f"URL error downloading {url}: {e}") from e


def _try_api_pred_url(uniprot_id: str, *, timeout_s: float = 20.0) -> str | None:
    """
    Best-effort: use AlphaFold DB API to find the current model PDB URL.
    Returns None if API is unavailable (timeout/blocked) or response is unexpected.
    """
    import json

    uid = str(uniprot_id).strip()
    api = f"https://alphafold.ebi.ac.uk/api/prediction/{uid}"
    req = urllib.request.Request(
        api,
        headers={"User-Agent": "Omega-Protein-Folding/afdb_download.py", "Accept": "application/json"},
        method="GET",
    )
    try:
        with urllib.request.urlopen(req, timeout=float(timeout_s)) as resp:
            raw = resp.read()
    except Exception:
        return None
    try:
        obj = json.loads(raw.decode("utf-8", "ignore"))
    except Exception:
        return None
    if isinstance(obj, list) and obj:
        # Historically: list[ { "pdbUrl": "...", "cifUrl": "...", ... } ]
        it = obj[0]
        if isinstance(it, dict):
            for k in ("pdbUrl", "pdb_url", "pdb"):
                v = it.get(k)
                if isinstance(v, str) and v.startswith("http"):
                    return v
    return None


def parse_pdb_ca_plddt(pdb_path: Path) -> Tuple[np.ndarray, np.ndarray, List[Tuple[str, int, str]]]:
    """
    Parse CA coordinates and pLDDT (stored in B-factor) from an AlphaFold DB PDB.

    Returns:
      - coords: (N,3) float64
      - plddt: (N,) float64
      - res_ids: list of (chain, resseq, icode)
    """
    coords: List[List[float]] = []
    plddt: List[float] = []
    res_ids: List[Tuple[str, int, str]] = []
    seen = set()

    # PDB fixed columns:
    # 13-16 atom name, 17 altLoc, 22 chainID, 23-26 resSeq, 27 iCode
    # 31-38 x, 39-46 y, 47-54 z
    # 61-66 B-factor
    for line in pdb_path.read_text(encoding="utf-8", errors="replace").splitlines():
        if not line.startswith("ATOM"):
            continue
        if len(line) < 66:
            continue

        atom = line[12:16].strip()
        if atom != "CA":
            continue

        altloc = line[16:17]
        if altloc not in {" ", "A"}:
            continue

        chain = line[21:22].strip() or "?"
        resseq_str = line[22:26].strip()
        icode = line[26:27].strip()
        try:
            resseq = int(resseq_str)
        except ValueError:
            continue

        rid = (chain, resseq, icode)
        if rid in seen:
            continue
        seen.add(rid)

        try:
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
            b = float(line[60:66])
        except ValueError:
            continue

        coords.append([x, y, z])
        plddt.append(b)
        res_ids.append(rid)

    if len(coords) < 2:
        raise ValueError(f"No CA trace parsed from {pdb_path}")

    return np.asarray(coords, dtype=np.float64), np.asarray(plddt, dtype=np.float64), res_ids


def write_plddt_csv(csv_path: Path, plddt: np.ndarray, res_ids: List[Tuple[str, int, str]]) -> None:
    ensure_dir(csv_path.parent)
    with csv_path.open("w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["i", "chain", "resseq", "icode", "pLDDT"])
        for i, ((chain, resseq, icode), p) in enumerate(zip(res_ids, plddt)):
            w.writerow([int(i), str(chain), int(resseq), str(icode), float(p)])


def write_summary_md(
    md_path: Path,
    *,
    uniprot_id: str,
    url: str,
    pdb_path: Path,
    plddt: np.ndarray,
    extra: Dict[str, str] | None = None,
) -> None:
    ensure_dir(md_path.parent)
    p = np.asarray(plddt, dtype=np.float64)
    frac_lt50 = float(np.mean(p < 50.0))
    frac_lt70 = float(np.mean(p < 70.0))
    lines = []
    lines.append(f"# AFDB download summary: {uniprot_id}\n")
    lines.append(f"- UniProt: `{uniprot_id}`\n")
    lines.append(f"- Source: `{url}`\n")
    lines.append(f"- Local PDB: `{pdb_path.as_posix()}`\n")
    lines.append(f"- SHA256: `{_sha256(pdb_path)}`\n")
    lines.append(f"- Downloaded at: `{time.strftime('%Y-%m-%d %H:%M:%S')}`\n")
    lines.append("\n## pLDDT statistics\n\n")
    lines.append(f"- N residues (CA parsed): **{len(p)}**\n")
    lines.append(f"- min / median / mean / max: **{p.min():.1f} / {np.median(p):.1f} / {p.mean():.1f} / {p.max():.1f}**\n")
    lines.append(f"- frac(pLDDT < 50): **{frac_lt50:.3f}**\n")
    lines.append(f"- frac(pLDDT < 70): **{frac_lt70:.3f}**\n")
    if extra:
        lines.append("\n## Notes\n\n")
        for k, v in extra.items():
            lines.append(f"- {k}: {v}\n")

    md_path.write_text("".join(lines), encoding="utf-8")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--uniprot", required=True, help="UniProt ID (e.g. P37840)")
    ap.add_argument(
        "--version",
        default="auto",
        help="AFDB model version for /files URL (v6/v5/v4) or 'auto' (default).",
    )
    ap.add_argument(
        "--out-dir",
        default="docs/runs/af_vs_omega",
        help="Output directory (relative to repo root unless absolute).",
    )
    ap.add_argument(
        "--cache-dir",
        default="data/raw/afdb_cache",
        help="Cache directory for downloaded PDBs (relative to repo root unless absolute).",
    )
    ap.add_argument("--force", action="store_true", help="Re-download even if cached exists.")
    args = ap.parse_args()

    root = Path(__file__).resolve().parents[2]
    out_dir = Path(args.out_dir)
    if not out_dir.is_absolute():
        out_dir = (root / out_dir).resolve()
    cache_dir = Path(args.cache_dir)
    if not cache_dir.is_absolute():
        cache_dir = (root / cache_dir).resolve()

    uid = str(args.uniprot).strip()
    vers = str(args.version).strip().lower()
    url: str | None = None
    if vers == "auto":
        url = _try_api_pred_url(uid)
        if url is None:
            # Fall back to /files pattern with recent versions first.
            for v in ("v6", "v5", "v4"):
                cand = afdb_pdb_url(uid, version=v)
                try:
                    cache_pdb = cache_dir / f"afdb_{uid}_model_{v}.pdb"
                    download(cand, cache_pdb, force=bool(args.force))
                    url = cand
                    vers = v
                    break
                except RuntimeError:
                    continue
            if url is None:
                raise RuntimeError(
                    f"Could not download AFDB model for {uid}. "
                    "Tried API and /files for v6/v5/v4."
                )
    else:
        url = afdb_pdb_url(uid, version=vers)
        cache_pdb = cache_dir / f"afdb_{uid}_model_{vers}.pdb"
        download(url, cache_pdb, force=bool(args.force))

    # If API resolved a URL directly, cache under a stable name.
    if vers == "auto":
        vers = "api"
    if "cache_pdb" not in locals():
        cache_pdb = cache_dir / f"afdb_{uid}_model_{vers}.pdb"
        download(str(url), cache_pdb, force=bool(args.force))

    # Copy into run folder for auditability (small, plain-text PDB).
    ensure_dir(out_dir)
    run_pdb = out_dir / f"afdb_{uid}_model.pdb"
    if (not run_pdb.exists()) or bool(args.force):
        run_pdb.write_bytes(cache_pdb.read_bytes())

    coords, plddt, res_ids = parse_pdb_ca_plddt(run_pdb)
    _ = coords  # kept for downstream scripts; we only write pLDDT artifacts here

    csv_path = out_dir / f"afdb_{uid}_plddt.csv"
    md_path = out_dir / f"afdb_{uid}_summary.md"
    write_plddt_csv(csv_path, plddt, res_ids)
    write_summary_md(
        md_path,
        uniprot_id=uid,
        url=url,
        pdb_path=run_pdb.relative_to(root) if run_pdb.is_absolute() else run_pdb,
        plddt=plddt,
        extra={"AFDB model version": str(vers)},
    )

    print(f"[ok] wrote: {csv_path}")
    print(f"[ok] wrote: {md_path}")


if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        print(f"[error] {e}", file=sys.stderr)
        raise

