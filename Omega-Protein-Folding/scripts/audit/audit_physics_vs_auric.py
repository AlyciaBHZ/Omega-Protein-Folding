#!/usr/bin/env python3
from __future__ import annotations

"""
Batch audit: compute CA-level physics proxies and Auric geometry certificate scores.

Default data source: toy decoy smoke PDBs (committed):
  docs/runs/toy_decoy_auric_smoke/pdbs/*.pdb

Outputs (committable):
  docs/runs/physics_audit/physics_vs_auric_toy_decoy.csv
"""

import argparse
import csv
import sys
from pathlib import Path
from typing import Dict, List

import numpy as np

from physics_metrics import calc_physics_metrics


def repo_root() -> Path:
    return Path(__file__).resolve().parents[2]


def ensure_dir(p: Path) -> None:
    p.mkdir(parents=True, exist_ok=True)


def load_ca_coords_fast(pdb_path: Path) -> np.ndarray:
    """
    Minimal CA-only loader (works for both full-atom and CA-only PDBs).
    """
    coords: List[List[float]] = []
    for ln in Path(pdb_path).read_text(encoding="utf-8", errors="ignore").splitlines():
        if not ln.startswith("ATOM"):
            continue
        if ln[12:16].strip() != "CA":
            continue
        try:
            x = float(ln[30:38])
            y = float(ln[38:46])
            z = float(ln[46:54])
        except Exception:
            continue
        coords.append([x, y, z])
    if len(coords) < 3:
        return np.zeros((0, 3), dtype=np.float64)
    return np.asarray(coords, dtype=np.float64)


def auric_geometry_scores_from_coords(coords: np.ndarray, *, alphabet: str, m: int) -> Dict[str, float]:
    """
    Compute Auric rhoA scores from y_perp built by oracle direction quantization.
    Returns type_entropy and smb_rate_hat for rhoA.
    """
    root = repo_root()
    sys.path.insert(0, str((root / "scripts" / "pdb_bench").resolve()))
    sys.path.insert(0, str((root / "scripts" / "auric").resolve()))

    from bench_utils import oracle_direction_reconstruct  # type: ignore  # noqa: E402
    from axis import shared_pca_axis  # type: ignore  # noqa: E402
    from fold_m import fold_sliding_windows_bits  # type: ignore  # noqa: E402
    from metrics import metrics_for_stream  # type: ignore  # noqa: E402
    from readout import rho_A_from_yperp  # type: ignore  # noqa: E402

    _, _n_path, y_perp, ph_rms, ph_max = oracle_direction_reconstruct(
        np.asarray(coords, dtype=np.float64),
        alphabet=str(alphabet),
    )
    u = shared_pca_axis(y_perp)
    bits_A = rho_A_from_yperp(y_perp, u=tuple(u.tolist()), u_mode="fixed", threshold="median")
    folded_A = fold_sliding_windows_bits(bits_A, m=int(m))
    met_A = metrics_for_stream(bits_A, folded_A)

    return {
        "ph_rms": float(ph_rms),
        "ph_max": float(ph_max),
        "auric_type_entropy_rhoA": float(met_A.get("type_entropy", float("nan"))),
        "auric_smb_rate_hat_rhoA": float(met_A.get("smb_rate_hat", float("nan"))),
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--pdb-glob",
        default="docs/runs/toy_decoy_auric_smoke/pdbs/*.pdb",
        help="Glob relative to repo root.",
    )
    ap.add_argument("--alphabet", default="triple232", choices=["axis12", "pair72", "triple232"])
    ap.add_argument("--m", type=int, default=8)
    ap.add_argument("--out-csv", default="docs/runs/physics_audit/physics_vs_auric_toy_decoy.csv")
    args = ap.parse_args()

    root = repo_root()
    g = Path(str(args.pdb_glob))
    glob_root = (root / g).parent if not g.is_absolute() else g.parent
    pat = g.name
    pdb_paths = sorted(glob_root.glob(pat))
    if not pdb_paths:
        raise SystemExit(f"No PDBs matched: {args.pdb_glob}")

    out_csv = Path(str(args.out_csv))
    if not out_csv.is_absolute():
        out_csv = (root / out_csv).resolve()
    ensure_dir(out_csv.parent)

    rows: List[Dict[str, object]] = []
    for p in pdb_paths:
        p = p.resolve()
        rel = p.relative_to(root).as_posix()
        name = p.name
        label = "native" if name.lower().startswith("native_") else "decoy"

        phys = calc_physics_metrics(p)
        coords = load_ca_coords_fast(p)
        aur = auric_geometry_scores_from_coords(coords, alphabet=str(args.alphabet), m=int(args.m)) if len(coords) else {}

        row: Dict[str, object] = {"pdb_path": rel, "filename": name, "label": label}
        row.update({k: phys.get(k, "") for k in ["N", "bond_mse", "clash_count", "clash_per_res", "rg"]})
        row.update(aur)
        rows.append(row)

    cols = [
        "pdb_path",
        "filename",
        "label",
        "N",
        "bond_mse",
        "clash_count",
        "clash_per_res",
        "rg",
        "ph_rms",
        "ph_max",
        "auric_type_entropy_rhoA",
        "auric_smb_rate_hat_rhoA",
    ]
    with out_csv.open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols)
        w.writeheader()
        for r in rows:
            w.writerow({c: r.get(c, "") for c in cols})

    print(f"[ok] wrote: {out_csv}")
    print(f"[ok] rows: {len(rows)}")


if __name__ == "__main__":
    main()

