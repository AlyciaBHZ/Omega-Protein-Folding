from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path
from typing import Dict, List

import numpy as np
import pandas as pd


def ensure_dir(p: Path) -> None:
    p.mkdir(parents=True, exist_ok=True)


def _summary_from_plddt_csv(path: Path) -> Dict[str, float]:
    df = pd.read_csv(path)
    p = df["pLDDT"].to_numpy(dtype=np.float64)
    return {
        "plddt_N": float(len(p)),
        "plddt_min": float(np.min(p)),
        "plddt_median": float(np.median(p)),
        "plddt_mean": float(np.mean(p)),
        "plddt_max": float(np.max(p)),
        "plddt_frac_lt50": float(np.mean(p < 50.0)),
        "plddt_frac_lt70": float(np.mean(p < 70.0)),
    }


def _parse_pdb_ca_coords(pdb_path: Path) -> np.ndarray:
    coords: List[List[float]] = []
    seen = set()
    for line in pdb_path.read_text(encoding="utf-8", errors="replace").splitlines():
        if not line.startswith("ATOM") or len(line) < 54:
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
        except ValueError:
            continue
        coords.append([x, y, z])
    if len(coords) < 2:
        raise ValueError(f"No CA coords parsed from {pdb_path}")
    return np.asarray(coords, dtype=np.float64)


def _auric_type_entropy_from_yperp(y_perp: np.ndarray, *, m: int) -> float:
    root = Path(__file__).resolve().parents[2]
    sys.path.insert(0, str((root / "scripts" / "auric").resolve()))
    from fold_m import fold_sliding_windows_bits  # noqa: E402
    from metrics import metrics_for_stream  # noqa: E402
    from readout import rho_A_from_yperp  # noqa: E402

    bits_A = rho_A_from_yperp(y_perp, u_mode="pca", threshold="median")
    folded = fold_sliding_windows_bits(bits_A, m=int(m))
    met = metrics_for_stream(bits_A, folded)
    return float(met.get("type_entropy", float("nan")))


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out-csv", default="docs/runs/af_vs_omega/case_table.csv")
    ap.add_argument("--m", type=int, default=8)
    args = ap.parse_args()

    root = Path(__file__).resolve().parents[2]
    out_csv = Path(args.out_csv)
    if not out_csv.is_absolute():
        out_csv = (root / out_csv).resolve()
    ensure_dir(out_csv.parent)

    # Inputs (convention: our run folder has these exact names)
    cases = [
        {
            "case": "AFDB_P06454",
            "pdb": root / "docs" / "runs" / "af_vs_omega" / "afdb_P06454_model.pdb",
            "plddt_csv": root / "docs" / "runs" / "af_vs_omega" / "afdb_P06454_plddt.csv",
        },
        {
            "case": "AFDB_P04637",
            "pdb": root / "docs" / "runs" / "af_vs_omega" / "afdb_P04637_model.pdb",
            "plddt_csv": root / "docs" / "runs" / "af_vs_omega" / "afdb_P04637_plddt.csv",
        },
        {
            "case": "OMEGA_1CRN",
            "pdb": root / "docs" / "runs" / "af_vs_omega" / "omega_pred_1CRN_v1.pdb",
            "plddt_csv": None,
        },
    ]

    # Reuse existing reconstruction helper.
    sys.path.insert(0, str((root / "scripts" / "pdb_bench").resolve()))
    from bench_utils import oracle_direction_reconstruct  # noqa: E402

    rows: List[Dict[str, float | str]] = []
    for c in cases:
        coords = _parse_pdb_ca_coords(Path(c["pdb"]))
        _, _, y_perp, ph_rms, ph_max = oracle_direction_reconstruct(coords, alphabet="triple232")
        H = _auric_type_entropy_from_yperp(y_perp, m=int(args.m))

        row: Dict[str, float | str] = {
            "case": str(c["case"]),
            "N_CA": float(coords.shape[0]),
            "ph_rms": float(ph_rms),
            "ph_max": float(ph_max),
            f"H_type_m{int(args.m)}_rhoA_pca_median": float(H),
        }
        if c["plddt_csv"] is not None:
            row.update(_summary_from_plddt_csv(Path(c["plddt_csv"])))
        rows.append(row)

    # Stable column order
    cols = [
        "case",
        "N_CA",
        "plddt_min",
        "plddt_median",
        "plddt_mean",
        "plddt_max",
        "plddt_frac_lt50",
        "plddt_frac_lt70",
        "ph_rms",
        "ph_max",
        f"H_type_m{int(args.m)}_rhoA_pca_median",
    ]
    # Write CSV with blanks for missing columns.
    with out_csv.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=cols)
        w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k, "") for k in cols})

    print(f"[ok] wrote: {out_csv}")


if __name__ == "__main__":
    main()

