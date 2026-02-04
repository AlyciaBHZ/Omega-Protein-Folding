#!/usr/bin/env python3
from __future__ import annotations

"""
Audit a decoy dataset via a manifest CSV (e.g., Decoys 'R' Us 4state_reduced).

Manifest format expected (as in docs/runs/decoys_4state_auric/decoy_manifest.csv):
  dataset,target_id,native_path,decoy_path
where paths are repo-root-relative (often into data/raw/decoy_cache/ which is local-only).

Outputs (committable):
  docs/runs/physics_audit/physics_vs_auric_<name>_{sampled|full}.csv

We also compute a simple normalized physics score:
  physics_score = (clash_count / N) + bond_mse * bond_scale
where bond_scale defaults to 1000.
"""

import argparse
import csv
import sys
from pathlib import Path
from typing import Dict, Iterable, List, Tuple

import numpy as np
import pandas as pd


def repo_root() -> Path:
    return Path(__file__).resolve().parents[2]


def ensure_dir(p: Path) -> None:
    p.mkdir(parents=True, exist_ok=True)


def load_ca_coords_fast(pdb_path: Path) -> np.ndarray:
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


def _radius_of_gyration(coords: np.ndarray) -> float:
    x = np.asarray(coords, dtype=np.float64)
    if x.ndim != 2 or x.shape[0] == 0:
        return float("nan")
    c = x.mean(axis=0, keepdims=True)
    return float(np.sqrt(np.mean(np.sum((x - c) ** 2, axis=1))))


def _bond_mse(coords: np.ndarray, *, ideal: float = 3.8) -> float:
    x = np.asarray(coords, dtype=np.float64)
    if x.ndim != 2 or x.shape[0] < 2:
        return float("nan")
    d = np.linalg.norm(x[1:] - x[:-1], axis=1)
    return float(np.mean((d - float(ideal)) ** 2))


def _clash_count(coords: np.ndarray, *, threshold: float = 3.5, min_sep: int = 1) -> float:
    x = np.asarray(coords, dtype=np.float64)
    N = int(x.shape[0])
    if x.ndim != 2 or N < 4:
        return 0.0

    diff = x[:, None, :] - x[None, :, :]
    dist2 = np.sum(diff * diff, axis=-1)

    idx = np.arange(N, dtype=np.int32)
    sep = np.abs(idx[:, None] - idx[None, :])
    mask = sep > int(min_sep)

    thr2 = float(threshold) ** 2
    clash_mat = (dist2 < thr2) & mask
    return float(np.sum(np.triu(clash_mat, k=1)))


def calc_physics_from_coords(coords: np.ndarray) -> Dict[str, float]:
    N = int(coords.shape[0])
    if N < 3:
        return {}
    c = _clash_count(coords)
    return {
        "N": float(N),
        "bond_mse": float(_bond_mse(coords)),
        "clash_count": float(c),
        "clash_per_res": float(c / max(1, N)),
        "rg": float(_radius_of_gyration(coords)),
    }


def _auric_imports() -> Tuple[object, object, object, object, object]:
    """
    Import non-packaged modules via sys.path inserts (repo scripts are not a Python package).
    Returns: oracle_direction_reconstruct, shared_pca_axis, rho_A_from_yperp, fold_sliding_windows_bits, metrics_for_stream
    """
    root = repo_root()
    sys.path.insert(0, str((root / "scripts" / "pdb_bench").resolve()))
    sys.path.insert(0, str((root / "scripts" / "auric").resolve()))

    from bench_utils import oracle_direction_reconstruct  # type: ignore  # noqa: E402
    from axis import shared_pca_axis  # type: ignore  # noqa: E402
    from fold_m import fold_sliding_windows_bits  # type: ignore  # noqa: E402
    from metrics import metrics_for_stream  # type: ignore  # noqa: E402
    from readout import rho_A_from_yperp  # type: ignore  # noqa: E402

    return oracle_direction_reconstruct, shared_pca_axis, rho_A_from_yperp, fold_sliding_windows_bits, metrics_for_stream


def _iter_sample(items: List[Path], *, k: int, seed: int) -> Iterable[Path]:
    if k <= 0 or k >= len(items):
        yield from items
        return
    rng = np.random.default_rng(int(seed))
    idx = rng.choice(len(items), size=int(k), replace=False)
    for i in sorted(idx.tolist()):
        yield items[int(i)]


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--manifest-csv", default="docs/runs/decoys_4state_auric/decoy_manifest.csv")
    ap.add_argument("--alphabet", default="triple232", choices=["axis12", "pair72", "triple232"])
    ap.add_argument("--m", type=int, default=8)
    ap.add_argument(
        "--max-decoys-per-target",
        type=int,
        default=200,
        help="Per-target cap for decoys. Use 0 (or negative) to audit ALL decoys.",
    )
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--bond-scale", type=float, default=1000.0)
    ap.add_argument("--out-csv", default="docs/runs/physics_audit/physics_vs_auric_4state_full.csv")
    args = ap.parse_args()

    root = repo_root()
    man = Path(str(args.manifest_csv))
    if not man.is_absolute():
        man = (root / man).resolve()
    if not man.exists():
        raise SystemExit(f"Manifest not found: {man}")

    df = pd.read_csv(man)
    needed = {"dataset", "target_id", "native_path", "decoy_path"}
    if not needed.issubset(set(df.columns)):
        raise SystemExit(f"Manifest missing columns. Need {sorted(needed)}; got {list(df.columns)}")

    oracle_direction_reconstruct, shared_pca_axis, rho_A_from_yperp, fold_sliding_windows_bits, metrics_for_stream = (
        _auric_imports()
    )

    out_csv = Path(str(args.out_csv))
    if not out_csv.is_absolute():
        out_csv = (root / out_csv).resolve()
    ensure_dir(out_csv.parent)

    rows: List[Dict[str, object]] = []
    grouped = df.groupby(["dataset", "target_id"], sort=True)

    for (dataset, target_id), g in grouped:
        native_rel = str(g["native_path"].iloc[0])
        native_abs = (root / native_rel).resolve()
        if not native_abs.exists():
            print(f"[skip] missing native: {dataset}/{target_id} -> {native_abs}")
            continue

        coords_nat = load_ca_coords_fast(native_abs)
        if len(coords_nat) < 3:
            print(f"[skip] native too short: {dataset}/{target_id}")
            continue

        # Shared axis derived from native y_perp (protocol-fixed; native-only).
        _, _n_path, y_perp_nat, _phrms, _phmax = oracle_direction_reconstruct(
            np.asarray(coords_nat, dtype=np.float64),
            alphabet=str(args.alphabet),
        )
        u = shared_pca_axis(y_perp_nat)

        def auric_from_coords(coords: np.ndarray) -> Dict[str, float]:
            _, _n_path2, y_perp, ph_rms, ph_max = oracle_direction_reconstruct(
                np.asarray(coords, dtype=np.float64),
                alphabet=str(args.alphabet),
            )
            bits_A = rho_A_from_yperp(y_perp, u=tuple(u.tolist()), u_mode="fixed", threshold="median")
            folded_A = fold_sliding_windows_bits(bits_A, m=int(args.m))
            met_A = metrics_for_stream(bits_A, folded_A)
            return {
                "ph_rms": float(ph_rms),
                "ph_max": float(ph_max),
                "auric_type_entropy_rhoA": float(met_A.get("type_entropy", float("nan"))),
                "auric_smb_rate_hat_rhoA": float(met_A.get("smb_rate_hat", float("nan"))),
            }

        # Native row
        phys_nat = calc_physics_from_coords(coords_nat)
        aur_nat = auric_from_coords(coords_nat)
        physics_score_nat = float(phys_nat.get("clash_per_res", float("nan"))) + float(args.bond_scale) * float(
            phys_nat.get("bond_mse", float("nan"))
        )
        rows.append(
            {
                "dataset": dataset,
                "target_id": target_id,
                "label": "native",
                "pdb_path": Path(native_rel).as_posix(),
                "filename": native_abs.name,
                **{k: phys_nat.get(k, "") for k in ["N", "bond_mse", "clash_count", "clash_per_res", "rg"]},
                **aur_nat,
                "physics_score": physics_score_nat,
            }
        )

        # Decoys (full or sampled)
        decoy_paths: List[Path] = []
        for rel in g["decoy_path"].astype(str).tolist():
            p = (root / rel).resolve()
            if p.exists():
                decoy_paths.append(p)
        if not decoy_paths:
            print(f"[warn] no decoys found on disk: {dataset}/{target_id}")
            continue

        k = int(args.max_decoys_per_target)
        for p in _iter_sample(decoy_paths, k=k, seed=int(args.seed) + int(abs(hash((dataset, target_id))) % 10_000)):
            coords = load_ca_coords_fast(p)
            if len(coords) < 3:
                continue
            phys = calc_physics_from_coords(coords)
            aur = auric_from_coords(coords)
            physics_score = float(phys.get("clash_per_res", float("nan"))) + float(args.bond_scale) * float(
                phys.get("bond_mse", float("nan"))
            )
            try:
                rel = p.relative_to(root).as_posix()
            except Exception:
                rel = str(p)
            rows.append(
                {
                    "dataset": dataset,
                    "target_id": target_id,
                    "label": "decoy",
                    "pdb_path": rel,
                    "filename": p.name,
                    **{k2: phys.get(k2, "") for k2 in ["N", "bond_mse", "clash_count", "clash_per_res", "rg"]},
                    **aur,
                    "physics_score": physics_score,
                }
            )

        n_dec = len(decoy_paths) if k <= 0 else min(k, len(decoy_paths))
        print(f"[ok] {dataset}/{target_id}: native + {n_dec} decoys")

    cols = [
        "dataset",
        "target_id",
        "label",
        "pdb_path",
        "filename",
        "N",
        "bond_mse",
        "clash_count",
        "clash_per_res",
        "rg",
        "physics_score",
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

