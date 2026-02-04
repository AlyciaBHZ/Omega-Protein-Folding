#!/usr/bin/env python3
from __future__ import annotations

"""
Lightweight CA-level physics audit metrics (NumPy-only).

This is intentionally a coarse-grained proxy layer (no Rosetta, no MD, no full stereochem).
It provides fast post-hoc checks that catch gross non-physical geometry.

Metrics (from CA trace):
- bond_mse: mean((||CA_i-CA_{i-1}|| - 3.8)^2)
- clash_count: number of non-neighbor CA pairs closer than a threshold
- clash_per_res: clash_count / N
- rg: radius of gyration
"""

import argparse
from pathlib import Path
from typing import Dict

import numpy as np
from Bio.PDB.PDBParser import PDBParser


def load_ca_coords_pdb(pdb_path: Path) -> np.ndarray:
    p = Path(pdb_path)
    parser = PDBParser(QUIET=True)
    s = parser.get_structure(p.stem, str(p))
    coords = []
    for atom in s.get_atoms():
        if atom.get_name() == "CA":
            coords.append(atom.get_coord())
    if len(coords) < 3:
        return np.zeros((0, 3), dtype=np.float64)
    return np.asarray(coords, dtype=np.float64)


def radius_of_gyration(coords: np.ndarray) -> float:
    x = np.asarray(coords, dtype=np.float64)
    if x.ndim != 2 or x.shape[0] == 0:
        return float("nan")
    c = x.mean(axis=0, keepdims=True)
    return float(np.sqrt(np.mean(np.sum((x - c) ** 2, axis=1))))


def bond_mse(coords: np.ndarray, *, ideal: float = 3.8) -> float:
    x = np.asarray(coords, dtype=np.float64)
    if x.ndim != 2 or x.shape[0] < 2:
        return float("nan")
    d = np.linalg.norm(x[1:] - x[:-1], axis=1)
    return float(np.mean((d - float(ideal)) ** 2))


def clash_count(coords: np.ndarray, *, threshold: float = 3.5, min_sep: int = 1) -> float:
    """
    Count CA-CA clashes for residue pairs with |i-j| > min_sep.
    Uses pure NumPy broadcasting; O(N^2), fine for N<=~350.
    """
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


def calc_physics_metrics(
    pdb_path: Path, *, ideal_bond: float = 3.8, clash_threshold: float = 3.5, clash_min_sep: int = 1
) -> Dict[str, float]:
    coords = load_ca_coords_pdb(Path(pdb_path))
    N = int(coords.shape[0])
    if N < 3:
        return {}

    b = bond_mse(coords, ideal=float(ideal_bond))
    c = clash_count(coords, threshold=float(clash_threshold), min_sep=int(clash_min_sep))
    rg = radius_of_gyration(coords)
    return {
        "N": float(N),
        "bond_mse": float(b),
        "clash_count": float(c),
        "clash_per_res": float(c / max(1, N)),
        "rg": float(rg),
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("pdb_path")
    ap.add_argument("--ideal-bond", type=float, default=3.8)
    ap.add_argument("--clash-threshold", type=float, default=3.5)
    ap.add_argument("--clash-min-sep", type=int, default=1)
    args = ap.parse_args()

    out = calc_physics_metrics(
        Path(args.pdb_path),
        ideal_bond=float(args.ideal_bond),
        clash_threshold=float(args.clash_threshold),
        clash_min_sep=int(args.clash_min_sep),
    )
    print(out)


if __name__ == "__main__":
    main()

