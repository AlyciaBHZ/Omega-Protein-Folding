#!/usr/bin/env python3
"""
Large-scale "oracle codec" benchmark on real PDB structures.

This implements the Step-2 idea described in
`artifacts/reports/step1_theta_star_step2_triple232_evidence_pack.md`:

  - Treat backbone as bond vectors d_t = x_{t+1} - x_t
  - For a given 6D alphabet (axis12/pair72/triple232), choose the closest 6D step
    (oracle) to match d_t after 6D->3D projection via an icosahedral basis.
  - Integrate steps to get a reconstructed backbone, then Kabsch-align and compute TM-score.

This is an expressivity *upper bound* (not a folding search).

Inputs:
  - data/raw/pdb_cache/*.cif  (downloaded by pdb_sample_and_download.py)
  - data/processed/pdb_bench/pdb_ids_*.txt

Outputs:
  - data/processed/pdb_bench/oracle_codec_bench_*.csv
  - artifacts/reports/pdb_oracle_codec_bench_*.md
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

import numpy as np
import pandas as pd
from Bio.PDB.MMCIFParser import MMCIFParser


PHI = (1.0 + 5.0**0.5) / 2.0


def ensure_dir(p: Path) -> None:
    p.mkdir(parents=True, exist_ok=True)


def tm_d0(L: int) -> float:
    # Standard TM-score d0 definition (Zhang & Skolnick).
    if L <= 0:
        return 1.0
    d0 = 1.24 * (L - 15) ** (1.0 / 3.0) - 1.8
    return float(max(d0, 0.5))


def kabsch(P: np.ndarray, Q: np.ndarray) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Return aligned P_aligned, rotation R, translation t such that:
      P_aligned ≈ P @ R + t  matches Q (least-squares).
    """
    assert P.shape == Q.shape and P.shape[1] == 3
    Pc = P - P.mean(axis=0, keepdims=True)
    Qc = Q - Q.mean(axis=0, keepdims=True)
    C = Pc.T @ Qc
    V, S, Wt = np.linalg.svd(C)
    d = np.sign(np.linalg.det(V @ Wt))
    D = np.diag([1.0, 1.0, d])
    R = V @ D @ Wt
    t = Q.mean(axis=0) - P.mean(axis=0) @ R
    P_aligned = P @ R + t
    return P_aligned, R, t


def tm_score(P: np.ndarray, Q: np.ndarray) -> float:
    """
    TM-score after optimal Kabsch alignment (single-pass).
    """
    P_aligned, _, _ = kabsch(P, Q)
    d = np.linalg.norm(P_aligned - Q, axis=1)
    d0 = tm_d0(len(d))
    return float(np.mean(1.0 / (1.0 + (d / d0) ** 2)))


def best_local_tm(P: np.ndarray, Q: np.ndarray, window: int = 100) -> float:
    if len(P) < window:
        return float("nan")
    best = -1.0
    for i in range(0, len(P) - window + 1):
        best = max(best, tm_score(P[i : i + window], Q[i : i + window]))
    return float(best)


def icosa_B(phi: float) -> np.ndarray:
    """
    Build a 3x6 icosahedral projection basis B (columns correspond to +e_i).
    """
    cols = [
        np.array([0.0, 1.0, phi]),
        np.array([0.0, 1.0, -phi]),
        np.array([1.0, phi, 0.0]),
        np.array([1.0, -phi, 0.0]),
        np.array([phi, 0.0, 1.0]),
        np.array([phi, 0.0, -1.0]),
    ]
    cols = [c / np.linalg.norm(c) for c in cols]
    return np.stack(cols, axis=1)  # (3,6)


def alphabet_steps(name: str) -> np.ndarray:
    """
    Return (K,6) int steps for axis12/pair72/triple232.
    """
    if name not in {"axis12", "pair72", "triple232"}:
        raise ValueError(f"Unknown alphabet: {name}")

    steps: List[np.ndarray] = []
    # axis steps
    for i in range(6):
        e = np.zeros(6, dtype=np.int32)
        e[i] = 1
        steps.append(e.copy())
        steps.append((-e).copy())

    if name in {"pair72", "triple232"}:
        for i in range(6):
            for j in range(i + 1, 6):
                for si in (-1, 1):
                    for sj in (-1, 1):
                        v = np.zeros(6, dtype=np.int32)
                        v[i] = si
                        v[j] = sj
                        steps.append(v)

    if name == "triple232":
        for i in range(6):
            for j in range(i + 1, 6):
                for k in range(j + 1, 6):
                    for si in (-1, 1):
                        for sj in (-1, 1):
                            for sk in (-1, 1):
                                v = np.zeros(6, dtype=np.int32)
                                v[i] = si
                                v[j] = sj
                                v[k] = sk
                                steps.append(v)

    out = np.stack(steps, axis=0)
    # Sanity counts
    expected = {"axis12": 12, "pair72": 72, "triple232": 232}[name]
    if out.shape[0] != expected:
        raise RuntimeError(f"Alphabet {name} size mismatch: got {out.shape[0]}, expected {expected}")
    return out


def load_ca_coords_longest_chain(cif_path: Path) -> Tuple[str, np.ndarray]:
    """
    Return (chain_id, coords) where coords is (N,3) CA trace for the longest chain.
    """
    parser = MMCIFParser(QUIET=True)
    structure_id = cif_path.stem
    s = parser.get_structure(structure_id, str(cif_path))

    # Use first model only
    model = next(iter(s.get_models()))

    best_chain = None
    best_coords = None
    for chain in model:
        coords = []
        for res in chain:
            if "CA" in res:
                coords.append(res["CA"].get_coord())
        if len(coords) >= 2:
            coords_arr = np.asarray(coords, dtype=np.float64)
            if best_coords is None or len(coords_arr) > len(best_coords):
                best_coords = coords_arr
                best_chain = chain.id

    if best_coords is None or best_chain is None:
        raise ValueError("No CA trace found")
    return str(best_chain), best_coords


def oracle_reconstruct(
    coords_native: np.ndarray,
    *,
    alphabet: str,
    phi: float = PHI,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, float]:
    """
    Return (coords_pred, n_path, y_perp, ph_rms) where:
      - coords_pred: (N,3) reconstructed 3D coords (before Kabsch)
      - n_path: (N,6) integer walk positions
      - y_perp: (N,3) perp coords
      - ph_rms: RMS(||y_perp - mean(y_perp)||)
    """
    N = coords_native.shape[0]
    d = coords_native[1:] - coords_native[:-1]  # (N-1,3)
    bond_lengths = np.linalg.norm(d, axis=1)
    median_len = float(np.median(bond_lengths))

    Bpar = icosa_B(phi) * median_len  # scale so unit step ≈ median CA-CA bond length
    Bperp = icosa_B(-1.0 / phi) * median_len

    steps = alphabet_steps(alphabet)  # (K,6)
    proj = (Bpar @ steps.T).T  # (K,3)
    proj_norm = np.linalg.norm(proj, axis=1)
    # Avoid divide-by-zero; should not happen for these alphabets.
    proj_dir = proj / np.clip(proj_norm[:, None], 1e-12, None)

    coords_pred = np.zeros((N, 3), dtype=np.float64)
    n_path = np.zeros((N, 6), dtype=np.int32)

    # Greedy oracle: per-bond choose closest *direction* from alphabet.
    # We scale the chosen direction to the observed bond length (CA bonds are ~constant),
    # so multi-axis steps add angular expressivity without forcing longer bond lengths.
    for t in range(N - 1):
        dt = d[t]
        dt_len = float(np.linalg.norm(dt))
        if dt_len < 1e-8:
            # Degenerate (missing coords); skip.
            k = 0
            step_vec = proj_dir[k] * 0.0
        else:
            dt_dir = dt / dt_len
            # Maximize cosine similarity
            cos = proj_dir @ dt_dir
            k = int(np.argmax(cos))
            step_vec = proj_dir[k] * dt_len
        step = steps[k]
        coords_pred[t + 1] = coords_pred[t] + step_vec
        n_path[t + 1] = n_path[t] + step

    y_perp = (Bperp @ n_path.T).T  # (N,3)
    w0 = y_perp.mean(axis=0, keepdims=True)
    ph_rms = float(np.sqrt(np.mean(np.sum((y_perp - w0) ** 2, axis=1))))
    return coords_pred, n_path, y_perp, ph_rms


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--ids", type=str, required=True, help="Path to pdb_ids_*.txt (one PDB ID per line).")
    ap.add_argument("--min-len", type=int, default=60)
    ap.add_argument("--max-len", type=int, default=350)
    ap.add_argument("--out-tag", type=str, default="v1")
    args = ap.parse_args()

    root = Path(__file__).resolve().parents[2]
    ids_path = Path(args.ids)
    if not ids_path.is_absolute():
        ids_path = (root / ids_path).resolve()

    cache_dir = root / "data" / "raw" / "pdb_cache"
    out_dir = root / "data" / "processed" / "pdb_bench"
    rep_dir = root / "artifacts" / "reports"
    ensure_dir(out_dir)
    ensure_dir(rep_dir)

    ids = [x.strip().upper() for x in ids_path.read_text(encoding="utf-8").splitlines() if x.strip()]
    rows: List[Dict[str, object]] = []

    for pdb_id in ids:
        cif = cache_dir / f"{pdb_id}.cif"
        if not cif.exists():
            continue
        try:
            chain_id, coords = load_ca_coords_longest_chain(cif)
        except Exception:
            continue

        N = len(coords)
        if N < args.min_len or N > args.max_len:
            continue

        recs = {}
        for alpha in ("axis12", "pair72", "triple232"):
            coords_pred, n_path, y_perp, ph_rms = oracle_reconstruct(coords, alphabet=alpha)
            tm = tm_score(coords_pred, coords)
            tm_local100 = best_local_tm(coords_pred, coords, window=100)
            recs[alpha] = (tm, tm_local100, ph_rms)

        rows.append(
            {
                "pdb_id": pdb_id,
                "chain": chain_id,
                "N": N,
                "tm_axis12": recs["axis12"][0],
                "tm_pair72": recs["pair72"][0],
                "tm_triple232": recs["triple232"][0],
                "tm_local100_axis12": recs["axis12"][1],
                "tm_local100_pair72": recs["pair72"][1],
                "tm_local100_triple232": recs["triple232"][1],
                "ph_rms_axis12": recs["axis12"][2],
                "ph_rms_pair72": recs["pair72"][2],
                "ph_rms_triple232": recs["triple232"][2],
            }
        )

    df = pd.DataFrame(rows)
    out_csv = out_dir / f"oracle_codec_bench_{args.out_tag}.csv"
    df.to_csv(out_csv, index=False)

    # Minimal report
    report_path = rep_dir / f"pdb_oracle_codec_bench_{args.out_tag}.md"
    ids_rel = str(ids_path.relative_to(root)).replace("\\", "/")
    summary = {
        "n_total_ids": len(ids),
        "n_used": int(len(df)),
        "min_len": args.min_len,
        "max_len": args.max_len,
        "csv": str(out_csv.relative_to(root)).replace("\\", "/"),
    }
    report_path.write_text(
        "\n".join(
            [
                "# PDB oracle codec benchmark (large-scale)",
                "",
                "This is an expressivity-only benchmark (oracle per-bond step choice + Kabsch + TM-score).",
                "",
                "## Dataset",
                f"- IDs file: `{ids_rel}`",
                f"- Used chains: {summary['n_used']} (filtered by N in [{args.min_len},{args.max_len}])",
                "",
                "## Outputs",
                f"- CSV: `{summary['csv']}`",
                f"- Suggested plots (run `scripts/pdb_bench/make_plots.py`):",
                f"  - `figures/pdb_bench/oracle_codec_tm_hist_{args.out_tag}.png`",
                f"  - `figures/pdb_bench/oracle_codec_tm_ecdf_{args.out_tag}.png`",
                "",
                "## Quick stats (median TM)",
                f"- axis12: {df['tm_axis12'].median():.3f}  (n={len(df)})",
                f"- pair72: {df['tm_pair72'].median():.3f}  (n={len(df)})",
                f"- triple232: {df['tm_triple232'].median():.3f}  (n={len(df)})",
                "",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    print(f"Wrote: {out_csv}")
    print(f"Wrote: {report_path}")


if __name__ == "__main__":
    main()

