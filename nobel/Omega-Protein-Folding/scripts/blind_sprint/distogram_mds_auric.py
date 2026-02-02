#!/usr/bin/env python3
"""
Blind Sprint (sparse-distogram MDS) + Auric reranking.

Target signal: native-derived *distogram* on a sampled subset of residue pairs.
This is stronger than pure contacts but still does not use oracle direction guidance.

We optimize a weighted stress objective over:
- chain edges (i,i+1) at bond length
- sampled distogram edges (i,j) at native distance d_true(i,j)

Then we compute Auric metrics on the resulting candidate and rerank across restarts.
"""

from __future__ import annotations

import argparse
import math
import time
from dataclasses import dataclass
from pathlib import Path
from typing import List

import numpy as np
import pandas as pd

import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "pdb_bench"))
sys.path.insert(0, str((Path(__file__).resolve().parents[1] / "auric").resolve()))

from bench_utils import PHI, ensure_dir, load_ca_coords_longest_chain, oracle_direction_reconstruct, tm_score  # noqa: E402
from fold_m import fold_sliding_windows_bits  # noqa: E402
from metrics import metrics_for_stream  # noqa: E402
from readout import rho_A_from_yperp, rho_B_from_npath  # noqa: E402


@dataclass
class EdgeSet:
    ii: np.ndarray
    jj: np.ndarray
    dij: np.ndarray
    wij: np.ndarray


def build_edges(
    d_true: np.ndarray,
    *,
    bond_len: float,
    n_pairs: int,
    min_sep: int,
    w_chain: float,
    w_dist: float,
    rng: np.random.Generator,
) -> EdgeSet:
    N = d_true.shape[0]
    ii: List[int] = []
    jj: List[int] = []
    dij: List[float] = []
    wij: List[float] = []

    # chain edges
    for i in range(N - 1):
        ii.append(i)
        jj.append(i + 1)
        dij.append(float(bond_len))
        wij.append(float(w_chain))

    # sampled distogram edges
    # sample uniformly among pairs with |i-j|>=min_sep+1
    pairs = []
    for _ in range(int(n_pairs) * 2):  # oversample then unique
        a = int(rng.integers(0, N))
        b = int(rng.integers(0, N))
        if a == b:
            continue
        if abs(a - b) <= int(min_sep):
            continue
        if a > b:
            a, b = b, a
        pairs.append((a, b))
        if len(pairs) >= int(n_pairs) * 5:
            break
    pairs = list(dict.fromkeys(pairs))[: int(n_pairs)]
    for a, b in pairs:
        ii.append(int(a))
        jj.append(int(b))
        dij.append(float(d_true[a, b]))
        wij.append(float(w_dist))

    return EdgeSet(ii=np.asarray(ii, np.int32), jj=np.asarray(jj, np.int32), dij=np.asarray(dij, np.float64), wij=np.asarray(wij, np.float64))


def optimize_mds(edges: EdgeSet, *, N: int, seed: int, steps: int, lr: float, clamp: float) -> np.ndarray:
    rng = np.random.default_rng(int(seed))
    X = rng.normal(size=(N, 3)).astype(np.float64)
    X -= X.mean(axis=0, keepdims=True)

    ii, jj, dij, wij = edges.ii, edges.jj, edges.dij, edges.wij

    for t in range(int(steps)):
        Xi = X[ii]
        Xj = X[jj]
        d = Xi - Xj
        dist = np.sqrt(np.sum(d * d, axis=1)) + 1e-12
        gcoef = 2.0 * wij * (dist - dij) / dist
        g = d * gcoef[:, None]

        G = np.zeros_like(X)
        np.add.at(G, ii, g)
        np.add.at(G, jj, -g)

        G = np.clip(G, -float(clamp), float(clamp))
        X -= float(lr) * G
        X -= X.mean(axis=0, keepdims=True)
        if (t + 1) % 300 == 0:
            lr *= 0.85
    return X


def auric_metrics(coords: np.ndarray, *, alphabet: str, m: int) -> dict:
    _, n_path, y_perp, _, _ = oracle_direction_reconstruct(coords, alphabet=alphabet, phi=PHI)
    # candidate-local PCA uvec
    Yc = y_perp - y_perp.mean(axis=0, keepdims=True)
    if np.allclose(Yc, 0.0):
        u = np.array([1.0, 0.0, 0.0], dtype=np.float64)
    else:
        _, _, Vt = np.linalg.svd(Yc, full_matrices=False)
        u = np.asarray(Vt[0], dtype=np.float64)
        nu = float(np.linalg.norm(u))
        u = (u / nu) if nu != 0.0 else np.array([1.0, 0.0, 0.0], dtype=np.float64)
        for k in range(3):
            if abs(u[k]) > 1e-12:
                if u[k] < 0:
                    u = -u
                break

    bits_A = rho_A_from_yperp(y_perp, u=tuple(u.tolist()), u_mode="fixed", threshold="median")
    bits_B = rho_B_from_npath(n_path)
    met_A = metrics_for_stream(bits_A, fold_sliding_windows_bits(bits_A, m=m))
    met_B = metrics_for_stream(bits_B, fold_sliding_windows_bits(bits_B, m=m))
    return {
        "type_entropy_A": float(met_A["type_entropy"]),
        "smb_rate_hat_B": float(met_B["smb_rate_hat"]),
        "uvec_x": float(u[0]),
        "uvec_y": float(u[1]),
        "uvec_z": float(u[2]),
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--pdb-id", required=True)
    ap.add_argument("--tag", default="blind_sprint_distogram_mds_auric")
    ap.add_argument("--seeds", default="0,1,2,3,4")
    ap.add_argument("--steps", type=int, default=1500)
    ap.add_argument("--lr", type=float, default=0.03)
    ap.add_argument("--clamp", type=float, default=5.0)
    ap.add_argument("--min-sep", type=int, default=3)
    ap.add_argument("--n-pairs", type=int, default=800, help="Number of sampled distogram pairs.")
    ap.add_argument("--w-chain", type=float, default=5.0)
    ap.add_argument("--w-dist", type=float, default=1.0)
    ap.add_argument("--alphabet", default="pair72", choices=["pair72", "triple232"])
    ap.add_argument("--auric-m", type=int, default=8)
    ap.add_argument("--rank-w-tm", type=float, default=0.0, help="(debug) optionally include TM in rank score (keep 0 in real blind runs).")
    ap.add_argument("--rank-w-auricA", type=float, default=0.2)
    ap.add_argument("--rank-w-auricB", type=float, default=0.2)
    args = ap.parse_args()

    root = Path(__file__).resolve().parents[2]
    cif = root / "data" / "raw" / "pdb_cache" / f"{str(args.pdb_id).upper()}.cif"
    if not cif.exists():
        raise FileNotFoundError(f"Missing mmCIF: {cif}")

    chain_id, coords_native = load_ca_coords_longest_chain(cif)
    coords_native = np.asarray(coords_native, dtype=np.float64)
    N = int(coords_native.shape[0])
    bond_len = float(np.median(np.linalg.norm(coords_native[1:] - coords_native[:-1], axis=1)))
    d_true = np.linalg.norm(coords_native[:, None, :] - coords_native[None, :, :], axis=2)

    out_dir = root / "data" / "processed" / "blind_sprint"
    ensure_dir(out_dir)
    run_dir = root / "docs" / "runs" / "blind_sprint_auric_small"
    ensure_dir(run_dir)

    seeds = [int(x.strip()) for x in str(args.seeds).split(",") if x.strip()]
    rows = []
    t0 = time.perf_counter()
    for s in seeds:
        rng = np.random.default_rng(int(s))
        edges = build_edges(
            d_true,
            bond_len=bond_len,
            n_pairs=int(args.n_pairs),
            min_sep=int(args.min_sep),
            w_chain=float(args.w_chain),
            w_dist=float(args.w_dist),
            rng=rng,
        )
        X = optimize_mds(edges, N=N, seed=s, steps=int(args.steps), lr=float(args.lr), clamp=float(args.clamp))
        tm = tm_score(X, coords_native)
        aur = auric_metrics(X, alphabet=str(args.alphabet), m=int(args.auric_m))
        rank = float(args.rank_w_tm) * float(tm) - float(args.rank_w_auricA) * float(aur["type_entropy_A"]) - float(args.rank_w_auricB) * float(aur["smb_rate_hat_B"])
        rows.append(
            {
                "pdb_id": str(args.pdb_id).upper(),
                "chain": str(chain_id),
                "N": N,
                "seed": int(s),
                "tm": float(tm),
                "rank_score": float(rank),
                "bond_len": float(bond_len),
                "n_pairs": int(args.n_pairs),
                **aur,
            }
        )
        print(f"{args.pdb_id}:{chain_id} seed={s} TM={tm:.3f} rank={rank:.3f} H_A={aur['type_entropy_A']:.3f}", flush=True)

    df = pd.DataFrame(rows).sort_values("tm", ascending=False)
    out_csv = out_dir / f"blind_sprint_distogram_mds_{args.tag}_{str(args.pdb_id).upper()}.csv"
    df.to_csv(out_csv, index=False)

    dt = time.perf_counter() - t0
    report = run_dir / f"{args.tag}_{str(args.pdb_id).upper()}.md"
    cols = ["seed", "tm", "type_entropy_A", "smb_rate_hat_B", "rank_score"]
    lines = []
    lines.append(f"# Blind Sprint (sparse-distogram MDS) + Auric: {str(args.pdb_id).upper()}:{chain_id}")
    lines.append("")
    lines.append(f"- N: {N}")
    lines.append(f"- distogram pairs: {int(args.n_pairs)} (min_sep={int(args.min_sep)})")
    lines.append(f"- MDS: steps={int(args.steps)}, lr={float(args.lr)}, clamp={float(args.clamp)}")
    lines.append(f"- weights: w_chain={float(args.w_chain)}, w_dist={float(args.w_dist)}")
    lines.append(f"- Auric: alphabet={args.alphabet}, m={int(args.auric_m)}")
    lines.append(f"- CSV (gitignored): `{out_csv.relative_to(root).as_posix()}`")
    lines.append(f"- Runtime: {dt:.1f}s")
    lines.append("")
    lines.append("## Results (sorted by TM for evaluation)")
    lines.append("")
    lines.append("| " + " | ".join(cols) + " |")
    lines.append("| " + " | ".join(["---"] * len(cols)) + " |")
    for _, r in df.iterrows():
        lines.append(
            "| "
            + " | ".join(
                [
                    str(int(r["seed"])),
                    f"{float(r['tm']):.3f}",
                    f"{float(r['type_entropy_A']):.3f}",
                    f"{float(r['smb_rate_hat_B']):.3f}",
                    f"{float(r['rank_score']):.3f}",
                ]
            )
            + " |"
        )
    lines.append("")
    report.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"Wrote: {out_csv}")
    print(f"Wrote: {report}")


if __name__ == "__main__":
    main()

