#!/usr/bin/env python3
"""
Blind Sprint (contact-MDS) + Auric reranking.

This runner is intentionally simple and robust:
- Uses only a native-derived contact map as the target signal (allowed).
- Generates candidate 3D chains by optimizing a weighted stress objective over:
  - chain edges (i,i+1) at bond length
  - contact edges (i,j) at a target contact distance
- Uses Auric metrics (computed from the candidate coords) as an additional certificate
  and as a reranking signal across random restarts.

This is not a full Omega engine; it is a minimal "blind sprint" feasibility demo.
"""

from __future__ import annotations

import argparse
import math
import time
from dataclasses import dataclass
from pathlib import Path
from typing import List, Tuple

import numpy as np
import pandas as pd

import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "pdb_bench"))
sys.path.insert(0, str((Path(__file__).resolve().parents[1] / "auric").resolve()))

from bench_utils import PHI, ensure_dir, load_ca_coords_longest_chain, oracle_direction_reconstruct, tm_score  # noqa: E402
from fold_m import fold_sliding_windows_bits  # noqa: E402
from metrics import metrics_for_stream  # noqa: E402
from readout import rho_A_from_yperp, rho_B_from_npath  # noqa: E402


def contact_matrix(coords: np.ndarray, *, cutoff: float, min_sep: int) -> np.ndarray:
    coords = np.asarray(coords, dtype=np.float64)
    N = coords.shape[0]
    d = coords[:, None, :] - coords[None, :, :]
    dist = np.sqrt(np.sum(d * d, axis=2))
    C = dist <= float(cutoff)
    for k in range(int(min_sep) + 1):
        idx = np.arange(N - k)
        C[idx, idx + k] = False
        C[idx + k, idx] = False
    return C


def contact_f1(C_true: np.ndarray, coords_pred: np.ndarray, *, cutoff: float, min_sep: int) -> float:
    C_pred = contact_matrix(coords_pred, cutoff=cutoff, min_sep=min_sep)
    # upper triangle pairs only
    triu = np.triu_indices(C_true.shape[0], k=min_sep + 1)
    t = C_true[triu]
    p = C_pred[triu]
    tp = int(np.sum(t & p))
    fp = int(np.sum((~t) & p))
    fn = int(np.sum(t & (~p)))
    if tp <= 0:
        return 0.0
    prec = tp / max(1, tp + fp)
    rec = tp / max(1, tp + fn)
    return float(2 * prec * rec / max(1e-12, (prec + rec)))


@dataclass
class EdgeSet:
    ii: np.ndarray  # (E,)
    jj: np.ndarray  # (E,)
    dij: np.ndarray  # (E,)
    wij: np.ndarray  # (E,)


def build_edges(C_true: np.ndarray, *, bond_len: float, contact_dist: float, w_chain: float, w_contact: float) -> EdgeSet:
    N = C_true.shape[0]
    ii = []
    jj = []
    dij = []
    wij = []
    # chain edges
    for i in range(N - 1):
        ii.append(i)
        jj.append(i + 1)
        dij.append(float(bond_len))
        wij.append(float(w_chain))
    # contact edges
    ci, cj = np.where(np.triu(C_true, k=1))
    for a, b in zip(ci.tolist(), cj.tolist()):
        ii.append(int(a))
        jj.append(int(b))
        dij.append(float(contact_dist))
        wij.append(float(w_contact))
    return EdgeSet(ii=np.asarray(ii, dtype=np.int32), jj=np.asarray(jj, dtype=np.int32), dij=np.asarray(dij, dtype=np.float64), wij=np.asarray(wij, dtype=np.float64))


def optimize_mds(
    edges: EdgeSet,
    *,
    N: int,
    seed: int,
    steps: int,
    lr: float,
    clamp: float,
) -> np.ndarray:
    rng = np.random.default_rng(int(seed))
    X = rng.normal(size=(N, 3)).astype(np.float64)
    X -= X.mean(axis=0, keepdims=True)

    ii = edges.ii
    jj = edges.jj
    dij = edges.dij
    wij = edges.wij

    for t in range(int(steps)):
        Xi = X[ii]
        Xj = X[jj]
        d = Xi - Xj
        dist = np.sqrt(np.sum(d * d, axis=1)) + 1e-12
        # gradient of wij*(dist-dij)^2 wrt Xi is 2*wij*(dist-dij)*d/dist
        gcoef = 2.0 * wij * (dist - dij) / dist
        g = d * gcoef[:, None]

        G = np.zeros_like(X)
        np.add.at(G, ii, g)
        np.add.at(G, jj, -g)

        # step
        G = np.clip(G, -float(clamp), float(clamp))
        X -= float(lr) * G
        X -= X.mean(axis=0, keepdims=True)

        # mild lr decay
        if (t + 1) % 200 == 0:
            lr *= 0.8

    return X


def auric_score(coords: np.ndarray, *, alphabet: str, m: int) -> dict:
    _, n_path, y_perp, _, _ = oracle_direction_reconstruct(coords, alphabet=alphabet, phi=PHI)
    # Shared-PCA uvec from this candidate (no oracle); deterministic SVD
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

    folded_A = fold_sliding_windows_bits(bits_A, m=m)
    met_A = metrics_for_stream(bits_A, folded_A)
    folded_B = fold_sliding_windows_bits(bits_B, m=m)
    met_B = metrics_for_stream(bits_B, folded_B)

    return {
        "type_entropy_A": float(met_A["type_entropy"]),
        "smb_rate_hat_A": float(met_A["smb_rate_hat"]),
        "type_entropy_B": float(met_B["type_entropy"]),
        "smb_rate_hat_B": float(met_B["smb_rate_hat"]),
        "uvec_x": float(u[0]),
        "uvec_y": float(u[1]),
        "uvec_z": float(u[2]),
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--pdb-id", required=True)
    ap.add_argument("--tag", default="blind_sprint_contact_mds")
    ap.add_argument("--seeds", default="0,1,2,3,4")
    ap.add_argument("--steps", type=int, default=1200)
    ap.add_argument("--lr", type=float, default=0.03)
    ap.add_argument("--clamp", type=float, default=5.0)
    ap.add_argument("--cutoff", type=float, default=8.0)
    ap.add_argument("--min-sep", type=int, default=3)
    ap.add_argument("--contact-dist", type=float, default=6.0)
    ap.add_argument("--w-chain", type=float, default=5.0)
    ap.add_argument("--w-contact", type=float, default=1.0)
    ap.add_argument("--alphabet", default="pair72", choices=["pair72", "triple232"])
    ap.add_argument("--auric-m", type=int, default=8)
    ap.add_argument("--rank-w-contact", type=float, default=1.0)
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
    C_true = contact_matrix(coords_native, cutoff=float(args.cutoff), min_sep=int(args.min_sep))
    edges = build_edges(C_true, bond_len=bond_len, contact_dist=float(args.contact_dist), w_chain=float(args.w_chain), w_contact=float(args.w_contact))

    out_dir = root / "data" / "processed" / "blind_sprint"
    ensure_dir(out_dir)
    run_dir = root / "docs" / "runs" / "blind_sprint_auric_small"
    ensure_dir(run_dir)

    seeds = [int(x.strip()) for x in str(args.seeds).split(",") if x.strip()]
    rows = []
    t0 = time.perf_counter()
    for s in seeds:
        X = optimize_mds(edges, N=N, seed=s, steps=int(args.steps), lr=float(args.lr), clamp=float(args.clamp))
        f1 = contact_f1(C_true, X, cutoff=float(args.cutoff), min_sep=int(args.min_sep))
        tm = tm_score(X, coords_native)
        aur = auric_score(X, alphabet=str(args.alphabet), m=int(args.auric_m))
        score = float(args.rank_w_contact) * float(f1) - float(args.rank_w_auricA) * float(aur["type_entropy_A"]) - float(args.rank_w_auricB) * float(aur["smb_rate_hat_B"])
        row = {
            "pdb_id": str(args.pdb_id).upper(),
            "chain": str(chain_id),
            "N": N,
            "seed": int(s),
            "contact_f1": float(f1),
            "tm": float(tm),
            "rank_score": float(score),
            "bond_len": float(bond_len),
            **aur,
        }
        rows.append(row)
        print(f"{args.pdb_id}:{chain_id} seed={s} TM={tm:.3f} F1={f1:.3f} rank={score:.3f} H_A={aur['type_entropy_A']:.3f}", flush=True)

    df = pd.DataFrame(rows).sort_values("rank_score", ascending=False)
    out_csv = out_dir / f"blind_sprint_contact_mds_{args.tag}_{str(args.pdb_id).upper()}.csv"
    df.to_csv(out_csv, index=False)

    dt = time.perf_counter() - t0
    report = run_dir / f"{args.tag}_{str(args.pdb_id).upper()}.md"
    # small manual table
    cols = ["seed", "tm", "contact_f1", "type_entropy_A", "smb_rate_hat_B", "rank_score"]
    lines = []
    lines.append(f"# Blind Sprint (contact-MDS) + Auric rerank: {str(args.pdb_id).upper()}:{chain_id}")
    lines.append("")
    lines.append(f"- N: {N}")
    lines.append(f"- seeds: {str(args.seeds)}")
    lines.append(f"- contact cutoff: {float(args.cutoff)} Å (min_sep={int(args.min_sep)})")
    lines.append(f"- MDS: steps={int(args.steps)}, lr={float(args.lr)}, clamp={float(args.clamp)}")
    lines.append(f"- weights: w_chain={float(args.w_chain)}, w_contact={float(args.w_contact)}, contact_dist={float(args.contact_dist)}")
    lines.append(f"- Auric: alphabet={args.alphabet}, m={int(args.auric_m)}")
    lines.append(f"- CSV (gitignored): `{out_csv.relative_to(root).as_posix()}`")
    lines.append(f"- Runtime: {dt:.1f}s")
    lines.append("")
    lines.append("## Results (sorted by rank_score)")
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
                    f"{float(r['contact_f1']):.3f}",
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

