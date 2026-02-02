#!/usr/bin/env python3
"""
Small Blind Sprint: search with Contact objective + Auric scoring (no oracle direction guidance).

This is a self-contained runner inside this repo (does NOT require an external engine).

Rules enforced:
- We do NOT use oracle direction quantization for proposing moves.
- We DO allow a native-derived contact map as the target signal (as confirmed).
- TM-score is computed for evaluation only at the end.

Search model:
- Build a chain of length N using a discrete 6D alphabet (pair72 or triple232),
  but with constant 3D bond length (scaled to the native median CA-CA bond length).
- Use a beam search over partial chains; the incremental objective is based on newly-formed
  contacts involving the newly appended residue.
- Optionally add Auric penalties every K steps from the current (y_perp, n_path).

Outputs:
- data/processed/blind_sprint/blind_sprint_auric_runs_<tag>.csv (gitignored)
- docs/runs/blind_sprint_auric_small/<tag>_report.md (committable)
"""

from __future__ import annotations

import argparse
import math
import time
from dataclasses import dataclass
from pathlib import Path
from typing import List, Optional, Tuple

import numpy as np
import pandas as pd

import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "pdb_bench"))
sys.path.insert(0, str((Path(__file__).resolve().parents[1] / "auric").resolve()))

from bench_utils import (  # noqa: E402
    PHI,
    alphabet_steps,
    ensure_dir,
    icosa_B,
    load_ca_coords_longest_chain,
    tm_score,
)

from fold_m import fold_sliding_windows_bits  # noqa: E402
from metrics import metrics_for_stream  # noqa: E402
from readout import rho_A_from_yperp, rho_B_from_npath  # noqa: E402


def contact_matrix(coords: np.ndarray, *, cutoff: float, min_sep: int) -> np.ndarray:
    coords = np.asarray(coords, dtype=np.float64)
    N = coords.shape[0]
    d = coords[:, None, :] - coords[None, :, :]
    dist = np.sqrt(np.sum(d * d, axis=2))
    C = dist <= float(cutoff)
    # remove diagonal and near-neighbors
    for k in range(int(min_sep) + 1):
        idx = np.arange(N - k)
        C[idx, idx + k] = False
        C[idx + k, idx] = False
    return C


def f1_from_counts(tp: int, fp: int, fn: int) -> float:
    tp = int(tp)
    fp = int(fp)
    fn = int(fn)
    if tp <= 0:
        return 0.0
    prec = tp / max(1, tp + fp)
    rec = tp / max(1, tp + fn)
    if prec + rec <= 0:
        return 0.0
    return float(2.0 * prec * rec / (prec + rec))


@dataclass
class State:
    # Geometry
    coords: np.ndarray  # (t,3)
    n_path: np.ndarray  # (t,6) int walk
    y_perp: np.ndarray  # (t,3)
    # incremental contact counts (over considered pairs among built residues)
    tp: int
    fp: int
    fn: int
    # cached score components
    auric_penalty: float


def auric_penalty_from_streams(
    y_perp: np.ndarray,
    n_path: np.ndarray,
    *,
    uvec: np.ndarray,
    m: int,
    wA: float,
    wB: float,
) -> float:
    """
    Smaller is better. Return weighted penalty.
    """
    bits_A = rho_A_from_yperp(y_perp, u=tuple(uvec.tolist()), u_mode="fixed", threshold="median")
    bits_B = rho_B_from_npath(n_path)

    folded_A = fold_sliding_windows_bits(bits_A, m=m)
    met_A = metrics_for_stream(bits_A, folded_A)
    folded_B = fold_sliding_windows_bits(bits_B, m=m)
    met_B = metrics_for_stream(bits_B, folded_B)

    # Both metrics are interpreted as "order" certificates: lower tends to be more ordered.
    # Penalize with type entropy for A and smb proxy for B (hybrid intuition).
    pen = 0.0
    pen += float(wA) * float(met_A["type_entropy"])
    pen += float(wB) * float(met_B["smb_rate_hat"])
    return float(pen)


def shared_pca_uvec_from_native(coords_native: np.ndarray, *, bond_len: float, alphabet: str) -> np.ndarray:
    """
    Compute u from native y_perp using oracle direction quantization only for u-definition.
    This uses native geometry (acceptable since target contact map is also native-derived).
    """
    # Build y_perp from oracle direction quantization (same helper used in benches)
    # We replicate the minimal transform here by calling existing oracle routine indirectly:
    # use alphabet steps + nearest direction on Bpar and accumulate n_path.
    steps = alphabet_steps(alphabet)  # (K,6)
    Bpar = icosa_B(PHI) * float(bond_len)
    Bperp = icosa_B(-1.0 / PHI) * float(bond_len)

    d = coords_native[1:] - coords_native[:-1]
    n = np.zeros(6, dtype=np.int32)
    n_path = [n.copy()]
    for t in range(d.shape[0]):
        dt = d[t]
        dt_len = float(np.linalg.norm(dt))
        if dt_len < 1e-8:
            k = 0
        else:
            dt_dir = dt / dt_len
            proj = (Bpar @ steps.T).T  # (K,3)
            proj_norm = np.linalg.norm(proj, axis=1)
            proj_dir = proj / np.clip(proj_norm[:, None], 1e-12, None)
            cos = proj_dir @ dt_dir
            k = int(np.argmax(cos))
        n = n + steps[k]
        n_path.append(n.copy())

    n_path = np.stack(n_path, axis=0)
    y_perp = (Bperp @ n_path.T).T
    Yc = y_perp - y_perp.mean(axis=0, keepdims=True)
    if np.allclose(Yc, 0.0):
        return np.array([1.0, 0.0, 0.0], dtype=np.float64)
    _, _, Vt = np.linalg.svd(Yc, full_matrices=False)
    u = np.asarray(Vt[0], dtype=np.float64)
    nu = float(np.linalg.norm(u))
    u = (u / nu) if nu != 0.0 else np.array([1.0, 0.0, 0.0], dtype=np.float64)
    for k in range(3):
        if abs(u[k]) > 1e-12:
            if u[k] < 0:
                u = -u
            break
    return u


def run_one(
    coords_native: np.ndarray,
    *,
    seed: int,
    codec: str,
    beam: int,
    K: int,
    cutoff: float,
    min_sep: int,
    auric_m: int,
    auric_every: int,
    w_contact: float,
    wA: float,
    wB: float,
) -> Tuple[State, dict]:
    rng = np.random.default_rng(int(seed))
    N = int(coords_native.shape[0])
    if N < 2:
        raise ValueError("Need N>=2")

    # Use native median bond length to set constant step size
    L = np.linalg.norm(coords_native[1:] - coords_native[:-1], axis=1)
    bond_len = float(np.median(L))
    Bpar = icosa_B(PHI) * bond_len
    Bperp = icosa_B(-1.0 / PHI) * bond_len

    steps = alphabet_steps(codec).astype(np.int32)  # (S,6)
    proj = (Bpar @ steps.T).T  # (S,3)
    proj_norm = np.linalg.norm(proj, axis=1)
    proj_dir = proj / np.clip(proj_norm[:, None], 1e-12, None)

    # Target contact matrix from native coords
    C_true = contact_matrix(coords_native, cutoff=cutoff, min_sep=min_sep)

    # Shared protocol axis uvec (derived once from native y_perp)
    uvec = shared_pca_uvec_from_native(coords_native, bond_len=bond_len, alphabet=codec)

    # init state at origin
    coords0 = np.zeros((1, 3), dtype=np.float64)
    n0 = np.zeros((1, 6), dtype=np.int32)
    y0 = (Bperp @ n0.T).T
    init = State(coords=coords0, n_path=n0, y_perp=y0, tp=0, fp=0, fn=0, auric_penalty=0.0)
    beam_states: List[State] = [init]

    # Precompute true-contact counts per column for fast fn updates? We'll do per-step scan.
    for t in range(1, N):
        cand_states: List[State] = []
        for st in beam_states:
            # propose K random steps
            idxs = rng.integers(0, steps.shape[0], size=int(K))
            for k in idxs:
                step6 = steps[int(k)]
                step3 = proj_dir[int(k)] * bond_len

                coords_new = np.vstack([st.coords, st.coords[-1] + step3])
                n_new_last = st.n_path[-1] + step6
                n_new = np.vstack([st.n_path, n_new_last[None, :]])
                y_new_last = (Bperp @ n_new_last.astype(np.float64)).astype(np.float64)
                y_new = np.vstack([st.y_perp, y_new_last[None, :]])

                tp, fp, fn = st.tp, st.fp, st.fn
                # update contact counts for pairs (i,t)
                # only consider i <= t-min_sep-1
                i_max = t - int(min_sep)
                if i_max > 0:
                    d = coords_new[:i_max] - coords_new[t]
                    dist = np.sqrt(np.sum(d * d, axis=1))
                    pred = dist <= float(cutoff)
                    true = C_true[:i_max, t]
                    tp += int(np.sum(pred & true))
                    fp += int(np.sum(pred & (~true)))
                    fn += int(np.sum((~pred) & true))

                aur_pen = st.auric_penalty
                if int(auric_every) > 0 and t >= int(auric_m) and (t % int(auric_every) == 0):
                    aur_pen = auric_penalty_from_streams(y_new, n_new, uvec=uvec, m=int(auric_m), wA=float(wA), wB=float(wB))

                cand_states.append(State(coords=coords_new, n_path=n_new, y_perp=y_new, tp=tp, fp=fp, fn=fn, auric_penalty=aur_pen))

        # select top beam by score
        scored = []
        for st in cand_states:
            f1 = f1_from_counts(st.tp, st.fp, st.fn)
            score = float(w_contact) * f1 - float(st.auric_penalty)
            scored.append((score, f1, st))
        scored.sort(key=lambda x: x[0], reverse=True)
        beam_states = [x[2] for x in scored[: int(beam)]]

    best = beam_states[0]
    tm = tm_score(best.coords, coords_native)
    f1_final = f1_from_counts(best.tp, best.fp, best.fn)
    meta = {
        "seed": int(seed),
        "codec": str(codec),
        "beam": int(beam),
        "K": int(K),
        "cutoff": float(cutoff),
        "min_sep": int(min_sep),
        "auric_m": int(auric_m),
        "auric_every": int(auric_every),
        "w_contact": float(w_contact),
        "wA": float(wA),
        "wB": float(wB),
        "tm": float(tm),
        "contact_f1": float(f1_final),
        "auric_penalty": float(best.auric_penalty),
        "bond_len": float(bond_len),
        "uvec_x": float(uvec[0]),
        "uvec_y": float(uvec[1]),
        "uvec_z": float(uvec[2]),
    }
    return best, meta


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--pdb-id", required=True, help="PDB ID (mmCIF must exist in data/raw/pdb_cache/).")
    ap.add_argument("--tag", default="blind_sprint_auric_small")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--seeds", default="", help="Optional comma list of seeds (overrides --seed).")
    ap.add_argument("--codec", default="pair72", choices=["pair72", "triple232"])
    ap.add_argument("--beam", type=int, default=64)
    ap.add_argument("--K", type=int, default=32)
    ap.add_argument("--cutoff", type=float, default=8.0)
    ap.add_argument("--min-sep", type=int, default=3)
    ap.add_argument("--auric-m", type=int, default=8)
    ap.add_argument("--auric-every", type=int, default=10)
    ap.add_argument("--w-contact", type=float, default=1.0)
    ap.add_argument("--wA", type=float, default=0.2, help="Auric A weight (type_entropy penalty).")
    ap.add_argument("--wB", type=float, default=0.2, help="Auric B weight (smb_rate_hat penalty).")
    args = ap.parse_args()

    root = Path(__file__).resolve().parents[2]
    cache = root / "data" / "raw" / "pdb_cache"
    cif = cache / f"{str(args.pdb_id).upper()}.cif"
    if not cif.exists():
        raise FileNotFoundError(f"Missing mmCIF: {cif}")

    chain_id, coords_native = load_ca_coords_longest_chain(cif)
    coords_native = np.asarray(coords_native, dtype=np.float64)

    seeds = [int(args.seed)]
    if str(args.seeds).strip():
        seeds = [int(x.strip()) for x in str(args.seeds).split(",") if x.strip()]

    out_dir = root / "data" / "processed" / "blind_sprint"
    ensure_dir(out_dir)
    run_dir = root / "docs" / "runs" / "blind_sprint_auric_small"
    ensure_dir(run_dir)

    rows = []
    t0 = time.perf_counter()
    for s in seeds:
        _, meta = run_one(
            coords_native,
            seed=s,
            codec=str(args.codec),
            beam=int(args.beam),
            K=int(args.K),
            cutoff=float(args.cutoff),
            min_sep=int(args.min_sep),
            auric_m=int(args.auric_m),
            auric_every=int(args.auric_every),
            w_contact=float(args.w_contact),
            wA=float(args.wA),
            wB=float(args.wB),
        )
        meta.update({"pdb_id": str(args.pdb_id).upper(), "chain": str(chain_id), "N": int(coords_native.shape[0])})
        rows.append(meta)
        print(f"{args.pdb_id}:{chain_id} seed={s} TM={meta['tm']:.3f} F1={meta['contact_f1']:.3f} auric_pen={meta['auric_penalty']:.3f}", flush=True)

    df = pd.DataFrame(rows)
    out_csv = out_dir / f"blind_sprint_auric_runs_{args.tag}_{str(args.pdb_id).upper()}.csv"
    df.to_csv(out_csv, index=False)

    dt = time.perf_counter() - t0
    report = run_dir / f"{args.tag}_{str(args.pdb_id).upper()}.md"
    csv_rel = out_csv.relative_to(root).as_posix()
    report.write_text(
        "\n".join(
            [
                f"# Blind Sprint + Auric (target={str(args.pdb_id).upper()}:{chain_id})",
                "",
                f"- N: {int(coords_native.shape[0])}",
                f"- codec: `{args.codec}`",
                f"- beam: {int(args.beam)}",
                f"- K per state: {int(args.K)}",
                f"- contact cutoff: {float(args.cutoff)} Å (min_sep={int(args.min_sep)})",
                f"- auric_m: {int(args.auric_m)} (every {int(args.auric_every)} steps)",
                f"- weights: w_contact={float(args.w_contact)}, wA={float(args.wA)}, wB={float(args.wB)}",
                "",
                f"- CSV (gitignored): `{csv_rel}`",
                f"- Runtime: {dt:.1f}s",
                "",
                "## Results",
                "",
                df.sort_values(\"tm\", ascending=False).to_markdown(index=False),
                "",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    print(f"Wrote: {out_csv}")
    print(f"Wrote: {report}")


if __name__ == "__main__":
    main()

