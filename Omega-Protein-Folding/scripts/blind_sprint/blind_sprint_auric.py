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
    phason_stats_from_yperp,
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


def df_to_markdown_table(df: pd.DataFrame) -> str:
    """
    Render a small dataframe as a markdown table without requiring `tabulate`.
    """
    if df is None or len(df) == 0:
        return "_(empty)_"
    cols = list(df.columns)
    lines = []
    lines.append("| " + " | ".join(str(c) for c in cols) + " |")
    lines.append("| " + " | ".join(["---"] * len(cols)) + " |")
    for _, row in df.iterrows():
        cells = []
        for c in cols:
            v = row[c]
            if isinstance(v, float):
                if math.isfinite(v):
                    cells.append(f"{v:.4f}")
                else:
                    cells.append("nan")
            else:
                cells.append(str(v))
        lines.append("| " + " | ".join(cells) + " |")
    return "\n".join(lines)


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
    # incremental distogram RMSE accumulator over considered pairs
    sum_sq: float
    cnt: int
    # cached score components
    auric_penalty: float


def auric_penalty_from_streams(
    y_perp: np.ndarray,
    n_path: np.ndarray,
    *,
    Bperp: np.ndarray,
    uvec: np.ndarray,
    m: int,
    wA: float,
    wB: float,
    rhoB_mode: str = "parity",
    rhoB_threshold: str = "median",
) -> float:
    """
    Smaller is better. Return weighted penalty.
    """
    bits_A = rho_A_from_yperp(y_perp, u=tuple(uvec.tolist()), u_mode="fixed", threshold="median")
    bits_B = rho_B_from_npath(
        n_path,
        mode=str(rhoB_mode),
        Bperp=Bperp,
        u=tuple(uvec.tolist()),
        u_mode="fixed",
        threshold=str(rhoB_threshold),
    )

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
    w_rmse: float,
    wA: float,
    wB: float,
    auric_rhoB_mode: str = "parity",
    auric_rhoB_threshold: str = "median",
    audit_out: Path | None = None,
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
    # Target distogram (rotation-invariant): pairwise distances
    d_true = np.linalg.norm(coords_native[:, None, :] - coords_native[None, :, :], axis=2)

    # Shared protocol axis uvec (derived once from native y_perp)
    uvec = shared_pca_uvec_from_native(coords_native, bond_len=bond_len, alphabet=codec)

    # init state at origin
    coords0 = np.zeros((1, 3), dtype=np.float64)
    n0 = np.zeros((1, 6), dtype=np.int32)
    y0 = (Bperp @ n0.T).T
    init = State(coords=coords0, n_path=n0, y_perp=y0, tp=0, fp=0, fn=0, sum_sq=0.0, cnt=0, auric_penalty=0.0)
    beam_states: List[State] = [init]

    # Precompute true-contact counts per column for fast fn updates? We'll do per-step scan.
    audit_rows: List[dict] = []
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
                sum_sq, cnt = float(st.sum_sq), int(st.cnt)
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

                    # distogram RMSE increment (on same considered pairs)
                    dt = d_true[:i_max, t]
                    err = dist - dt
                    sum_sq += float(np.sum(err * err))
                    cnt += int(err.shape[0])

                aur_pen = st.auric_penalty
                if int(auric_every) > 0 and t >= int(auric_m) and (t % int(auric_every) == 0):
                    aur_pen = auric_penalty_from_streams(
                        y_new,
                        n_new,
                        Bperp=Bperp,
                        uvec=uvec,
                        m=int(auric_m),
                        wA=float(wA),
                        wB=float(wB),
                        rhoB_mode=str(auric_rhoB_mode),
                        rhoB_threshold=str(auric_rhoB_threshold),
                    )

                cand_states.append(
                    State(
                        coords=coords_new,
                        n_path=n_new,
                        y_perp=y_new,
                        tp=tp,
                        fp=fp,
                        fn=fn,
                        sum_sq=sum_sq,
                        cnt=cnt,
                        auric_penalty=aur_pen,
                    )
                )

        # select top beam by score
        scored = []
        for st in cand_states:
            f1 = f1_from_counts(st.tp, st.fp, st.fn)
            rmse = math.sqrt(float(st.sum_sq) / max(1, int(st.cnt)))
            score = float(w_contact) * f1 - float(w_rmse) * rmse - float(st.auric_penalty)
            scored.append((score, f1, st))
        scored.sort(key=lambda x: x[0], reverse=True)
        beam_states = [x[2] for x in scored[: int(beam)]]

        if audit_out is not None and scored:
            best_score, best_f1, best_state = scored[0]
            best_rmse = math.sqrt(float(best_state.sum_sq) / max(1, int(best_state.cnt)))
            ph_rms, ph_max = phason_stats_from_yperp(best_state.y_perp)
            audit_rows.append(
                {
                    "t": int(t),
                    "score": float(best_score),
                    "contact_f1": float(best_f1),
                    "dist_rmse": float(best_rmse),
                    "auric_penalty": float(best_state.auric_penalty),
                    "ph_rms": float(ph_rms),
                    "ph_max": float(ph_max),
                    "tp": int(best_state.tp),
                    "fp": int(best_state.fp),
                    "fn": int(best_state.fn),
                }
            )

    best = beam_states[0]
    tm = tm_score(best.coords, coords_native)
    f1_final = f1_from_counts(best.tp, best.fp, best.fn)
    rmse_final = math.sqrt(float(best.sum_sq) / max(1, int(best.cnt)))
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
        "w_rmse": float(w_rmse),
        "wA": float(wA),
        "wB": float(wB),
        "tm": float(tm),
        "contact_f1": float(f1_final),
        "dist_rmse": float(rmse_final),
        "auric_penalty": float(best.auric_penalty),
        "bond_len": float(bond_len),
        "uvec_x": float(uvec[0]),
        "uvec_y": float(uvec[1]),
        "uvec_z": float(uvec[2]),
    }

    if audit_out is not None and audit_rows:
        ensure_dir(audit_out.parent)
        df_a = pd.DataFrame(audit_rows)
        # attach run-level context for convenience (repeated columns)
        df_a["seed"] = int(seed)
        df_a["codec"] = str(codec)
        df_a["beam"] = int(beam)
        df_a["K"] = int(K)
        df_a["cutoff"] = float(cutoff)
        df_a["min_sep"] = int(min_sep)
        df_a["auric_m"] = int(auric_m)
        df_a["auric_every"] = int(auric_every)
        df_a["w_contact"] = float(w_contact)
        df_a["w_rmse"] = float(w_rmse)
        df_a["wA"] = float(wA)
        df_a["wB"] = float(wB)
        df_a.to_csv(audit_out, index=False)

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
    ap.add_argument("--auric-rhoB-mode", default="parity", choices=["parity", "vel"], help="rhoB readout mode.")
    ap.add_argument(
        "--auric-rhoB-threshold",
        default="median",
        choices=["median", "zero"],
        help="Threshold for rhoB binarization (used for vel mode).",
    )
    ap.add_argument("--w-contact", type=float, default=1.0)
    ap.add_argument("--w-rmse", type=float, default=1.0, help="Distogram RMSE weight (penalty).")
    ap.add_argument("--wA", type=float, default=0.2, help="Auric A weight (type_entropy penalty).")
    ap.add_argument("--wB", type=float, default=0.2, help="Auric B weight (smb_rate_hat penalty).")
    ap.add_argument(
        "--run-dir",
        default="",
        help="Optional directory for markdown reports (absolute or relative to repo root).",
    )
    ap.add_argument(
        "--out-csv",
        default="",
        help="Optional explicit output CSV path (absolute or relative to repo root).",
    )
    ap.add_argument("--no-report", action="store_true", help="Do not write the markdown report.")
    ap.add_argument("--quiet", action="store_true", help="Suppress console prints.")
    ap.add_argument(
        "--audit-out",
        default="",
        help="Optional CSV path to write a per-step audit log (best beam state per t).",
    )
    ap.add_argument(
        "--pred-pdb-out",
        default="",
        help="Optional PDB path to write the final predicted CA trace (single-seed only).",
    )
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
    if str(args.run_dir).strip():
        run_dir = Path(str(args.run_dir))
        if not run_dir.is_absolute():
            run_dir = (root / run_dir).resolve()
    ensure_dir(run_dir)

    rows = []
    t0 = time.perf_counter()
    for s in seeds:
        audit_out: Path | None = None
        if str(args.audit_out).strip() and len(seeds) == 1:
            audit_out = Path(str(args.audit_out))
            if not audit_out.is_absolute():
                audit_out = (root / audit_out).resolve()
        best_state, meta = run_one(
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
            w_rmse=float(args.w_rmse),
            wA=float(args.wA),
            wB=float(args.wB),
            auric_rhoB_mode=str(args.auric_rhoB_mode),
            auric_rhoB_threshold=str(args.auric_rhoB_threshold),
            audit_out=audit_out,
        )
        if str(args.pred_pdb_out).strip() and len(seeds) == 1:
            out_pdb = Path(str(args.pred_pdb_out))
            if not out_pdb.is_absolute():
                out_pdb = (root / out_pdb).resolve()
            ensure_dir(out_pdb.parent)
            # Write a minimal CA-only PDB for visualization/debugging.
            lines = []
            for i, (x, y, z) in enumerate(np.asarray(best_state.coords, dtype=np.float64), start=1):
                # PDB fixed columns; keep it simple.
                lines.append(
                    f"ATOM  {i:5d}  CA  ALA A{i:4d}    "
                    f"{x:8.3f}{y:8.3f}{z:8.3f}"
                    f"{1.00:6.2f}{0.00:6.2f}           C"
                )
            lines.append("END")
            out_pdb.write_text("\n".join(lines) + "\n", encoding="utf-8")
        meta.update({"pdb_id": str(args.pdb_id).upper(), "chain": str(chain_id), "N": int(coords_native.shape[0])})
        rows.append(meta)
        if not bool(args.quiet):
            print(
                f"{args.pdb_id}:{chain_id} seed={s} TM={meta['tm']:.3f} F1={meta['contact_f1']:.3f} "
                f"RMSE={meta['dist_rmse']:.3f} auric_pen={meta['auric_penalty']:.3f}",
                flush=True,
            )

    df = pd.DataFrame(rows)
    out_csv = out_dir / f"blind_sprint_auric_runs_{args.tag}_{str(args.pdb_id).upper()}.csv"
    if str(args.out_csv).strip():
        out_csv = Path(str(args.out_csv))
        if not out_csv.is_absolute():
            out_csv = (root / out_csv).resolve()
        ensure_dir(out_csv.parent)
    df.to_csv(out_csv, index=False)

    dt = time.perf_counter() - t0
    if not bool(args.no_report):
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
                    f"- auric_rhoB_mode: `{args.auric_rhoB_mode}` (threshold={args.auric_rhoB_threshold})",
                    f"- weights: w_contact={float(args.w_contact)}, wA={float(args.wA)}, wB={float(args.wB)}",
                    "",
                    f"- CSV (gitignored): `{csv_rel}`",
                    f"- Runtime: {dt:.1f}s",
                    "",
                    "## Results",
                    "",
                    df_to_markdown_table(df.sort_values("tm", ascending=False)),
                    "",
                ]
            )
            + "\n",
            encoding="utf-8",
        )
        if not bool(args.quiet):
            print(f"Wrote: {report}")
    if not bool(args.quiet):
        print(f"Wrote: {out_csv}")


if __name__ == "__main__":
    main()

