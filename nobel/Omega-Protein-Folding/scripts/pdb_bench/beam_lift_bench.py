#!/usr/bin/env python3
"""
Chain-consistent 6D beam-lift benchmark on real PDB backbones vs random controls.

This is a scalable version of the idea described in:
  artifacts/reports/chain_lift_extalpha_phason_random_controls.md

Pipeline:
  - load longest-chain CA coords from mmCIF
  - represent native backbone via bond directions d_t
  - run beam search over a 6D integer walk n_t using an alphabet (axis12/pair72)
  - objective: angular mismatch + optional perp-space regularizer
  - compute ph_rms/ph_max from y_perp=B_perp n_t (w0=mean)
  - compare ph_rms against random-chain controls with same bond-length sequence

Outputs:
  - data/processed/pdb_bench/beam_lift_bench_summary_{tag}.csv
  - artifacts/reports/pdb_beam_lift_bench_{tag}.md
  - (optional) long-form samples CSV
"""

from __future__ import annotations

import argparse
import time
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd

import sys

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str((Path(__file__).resolve().parents[1] / "auric").resolve()))

from bench_utils import (  # noqa: E402
    PHI,
    alphabet_steps,
    cliffs_delta_one_vs_many,
    contact_density,
    ensure_dir,
    icosa_B,
    load_ca_coords_longest_chain,
    phason_stats_from_yperp,
    phason_stats_piecewise,
    random_chain_from_lengths,
    radius_of_gyration,
    tm_score,
)

try:
    from fold_m import fold_sliding_windows_bits  # type: ignore
    from metrics import metrics_for_stream  # type: ignore
    from readout import rho_A_from_yperp, rho_B_from_npath  # type: ignore
except Exception:  # pragma: no cover
    fold_sliding_windows_bits = None  # type: ignore
    metrics_for_stream = None  # type: ignore
    rho_A_from_yperp = None  # type: ignore
    rho_B_from_npath = None  # type: ignore


def beam_lift(
    coords_native: np.ndarray,
    *,
    alphabet: str,
    beam_width: int = 16,
    topk: int = 6,
    w_perp: float = 0.01,
    piecewise_seg_len: int = 0,
    phi: float = PHI,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, float, float, float]:
    """
    Return (coords_pred, n_path, y_perp, tm, ph_rms, ph_max).
    """
    N = coords_native.shape[0]
    d = coords_native[1:] - coords_native[:-1]  # (N-1,3)
    L = np.linalg.norm(d, axis=1)
    median_len = float(np.median(L))

    Bpar = icosa_B(phi) * median_len
    Bperp = icosa_B(-1.0 / phi) * median_len

    steps = alphabet_steps(alphabet)  # (K,6)
    proj = (Bpar @ steps.T).T  # (K,3)
    proj_dir = proj / np.clip(np.linalg.norm(proj, axis=1)[:, None], 1e-12, None)

    # Precompute candidate step indices per bond direction
    dt_dir = d / np.clip(L[:, None], 1e-12, None)
    cand = []
    for t in range(N - 1):
        cos = proj_dir @ dt_dir[t]
        idx = np.argsort(-cos)[:topk]
        cand.append(idx.astype(np.int32))

    B = int(max(1, beam_width))
    K = int(max(1, topk))

    # Beam arrays
    n_beam = np.zeros((1, 6), dtype=np.int32)
    coord_beam = np.zeros((1, 3), dtype=np.float64)
    cost_beam = np.zeros((1,), dtype=np.float64)
    # For piecewise phason regularization: store y_perp at current segment start per beam.
    yseg0_beam = np.zeros((1, 3), dtype=np.float64)

    # Backpointers per step
    parents: List[np.ndarray] = []
    step_ids: List[np.ndarray] = []

    for t in range(N - 1):
        idx = cand[t]  # (K,)
        # Expand all beams
        nb = n_beam.shape[0]
        # broadcast (nb,K,6)
        n_new = n_beam[:, None, :] + steps[idx][None, :, :]
        # y_perp for new nodes: (nb,K,3)
        y_new = (Bperp @ n_new.reshape(-1, 6).T).T.reshape(nb, K, 3)
        # cost: angular mismatch + perp regularizer
        cos = (proj_dir[idx] @ dt_dir[t]).astype(np.float64)  # (K,)
        ang_cost = (1.0 - cos)[None, :]  # (1,K)
        if piecewise_seg_len and piecewise_seg_len > 0:
            # reset segment anchor at boundaries (approximate piecewise slice drift)
            if (t % int(piecewise_seg_len)) == 0:
                yseg0_beam = coord_beam * 0.0 + (Bperp @ n_beam.T).T  # current y_perp
            perp_cost = np.sum((y_new - yseg0_beam[:, None, :]) ** 2, axis=2)
        else:
            perp_cost = np.sum(y_new**2, axis=2)  # (nb,K)
        cost_new = cost_beam[:, None] + ang_cost + w_perp * perp_cost

        # coords update: step_dir * observed bond length
        step_vec = proj_dir[idx] * float(L[t])  # (K,3)
        coord_new = coord_beam[:, None, :] + step_vec[None, :, :]

        # Select best B across all nb*K
        flat_cost = cost_new.reshape(-1)
        take = np.argsort(flat_cost)[:B]
        parent = (take // K).astype(np.int32)
        choice = (take % K).astype(np.int32)

        n_beam = n_new.reshape(-1, 6)[take]
        coord_beam = coord_new.reshape(-1, 3)[take]
        cost_beam = flat_cost[take]
        if piecewise_seg_len and piecewise_seg_len > 0:
            yseg0_beam = yseg0_beam[parent]

        parents.append(parent)
        step_ids.append(idx[choice])

    # Pick best final
    best = int(np.argmin(cost_beam))
    # Backtrack step indices to recover n_path
    n_path = np.zeros((N, 6), dtype=np.int32)
    coords_pred = np.zeros((N, 3), dtype=np.float64)
    cur = best
    chosen_steps = []
    for t in range(N - 2, -1, -1):
        chosen_steps.append(int(step_ids[t][cur]))
        cur = int(parents[t][cur])
    chosen_steps = chosen_steps[::-1]

    # Rebuild path deterministically
    for t, sid in enumerate(chosen_steps):
        dt_len = float(L[t])
        n_path[t + 1] = n_path[t] + steps[sid]
        coords_pred[t + 1] = coords_pred[t] + proj_dir[sid] * dt_len

    y_perp = (Bperp @ n_path.T).T
    ph_rms, ph_max = phason_stats_from_yperp(y_perp)
    if piecewise_seg_len and piecewise_seg_len > 0:
        ph_rms_pw, ph_max_pw = phason_stats_piecewise(y_perp, segment_len=int(piecewise_seg_len))
        # Report piecewise metrics (override) to match "piecewise slice" intent.
        ph_rms, ph_max = ph_rms_pw, ph_max_pw
    tm = tm_score(coords_pred, coords_native)
    return coords_pred, n_path, y_perp, tm, ph_rms, ph_max


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--ids", required=True)
    ap.add_argument("--tag", default="n1000_seed0_pair72")
    ap.add_argument("--alphabet", default="pair72", choices=["axis12", "pair72", "triple232"])
    ap.add_argument("--beam-width", type=int, default=12)
    ap.add_argument("--topk", type=int, default=6)
    ap.add_argument("--w-perp", type=float, default=0.01)
    ap.add_argument("--piecewise-seg-len", type=int, default=0, help="If >0, use piecewise phason regularization/metric with this segment length.")
    ap.add_argument("--random-reps", type=int, default=2)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--min-len", type=int, default=60)
    ap.add_argument("--max-len", type=int, default=350)
    ap.add_argument("--max-proteins", type=int, default=0, help="0 means no limit.")
    ap.add_argument("--auric", action="store_true", help="If set, compute auric Fold_m certificate metrics.")
    ap.add_argument("--auric-m", default="6,8,10,12", help="Comma list of m values for Fold_m (e.g. 6,8,10,12).")
    ap.add_argument("--auric-readouts", default="all", choices=["A", "B", "all"], help="Which ρ readouts to use.")
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
    auric_rows: List[Dict[str, object]] = []
    rng_global = np.random.default_rng(args.seed)
    auric_ms: List[int] = []
    if args.auric:
        if fold_sliding_windows_bits is None or metrics_for_stream is None or rho_A_from_yperp is None or rho_B_from_npath is None:
            raise RuntimeError("Auric is enabled but auric modules failed to import (scripts/auric/*).")
        auric_ms = [int(x.strip()) for x in str(args.auric_m).split(",") if x.strip()]
        if not auric_ms:
            raise ValueError("--auric-m must contain at least one integer m value")

    t0 = time.perf_counter()
    used = 0

    for pdb_id in ids:
        if args.max_proteins and used >= args.max_proteins:
            break
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

        d = coords[1:] - coords[:-1]
        bond_lengths = np.linalg.norm(d, axis=1)

        seed_i = (abs(hash((args.seed, pdb_id, "beam"))) % (2**32 - 1)) + 1
        rng = np.random.default_rng(seed_i)

        t1 = time.perf_counter()
        _, n_path, y_perp, tm, ph_rms, ph_max = beam_lift(
            coords,
            alphabet=args.alphabet,
            beam_width=args.beam_width,
            topk=args.topk,
            w_perp=args.w_perp,
            piecewise_seg_len=args.piecewise_seg_len,
        )
        dt_real = time.perf_counter() - t1

        # Random controls (mean over reps)
        ph_rms_rand = []
        auric_ent_random: Dict[Tuple[str, int], List[float]] = {}
        for r in range(args.random_reps):
            coords_r = random_chain_from_lengths(rng, bond_lengths)
            _, npr, ypr, _, pr, _ = beam_lift(
                coords_r,
                alphabet=args.alphabet,
                beam_width=args.beam_width,
                topk=args.topk,
                w_perp=args.w_perp,
                piecewise_seg_len=args.piecewise_seg_len,
            )
            ph_rms_rand.append(pr)
            if args.auric:
                bits_A = rho_A_from_yperp(ypr)
                bits_B = rho_B_from_npath(npr)
                for m in auric_ms:
                    if args.auric_readouts in {"A", "all"}:
                        folded = fold_sliding_windows_bits(bits_A, m=m)
                        met = metrics_for_stream(bits_A, folded)
                        auric_ent_random.setdefault(("A", m), []).append(float(met["type_entropy"]))
                    if args.auric_readouts in {"B", "all"}:
                        folded = fold_sliding_windows_bits(bits_B, m=m)
                        met = metrics_for_stream(bits_B, folded)
                        auric_ent_random.setdefault(("B", m), []).append(float(met["type_entropy"]))

        rg = radius_of_gyration(coords)
        cd = contact_density(coords, cutoff=8.0, min_sep=3)
        rows.append(
            {
                "pdb_id": pdb_id,
                "chain": chain_id,
                "N": N,
                "alphabet": args.alphabet,
                "beam_width": args.beam_width,
                "topk": args.topk,
                "w_perp": args.w_perp,
                "piecewise_seg_len": args.piecewise_seg_len,
                "tm_beam": tm,
                "ph_rms_real": ph_rms,
                "ph_max_real": ph_max,
                "ph_rms_random_mean": float(np.mean(ph_rms_rand)) if len(ph_rms_rand) else float("nan"),
                "rg_real": rg,
                "contact_density_real": cd,
                "runtime_s_real": dt_real,
            }
        )

        if args.auric:
            bits_A_real = rho_A_from_yperp(y_perp)
            bits_B_real = rho_B_from_npath(n_path)
            for m in auric_ms:
                if args.auric_readouts in {"A", "all"}:
                    folded = fold_sliding_windows_bits(bits_A_real, m=m)
                    met = metrics_for_stream(bits_A_real, folded)
                    ent_rand = np.asarray(auric_ent_random.get(("A", m), []), dtype=np.float64)
                    auric_rows.append(
                        {
                            "pdb_id": pdb_id,
                            "chain": chain_id,
                            "N": N,
                            "alphabet": args.alphabet,
                            "readout": "A",
                            "m": m,
                            "tm_beam": tm,
                            "ph_rms_real": ph_rms,
                            "ph_max_real": ph_max,
                            "type_entropy_real": float(met["type_entropy"]),
                            "type_entropy_random_mean": float(np.mean(ent_rand)) if len(ent_rand) else float("nan"),
                            "delta_type_entropy_real_vs_random": cliffs_delta_one_vs_many(float(met["type_entropy"]), ent_rand) if len(ent_rand) else float("nan"),
                            "type_support_real": float(met["type_support"]),
                            "run1_mean_real": float(met["run1_mean"]),
                            "run1_max_real": float(met["run1_max"]),
                            "smb_rate_hat_real": float(met["smb_rate_hat"]),
                        }
                    )
                if args.auric_readouts in {"B", "all"}:
                    folded = fold_sliding_windows_bits(bits_B_real, m=m)
                    met = metrics_for_stream(bits_B_real, folded)
                    ent_rand = np.asarray(auric_ent_random.get(("B", m), []), dtype=np.float64)
                    auric_rows.append(
                        {
                            "pdb_id": pdb_id,
                            "chain": chain_id,
                            "N": N,
                            "alphabet": args.alphabet,
                            "readout": "B",
                            "m": m,
                            "tm_beam": tm,
                            "ph_rms_real": ph_rms,
                            "ph_max_real": ph_max,
                            "type_entropy_real": float(met["type_entropy"]),
                            "type_entropy_random_mean": float(np.mean(ent_rand)) if len(ent_rand) else float("nan"),
                            "delta_type_entropy_real_vs_random": cliffs_delta_one_vs_many(float(met["type_entropy"]), ent_rand) if len(ent_rand) else float("nan"),
                            "type_support_real": float(met["type_support"]),
                            "run1_mean_real": float(met["run1_mean"]),
                            "run1_max_real": float(met["run1_max"]),
                            "smb_rate_hat_real": float(met["smb_rate_hat"]),
                        }
                    )
        used += 1

    df = pd.DataFrame(rows)
    out_csv = out_dir / f"beam_lift_bench_summary_{args.tag}.csv"
    df.to_csv(out_csv, index=False)
    auric_out_csv = out_dir / f"beam_lift_auric_summary_{args.tag}.csv"
    if args.auric:
        pd.DataFrame(auric_rows).to_csv(auric_out_csv, index=False)

    ids_rel = str(ids_path.relative_to(root)).replace("\\", "/")
    out_rel = str(out_csv.relative_to(root)).replace("\\", "/")
    report_path = rep_dir / f"pdb_beam_lift_bench_{args.tag}.md"
    lines = [
        "# PDB beam-lift benchmark (large-scale)",
        "",
        "Chain-consistent 6D walk via beam search over the alphabet; reports phason proxy spread and TM-score (after Kabsch).",
        "",
        "## Dataset",
        f"- IDs: `{ids_rel}`",
        f"- Used chains: {len(df)}",
        f"- Length filter: N in [{args.min_len},{args.max_len}]",
        "",
        "## Settings",
        f"- alphabet: {args.alphabet}",
        f"- beam_width: {args.beam_width}",
        f"- topk per step: {args.topk}",
        f"- w_perp: {args.w_perp}",
        f"- random reps: {args.random_reps}",
        "",
        "## Outputs",
        f"- CSV: `{out_rel}`",
        f"- Auric CSV: `{str(auric_out_csv.relative_to(root)).replace('\\\\', '/')}`" if args.auric else "",
        "",
        "## Quick stats (median)",
        f"- tm_beam: {df['tm_beam'].median():.3f}",
        f"- ph_rms_real: {df['ph_rms_real'].median():.3f}",
        f"- ph_rms_random_mean: {df['ph_rms_random_mean'].median():.3f}",
        "",
    ]
    report_path.write_text("\n".join(lines) + "\n", encoding="utf-8")

    dt = time.perf_counter() - t0
    print(f"Wrote: {out_csv}")
    if args.auric:
        print(f"Wrote: {auric_out_csv}")
    print(f"Wrote: {report_path}")
    print(f"Runtime: {dt:.1f}s")


if __name__ == "__main__":
    main()

