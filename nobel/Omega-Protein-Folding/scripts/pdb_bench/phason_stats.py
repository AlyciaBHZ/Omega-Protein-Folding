#!/usr/bin/env python3
"""
Large-scale phason proxy statistics on real PDB chains vs controls.

We use an operational proxy compatible with the evidence pack:
  - recover a 6D integer walk n_i via oracle direction quantization under a given alphabet
  - compute y_i^⊥ = B_perp n_i, set w0 = mean(y^⊥)
  - report ph_rms, ph_max = max ||y^⊥-w0||

Controls:
  - random chains with same per-bond length sequence
  - perturbed-native chains (direction noise, bond lengths preserved)

Outputs:
  - data/processed/pdb_bench/phason_stats_samples_{tag}.csv  (long-form per-replicate)
  - data/processed/pdb_bench/phason_stats_summary_{tag}.csv  (per-protein summary + effect sizes)
  - artifacts/reports/pdb_phason_stats_{tag}.md  (dataset-level summary)
"""

from __future__ import annotations

import argparse
import csv
import time
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd

import sys

# Allow running as a script without installing a package.
sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str((Path(__file__).resolve().parents[1] / "auric").resolve()))

from bench_utils import (  # noqa: E402
    PHI,
    cliffs_delta_one_vs_many,
    contact_density,
    ensure_dir,
    load_ca_coords_longest_chain,
    oracle_direction_reconstruct,
    perturb_chain_directions,
    phason_stats_piecewise,
    phason_stats_from_yperp,
    phason_stats_linear_drift,
    radius_of_gyration,
    random_chain_from_lengths,
    tm_score,
    w0_piecewise_jumps,
)

from controls import PseudoMoltenParams, generate_pseudo_molten_globule  # noqa: E402

try:
    # Optional auric integration (kept lightweight and script-friendly).
    from fold_m import fold_sliding_windows_bits  # type: ignore
    from metrics import metrics_for_stream  # type: ignore
    from readout import rho_A_from_yperp, rho_B_from_npath  # type: ignore
except Exception:  # pragma: no cover
    fold_sliding_windows_bits = None  # type: ignore
    metrics_for_stream = None  # type: ignore
    rho_A_from_yperp = None  # type: ignore
    rho_B_from_npath = None  # type: ignore


def bootstrap_ci_mean(xs: np.ndarray, rng: np.random.Generator, iters: int = 2000, alpha: float = 0.05) -> Tuple[float, float]:
    xs = np.asarray(xs, dtype=np.float64)
    n = len(xs)
    if n == 0:
        return float("nan"), float("nan")
    means = []
    for _ in range(iters):
        samp = xs[rng.integers(0, n, size=n)]
        means.append(float(np.mean(samp)))
    means.sort()
    lo = means[int((alpha / 2) * iters)]
    hi = means[int((1 - alpha / 2) * iters) - 1]
    return lo, hi


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--ids", required=True, help="IDs file (one PDB ID per line).")
    ap.add_argument("--tag", default="n1000_seed0")
    ap.add_argument("--alphabet", default="triple232", choices=["axis12", "pair72", "triple232", "all"])
    ap.add_argument("--random-reps", type=int, default=5)
    ap.add_argument("--perturb-reps", type=int, default=3)
    ap.add_argument("--perturb-noise", type=float, default=0.15, help="Direction noise std (dimensionless).")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--min-len", type=int, default=60)
    ap.add_argument("--max-len", type=int, default=350)
    ap.add_argument("--piecewise-seg-len", type=int, default=0, help="If >0, compute piecewise phason proxy with per-segment w0.")
    ap.add_argument("--linear-drift", action="store_true", help="If set, compute linear-drift-corrected phason proxy.")
    ap.add_argument("--random-model", default="gaussian", choices=["gaussian", "pseudo_molten"], help="Which random control generator to use.")
    ap.add_argument("--pmg-persistence", type=float, default=0.7)
    ap.add_argument("--pmg-min-dist", type=float, default=3.5)
    ap.add_argument("--pmg-conf-k", type=float, default=1.0)
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
    alphabets = ["axis12", "pair72", "triple232"] if args.alphabet == "all" else [args.alphabet]

    rng_global = np.random.default_rng(args.seed)

    samples_rows: List[Dict[str, object]] = []
    summary_rows: List[Dict[str, object]] = []
    auric_samples_rows: List[Dict[str, object]] = []
    auric_summary_rows: List[Dict[str, object]] = []

    t0 = time.perf_counter()
    used = 0

    auric_ms: List[int] = []
    if args.auric:
        if fold_sliding_windows_bits is None or metrics_for_stream is None or rho_A_from_yperp is None or rho_B_from_npath is None:
            raise RuntimeError("Auric is enabled but auric modules failed to import (scripts/auric/*).")
        auric_ms = [int(x.strip()) for x in str(args.auric_m).split(",") if x.strip()]
        if not auric_ms:
            raise ValueError("--auric-m must contain at least one integer m value")
        for m in auric_ms:
            if m <= 0 or m > 24:
                raise ValueError(f"Unsupported m={m}. Suggested range is 1..24 for runtime/memory safety.")

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

        # Deterministic per-protein RNG
        seed_i = (abs(hash((args.seed, pdb_id))) % (2**32 - 1)) + 1
        rng = np.random.default_rng(seed_i)

        for alpha_name in alphabets:
            # Real
            pred, n_path, y_perp, ph_rms_real, ph_max_real = oracle_direction_reconstruct(coords, alphabet=alpha_name, phi=PHI)
            if args.piecewise_seg_len and args.piecewise_seg_len > 0:
                ph_rms_pw, ph_max_pw = phason_stats_piecewise(y_perp, segment_len=int(args.piecewise_seg_len))
            else:
                ph_rms_pw, ph_max_pw = float("nan"), float("nan")
            if args.linear_drift:
                ph_rms_lin, ph_max_lin, drift_speed = phason_stats_linear_drift(y_perp)
            else:
                ph_rms_lin, ph_max_lin, drift_speed = float("nan"), float("nan"), float("nan")
            if args.piecewise_seg_len and args.piecewise_seg_len > 0:
                w0_jump_mean, w0_jump_max = w0_piecewise_jumps(y_perp, segment_len=int(args.piecewise_seg_len))
            else:
                w0_jump_mean, w0_jump_max = float("nan"), float("nan")
            tm = tm_score(pred, coords)
            rg = radius_of_gyration(coords)
            cd = contact_density(coords, cutoff=8.0, min_sep=3)
            samples_rows.append(
                {
                    "pdb_id": pdb_id,
                    "chain": chain_id,
                    "N": N,
                    "alphabet": alpha_name,
                    "group": "real",
                    "rep": 0,
                    "ph_rms": ph_rms_real,
                    "ph_max": ph_max_real,
                    "ph_rms_piecewise": ph_rms_pw,
                    "ph_max_piecewise": ph_max_pw,
                    "ph_rms_lin": ph_rms_lin,
                    "ph_max_lin": ph_max_lin,
                    "drift_speed": drift_speed,
                    "w0_jump_mean": w0_jump_mean,
                    "w0_jump_max": w0_jump_max,
                    "tm_oracle": tm,
                    "rg": rg,
                    "contact_density": cd,
                }
            )

            # Auric (real)
            if args.auric:
                bits_A_real = rho_A_from_yperp(y_perp)
                bits_B_real = rho_B_from_npath(n_path)

            # Random controls
            ph_rms_rand = []
            ph_max_rand = []
            ph_rms_rand_pw = []
            ph_rms_rand_lin = []
            # Auric controls (store entropy streams for effect sizes)
            auric_ent_random: Dict[Tuple[str, int], List[float]] = {}
            auric_ent_pert: Dict[Tuple[str, int], List[float]] = {}

            for r in range(args.random_reps):
                if args.random_model == "pseudo_molten":
                    pmg = PseudoMoltenParams(
                        persistence=float(args.pmg_persistence),
                        min_dist=float(args.pmg_min_dist),
                        conf_k=float(args.pmg_conf_k),
                    )
                    coords_r = generate_pseudo_molten_globule(rng, bond_lengths, target_rg=rg, params=pmg)
                else:
                    coords_r = random_chain_from_lengths(rng, bond_lengths)
                _, npr, ypr, pr, pm = oracle_direction_reconstruct(coords_r, alphabet=alpha_name, phi=PHI)
                if args.piecewise_seg_len and args.piecewise_seg_len > 0:
                    pr_pw, pm_pw = phason_stats_piecewise(ypr, segment_len=int(args.piecewise_seg_len))
                else:
                    pr_pw, pm_pw = float("nan"), float("nan")
                if args.linear_drift:
                    pr_lin, pm_lin, drift_lin = phason_stats_linear_drift(ypr)
                else:
                    pr_lin, pm_lin, drift_lin = float("nan"), float("nan"), float("nan")
                if args.piecewise_seg_len and args.piecewise_seg_len > 0:
                    jump_m, jump_M = w0_piecewise_jumps(ypr, segment_len=int(args.piecewise_seg_len))
                else:
                    jump_m, jump_M = float("nan"), float("nan")
                ph_rms_rand.append(pr)
                ph_max_rand.append(pm)
                ph_rms_rand_pw.append(pr_pw)
                ph_rms_rand_lin.append(pr_lin)
                samples_rows.append(
                    {
                        "pdb_id": pdb_id,
                        "chain": chain_id,
                        "N": N,
                        "alphabet": alpha_name,
                        "group": "random",
                        "rep": r,
                        "ph_rms": pr,
                        "ph_max": pm,
                        "ph_rms_piecewise": pr_pw,
                        "ph_max_piecewise": pm_pw,
                        "ph_rms_lin": pr_lin,
                        "ph_max_lin": pm_lin,
                        "drift_speed": drift_lin,
                        "w0_jump_mean": jump_m,
                        "w0_jump_max": jump_M,
                        "tm_oracle": float("nan"),
                        "rg": radius_of_gyration(coords_r),
                        "contact_density": contact_density(coords_r, cutoff=8.0, min_sep=3),
                    }
                )

                if args.auric:
                    bits_A = rho_A_from_yperp(ypr)
                    bits_B = rho_B_from_npath(npr)
                    for m in auric_ms:
                        if args.auric_readouts in {"A", "all"}:
                            folded = fold_sliding_windows_bits(bits_A, m=m)
                            met = metrics_for_stream(bits_A, folded)
                            auric_samples_rows.append(
                                {
                                    "pdb_id": pdb_id,
                                    "chain": chain_id,
                                    "N": N,
                                    "alphabet": alpha_name,
                                    "group": "random",
                                    "rep": r,
                                    "readout": "A",
                                    "m": m,
                                    "bits_len": int(bits_A.shape[0]),
                                    "types_len": int(folded.shape[0]),
                                    **met,
                                }
                            )
                            auric_ent_random.setdefault(("A", m), []).append(float(met["type_entropy"]))
                        if args.auric_readouts in {"B", "all"}:
                            folded = fold_sliding_windows_bits(bits_B, m=m)
                            met = metrics_for_stream(bits_B, folded)
                            auric_samples_rows.append(
                                {
                                    "pdb_id": pdb_id,
                                    "chain": chain_id,
                                    "N": N,
                                    "alphabet": alpha_name,
                                    "group": "random",
                                    "rep": r,
                                    "readout": "B",
                                    "m": m,
                                    "bits_len": int(bits_B.shape[0]),
                                    "types_len": int(folded.shape[0]),
                                    **met,
                                }
                            )
                            auric_ent_random.setdefault(("B", m), []).append(float(met["type_entropy"]))

            # Perturbed native controls
            ph_rms_pert = []
            ph_max_pert = []
            ph_rms_pert_pw = []
            ph_rms_pert_lin = []
            for r in range(args.perturb_reps):
                coords_p = perturb_chain_directions(rng, coords, noise=args.perturb_noise)
                _, npp, ypp, pr, pm = oracle_direction_reconstruct(coords_p, alphabet=alpha_name, phi=PHI)
                if args.piecewise_seg_len and args.piecewise_seg_len > 0:
                    pr_pw, pm_pw = phason_stats_piecewise(ypp, segment_len=int(args.piecewise_seg_len))
                else:
                    pr_pw, pm_pw = float("nan"), float("nan")
                if args.linear_drift:
                    pr_lin, pm_lin, drift_lin = phason_stats_linear_drift(ypp)
                else:
                    pr_lin, pm_lin, drift_lin = float("nan"), float("nan"), float("nan")
                if args.piecewise_seg_len and args.piecewise_seg_len > 0:
                    jump_m, jump_M = w0_piecewise_jumps(ypp, segment_len=int(args.piecewise_seg_len))
                else:
                    jump_m, jump_M = float("nan"), float("nan")
                ph_rms_pert.append(pr)
                ph_max_pert.append(pm)
                ph_rms_pert_pw.append(pr_pw)
                ph_rms_pert_lin.append(pr_lin)
                samples_rows.append(
                    {
                        "pdb_id": pdb_id,
                        "chain": chain_id,
                        "N": N,
                        "alphabet": alpha_name,
                        "group": "perturbed",
                        "rep": r,
                        "ph_rms": pr,
                        "ph_max": pm,
                        "ph_rms_piecewise": pr_pw,
                        "ph_max_piecewise": pm_pw,
                        "ph_rms_lin": pr_lin,
                        "ph_max_lin": pm_lin,
                        "drift_speed": drift_lin,
                        "w0_jump_mean": jump_m,
                        "w0_jump_max": jump_M,
                        "tm_oracle": float("nan"),
                        "rg": radius_of_gyration(coords_p),
                        "contact_density": contact_density(coords_p, cutoff=8.0, min_sep=3),
                    }
                )

                if args.auric:
                    bits_A = rho_A_from_yperp(ypp)
                    bits_B = rho_B_from_npath(npp)
                    for m in auric_ms:
                        if args.auric_readouts in {"A", "all"}:
                            folded = fold_sliding_windows_bits(bits_A, m=m)
                            met = metrics_for_stream(bits_A, folded)
                            auric_samples_rows.append(
                                {
                                    "pdb_id": pdb_id,
                                    "chain": chain_id,
                                    "N": N,
                                    "alphabet": alpha_name,
                                    "group": "perturbed",
                                    "rep": r,
                                    "readout": "A",
                                    "m": m,
                                    "bits_len": int(bits_A.shape[0]),
                                    "types_len": int(folded.shape[0]),
                                    **met,
                                }
                            )
                            auric_ent_pert.setdefault(("A", m), []).append(float(met["type_entropy"]))
                        if args.auric_readouts in {"B", "all"}:
                            folded = fold_sliding_windows_bits(bits_B, m=m)
                            met = metrics_for_stream(bits_B, folded)
                            auric_samples_rows.append(
                                {
                                    "pdb_id": pdb_id,
                                    "chain": chain_id,
                                    "N": N,
                                    "alphabet": alpha_name,
                                    "group": "perturbed",
                                    "rep": r,
                                    "readout": "B",
                                    "m": m,
                                    "bits_len": int(bits_B.shape[0]),
                                    "types_len": int(folded.shape[0]),
                                    **met,
                                }
                            )
                            auric_ent_pert.setdefault(("B", m), []).append(float(met["type_entropy"]))

            ph_rms_rand = np.asarray(ph_rms_rand, dtype=np.float64)
            ph_rms_pert = np.asarray(ph_rms_pert, dtype=np.float64)
            ph_rms_rand_pw = np.asarray(ph_rms_rand_pw, dtype=np.float64)
            ph_rms_pert_pw = np.asarray(ph_rms_pert_pw, dtype=np.float64)
            ph_rms_rand_lin = np.asarray(ph_rms_rand_lin, dtype=np.float64)
            ph_rms_pert_lin = np.asarray(ph_rms_pert_lin, dtype=np.float64)

            summary_rows.append(
                {
                    "pdb_id": pdb_id,
                    "chain": chain_id,
                    "N": N,
                    "alphabet": alpha_name,
                    "ph_rms_real": ph_rms_real,
                    "ph_max_real": ph_max_real,
                    "ph_rms_piecewise_real": ph_rms_pw,
                    "ph_max_piecewise_real": ph_max_pw,
                    "ph_rms_lin_real": ph_rms_lin,
                    "ph_max_lin_real": ph_max_lin,
                    "drift_speed_real": drift_speed,
                    "w0_jump_mean_real": w0_jump_mean,
                    "w0_jump_max_real": w0_jump_max,
                    "tm_oracle": tm,
                    "rg_real": rg,
                    "contact_density_real": cd,
                    "ph_rms_random_mean": float(np.mean(ph_rms_rand)) if len(ph_rms_rand) else float("nan"),
                    "ph_rms_random_median": float(np.median(ph_rms_rand)) if len(ph_rms_rand) else float("nan"),
                    "ph_rms_random_std": float(np.std(ph_rms_rand)) if len(ph_rms_rand) else float("nan"),
                    "ph_rms_pert_mean": float(np.mean(ph_rms_pert)) if len(ph_rms_pert) else float("nan"),
                    "ph_rms_pert_median": float(np.median(ph_rms_pert)) if len(ph_rms_pert) else float("nan"),
                    "ph_rms_pert_std": float(np.std(ph_rms_pert)) if len(ph_rms_pert) else float("nan"),
                    "cliffs_delta_real_vs_random": cliffs_delta_one_vs_many(ph_rms_real, ph_rms_rand) if len(ph_rms_rand) else float("nan"),
                    "cliffs_delta_real_vs_perturbed": cliffs_delta_one_vs_many(ph_rms_real, ph_rms_pert) if len(ph_rms_pert) else float("nan"),
                    "cliffs_delta_piecewise_vs_random": cliffs_delta_one_vs_many(ph_rms_pw, ph_rms_rand_pw) if len(ph_rms_rand_pw) else float("nan"),
                    "cliffs_delta_piecewise_vs_perturbed": cliffs_delta_one_vs_many(ph_rms_pw, ph_rms_pert_pw) if len(ph_rms_pert_pw) else float("nan"),
                    "cliffs_delta_lin_vs_random": cliffs_delta_one_vs_many(ph_rms_lin, ph_rms_rand_lin) if len(ph_rms_rand_lin) else float("nan"),
                    "cliffs_delta_lin_vs_perturbed": cliffs_delta_one_vs_many(ph_rms_lin, ph_rms_pert_lin) if len(ph_rms_pert_lin) else float("nan"),
                }
            )

            # Auric summary (real vs controls) in long-form rows
            if args.auric:
                for m in auric_ms:
                    if args.auric_readouts in {"A", "all"}:
                        folded = fold_sliding_windows_bits(bits_A_real, m=m)
                        met_real = metrics_for_stream(bits_A_real, folded)
                        ent_rand = np.asarray(auric_ent_random.get(("A", m), []), dtype=np.float64)
                        ent_pert = np.asarray(auric_ent_pert.get(("A", m), []), dtype=np.float64)
                        auric_summary_rows.append(
                            {
                                "pdb_id": pdb_id,
                                "chain": chain_id,
                                "N": N,
                                "alphabet": alpha_name,
                                "readout": "A",
                                "m": m,
                                "type_entropy_real": float(met_real["type_entropy"]),
                                "type_entropy_random_mean": float(np.mean(ent_rand)) if len(ent_rand) else float("nan"),
                                "type_entropy_pert_mean": float(np.mean(ent_pert)) if len(ent_pert) else float("nan"),
                                "delta_type_entropy_real_vs_random": cliffs_delta_one_vs_many(float(met_real["type_entropy"]), ent_rand) if len(ent_rand) else float("nan"),
                                "delta_type_entropy_real_vs_perturbed": cliffs_delta_one_vs_many(float(met_real["type_entropy"]), ent_pert) if len(ent_pert) else float("nan"),
                                "type_support_real": float(met_real["type_support"]),
                                "run1_mean_real": float(met_real["run1_mean"]),
                                "run1_max_real": float(met_real["run1_max"]),
                                "smb_rate_hat_real": float(met_real["smb_rate_hat"]),
                            }
                        )
                        auric_samples_rows.append(
                            {
                                "pdb_id": pdb_id,
                                "chain": chain_id,
                                "N": N,
                                "alphabet": alpha_name,
                                "group": "real",
                                "rep": 0,
                                "readout": "A",
                                "m": m,
                                "bits_len": int(bits_A_real.shape[0]),
                                "types_len": int(folded.shape[0]),
                                **met_real,
                            }
                        )
                    if args.auric_readouts in {"B", "all"}:
                        folded = fold_sliding_windows_bits(bits_B_real, m=m)
                        met_real = metrics_for_stream(bits_B_real, folded)
                        ent_rand = np.asarray(auric_ent_random.get(("B", m), []), dtype=np.float64)
                        ent_pert = np.asarray(auric_ent_pert.get(("B", m), []), dtype=np.float64)
                        auric_summary_rows.append(
                            {
                                "pdb_id": pdb_id,
                                "chain": chain_id,
                                "N": N,
                                "alphabet": alpha_name,
                                "readout": "B",
                                "m": m,
                                "type_entropy_real": float(met_real["type_entropy"]),
                                "type_entropy_random_mean": float(np.mean(ent_rand)) if len(ent_rand) else float("nan"),
                                "type_entropy_pert_mean": float(np.mean(ent_pert)) if len(ent_pert) else float("nan"),
                                "delta_type_entropy_real_vs_random": cliffs_delta_one_vs_many(float(met_real["type_entropy"]), ent_rand) if len(ent_rand) else float("nan"),
                                "delta_type_entropy_real_vs_perturbed": cliffs_delta_one_vs_many(float(met_real["type_entropy"]), ent_pert) if len(ent_pert) else float("nan"),
                                "type_support_real": float(met_real["type_support"]),
                                "run1_mean_real": float(met_real["run1_mean"]),
                                "run1_max_real": float(met_real["run1_max"]),
                                "smb_rate_hat_real": float(met_real["smb_rate_hat"]),
                            }
                        )
                        auric_samples_rows.append(
                            {
                                "pdb_id": pdb_id,
                                "chain": chain_id,
                                "N": N,
                                "alphabet": alpha_name,
                                "group": "real",
                                "rep": 0,
                                "readout": "B",
                                "m": m,
                                "bits_len": int(bits_B_real.shape[0]),
                                "types_len": int(folded.shape[0]),
                                **met_real,
                            }
                        )

        used += 1

    samples_df = pd.DataFrame(samples_rows)
    summary_df = pd.DataFrame(summary_rows)

    samples_csv = out_dir / f"phason_stats_samples_{args.tag}.csv"
    summary_csv = out_dir / f"phason_stats_summary_{args.tag}.csv"
    samples_df.to_csv(samples_csv, index=False)
    summary_df.to_csv(summary_csv, index=False)

    if args.auric:
        auric_samples_df = pd.DataFrame(auric_samples_rows)
        auric_summary_df = pd.DataFrame(auric_summary_rows)
        auric_samples_csv = out_dir / f"auric_stats_samples_{args.tag}.csv"
        auric_summary_csv = out_dir / f"auric_stats_summary_{args.tag}.csv"
        auric_samples_df.to_csv(auric_samples_csv, index=False)
        auric_summary_df.to_csv(auric_summary_csv, index=False)

    # Dataset-level report
    report_path = rep_dir / f"pdb_phason_stats_{args.tag}.md"
    ids_rel = str(ids_path.relative_to(root)).replace("\\", "/")
    samples_rel = str(samples_csv.relative_to(root)).replace("\\", "/")
    summary_rel = str(summary_csv.relative_to(root)).replace("\\", "/")
    proxy_desc = []
    proxy_desc.append("Operational proxy: oracle direction quantization → 6D walk → perp-space spread.")
    proxy_desc.append("Global: w0 = mean(y^⊥).")
    if args.piecewise_seg_len and args.piecewise_seg_len > 0:
        proxy_desc.append(f"Piecewise: per-segment w0 (seg_len={int(args.piecewise_seg_len)}).")
    if args.linear_drift:
        proxy_desc.append("Linear drift: fit w0(i)=b+a·i and measure residual.")

    lines = [
        "# PDB phason proxy statistics (large-scale)",
        "",
        " ".join(proxy_desc),
        "",
        "## Dataset",
        f"- IDs: `{ids_rel}`",
        f"- Used chains: {used}",
        f"- Length filter: N in [{args.min_len},{args.max_len}]",
        f"- Random reps per protein: {args.random_reps}",
        f"- Perturbed reps per protein: {args.perturb_reps} (noise={args.perturb_noise})",
        "",
        "## Outputs",
        f"- Samples CSV: `{samples_rel}`",
        f"- Summary CSV: `{summary_rel}`",
        "",
        "## Quick stats (median ph_rms_real)",
    ]
    for alpha_name in alphabets:
        sub = summary_df[summary_df["alphabet"] == alpha_name]
        if len(sub) == 0:
            continue
        lines.append(f"- {alpha_name}: {sub['ph_rms_real'].median():.3f} (n={len(sub)})")

    # Dataset-level mean Cliff's delta (bootstrap CI)
    lines += ["", "## Effect size (mean Cliff's delta across proteins; bootstrap CI)", ""]
    rng_ci = np.random.default_rng(args.seed + 123)
    for alpha_name in alphabets:
        sub = summary_df[summary_df["alphabet"] == alpha_name]
        if len(sub) == 0:
            continue
        def emit(name: str, col_r: str, col_p: str) -> None:
            d_rand = sub[col_r].dropna().to_numpy()
            d_pert = sub[col_p].dropna().to_numpy()
            lo_r, hi_r = bootstrap_ci_mean(d_rand, rng_ci)
            lo_p, hi_p = bootstrap_ci_mean(d_pert, rng_ci)
            lines.append(f"- {alpha_name} {name} real-vs-random: mean={np.mean(d_rand):.3f}, CI=[{lo_r:.3f},{hi_r:.3f}] (n={len(d_rand)})")
            lines.append(f"- {alpha_name} {name} real-vs-perturbed: mean={np.mean(d_pert):.3f}, CI=[{lo_p:.3f},{hi_p:.3f}] (n={len(d_pert)})")

        emit("global(ph_rms)", "cliffs_delta_real_vs_random", "cliffs_delta_real_vs_perturbed")
        if args.piecewise_seg_len and args.piecewise_seg_len > 0:
            emit("piecewise(ph_rms_pw)", "cliffs_delta_piecewise_vs_random", "cliffs_delta_piecewise_vs_perturbed")
        if args.linear_drift:
            emit("lin(ph_rms_lin)", "cliffs_delta_lin_vs_random", "cliffs_delta_lin_vs_perturbed")

    dt = time.perf_counter() - t0
    lines += ["", f"Runtime: {dt:.1f}s", ""]
    report_path.write_text("\n".join(lines) + "\n", encoding="utf-8")

    if args.auric:
        auric_report_path = rep_dir / f"pdb_auric_stats_{args.tag}.md"
        auric_samples_rel = str((out_dir / f"auric_stats_samples_{args.tag}.csv").relative_to(root)).replace("\\", "/")
        auric_summary_rel = str((out_dir / f"auric_stats_summary_{args.tag}.csv").relative_to(root)).replace("\\", "/")
        lines2 = [
            "# PDB auric Fold_m certificate metrics (large-scale)",
            "",
            "Auric protocol: ρ(t) → sliding windows → Fold_m → type stream metrics.",
            "",
            "## Dataset",
            f"- IDs: `{ids_rel}`",
            f"- Used chains: {used}",
            f"- Length filter: N in [{args.min_len},{args.max_len}]",
            f"- Alphabets: {', '.join(alphabets)}",
            f"- m list: {', '.join(str(m) for m in auric_ms)}",
            f"- readouts: {args.auric_readouts}",
            f"- Random reps per protein: {args.random_reps}",
            f"- Perturbed reps per protein: {args.perturb_reps} (noise={args.perturb_noise})",
            "",
            "## Outputs",
            f"- Samples CSV: `{auric_samples_rel}`",
            f"- Summary CSV: `{auric_summary_rel}`",
            "",
            "## Quick stats (median type_entropy_real; mean delta across proteins)",
        ]
        if len(auric_summary_df):
            for alpha_name in alphabets:
                suba = auric_summary_df[auric_summary_df["alphabet"] == alpha_name]
                for ro in (["A", "B"] if args.auric_readouts == "all" else [args.auric_readouts]):
                    subr = suba[suba["readout"] == ro]
                    for m in auric_ms:
                        subm = subr[subr["m"] == m]
                        if len(subm) == 0:
                            continue
                        med = float(np.nanmedian(subm["type_entropy_real"].to_numpy(dtype=np.float64)))
                        d_r = float(np.nanmean(subm["delta_type_entropy_real_vs_random"].to_numpy(dtype=np.float64)))
                        d_p = float(np.nanmean(subm["delta_type_entropy_real_vs_perturbed"].to_numpy(dtype=np.float64)))
                        lines2.append(f"- {alpha_name} ρ{ro} m={m}: median(H_type_real)={med:.3f}, mean δ(real,random)={d_r:.3f}, mean δ(real,pert)={d_p:.3f}")
        lines2.append("")
        auric_report_path.write_text("\n".join(lines2), encoding="utf-8")

    print(f"Wrote: {samples_csv}")
    print(f"Wrote: {summary_csv}")
    print(f"Wrote: {report_path}")
    if args.auric:
        print(f"Wrote: {auric_samples_csv}")  # type: ignore[name-defined]
        print(f"Wrote: {auric_summary_csv}")  # type: ignore[name-defined]
        print(f"Wrote: {auric_report_path}")  # type: ignore[name-defined]


if __name__ == "__main__":
    main()

