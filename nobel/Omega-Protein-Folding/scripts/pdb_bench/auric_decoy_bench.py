#!/usr/bin/env python3
"""
Auric Native-vs-Decoy benchmark (hard-negative control).

Goal:
- Compare native structures to realistic "near-native but wrong" decoys (e.g. Decoys 'R' Us)
  using the same Auric hybrid protocol as the full1k run:
    - rhoA: shared-PCA axis u computed from the native y_perp, threshold=median
    - rhoB: parity from the 6D step stream (n_path)
    - m list: 6,8,10

Inputs:
- A manifest CSV with columns: dataset,target_id,native_path,decoy_path
  (produced by scripts/pdb_bench/decoys_r_us_download.py)

Outputs:
- data/processed/decoy_bench/auric_decoys_samples_{tag}.csv (long-form; gitignored)
- data/processed/decoy_bench/auric_decoys_summary_{tag}.csv (per-target; gitignored)
- docs/runs/<run-name>/decoy_auric_report_{tag}.md (committable report)
- docs/runs/<run-name>/plots/*.png (committable plots)
"""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Tuple

import sys

import matplotlib
matplotlib.use("Agg")  # Force non-interactive backend (avoid GUI stalls on Windows).
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# Allow running as a script without installing a package.
sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str((Path(__file__).resolve().parents[1] / "auric").resolve()))

from bench_utils import PHI, cliffs_delta_one_vs_many, ensure_dir, oracle_direction_reconstruct
from bench_utils import geometry_qc_ca_trace, icosa_B
from coords_io import load_ca_coords

from axis import best_native_axis, shared_pca_axis
from fold_m import fold_sliding_windows_bits
from metrics import metrics_for_stream
from readout import rho_A_from_yperp, rho_B_from_npath


def _resolve_under_root(root: Path, p: str | Path) -> Path:
    pp = Path(str(p))
    return pp if pp.is_absolute() else (root / pp).resolve()


def _relpath_posix(p: Path, root: Path) -> str:
    try:
        return p.relative_to(root).as_posix()
    except Exception:
        return str(p)


@dataclass(frozen=True)
class AuricCfg:
    alphabet: str
    ms: Tuple[int, ...]
    rhoA_threshold: str = "median"


def _lift_from_coords(
    coords: np.ndarray,
    *,
    alphabet: str,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, float, float]:
    """
    Lift a CA trace to (n_path, y_perp), plus Bperp (for rhoB=vel) and phason proxies.
    Returns: (n_path, y_perp, Bperp, ph_rms, ph_max)
    """
    coords = np.asarray(coords, dtype=np.float64)
    _, n_path, y_perp, ph_rms, ph_max = oracle_direction_reconstruct(coords, alphabet=str(alphabet), phi=PHI)

    d = coords[1:] - coords[:-1]
    bond_lengths = np.linalg.norm(d, axis=1)
    median_len = float(np.median(bond_lengths)) if bond_lengths.size else 1.0
    Bperp = icosa_B(-1.0 / PHI) * median_len  # (3,6)
    return n_path, y_perp, Bperp, float(ph_rms), float(ph_max)


def compute_auric_from_lift(
    n_path: np.ndarray,
    y_perp: np.ndarray,
    Bperp: np.ndarray,
    *,
    cfg: AuricCfg,
    uvec: np.ndarray,
) -> Dict[Tuple[str, int, str], float]:
    """
    Return metric values keyed by (readout, m, metric_name) for readout in {A,B}.
    """
    bits_A = rho_A_from_yperp(y_perp, u=tuple(uvec.tolist()), u_mode="fixed", threshold=cfg.rhoA_threshold)
    bits_B = rho_B_from_npath(n_path, mode="parity")
    bits_Bvel = rho_B_from_npath(
        n_path,
        mode="vel",
        Bperp=Bperp,
        u=tuple(uvec.tolist()),
        u_mode="fixed",
        threshold="median",
    )

    out: Dict[Tuple[str, int, str], float] = {}
    for m in cfg.ms:
        folded_A = fold_sliding_windows_bits(bits_A, m=m)
        met_A = metrics_for_stream(bits_A, folded_A)
        out[("A", m, "type_entropy")] = float(met_A["type_entropy"])
        out[("A", m, "smb_rate_hat")] = float(met_A["smb_rate_hat"])

        folded_B = fold_sliding_windows_bits(bits_B, m=m)
        met_B = metrics_for_stream(bits_B, folded_B)
        out[("B", m, "type_entropy")] = float(met_B["type_entropy"])
        out[("B", m, "smb_rate_hat")] = float(met_B["smb_rate_hat"])

        folded_Bv = fold_sliding_windows_bits(bits_Bvel, m=m)
        met_Bv = metrics_for_stream(bits_Bvel, folded_Bv)
        out[("Bvel", m, "type_entropy")] = float(met_Bv["type_entropy"])
        out[("Bvel", m, "smb_rate_hat")] = float(met_Bv["smb_rate_hat"])
    return out


def percentile_among_decoys(x_native: float, xs_decoy: np.ndarray) -> float:
    xs = np.asarray(xs_decoy, dtype=np.float64)
    xs = xs[np.isfinite(xs)]
    if xs.size == 0 or not np.isfinite(x_native):
        return float("nan")
    return float(np.mean(xs <= float(x_native)))


def _iter_manifest_rows(path: Path) -> Iterable[dict]:
    with path.open("r", encoding="utf-8", newline="") as f:
        r = csv.DictReader(f)
        for row in r:
            yield row


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--manifest", required=True, help="decoy_manifest.csv from decoys_r_us_download.py")
    ap.add_argument("--run-name", default="decoys_4state_auric", help="docs/runs/<run-name>/ output folder.")
    ap.add_argument("--tag", default="decoys_4state_reduced_triple232_hybrid", help="Output tag suffix.")
    ap.add_argument("--alphabet", default="triple232", choices=["axis12", "pair72", "triple232"])
    ap.add_argument("--m", default="6,8,10", help="Comma list of m values.")
    ap.add_argument(
        "--axis-mode",
        default="shared_pca",
        choices=["shared_pca", "best_native_axis"],
        help="Axis selection mode for rhoA: shared_pca (default) or strict native-only best_native_axis.",
    )
    ap.add_argument(
        "--ph-fail-rms",
        type=float,
        default=100.0,
        help="Flag/skip entries where ph_rms exceeds this (FailedGeometry via phason proxy).",
    )
    ap.add_argument("--min-len", type=int, default=40)
    ap.add_argument("--max-len", type=int, default=400)
    ap.add_argument("--max-targets", type=int, default=0, help="If >0, limit number of targets (debug/quick runs).")
    ap.add_argument("--max-decoys-per-target", type=int, default=0, help="If >0, cap decoys per target (debug/quick runs).")
    ap.add_argument("--no-plots", action="store_true", help="If set, skip plot generation (fast/headless).")
    ap.add_argument("--progress-every", type=int, default=200, help="Progress print frequency in decoy loop.")
    ap.add_argument("--qc-bond-nominal", type=float, default=3.8, help="Nominal CA-CA bond length (Å) for QC.")
    ap.add_argument("--qc-bond-tol", type=float, default=0.5, help="Outlier tolerance for |d-3.8| (Å) in QC.")
    ap.add_argument("--qc-ca-ca-max-fail", type=float, default=10.0, help="Fail if max consecutive CA-CA distance exceeds this (Å).")
    ap.add_argument(
        "--qc-outlier-frac-fail",
        type=float,
        default=0.20,
        help="Fail if fraction of consecutive bonds outside tolerance exceeds this.",
    )
    args = ap.parse_args()

    root = Path(__file__).resolve().parents[2]
    man = Path(args.manifest)
    if not man.is_absolute():
        man = (root / man).resolve()

    ms = tuple(int(x.strip()) for x in str(args.m).split(",") if x.strip())
    cfg = AuricCfg(alphabet=str(args.alphabet), ms=ms)

    out_dir = root / "data" / "processed" / "decoy_bench"
    ensure_dir(out_dir)
    run_dir = root / "docs" / "runs" / str(args.run_name)
    plots_dir = run_dir / "plots"
    ensure_dir(run_dir)
    ensure_dir(plots_dir)

    samples_csv = out_dir / f"auric_decoys_samples_{args.tag}.csv"
    summary_csv = out_dir / f"auric_decoys_summary_{args.tag}.csv"
    report_md = run_dir / f"decoy_auric_report_{args.tag}.md"

    # Group manifest rows by target
    rows_by_target: Dict[str, List[dict]] = {}
    for row in _iter_manifest_rows(man):
        tid = str(row.get("target_id", "")).strip()
        if not tid:
            continue
        rows_by_target.setdefault(tid, []).append(row)

    sample_rows: List[Dict[str, object]] = []
    summary_rows: List[Dict[str, object]] = []

    items = sorted(rows_by_target.items())
    if int(args.max_targets) > 0:
        items = items[: int(args.max_targets)]

    print(f"Loaded manifest targets={len(rows_by_target)}; running targets={len(items)}", flush=True)

    for tid, rows in items:
        # find a native path from the first row
        native_path = _resolve_under_root(root, str(rows[0].get("native_path", "")).strip())
        if not native_path.exists():
            continue
        try:
            native_chain, native_coords = load_ca_coords(native_path)
        except Exception:
            continue
        Nn = int(native_coords.shape[0])
        if Nn < int(args.min_len) or Nn > int(args.max_len):
            continue

        qc_native = geometry_qc_ca_trace(
            native_coords,
            bond_nominal=float(args.qc_bond_nominal),
            bond_tol=float(args.qc_bond_tol),
            ca_ca_max_fail=float(args.qc_ca_ca_max_fail),
            outlier_frac_fail=float(args.qc_outlier_frac_fail),
        )
        native_failed_geometry = bool(qc_native.get("failed_geometry", False))

        # Lift native once (also supplies y_perp for axis selection + phason proxy QC).
        try:
            n_path_native, y_perp_native, Bperp_native, ph_rms_native, ph_max_native = _lift_from_coords(
                native_coords, alphabet=cfg.alphabet
            )
        except Exception:
            continue

        native_failed_ph = bool(np.isfinite(ph_rms_native) and (ph_rms_native > float(args.ph_fail_rms)))
        if native_failed_ph:
            # If the native itself fails the phason-geometry check, skip this target:
            # axis selection and separation claims are not meaningful for a broken lift.
            continue

        # Axis selection from NATIVE only (protocol-safe; no decoy leakage).
        m_axis = 8 if 8 in cfg.ms else cfg.ms[0]
        if str(args.axis_mode) == "best_native_axis":
            uvec = best_native_axis(y_perp_native, m=int(m_axis), threshold=str(cfg.rhoA_threshold), metric="type_entropy")
        else:
            uvec = shared_pca_axis(y_perp_native)

        # native metrics
        met_native = compute_auric_from_lift(n_path_native, y_perp_native, Bperp_native, cfg=cfg, uvec=uvec)

        # decoys metrics
        met_decoys: Dict[Tuple[str, int, str], List[float]] = {}
        n_decoys_used = 0
        n_decoys_failed_geometry = 0
        n_decoys_failed_ph = 0
        rows_dec = rows
        if int(args.max_decoys_per_target) > 0:
            rows_dec = rows_dec[: int(args.max_decoys_per_target)]

        for j, r in enumerate(rows_dec, start=1):
            decoy_path = _resolve_under_root(root, str(r.get("decoy_path", "")).strip())
            if not decoy_path.exists():
                continue
            try:
                _, decoy_coords = load_ca_coords(decoy_path)
            except Exception:
                continue
            Nd = int(decoy_coords.shape[0])
            if Nd < int(args.min_len) or Nd > int(args.max_len):
                continue

            qc_decoy = geometry_qc_ca_trace(
                decoy_coords,
                bond_nominal=float(args.qc_bond_nominal),
                bond_tol=float(args.qc_bond_tol),
                ca_ca_max_fail=float(args.qc_ca_ca_max_fail),
                outlier_frac_fail=float(args.qc_outlier_frac_fail),
            )
            if bool(qc_decoy.get("failed_geometry", False)):
                n_decoys_failed_geometry += 1
                # Keep a lightweight record (QC-only row) for audit, but exclude from Auric stats.
                sample_rows.append(
                    {
                        "dataset": str(r.get("dataset", "(unknown)")),
                        "target_id": tid,
                        "native_chain": native_chain,
                        "N_native": Nn,
                        "decoy_path": _relpath_posix(decoy_path, root),
                        "N_decoy": Nd,
                        "readout": "QC",
                        "m": -1,
                        "metric": "failed_geometry",
                        "value": 1.0,
                        "ph_rms": float("nan"),
                        "ph_max": float("nan"),
                        "ca_ca_median": float(qc_decoy.get("ca_ca_median", float("nan"))),
                        "ca_ca_max": float(qc_decoy.get("ca_ca_max", float("nan"))),
                        "ca_ca_outlier_frac": float(qc_decoy.get("ca_ca_outlier_frac", float("nan"))),
                    }
                )
                continue

            try:
                n_path_d, y_perp_d, Bperp_d, ph_rms_d, ph_max_d = _lift_from_coords(decoy_coords, alphabet=cfg.alphabet)
            except Exception:
                continue
            if np.isfinite(ph_rms_d) and (float(ph_rms_d) > float(args.ph_fail_rms)):
                n_decoys_failed_ph += 1
                sample_rows.append(
                    {
                        "dataset": str(r.get("dataset", "(unknown)")),
                        "target_id": tid,
                        "native_chain": native_chain,
                        "N_native": Nn,
                        "decoy_path": _relpath_posix(decoy_path, root),
                        "N_decoy": Nd,
                        "readout": "QC",
                        "m": -1,
                        "metric": "failed_ph_rms",
                        "value": 1.0,
                        "ph_rms": float(ph_rms_d),
                        "ph_max": float(ph_max_d),
                        "ca_ca_median": float(qc_decoy.get("ca_ca_median", float("nan"))),
                        "ca_ca_max": float(qc_decoy.get("ca_ca_max", float("nan"))),
                        "ca_ca_outlier_frac": float(qc_decoy.get("ca_ca_outlier_frac", float("nan"))),
                    }
                )
                continue

            try:
                met_d = compute_auric_from_lift(n_path_d, y_perp_d, Bperp_d, cfg=cfg, uvec=uvec)
            except Exception:
                continue
            n_decoys_used += 1
            for k, v in met_d.items():
                met_decoys.setdefault(k, []).append(float(v))

            for (ro, m, metric_name), v in met_d.items():
                sample_rows.append(
                    {
                        "dataset": str(r.get("dataset", "(unknown)")),
                        "target_id": tid,
                        "native_chain": native_chain,
                        "N_native": Nn,
                        "decoy_path": _relpath_posix(decoy_path, root),
                        "N_decoy": Nd,
                        "readout": ro,
                        "m": int(m),
                        "metric": metric_name,
                        "value": float(v),
                        "ph_rms": float(ph_rms_d),
                        "ph_max": float(ph_max_d),
                        "ca_ca_median": float(qc_decoy.get("ca_ca_median", float("nan"))),
                        "ca_ca_max": float(qc_decoy.get("ca_ca_max", float("nan"))),
                        "ca_ca_outlier_frac": float(qc_decoy.get("ca_ca_outlier_frac", float("nan"))),
                    }
                )

            if int(args.progress_every) > 0 and (j % int(args.progress_every) == 0):
                print(f"[{tid}] processed decoys: {j}/{len(rows_dec)} (used={n_decoys_used})", flush=True)

        # per-target summary (native vs decoys)
        out: Dict[str, object] = {
            "target_id": tid,
            "N_native": Nn,
            "native_failed_geometry": int(native_failed_geometry),
            "n_decoys_used": n_decoys_used,
            "n_decoys_failed_geometry": int(n_decoys_failed_geometry),
            "n_decoys_failed_ph_rms": int(n_decoys_failed_ph),
            "native_ca_ca_median": float(qc_native.get("ca_ca_median", float("nan"))),
            "native_ca_ca_max": float(qc_native.get("ca_ca_max", float("nan"))),
            "native_ca_ca_outlier_frac": float(qc_native.get("ca_ca_outlier_frac", float("nan"))),
            "native_ph_rms": float(ph_rms_native),
            "native_ph_max": float(ph_max_native),
            "axis_mode": str(args.axis_mode),
            "uvec_x": float(uvec[0]),
            "uvec_y": float(uvec[1]),
            "uvec_z": float(uvec[2]),
        }
        for (ro, m, metric_name), xnat in met_native.items():
            ys = np.asarray(met_decoys.get((ro, m, metric_name), []), dtype=np.float64)
            out[f"{metric_name}_native_rho{ro}_m{m}"] = float(xnat)
            out[f"{metric_name}_decoy_mean_rho{ro}_m{m}"] = float(np.nanmean(ys)) if ys.size else float("nan")
            out[f"delta_native_vs_decoys_{metric_name}_rho{ro}_m{m}"] = float(cliffs_delta_one_vs_many(float(xnat), ys)) if ys.size else float("nan")
            out[f"pct_native_among_decoys_{metric_name}_rho{ro}_m{m}"] = percentile_among_decoys(float(xnat), ys)
        summary_rows.append(out)

        if not args.no_plots:
            # Plots: per-target decoy distributions with native marker (m=8 default)
            for metric_name in ("type_entropy", "smb_rate_hat"):
                for ro in ("A", "B", "Bvel"):
                    m = 8 if 8 in cfg.ms else cfg.ms[0]
                    ys = np.asarray(met_decoys.get((ro, m, metric_name), []), dtype=np.float64)
                    ys = ys[np.isfinite(ys)]
                    if ys.size == 0:
                        continue
                    xnat = float(met_native[(ro, m, metric_name)])
                    fig, ax = plt.subplots(figsize=(6.8, 3.6))
                    ax.hist(ys, bins=40, alpha=0.85, label="decoys")
                    ax.axvline(xnat, color="k", linestyle="--", linewidth=2, label="native")
                    ax.set_title(f"{tid}: {metric_name} (rho{ro}, m={m})")
                    ax.set_xlabel(metric_name)
                    ax.set_ylabel("count")
                    ax.legend(frameon=False)
                    fig.tight_layout()
                    outp = plots_dir / f"decoys_hist_{metric_name}_{tid}_rho{ro}_m{m}_{args.tag}.png"
                    fig.savefig(outp, dpi=200)
                    plt.close(fig)

    # Write CSVs (gitignored)
    pd.DataFrame(sample_rows).to_csv(samples_csv, index=False)
    sdf = pd.DataFrame(summary_rows)
    sdf.to_csv(summary_csv, index=False)

    # Write report (committable)
    lines: List[str] = []
    lines.append(f"# Auric decoy benchmark ({args.tag})")
    lines.append("")
    lines.append(f"- Manifest: `{man.relative_to(root).as_posix()}`")
    lines.append(f"- Alphabet: `{cfg.alphabet}`")
    lines.append(f"- m list: {', '.join(str(m) for m in cfg.ms)}")
    lines.append(f"- Axis mode: `{str(args.axis_mode)}`")
    lines.append(f"- FailedGeometry (phason): ph_rms <= {float(args.ph_fail_rms):.1f}")
    lines.append(f"- Targets used: {int(sdf.shape[0])}")
    lines.append(f"- Samples CSV (gitignored): `{samples_csv.relative_to(root).as_posix()}`")
    lines.append(f"- Summary CSV (gitignored): `{summary_csv.relative_to(root).as_posix()}`")
    lines.append(f"- Plots: `docs/runs/{args.run_name}/plots/`")
    lines.append(f"- Geometry QC: max CA-CA <= {float(args.qc_ca_ca_max_fail):.1f} Å and outlier_frac <= {float(args.qc_outlier_frac_fail):.2f}")
    lines.append("")
    if len(sdf) == 0:
        lines.append("No targets produced usable decoy metrics (check parsing / length filters).")
    else:
        # A compact table: m=8 (or first m)
        m0 = 8 if 8 in cfg.ms else cfg.ms[0]
        if "native_failed_geometry" in sdf.columns:
            n_fail_nat = int(np.sum(sdf["native_failed_geometry"].to_numpy(dtype=np.int64) > 0))
            lines.append(f"- Targets with native FailedGeometry under this QC: {n_fail_nat}/{int(sdf.shape[0])}")
            lines.append("")
        lines.append("## Native vs Decoy separation (summary)")
        lines.append("")
        lines.append(f"Using m={m0}. δ<0 means native has lower metric (more ordered) than decoys.")
        lines.append("")
        hdr = [
            "target_id",
            "n_decoys_used",
            f"pct_native_type_entropy_rhoA_m{m0}",
            f"delta_native_vs_decoys_type_entropy_rhoA_m{m0}",
            f"pct_native_smb_rate_hat_rhoA_m{m0}",
            f"delta_native_vs_decoys_smb_rate_hat_rhoA_m{m0}",
            f"pct_native_type_entropy_rhoB_m{m0}",
            f"pct_native_type_entropy_rhoBvel_m{m0}",
            f"pct_combined_min(A,B)_type_entropy_m{m0}",
        ]
        lines.append("| " + " | ".join(hdr) + " |")
        lines.append("| " + " | ".join(["---"] * len(hdr)) + " |")
        for _, r in sdf.sort_values("target_id").iterrows():
            pA = float(r.get(f"pct_native_among_decoys_type_entropy_rhoA_m{m0}", float("nan")))
            pB = float(r.get(f"pct_native_among_decoys_type_entropy_rhoB_m{m0}", float("nan")))
            pBv = float(r.get(f"pct_native_among_decoys_type_entropy_rhoBvel_m{m0}", float("nan")))
            pComb = float(np.nanmin(np.asarray([pA, pB], dtype=np.float64)))
            row = [
                str(r["target_id"]),
                f"{int(r.get('n_decoys_used', 0))} (failQC={int(r.get('n_decoys_failed_geometry', 0))})",
                f"{float(r.get(f'pct_native_among_decoys_type_entropy_rhoA_m{m0}', float('nan'))):.3f}",
                f"{float(r.get(f'delta_native_vs_decoys_type_entropy_rhoA_m{m0}', float('nan'))):.3f}",
                f"{float(r.get(f'pct_native_among_decoys_smb_rate_hat_rhoA_m{m0}', float('nan'))):.3f}",
                f"{float(r.get(f'delta_native_vs_decoys_smb_rate_hat_rhoA_m{m0}', float('nan'))):.3f}",
                f"{pB:.3f}",
                f"{pBv:.3f}",
                f"{pComb:.3f}",
            ]
            lines.append("| " + " | ".join(row) + " |")
        lines.append("")

        # Aggregate: how often native is in lowest 10% of decoys
        pAcol = f"pct_native_among_decoys_type_entropy_rhoA_m{m0}"
        pBcol = f"pct_native_among_decoys_type_entropy_rhoB_m{m0}"
        pBvcol = f"pct_native_among_decoys_type_entropy_rhoBvel_m{m0}"
        if pAcol in sdf.columns:
            frac10A = float(np.mean(sdf[pAcol].to_numpy(dtype=np.float64) <= 0.10))
            lines.append(f"- Fraction targets with native in lowest 10% of decoys (rhoA type_entropy, m={m0}): {frac10A:.3f}")
        if pBcol in sdf.columns:
            frac10B = float(np.mean(sdf[pBcol].to_numpy(dtype=np.float64) <= 0.10))
            lines.append(f"- Fraction targets with native in lowest 10% of decoys (rhoB=parity type_entropy, m={m0}): {frac10B:.3f}")
        if pBvcol in sdf.columns:
            frac10Bv = float(np.mean(sdf[pBvcol].to_numpy(dtype=np.float64) <= 0.10))
            lines.append(f"- Fraction targets with native in lowest 10% of decoys (rhoB=vel type_entropy, m={m0}): {frac10Bv:.3f}")
        if pAcol in sdf.columns and pBcol in sdf.columns:
            pComb = np.nanmin(np.stack([sdf[pAcol].to_numpy(dtype=np.float64), sdf[pBcol].to_numpy(dtype=np.float64)], axis=1), axis=1)
            frac10Comb = float(np.mean(pComb <= 0.10))
            lines.append(f"- Fraction targets with native in lowest 10% by min(A,Bparity) (type_entropy, m={m0}): {frac10Comb:.3f}")
        lines.append("")

    report_md.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"Wrote: {samples_csv}")
    print(f"Wrote: {summary_csv}")
    print(f"Wrote: {report_md}")


if __name__ == "__main__":
    main()

