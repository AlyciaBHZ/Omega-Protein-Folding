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

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from Bio.PDB.PDBParser import PDBParser

# Allow running as a script without installing a package.
sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str((Path(__file__).resolve().parents[1] / "auric").resolve()))

from bench_utils import PHI, cliffs_delta_one_vs_many, ensure_dir, oracle_direction_reconstruct

from fold_m import fold_sliding_windows_bits
from metrics import metrics_for_stream
from readout import rho_A_from_yperp, rho_B_from_npath


def _pca_uvec(Y: np.ndarray) -> np.ndarray:
    Y = np.asarray(Y, dtype=np.float64)
    Yc = Y - Y.mean(axis=0, keepdims=True)
    if np.allclose(Yc, 0.0):
        return np.array([1.0, 0.0, 0.0], dtype=np.float64)
    _, _, Vt = np.linalg.svd(Yc, full_matrices=False)
    u = np.asarray(Vt[0], dtype=np.float64)
    nu = float(np.linalg.norm(u))
    u = (u / nu) if nu != 0.0 else np.array([1.0, 0.0, 0.0], dtype=np.float64)
    # deterministic sign
    for k in range(3):
        if abs(u[k]) > 1e-12:
            if u[k] < 0:
                u = -u
            break
    return u


def load_ca_coords_longest_chain_pdb(pdb_path: Path) -> Tuple[str, np.ndarray]:
    parser = PDBParser(QUIET=True)
    s = parser.get_structure(pdb_path.stem, str(pdb_path))
    model = next(iter(s.get_models()))

    best_chain = None
    best_coords = None
    for chain in model:
        coords = []
        for res in chain:
            if "CA" in res:
                coords.append(res["CA"].get_coord())
        if len(coords) >= 2:
            arr = np.asarray(coords, dtype=np.float64)
            if best_coords is None or len(arr) > len(best_coords):
                best_coords = arr
                best_chain = chain.id
    if best_chain is None or best_coords is None:
        raise ValueError("No CA trace found")
    return str(best_chain), best_coords


@dataclass(frozen=True)
class AuricCfg:
    alphabet: str
    ms: Tuple[int, ...]
    rhoA_threshold: str = "median"


def compute_auric_from_coords(coords: np.ndarray, *, cfg: AuricCfg, uvec: np.ndarray) -> Dict[Tuple[str, int, str], float]:
    """
    Return metric values keyed by (readout, m, metric_name) for readout in {A,B}.
    """
    _, n_path, y_perp, _, _ = oracle_direction_reconstruct(coords, alphabet=cfg.alphabet, phi=PHI)
    bits_A = rho_A_from_yperp(y_perp, u=tuple(uvec.tolist()), u_mode="fixed", threshold=cfg.rhoA_threshold)
    bits_B = rho_B_from_npath(n_path)

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
    ap.add_argument("--min-len", type=int, default=40)
    ap.add_argument("--max-len", type=int, default=400)
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

    for tid, rows in sorted(rows_by_target.items()):
        # find a native path from the first row
        native_path = (root / str(rows[0]["native_path"])).resolve()
        if not native_path.exists():
            continue
        try:
            native_chain, native_coords = load_ca_coords_longest_chain_pdb(native_path)
        except Exception:
            continue
        Nn = int(native_coords.shape[0])
        if Nn < int(args.min_len) or Nn > int(args.max_len):
            continue

        # shared PCA axis from native y_perp
        _, _, y_perp_native, _, _ = oracle_direction_reconstruct(native_coords, alphabet=cfg.alphabet, phi=PHI)
        uvec = _pca_uvec(y_perp_native)

        # native metrics
        met_native = compute_auric_from_coords(native_coords, cfg=cfg, uvec=uvec)

        # decoys metrics
        met_decoys: Dict[Tuple[str, int, str], List[float]] = {}
        n_decoys_used = 0
        for r in rows:
            decoy_path = (root / str(r["decoy_path"])).resolve()
            if not decoy_path.exists():
                continue
            try:
                _, decoy_coords = load_ca_coords_longest_chain_pdb(decoy_path)
            except Exception:
                continue
            Nd = int(decoy_coords.shape[0])
            if Nd < int(args.min_len) or Nd > int(args.max_len):
                continue
            try:
                met_d = compute_auric_from_coords(decoy_coords, cfg=cfg, uvec=uvec)
            except Exception:
                continue
            n_decoys_used += 1
            for k, v in met_d.items():
                met_decoys.setdefault(k, []).append(float(v))

            for (ro, m, metric_name), v in met_d.items():
                sample_rows.append(
                    {
                        "dataset": str(r.get("dataset", "")),
                        "target_id": tid,
                        "native_chain": native_chain,
                        "N_native": Nn,
                        "decoy_path": str(decoy_path.relative_to(root)).replace("\\", "/"),
                        "N_decoy": Nd,
                        "readout": ro,
                        "m": int(m),
                        "metric": metric_name,
                        "value": float(v),
                    }
                )

        # per-target summary (native vs decoys)
        out: Dict[str, object] = {
            "target_id": tid,
            "N_native": Nn,
            "n_decoys_used": n_decoys_used,
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

        # Plots: per-target decoy distributions with native marker (m=8 default)
        for metric_name in ("type_entropy", "smb_rate_hat"):
            for ro in ("A", "B"):
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
    lines.append(f"- Targets used: {int(sdf.shape[0])}")
    lines.append(f"- Samples CSV (gitignored): `{samples_csv.relative_to(root).as_posix()}`")
    lines.append(f"- Summary CSV (gitignored): `{summary_csv.relative_to(root).as_posix()}`")
    lines.append(f"- Plots: `docs/runs/{args.run_name}/plots/`")
    lines.append("")
    if len(sdf) == 0:
        lines.append("No targets produced usable decoy metrics (check parsing / length filters).")
    else:
        # A compact table: m=8 (or first m)
        m0 = 8 if 8 in cfg.ms else cfg.ms[0]
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
        ]
        lines.append("| " + " | ".join(hdr) + " |")
        lines.append("| " + " | ".join(["---"] * len(hdr)) + " |")
        for _, r in sdf.sort_values("target_id").iterrows():
            row = [
                str(r["target_id"]),
                str(int(r.get("n_decoys_used", 0))),
                f"{float(r.get(f'pct_native_among_decoys_type_entropy_rhoA_m{m0}', float('nan'))):.3f}",
                f"{float(r.get(f'delta_native_vs_decoys_type_entropy_rhoA_m{m0}', float('nan'))):.3f}",
                f"{float(r.get(f'pct_native_among_decoys_smb_rate_hat_rhoA_m{m0}', float('nan'))):.3f}",
                f"{float(r.get(f'delta_native_vs_decoys_smb_rate_hat_rhoA_m{m0}', float('nan'))):.3f}",
            ]
            lines.append("| " + " | ".join(row) + " |")
        lines.append("")

        # Aggregate: how often native is in lowest 10% of decoys
        pcol = f"pct_native_among_decoys_type_entropy_rhoA_m{m0}"
        if pcol in sdf.columns:
            frac10 = float(np.mean(sdf[pcol].to_numpy(dtype=np.float64) <= 0.10))
            lines.append(f"- Fraction targets with native in lowest 10% of decoys (rhoA type_entropy, m={m0}): {frac10:.3f}")
            lines.append("")

    report_md.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"Wrote: {samples_csv}")
    print(f"Wrote: {summary_csv}")
    print(f"Wrote: {report_md}")


if __name__ == "__main__":
    main()

