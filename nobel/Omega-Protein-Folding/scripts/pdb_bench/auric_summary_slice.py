#!/usr/bin/env python3
"""
Paper-facing analysis slice for auric_stats_summary_*.csv.

This summarizes (by m, readout) the effect sizes vs strong nulls and reports
cross-readout agreement corr(z_A, z_B) across proteins.

Intended output location: docs/runs/<run-name>/auric_summary_slice_<tag>.md
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


def _finite(x: pd.Series) -> np.ndarray:
    v = x.to_numpy(dtype=np.float64)
    return v[np.isfinite(v)]


def _mean_med(v: np.ndarray) -> tuple[float, float]:
    if v.size == 0:
        return float("nan"), float("nan")
    return float(np.mean(v)), float(np.median(v))


def _corr(a: np.ndarray, b: np.ndarray) -> float:
    if a.size == 0 or b.size == 0:
        return float("nan")
    m = np.isfinite(a) & np.isfinite(b)
    if int(np.sum(m)) < 3:
        return float("nan")
    aa = a[m]
    bb = b[m]
    sa = float(np.std(aa))
    sb = float(np.std(bb))
    if sa < 1e-12 or sb < 1e-12:
        return float("nan")
    return float(np.corrcoef(aa, bb)[0, 1])


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--summary", required=True, help="auric_stats_summary_*.csv (from phason_stats.py --auric)")
    ap.add_argument("--tag", required=True, help="Tag used for naming the output markdown.")
    ap.add_argument("--run-name", default="full1k_hybrid_auric", help="docs/runs/<run-name>/ output folder.")
    ap.add_argument("--alphabet", default="triple232")
    args = ap.parse_args()

    root = Path(__file__).resolve().parents[2]
    summary_path = Path(args.summary)
    if not summary_path.is_absolute():
        summary_path = (root / summary_path).resolve()

    out_dir = root / "docs" / "runs" / str(args.run_name)
    out_dir.mkdir(parents=True, exist_ok=True)
    out_md = out_dir / f"auric_summary_slice_{args.tag}.md"

    df = pd.read_csv(summary_path)
    if len(df) == 0:
        raise ValueError("Empty auric summary CSV")

    df = df[df["alphabet"].astype(str).str.lower() == str(args.alphabet).lower()].copy()
    if len(df) == 0:
        raise ValueError(f"No rows for alphabet={args.alphabet!r}")

    metrics = ["type_entropy", "smb_rate_hat"]
    ms = sorted(df["m"].unique().tolist())

    # Detect available blockshuffle families (k values) per metric.
    blk_ks: dict[str, list[int]] = {}
    for metric in metrics:
        ks = set()
        for c in df.columns:
            pref = f"delta_{metric}_real_vs_blockshuffle_k"
            if c.startswith(pref):
                # c like delta_type_entropy_real_vs_blockshuffle_k8
                try:
                    ks.add(int(c[len(pref) :]))
                except Exception:
                    pass
        blk_ks[metric] = sorted(ks)

    lines: list[str] = []
    lines.append(f"# Auric summary slice ({args.tag})")
    lines.append("")
    lines.append(f"- Source: `{summary_path.relative_to(root).as_posix()}`")
    lines.append(f"- Alphabet: `{args.alphabet}`")
    lines.append(f"- n proteins (unique pdb_id): {df['pdb_id'].astype(str).str.upper().nunique()}")
    lines.append("")

    def fmt(x: float) -> str:
        return "nan" if not np.isfinite(x) else f"{x:.3f}"

    # Per-m/readout stats vs shuffle and blockshuffle(k)
    for metric in metrics:
        lines.append(f"## Metric: `{metric}`")
        lines.append("")

        # Table header
        hdr = ["m", "readout", "mean δ(real,shuffle)", "median δ(real,shuffle)", "mean pct(real,shuffle)", "median pct(real,shuffle)"]
        for k in blk_ks[metric]:
            hdr += [
                f"mean δ(real,blk{k})",
                f"median δ(real,blk{k})",
                f"mean pct(real,blk{k})",
                f"median pct(real,blk{k})",
            ]
        lines.append("| " + " | ".join(hdr) + " |")
        lines.append("| " + " | ".join(["---"] * len(hdr)) + " |")

        for m in ms:
            for ro in ["A", "B"]:
                sub = df[(df["m"] == m) & (df["readout"] == ro)]
                if len(sub) == 0:
                    continue

                d_shuf = _finite(sub.get(f"delta_{metric}_real_vs_shuffle", pd.Series(dtype=float)))
                p_shuf = _finite(sub.get(f"pct_{metric}_real_vs_shuffle", pd.Series(dtype=float)))
                d_shuf_mu, d_shuf_med = _mean_med(d_shuf)
                p_shuf_mu, p_shuf_med = _mean_med(p_shuf)

                row = [str(int(m)), ro, fmt(d_shuf_mu), fmt(d_shuf_med), fmt(p_shuf_mu), fmt(p_shuf_med)]

                for k in blk_ks[metric]:
                    dcol = f"delta_{metric}_real_vs_blockshuffle_k{k}"
                    pcol = f"pct_{metric}_real_vs_blockshuffle_k{k}"
                    dv = _finite(sub.get(dcol, pd.Series(dtype=float)))
                    pv = _finite(sub.get(pcol, pd.Series(dtype=float)))
                    dmu, dmed = _mean_med(dv)
                    pmu, pmed = _mean_med(pv)
                    row += [fmt(dmu), fmt(dmed), fmt(pmu), fmt(pmed)]

                lines.append("| " + " | ".join(row) + " |")

        lines.append("")

        # Cross-readout agreement corr(z_A, z_B) across proteins, per m
        zcol = f"z_{metric}_real_vs_shuffle"
        if zcol in df.columns:
            lines.append("### Cross-readout agreement (shuffle z-scores)")
            lines.append("")
            lines.append("| m | corr(z_A, z_B) | n pairs |")
            lines.append("| --- | --- | --- |")
            for m in ms:
                sub = df[df["m"] == m]
                a = sub[sub["readout"] == "A"][["pdb_id", zcol]].rename(columns={zcol: "zA"})
                b = sub[sub["readout"] == "B"][["pdb_id", zcol]].rename(columns={zcol: "zB"})
                j = a.merge(b, on="pdb_id", how="inner")
                zA = j["zA"].to_numpy(dtype=np.float64)
                zB = j["zB"].to_numpy(dtype=np.float64)
                c = _corr(zA, zB)
                n = int(np.sum(np.isfinite(zA) & np.isfinite(zB)))
                lines.append(f"| {int(m)} | {fmt(c)} | {n} |")
            lines.append("")

    out_md.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"Wrote: {out_md}")


if __name__ == "__main__":
    main()

