#!/usr/bin/env python3
"""
Blind Sprint ablation suite (baseline / +A / +B / +A+B; multi-seed) with seed-level parallelism.

This suite shells out to:
  scripts/blind_sprint/blind_sprint_auric.py

Outputs:
- data/processed/blind_sprint/blind_sprint_ablations_<tag>.csv   (gitignored)
- docs/runs/<run-dir>/<tag>.md                                  (committable)
- docs/runs/<run-dir>/plots/*.png                               (committable)
"""

from __future__ import annotations

import argparse
import concurrent.futures as cf
import json
import math
import subprocess
import time
from pathlib import Path
from typing import List, Tuple

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def _run_one_seed(
    pdb_id: str,
    *,
    seed: int,
    codec: str,
    beam: int,
    K: int,
    wA: float,
    wB: float,
    rhoB_mode: str,
    rhoB_threshold: str,
    run_dir: Path,
    out_csv: Path,
    tag: str,
) -> pd.DataFrame:
    cmd = [
        "python",
        "scripts/blind_sprint/blind_sprint_auric.py",
        "--pdb-id",
        str(pdb_id),
        "--seed",
        str(int(seed)),
        "--codec",
        str(codec),
        "--beam",
        str(int(beam)),
        "--K",
        str(int(K)),
        "--wA",
        str(float(wA)),
        "--wB",
        str(float(wB)),
        "--auric-rhoB-mode",
        str(rhoB_mode),
        "--auric-rhoB-threshold",
        str(rhoB_threshold),
        "--run-dir",
        str(run_dir.as_posix()),
        "--out-csv",
        str(out_csv.as_posix()),
        "--no-report",
        "--quiet",
        "--tag",
        str(tag),
    ]
    subprocess.run(cmd, check=True)
    return pd.read_csv(out_csv)


def _boxplot(df: pd.DataFrame, *, out_path: Path, title: str) -> None:
    order = ["baseline", "A", "Bparity", "A+Bparity", "Bvel", "A+Bvel"]
    present = [g for g in order if g in set(df["group"])]
    data = [df[df["group"] == g]["tm"].to_numpy(dtype=np.float64) for g in present]
    fig, ax = plt.subplots(figsize=(7.8, 3.8))
    ax.boxplot(data, labels=present, showfliers=False)
    ax.set_ylabel("TM-score (eval only)")
    ax.set_title(title)
    ax.set_ylim(0.0, 1.02)
    fig.tight_layout()
    fig.savefig(out_path, dpi=220)
    plt.close(fig)


def _ecdf(df: pd.DataFrame, *, out_path: Path, title: str) -> None:
    fig, ax = plt.subplots(figsize=(7.0, 4.0))
    for g, sub in df.groupby("group"):
        x = np.sort(sub["tm"].to_numpy(dtype=np.float64))
        y = np.arange(1, len(x) + 1) / max(1, len(x))
        ax.step(x, y, where="post", label=str(g))
    ax.set_xlabel("TM-score (eval only)")
    ax.set_ylabel("ECDF")
    ax.set_title(title)
    ax.set_xlim(0.0, 1.0)
    ax.set_ylim(0.0, 1.0)
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(out_path, dpi=220)
    plt.close(fig)


def _df_to_md(df: pd.DataFrame) -> str:
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
                cells.append(f"{v:.3f}" if math.isfinite(v) else "nan")
            else:
                cells.append(str(v))
        lines.append("| " + " | ".join(cells) + " |")
    return "\n".join(lines)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--ids", default="1CTF,1R69,1SN3,2CRO,3ICB,4PTI,4RXN", help="Comma PDB IDs.")
    ap.add_argument("--seeds", default="0,1,2,3,4", help="Comma seeds.")
    ap.add_argument("--codec", default="pair72", choices=["pair72", "triple232"])
    ap.add_argument("--beam", type=int, default=128)
    ap.add_argument("--K", type=int, default=64)
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--tag", default="blind_sprint_ablations_4state_v2")
    ap.add_argument("--run-dir", default="docs/runs/blind_sprint_ablations_4state_v2")
    ap.add_argument("--wA", type=float, default=0.2)
    ap.add_argument("--wB", type=float, default=0.2)
    args = ap.parse_args()

    ids = [x.strip().upper() for x in str(args.ids).split(",") if x.strip()]
    seeds = [int(x.strip()) for x in str(args.seeds).split(",") if x.strip()]

    root = Path(__file__).resolve().parents[2]
    run_dir = Path(str(args.run_dir))
    if not run_dir.is_absolute():
        run_dir = (root / run_dir).resolve()
    plots_dir = run_dir / "plots"
    run_dir.mkdir(parents=True, exist_ok=True)
    plots_dir.mkdir(parents=True, exist_ok=True)

    groups: List[Tuple[str, float, float, str, str]] = [
        ("baseline", 0.0, 0.0, "parity", "median"),
        ("A", float(args.wA), 0.0, "parity", "median"),
        ("Bparity", 0.0, float(args.wB), "parity", "median"),
        ("A+Bparity", float(args.wA), float(args.wB), "parity", "median"),
        ("Bvel", 0.0, float(args.wB), "vel", "median"),
        ("A+Bvel", float(args.wA), float(args.wB), "vel", "median"),
    ]

    tmp_dir = root / "data" / "processed" / "blind_sprint" / f"_tmp_ablations_{args.tag}"
    tmp_dir.mkdir(parents=True, exist_ok=True)

    jobs: List[Tuple[str, str, int, float, float, str, str]] = []
    for pdb_id in ids:
        for gname, wA, wB, rhoB_mode, rhoB_thr in groups:
            for seed in seeds:
                jobs.append((pdb_id, gname, int(seed), float(wA), float(wB), str(rhoB_mode), str(rhoB_thr)))

    t0 = time.perf_counter()
    rows_all: List[pd.DataFrame] = []
    max_workers = max(1, int(args.workers))
    with cf.ThreadPoolExecutor(max_workers=max_workers) as ex:
        fut_to_job = {}
        for pdb_id, gname, seed, wA, wB, rhoB_mode, rhoB_thr in jobs:
            tag2 = f"{args.tag}__{gname}"
            out_csv = tmp_dir / f"{pdb_id}_{gname}_seed{seed}_{rhoB_mode}.csv"
            fut = ex.submit(
                _run_one_seed,
                pdb_id,
                seed=int(seed),
                codec=str(args.codec),
                beam=int(args.beam),
                K=int(args.K),
                wA=float(wA),
                wB=float(wB),
                rhoB_mode=str(rhoB_mode),
                rhoB_threshold=str(rhoB_thr),
                run_dir=run_dir,
                out_csv=out_csv,
                tag=tag2,
            )
            fut_to_job[fut] = (pdb_id, gname, seed, rhoB_mode)

        for fut in cf.as_completed(list(fut_to_job.keys())):
            pdb_id, gname, seed, rhoB_mode = fut_to_job[fut]
            df = fut.result()
            df["group"] = str(gname)
            df["pdb_id"] = str(pdb_id)
            df["seed"] = int(seed)
            df["rhoB_mode"] = str(rhoB_mode)
            rows_all.append(df)
            print(f"[done] {pdb_id} group={gname} seed={seed} rhoB={rhoB_mode}", flush=True)

    df_all = pd.concat(rows_all, axis=0, ignore_index=True)
    out_csv = root / "data" / "processed" / "blind_sprint" / f"blind_sprint_ablations_{args.tag}.csv"
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    df_all.to_csv(out_csv, index=False)

    best = df_all.groupby(["pdb_id", "group"], as_index=False)["tm"].max().rename(columns={"tm": "tm_best"})
    best_pivot = best.pivot(index="pdb_id", columns="group", values="tm_best")

    frac_rows = []
    for g in sorted(best["group"].unique()):
        sub = best[best["group"] == g]
        frac = float(np.mean(sub["tm_best"].to_numpy(dtype=np.float64) > 0.6))
        frac_rows.append({"group": g, "frac_best_tm_gt_0p6": frac, "n_targets": int(sub.shape[0])})
    frac_df = pd.DataFrame(frac_rows)

    _boxplot(df_all, out_path=plots_dir / f"tm_box_{args.tag}.png", title=f"Blind Sprint ablations ({args.tag})")
    _ecdf(df_all, out_path=plots_dir / f"tm_ecdf_{args.tag}.png", title=f"Blind Sprint ablations ({args.tag})")

    report = run_dir / f"{args.tag}.md"
    meta = {
        "ids": ids,
        "seeds": seeds,
        "codec": str(args.codec),
        "beam": int(args.beam),
        "K": int(args.K),
        "workers": int(args.workers),
        "groups": [{"name": g, "wA": wA, "wB": wB, "rhoB_mode": m, "rhoB_threshold": thr} for (g, wA, wB, m, thr) in groups],
    }
    report.write_text(
        "\n".join(
            [
                f"# Blind Sprint ablations (tag={args.tag})",
                "",
                "## Protocol",
                "",
                f"- targets: {', '.join(ids)}",
                f"- seeds: {', '.join(str(s) for s in seeds)}",
                f"- codec: `{args.codec}`",
                f"- beam: {int(args.beam)}",
                f"- K: {int(args.K)}",
                f"- workers: {int(args.workers)}",
                "",
                "## Best-of-seeds TM per target",
                "",
                "(TM is evaluation-only; not used for acceptance.)",
                "",
                _df_to_md(best_pivot.reset_index().fillna(np.nan)),
                "",
                "## Acceptance proxy (best-of-seeds TM>0.6 fraction)",
                "",
                _df_to_md(frac_df),
                "",
                "## Artifacts",
                "",
                f"- CSV (gitignored): `{out_csv.relative_to(root).as_posix()}`",
                f"- Plots: `{plots_dir.relative_to(root).as_posix()}/`",
                "",
                "## Meta",
                "",
                "```json",
                json.dumps(meta, indent=2),
                "```",
                "",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    dt = time.perf_counter() - t0
    print(f"Wrote: {report}")
    print(f"Done in {dt:.1f}s")


if __name__ == "__main__":
    main()

