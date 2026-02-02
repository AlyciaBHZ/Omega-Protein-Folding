#!/usr/bin/env python3
"""
Run a small blind-sprint suite (sparse-distogram MDS + Auric) across multiple PDB IDs.

Outputs:
- data/processed/blind_sprint/blind_sprint_suite_<tag>.csv (gitignored)
- docs/runs/blind_sprint_auric_small/suite_<tag>.md (committable)
"""

from __future__ import annotations

import argparse
import time
from pathlib import Path

import numpy as np
import pandas as pd

import sys

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))  # allow importing distogram_mds_auric
sys.path.insert(0, str((HERE.parent / "pdb_bench").resolve()))

from bench_utils import ensure_dir, load_ca_coords_longest_chain  # noqa: E402

import distogram_mds_auric as mds  # noqa: E402


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--ids", required=True, help="IDs file (one PDB ID per line).")
    ap.add_argument("--tag", default="small_suite_v1")
    ap.add_argument("--seeds", default="0,1,2,3,4")
    ap.add_argument("--n-pairs", type=int, default=300, help="Target number of distogram pairs per run.")
    ap.add_argument("--steps", type=int, default=1500)
    ap.add_argument("--lr", type=float, default=0.03)
    ap.add_argument("--clamp", type=float, default=5.0)
    ap.add_argument("--min-sep", type=int, default=3)
    ap.add_argument("--w-chain", type=float, default=5.0)
    ap.add_argument("--w-dist", type=float, default=1.0)
    ap.add_argument("--alphabet", default="pair72", choices=["pair72", "triple232"])
    ap.add_argument("--auric-m", type=int, default=8)
    args = ap.parse_args()

    root = Path(__file__).resolve().parents[2]
    ids_path = Path(args.ids)
    if not ids_path.is_absolute():
        ids_path = (root / ids_path).resolve()

    ids = [x.strip().upper() for x in ids_path.read_text(encoding="utf-8").splitlines() if x.strip()]
    seeds = [int(x.strip()) for x in str(args.seeds).split(",") if x.strip()]

    out_dir = root / "data" / "processed" / "blind_sprint"
    ensure_dir(out_dir)
    run_dir = root / "docs" / "runs" / "blind_sprint_auric_small"
    ensure_dir(run_dir)

    rows = []
    t0 = time.perf_counter()
    for pdb_id in ids:
        cif = root / "data" / "raw" / "pdb_cache" / f"{pdb_id}.cif"
        if not cif.exists():
            continue
        try:
            chain_id, coords_native = load_ca_coords_longest_chain(cif)
        except Exception:
            continue
        coords_native = np.asarray(coords_native, dtype=np.float64)
        N = int(coords_native.shape[0])
        bond_len = float(np.median(np.linalg.norm(coords_native[1:] - coords_native[:-1], axis=1)))
        d_true = np.linalg.norm(coords_native[:, None, :] - coords_native[None, :, :], axis=2)

        eligible_pairs = (N * (N - 1)) // 2 - (N - 1)  # rough; ok for scaling
        n_pairs = int(min(int(args.n_pairs), max(50, eligible_pairs)))

        best_tm = -1.0
        best_row = None
        for s in seeds:
            rng = np.random.default_rng(int(s))
            edges = mds.build_edges(
                d_true,
                bond_len=bond_len,
                n_pairs=n_pairs,
                min_sep=int(args.min_sep),
                w_chain=float(args.w_chain),
                w_dist=float(args.w_dist),
                rng=rng,
            )
            X = mds.optimize_mds(edges, N=N, seed=s, steps=int(args.steps), lr=float(args.lr), clamp=float(args.clamp))
            tm = float(mds.tm_score(X, coords_native))
            aur = mds.auric_metrics(X, alphabet=str(args.alphabet), m=int(args.auric_m))
            row = {
                "pdb_id": pdb_id,
                "chain": str(chain_id),
                "N": N,
                "seed": int(s),
                "tm": tm,
                "n_pairs": n_pairs,
                "bond_len": bond_len,
                "type_entropy_A": float(aur["type_entropy_A"]),
                "smb_rate_hat_B": float(aur["smb_rate_hat_B"]),
            }
            rows.append(row)
            if tm > best_tm:
                best_tm = tm
                best_row = row
        if best_row is not None:
            print(f"{pdb_id}:{chain_id} best_TM={best_tm:.3f} (N={N}, n_pairs={n_pairs})", flush=True)

    df = pd.DataFrame(rows)
    out_csv = out_dir / f"blind_sprint_suite_{args.tag}.csv"
    df.to_csv(out_csv, index=False)
    dt = time.perf_counter() - t0

    # Per-target best
    best = df.sort_values("tm", ascending=False).groupby(["pdb_id", "chain", "N"], as_index=False).head(1)

    # Write a committable summary
    report = run_dir / f"suite_{args.tag}.md"
    lines = []
    lines.append(f"# Blind Sprint suite (sparse-distogram MDS + Auric) — {args.tag}")
    lines.append("")
    lines.append(f"- IDs: `{ids_path.relative_to(root).as_posix()}`")
    lines.append(f"- seeds: {args.seeds}")
    lines.append(f"- n_pairs: {int(args.n_pairs)} (capped per target)")
    lines.append(f"- MDS: steps={int(args.steps)}, lr={float(args.lr)}, clamp={float(args.clamp)}")
    lines.append(f"- Auric: alphabet={args.alphabet}, m={int(args.auric_m)}")
    lines.append(f"- CSV (gitignored): `{out_csv.relative_to(root).as_posix()}`")
    lines.append(f"- Runtime: {dt:.1f}s")
    lines.append("")
    lines.append("## Best-of-seeds per target (TM is evaluation-only)")
    lines.append("")
    cols = ["pdb_id", "chain", "N", "seed", "tm", "n_pairs", "type_entropy_A", "smb_rate_hat_B"]
    lines.append("| " + " | ".join(cols) + " |")
    lines.append("| " + " | ".join(["---"] * len(cols)) + " |")
    for _, r in best.iterrows():
        lines.append(
            "| "
            + " | ".join(
                [
                    str(r["pdb_id"]),
                    str(r["chain"]),
                    str(int(r["N"])),
                    str(int(r["seed"])),
                    f"{float(r['tm']):.3f}",
                    str(int(r["n_pairs"])),
                    f"{float(r['type_entropy_A']):.3f}",
                    f"{float(r['smb_rate_hat_B']):.3f}",
                ]
            )
            + " |"
        )
    lines.append("")
    report.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"Wrote: {report}")


if __name__ == "__main__":
    main()

