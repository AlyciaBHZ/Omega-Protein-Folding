#!/usr/bin/env python3
"""
Attach an Omega engine implementation and run an end-to-end benchmark on PDB targets.

This repository includes an evidence pack and an interface contract:
  - scripts/evidence_v1/specs/omega_engine_spec.md
  - scripts/evidence_v1/specs/omega_minimal_api.py

The full engine implementation may live elsewhere. This runner expects a python module
that provides:
  - Target, RunConfig dataclasses (or compatible)
  - run_omega(target: Target, cfg: RunConfig) -> dict with at least metrics keys

Example:
  python scripts/pdb_bench/engine_bench_runner.py ^
    --engine-module my_omega_engine.api ^
    --ids data/processed/pdb_bench/pdb_ids_n200_seed0.txt ^
    --tag n200_seed0_engine_v1
"""

from __future__ import annotations

import argparse
import importlib
import json
import time
from pathlib import Path
from typing import Any, Dict, List

import numpy as np
import pandas as pd

import sys

sys.path.insert(0, str(Path(__file__).resolve().parent))

from bench_utils import ensure_dir, load_ca_coords_longest_chain  # noqa: E402


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--engine-module", required=True, help="Python module path that exposes Target/RunConfig/run_omega.")
    ap.add_argument("--ids", required=True, help="PDB IDs file (one per line).")
    ap.add_argument("--tag", default="engine_v1")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--steps", type=int, default=200)
    ap.add_argument("--K", type=int, default=32)
    ap.add_argument("--codec", type=str, default="pair72")
    ap.add_argument("--modeA", type=str, default="A_baseline")
    ap.add_argument("--modeB", type=str, default="B_challenger")
    ap.add_argument("--min-len", type=int, default=60)
    ap.add_argument("--max-len", type=int, default=350)
    ap.add_argument("--max-proteins", type=int, default=0)
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

    mod = importlib.import_module(args.engine_module)
    Target = getattr(mod, "Target")
    RunConfig = getattr(mod, "RunConfig")
    run_omega = getattr(mod, "run_omega")

    ids = [x.strip().upper() for x in ids_path.read_text(encoding="utf-8").splitlines() if x.strip()]
    rows: List[Dict[str, Any]] = []
    t0 = time.perf_counter()
    used = 0

    def run_one(pdb_id: str, mode: str, seed_offset: int) -> Dict[str, Any]:
        cif = cache_dir / f"{pdb_id}.cif"
        chain_id, coords = load_ca_coords_longest_chain(cif)
        target = Target(target_id=f"{pdb_id}:{chain_id}", coords_native=np.asarray(coords, dtype=np.float64))
        cfg = RunConfig(seed=int(args.seed + seed_offset), steps=int(args.steps), K=int(args.K), codec=str(args.codec), mode=str(mode))
        out = run_omega(target, cfg)
        metrics = out.get("metrics", out)  # allow direct metrics dict
        flat = {"pdb_id": pdb_id, "chain": chain_id, "N": int(len(coords)), "mode": mode}
        # Copy common keys if present
        for k in ["tm", "rmse", "f1", "ph_max", "theta", "ph_ratio", "rejects_total", "projector_calls", "w0_delta_norm"]:
            if k in metrics:
                flat[k] = metrics[k]
        return flat

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

        try:
            rows.append(run_one(pdb_id, args.modeA, seed_offset=0))
            rows.append(run_one(pdb_id, args.modeB, seed_offset=100000))
        except Exception as e:
            rows.append({"pdb_id": pdb_id, "chain": chain_id, "N": int(N), "mode": "error", "error": f"{type(e).__name__}: {e}"})

        used += 1

    df = pd.DataFrame(rows)
    out_csv = out_dir / f"engine_bench_{args.tag}.csv"
    df.to_csv(out_csv, index=False)

    report_path = rep_dir / f"pdb_engine_bench_{args.tag}.md"
    report_path.write_text(
        "\n".join(
            [
                "# PDB Omega engine benchmark",
                "",
                "This run requires an external Omega engine implementation attached via `--engine-module`.",
                "",
                "## Dataset",
                f"- IDs: `{str(ids_path.relative_to(root)).replace('\\\\','/')}`",
                f"- Used chains: {used}",
                "",
                "## Settings",
                f"- engine_module: `{args.engine_module}`",
                f"- steps: {args.steps}",
                f"- K: {args.K}",
                f"- codec: {args.codec}",
                f"- modes: {args.modeA}, {args.modeB}",
                "",
                "## Outputs",
                f"- CSV: `{str(out_csv.relative_to(root)).replace('\\\\','/')}`",
                "",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    dt = time.perf_counter() - t0
    print(f"Wrote: {out_csv}")
    print(f"Wrote: {report_path}")
    print(f"Runtime: {dt:.1f}s")


if __name__ == "__main__":
    main()

