#!/usr/bin/env python3
from __future__ import annotations

"""
Generate Omega-style CA-only predicted traces (and per-step audit CSVs) for the
QUARK/I-TASSER homology-ablation targets, as inputs for Rosetta/PyRosetta validation.

Why this exists:
- Our in-repo Omega example `docs/runs/af_vs_omega/omega_pred_1CRN_v1.pdb` is a CA-only trace.
- For Rosetta relax validation, we need Omega outputs for the same 6 targets
  (`1R69`, `2CRO`, `4PTI`, `1CTF`, `1DTK`, `1SHF-A`).
- We generate them using the existing in-repo runner `scripts/blind_sprint/blind_sprint_auric.py`
  (native-derived contacts/distances; no oracle directions).

Outputs (committable):
- docs/runs/quark_itasser_homology_ablation/omega_models/omega_pred_<target_id>_seed0.pdb
- docs/runs/quark_itasser_homology_ablation/omega_models/omega_audit_<target_id>_seed0.csv
- docs/runs/quark_itasser_homology_ablation/omega_models_manifest.csv
"""

import csv
import runpy
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Tuple

import numpy as np


TARGETS_6: Tuple[str, ...] = ("1R69", "2CRO", "4PTI", "1CTF", "1DTK", "1SHF-A")


def repo_root() -> Path:
    return Path(__file__).resolve().parents[2]


def ensure_dir(p: Path) -> None:
    p.mkdir(parents=True, exist_ok=True)


def _iter_manifest_rows(path: Path) -> Iterable[dict]:
    with path.open("r", encoding="utf-8", newline="") as f:
        r = csv.DictReader(f)
        for row in r:
            yield row


def _native_paths_from_manifests(root: Path) -> Dict[str, Tuple[str, Path]]:
    """
    Return map target_id -> (dataset, native_path) using the existing decoy manifests.
    """
    out: Dict[str, Tuple[str, Path]] = {}
    man_paths = [
        root / "docs" / "runs" / "decoys_4state_auric" / "decoy_manifest.csv",
        root / "docs" / "runs" / "decoys_lmds_auric" / "decoy_manifest.csv",
    ]
    for mp in man_paths:
        if not mp.exists():
            continue
        for row in _iter_manifest_rows(mp):
            tid = str(row.get("target_id", "")).strip()
            ds = str(row.get("dataset", "(unknown)")).strip()
            npth = str(row.get("native_path", "")).strip()
            if tid and npth and tid not in out:
                out[tid] = (ds, (root / npth).resolve())
    return out


def _resolve_tid(native_map: Dict[str, Tuple[str, Path]], t: str) -> str:
    for k in native_map.keys():
        if k.lower() == t.lower():
            return k
    return t


def write_ca_only_pdb(coords: np.ndarray, out_pdb: Path, *, chain_id: str = "A") -> None:
    x = np.asarray(coords, dtype=np.float64)
    lines: List[str] = []
    for i, (xx, yy, zz) in enumerate(x, start=1):
        lines.append(
            f"ATOM  {i:5d}  CA  ALA {chain_id}{i:4d}    "
            f"{xx:8.3f}{yy:8.3f}{zz:8.3f}"
            f"{1.00:6.2f}{0.00:6.2f}           C"
        )
    lines.append("END")
    out_pdb.write_text("\n".join(lines) + "\n", encoding="utf-8")


@dataclass(frozen=True)
class OmegaModelRow:
    target_id: str
    dataset: str
    native_pdb: Path
    omega_model_path: Path
    omega_audit_path: Path
    tm_eval: float


def main() -> None:
    root = repo_root()

    # Robust CA coordinate loader (PDB/mmCIF)
    import sys

    sys.path.insert(0, str((root / "scripts" / "pdb_bench").resolve()))
    from coords_io import load_ca_coords  # type: ignore  # noqa: E402

    # Load blind sprint runner function via runpy.
    bs_path = (root / "scripts" / "blind_sprint" / "blind_sprint_auric.py").resolve()
    mod = runpy.run_path(str(bs_path))
    run_one = mod["run_one"]

    native_map = _native_paths_from_manifests(root)
    missing = [t for t in TARGETS_6 if t.lower() not in {k.lower() for k in native_map.keys()}]
    if missing:
        raise SystemExit(f"Missing native paths for targets in manifests: {missing}")

    run_dir = root / "docs" / "runs" / "quark_itasser_homology_ablation"
    omega_dir = run_dir / "omega_models"
    ensure_dir(omega_dir)

    rows: List[OmegaModelRow] = []
    for t in TARGETS_6:
        tid = _resolve_tid(native_map, t)
        ds, native_pdb = native_map[tid]
        if not native_pdb.exists():
            raise SystemExit(f"Missing native PDB: {native_pdb}")

        chain_id, coords_native = load_ca_coords(native_pdb)
        coords_native = np.asarray(coords_native, dtype=np.float64)

        out_pdb = omega_dir / f"omega_pred_{tid}_seed0.pdb"
        audit_csv = omega_dir / f"omega_audit_{tid}_seed0.csv"

        # Deterministic Omega-style blind sprint.
        best_state, meta = run_one(
            coords_native,
            seed=0,
            codec="pair72",
            beam=64,
            K=32,
            cutoff=8.0,
            min_sep=3,
            auric_m=8,
            auric_every=10,
            w_contact=1.0,
            w_rmse=1.0,
            wA=0.2,
            wB=0.2,
            audit_out=audit_csv,
        )

        write_ca_only_pdb(best_state.coords, out_pdb, chain_id="A")
        if not out_pdb.exists():
            raise RuntimeError(f"Failed to write omega_pred PDB: {out_pdb}")
        if not audit_csv.exists():
            raise RuntimeError(f"Failed to write omega audit CSV: {audit_csv}")

        tm = float(meta.get("tm", float("nan")))
        rows.append(
            OmegaModelRow(
                target_id=tid,
                dataset=ds,
                native_pdb=native_pdb,
                omega_model_path=out_pdb,
                omega_audit_path=audit_csv,
                tm_eval=tm,
            )
        )

        print(f"[ok] {tid} ({ds}) omega_pred written; TM(eval)={tm:.3f}", flush=True)

    # Write committable manifest used by the PyRosetta relax script.
    man_path = run_dir / "omega_models_manifest.csv"
    with man_path.open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(
            f,
            fieldnames=[
                "target_id",
                "dataset",
                "native_pdb",
                "omega_model_path",
                "omega_audit_path",
                "tm_eval_seed0",
            ],
        )
        w.writeheader()
        for r in rows:
            w.writerow(
                {
                    "target_id": r.target_id,
                    "dataset": r.dataset,
                    "native_pdb": r.native_pdb.relative_to(root).as_posix(),
                    "omega_model_path": r.omega_model_path.relative_to(root).as_posix(),
                    "omega_audit_path": r.omega_audit_path.relative_to(root).as_posix(),
                    "tm_eval_seed0": f"{r.tm_eval:.4f}" if np.isfinite(r.tm_eval) else "",
                }
            )

    print(f"[ok] wrote manifest: {man_path}", flush=True)


if __name__ == "__main__":
    main()

