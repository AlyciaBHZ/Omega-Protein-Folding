#!/usr/bin/env python3
from __future__ import annotations

"""
PyRosetta coordinate-constrained FastRelax validation for Omega/QUARK/I-TASSER models.

This implements the plan's "Rosetta as physics/geometry verifier" role with the
smallest actionable protocol:
- Build a full-atom pose from the target FASTA.
- Add CA coordinate constraints to match an input CA trace (Omega CA-only PDB).
- Run FastRelax for nstruct replicates, select best total_score.

Outputs (committable):
- docs/runs/quark_itasser_homology_ablation/rosetta_relax_summary.csv
- docs/runs/quark_itasser_homology_ablation/rosetta_relax_report.md (filled manually/iteratively)

Heavy outputs (gitignored):
- data/raw/rosetta_cache/quark_itasser_homology_ablation/<method>/<target_id>/.../*.pdb
"""

import argparse
import csv
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Tuple

import numpy as np


def repo_root() -> Path:
    return Path(__file__).resolve().parents[2]


def ensure_dir(p: Path) -> None:
    p.mkdir(parents=True, exist_ok=True)


def read_fasta(path: Path) -> Dict[str, str]:
    """
    Return mapping {target_id: sequence}.

    Accept headers like:
      >1r69|chain=A|len=63|dataset=4state_reduced
    """
    txt = path.read_text(encoding="utf-8", errors="ignore").splitlines()
    out: Dict[str, str] = {}
    cur_id = None
    buf: List[str] = []
    for ln in txt:
        ln = ln.strip()
        if not ln:
            continue
        if ln.startswith(">"):
            if cur_id is not None:
                out[cur_id] = "".join(buf).replace(" ", "").upper()
            hdr = ln[1:]
            cur_id = hdr.split("|", 1)[0].strip()
            buf = []
        else:
            buf.append(ln.strip())
    if cur_id is not None:
        out[cur_id] = "".join(buf).replace(" ", "").upper()
    return out


def read_csv_rows(path: Path) -> List[dict]:
    with path.open("r", encoding="utf-8", newline="") as f:
        r = csv.DictReader(f)
        return [row for row in r]


def ca_rmsd_kabsch(P: np.ndarray, Q: np.ndarray) -> float:
    P = np.asarray(P, dtype=np.float64)
    Q = np.asarray(Q, dtype=np.float64)
    if P.shape != Q.shape or P.ndim != 2 or P.shape[1] != 3:
        return float("nan")
    # reuse our in-repo Kabsch
    root = repo_root()
    import sys

    sys.path.insert(0, str((root / "scripts" / "pdb_bench").resolve()))
    from bench_utils import kabsch  # type: ignore  # noqa: E402

    P_aligned, _, _ = kabsch(P, Q)
    d2 = np.sum((P_aligned - Q) ** 2, axis=1)
    return float(np.sqrt(np.mean(d2)))


def _import_pyrosetta() -> object:
    try:
        import pyrosetta  # type: ignore

        return pyrosetta
    except Exception as e:  # pragma: no cover
        raise RuntimeError(
            "PyRosetta is not installed in this Python environment. "
            "Install PyRosetta (licensed) and rerun this script."
        ) from e


def _add_ca_coord_constraints(
    pose,
    ca_target: np.ndarray,
    *,
    sd: float,
) -> None:
    """
    Add coordinate constraints for CA atoms.

    Default protocol: constrain each CA_i relative to CA_1 with a harmonic well:
      || (CA_i - CA_1) - (target_i - target_1) || ~ Harmonic(0, sd)
    This avoids needing a virtual root residue and is sufficient for stability checks.
    """
    import pyrosetta  # type: ignore

    AtomID = pyrosetta.rosetta.core.id.AtomID
    HarmonicFunc = pyrosetta.rosetta.core.scoring.func.HarmonicFunc
    CoordinateConstraint = pyrosetta.rosetta.core.scoring.constraints.CoordinateConstraint
    Vec = pyrosetta.rosetta.numeric.xyzVector_double_t

    N = int(ca_target.shape[0])
    if N <= 0:
        return
    # reference: CA of residue 1
    ref = AtomID(pose.residue(1).atom_index("CA"), 1)
    origin = np.asarray(ca_target[0], dtype=np.float64).reshape(3)
    func = HarmonicFunc(0.0, float(sd))

    for i in range(1, N + 1):
        aid = AtomID(pose.residue(i).atom_index("CA"), i)
        tgt = np.asarray(ca_target[i - 1], dtype=np.float64).reshape(3) - origin
        v = Vec(float(tgt[0]), float(tgt[1]), float(tgt[2]))
        try:
            cst = CoordinateConstraint(aid, ref, v, func)
        except TypeError:
            # Fallback for alternate bindings (rare): absolute coordinate constraint
            v_abs = Vec(float(ca_target[i - 1, 0]), float(ca_target[i - 1, 1]), float(ca_target[i - 1, 2]))
            cst = CoordinateConstraint(aid, v_abs, func)
        pose.add_constraint(cst)


def _pose_ca_coords(pose) -> np.ndarray:
    import pyrosetta  # type: ignore

    AtomID = pyrosetta.rosetta.core.id.AtomID
    N = int(pose.total_residue())
    out = np.zeros((N, 3), dtype=np.float64)
    for i in range(1, N + 1):
        aid = AtomID(pose.residue(i).atom_index("CA"), i)
        xyz = pose.xyz(aid)
        out[i - 1, 0] = float(xyz.x)
        out[i - 1, 1] = float(xyz.y)
        out[i - 1, 2] = float(xyz.z)
    return out


@dataclass(frozen=True)
class RelaxResult:
    best_score: float
    median_score: float
    best_rmsd_ca: float
    best_pdb_path: Path


def relax_pose_replicates(
    *,
    seq: str,
    ca_target: np.ndarray,
    nstruct: int,
    coord_sd: float,
    coord_weight: float,
    out_dir: Path,
    tag: str,
) -> RelaxResult:
    pyrosetta = _import_pyrosetta()
    # init once per process (no-op if already)
    pyrosetta.init("-mute all")

    from pyrosetta.rosetta.core.scoring import get_score_function  # type: ignore
    from pyrosetta.rosetta.protocols.relax import FastRelax  # type: ignore
    from pyrosetta.rosetta.core.scoring import ScoreType  # type: ignore

    ensure_dir(out_dir)

    scorefxn = get_score_function()
    scorefxn.set_weight(ScoreType.coordinate_constraint, float(coord_weight))

    # base pose (PyRosetta API differs across releases; prefer the top-level helper)
    try:
        pose0 = pyrosetta.pose_from_sequence(str(seq), "fa_standard", auto_termini=True)
    except Exception:
        pose0 = pyrosetta.Pose()
        try:
            from pyrosetta.rosetta.core.pose import make_pose_from_sequence  # type: ignore

            make_pose_from_sequence(pose0, str(seq), "fa_standard", True)
        except Exception as e:
            raise RuntimeError("Failed to construct a pose from sequence (PyRosetta API mismatch).") from e
    _add_ca_coord_constraints(pose0, np.asarray(ca_target, dtype=np.float64), sd=float(coord_sd))

    scores: List[float] = []
    best_score = float("inf")
    best_pose = None
    best_rmsd = float("inf")

    relax = FastRelax()
    relax.set_scorefxn(scorefxn)

    for k in range(int(nstruct)):
        pose = pose0.clone()
        relax.apply(pose)
        s = float(scorefxn(pose))
        scores.append(s)
        ca = _pose_ca_coords(pose)
        rmsd = ca_rmsd_kabsch(ca, ca_target)
        if s < best_score - 1e-8:
            best_score = s
            best_pose = pose
            best_rmsd = float(rmsd)

    med = float(np.median(np.asarray(scores, dtype=np.float64))) if scores else float("nan")
    if best_pose is None:
        raise RuntimeError("FastRelax produced no poses")

    out_pdb = out_dir / f"{tag}_best.pdb"
    best_pose.dump_pdb(str(out_pdb))
    return RelaxResult(best_score=best_score, median_score=med, best_rmsd_ca=best_rmsd, best_pdb_path=out_pdb)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--targets-fasta",
        default="docs/runs/quark_itasser_homology_ablation/targets_6.fasta",
        help="FASTA with target sequences.",
    )
    ap.add_argument(
        "--omega-manifest",
        default="docs/runs/quark_itasser_homology_ablation/omega_models_manifest.csv",
        help="Manifest mapping target_id to Omega CA-only model paths.",
    )
    ap.add_argument(
        "--out-csv",
        default="docs/runs/quark_itasser_homology_ablation/rosetta_relax_summary.csv",
        help="Committable summary CSV path.",
    )
    ap.add_argument(
        "--cache-dir",
        default="data/raw/rosetta_cache/quark_itasser_homology_ablation",
        help="Untracked heavy output cache directory.",
    )
    ap.add_argument(
        "--targets",
        default="",
        help="Comma-separated target_ids to run (case-insensitive). If empty, run all in the manifest.",
    )
    ap.add_argument("--nstruct", type=int, default=20)
    ap.add_argument("--coord-sd", type=float, default=1.0)
    ap.add_argument("--coord-weight", type=float, default=1.0)
    ap.add_argument("--dry-run", action="store_true", help="Validate inputs and write headers only; no PyRosetta run.")
    args = ap.parse_args()

    root = repo_root()
    fasta_path = Path(str(args.targets_fasta))
    if not fasta_path.is_absolute():
        fasta_path = (root / fasta_path).resolve()
    man_path = Path(str(args.omega_manifest))
    if not man_path.is_absolute():
        man_path = (root / man_path).resolve()

    out_csv = Path(str(args.out_csv))
    if not out_csv.is_absolute():
        out_csv = (root / out_csv).resolve()
    ensure_dir(out_csv.parent)

    cache_dir = Path(str(args.cache_dir))
    if not cache_dir.is_absolute():
        cache_dir = (root / cache_dir).resolve()
    ensure_dir(cache_dir)

    seqs = read_fasta(fasta_path)
    rows = read_csv_rows(man_path)
    wanted = {t.strip().lower() for t in str(args.targets).split(",") if t.strip()}

    # CA coords loader
    import sys

    sys.path.insert(0, str((root / "scripts" / "pdb_bench").resolve()))
    from coords_io import load_ca_coords  # type: ignore  # noqa: E402

    out_rows: List[Dict[str, object]] = []
    for row in rows:
        tid = str(row["target_id"]).strip()
        if wanted and tid.lower() not in wanted:
            continue
        omega_pdb = Path(str(row["omega_model_path"]).strip())
        if not omega_pdb.is_absolute():
            omega_pdb = (root / omega_pdb).resolve()
        if tid not in seqs:
            # try case-insensitive match
            match = None
            for k in seqs.keys():
                if k.lower() == tid.lower():
                    match = k
                    break
            if match is None:
                raise RuntimeError(f"Missing FASTA entry for target_id={tid}")
            seq = seqs[match]
        else:
            seq = seqs[tid]

        _, ca = load_ca_coords(omega_pdb)
        ca = np.asarray(ca, dtype=np.float64)
        if len(seq) != int(ca.shape[0]):
            raise RuntimeError(f"Length mismatch target={tid}: fasta={len(seq)} ca_trace={ca.shape[0]} path={omega_pdb}")

        if bool(args.dry_run):
            out_rows.append(
                {
                    "method": "omega",
                    "target_id": tid,
                    "input_model": omega_pdb.relative_to(root).as_posix(),
                    "nstruct": int(args.nstruct),
                    "coord_sd": float(args.coord_sd),
                    "coord_weight": float(args.coord_weight),
                    "best_total_score": "",
                    "median_total_score": "",
                    "best_rmsd_ca": "",
                    "best_pdb_path": "",
                    "notes": "dry_run",
                }
            )
            continue

        out_dir = cache_dir / "omega" / tid
        res = relax_pose_replicates(
            seq=seq,
            ca_target=ca,
            nstruct=int(args.nstruct),
            coord_sd=float(args.coord_sd),
            coord_weight=float(args.coord_weight),
            out_dir=out_dir,
            tag=f"relax_{tid}",
        )

        out_rows.append(
            {
                "method": "omega",
                "target_id": tid,
                "input_model": omega_pdb.relative_to(root).as_posix(),
                "nstruct": int(args.nstruct),
                "coord_sd": float(args.coord_sd),
                "coord_weight": float(args.coord_weight),
                "best_total_score": f"{res.best_score:.4f}" if math.isfinite(res.best_score) else "",
                "median_total_score": f"{res.median_score:.4f}" if math.isfinite(res.median_score) else "",
                "best_rmsd_ca": f"{res.best_rmsd_ca:.4f}" if math.isfinite(res.best_rmsd_ca) else "",
                "best_pdb_path": res.best_pdb_path.relative_to(root).as_posix(),
                "notes": "",
            }
        )

        print(f"[ok] omega {tid}: best={res.best_score:.3f} rmsd_ca={res.best_rmsd_ca:.3f}", flush=True)

    # Write summary CSV
    with out_csv.open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(
            f,
            fieldnames=[
                "method",
                "target_id",
                "input_model",
                "nstruct",
                "coord_sd",
                "coord_weight",
                "best_total_score",
                "median_total_score",
                "best_rmsd_ca",
                "best_pdb_path",
                "notes",
            ],
        )
        w.writeheader()
        for r in out_rows:
            w.writerow(r)

    print(f"[ok] wrote: {out_csv}", flush=True)


if __name__ == "__main__":
    main()

