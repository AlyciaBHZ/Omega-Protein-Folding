# PyRosetta FastRelax validation (coordinate-constrained)

This report treats Rosetta/PyRosetta as a **physics/geometry validator** rather than a full predictor.

## Protocol (minimal)

- Input: CA-only Omega predicted traces under `omega_models/` (seed0) + the corresponding FASTA.
- Build a full-atom pose from FASTA.
- Add CA coordinate constraints (harmonic, SD=\(\sigma\)) to match the CA trace.
- Run `FastRelax` with `nstruct` repeats; pick best `total_score`.
- Record: best/median score and the best CA RMSD (Kabsch) to the target CA trace.

## Artifacts

- Omega inputs:
  - `omega_models_manifest.csv`
  - `omega_models/omega_pred_<target>_seed0.pdb`
  - `omega_models/omega_audit_<target>_seed0.csv`
- Relax summary CSV (committable):
  - `rosetta_relax_summary.csv`
- Heavy outputs (not committed):
  - `data/raw/rosetta_cache/quark_itasser_homology_ablation/omega/<target_id>/relax_<target_id>_best.pdb`

## How to run

Dry-run (no PyRosetta required; validates paths/lengths):

```bash
python scripts/benchmarks/pyrosetta_relax_validate.py --dry-run
```

Actual run (requires PyRosetta installed):

```bash
python scripts/benchmarks/pyrosetta_relax_validate.py --nstruct 20 --coord-sd 1.0 --coord-weight 1.0
```

## Results (fill after run)

Summaries live in `rosetta_relax_summary.csv`. Interpretations to emphasize:

- **Stable basin**: small CA drift + consistently low energies across repeats.\n- **Unstable**: large drift or pathological energies, suggesting the CA trace is not near a Rosetta basin.

