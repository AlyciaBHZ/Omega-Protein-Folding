# AF vs Omega (visual transparency benchmark) — run artifacts

This folder contains **auditable, paper-facing artifacts** used to build the
AlphaFold-vs-Omega transparency figure.

## AlphaFold DB examples (AF baseline)

We download an AlphaFold DB model PDB for a UniProt accession, then extract
**pLDDT** from the PDB **B-factor** field (per-residue; CA atoms).

### Chosen low-confidence example (used in the paper figure)

- UniProt: `P06454` (prothymosin alpha; intrinsically-disordered)
- Artifacts:
  - `afdb_P06454_model.pdb`
  - `afdb_P06454_plddt.csv`
  - `afdb_P06454_summary.md`

### Other probed examples (kept for traceability; not referenced by the paper)

- `P37840`: `afdb_P37840_summary.md`
- `P62328`: `afdb_P62328_summary.md`

## Omega example (white-box baseline)

We export a per-step audit log from an in-repo Omega runner (`blind_sprint_auric.py`)
and plot it as a trajectory panel.

- Target: `1CRN`
- Artifacts:
  - `omega_trajectory_1CRN.csv`
  - `omega_trajectory_1CRN.png`

