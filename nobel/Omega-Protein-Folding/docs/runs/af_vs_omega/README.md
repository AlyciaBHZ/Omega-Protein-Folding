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
- `P04637` (p53; mixed confidence): `afdb_P04637_summary.md`, `debugger_afdb_P04637.png`, `shadow_afdb_P04637.png`

## Omega example (white-box baseline)

We export a per-step audit log from an in-repo Omega runner (`blind_sprint_auric.py`)
and plot it as a trajectory panel.

- Target: `1CRN`
- Artifacts:
  - `omega_trajectory_1CRN.csv`
  - `omega_trajectory_1CRN.png`
  - `omega_trajectory_1CRN_v1.csv` (improved run; TM~0.74)
  - `omega_trajectory_1CRN_v1.png`
  - `omega_pred_1CRN_v1.pdb` (CA-only trace for plotting)

## 6D shadow projections (Auric-space)

- AlphaFold DB: `shadow_afdb_P06454.png`
- Omega: `shadow_omega_1CRN.png`

## Debugger views

- AlphaFold DB pLDDT: `debugger_afdb_P06454.png`
- Omega audit debugger: `debugger_omega_1CRN.png`

## Quick numeric sanity table (for reviewer-facing robustness)

- Table: `case_table.csv`

This table pins down the two AFDB examples and the Omega example with a few
protocol-fixed numbers (pLDDT stats where available, plus the same lift-derived
phason proxies and an Auric \(H(\mathrm{type})\) scalar from \(\rho_A\)).

## Paper asset

- Composite figure: `paper/assets/figures/fig_af_vs_omega.png`

