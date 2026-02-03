# QUARK / I-TASSER homology-ablation targets (N=6)

This folder supports the final benchmark described in the plan: QUARK vs I-TASSER,
and each with a homology/fragment exclusion ablation to assess template dependence.

## Target list

| target_id | dataset | chain | length | native_pdb | QUARK_A | QUARK_B | ITASSER_A | ITASSER_B | notes |
| --- | --- | --- | ---: | --- | --- | --- | --- | --- | --- |
| 1r69 | 4state_reduced | A | 63 | `data/raw/decoy_cache/decoys_r_us/4state_reduced/dd/multiple/4state_reduced/doc/pdb_orig/1r69.pdb` | QA16125 |  | S820065 |  | QUARK-A + I-TASSER-A results ingested |
| 2cro | 4state_reduced | A | 65 | `data/raw/decoy_cache/decoys_r_us/4state_reduced/dd/multiple/4state_reduced/doc/pdb_orig/2cro.pdb` |  |  |  |  |  |
| 4pti | 4state_reduced | A | 58 | `data/raw/decoy_cache/decoys_r_us/4state_reduced/dd/multiple/4state_reduced/doc/pdb_orig/4pti.pdb` |  |  |  |  |  |
| 1ctf | 4state_reduced | A | 68 | `data/raw/decoy_cache/decoys_r_us/4state_reduced/dd/multiple/4state_reduced/doc/pdb_orig/1ctf.pdb` |  |  |  |  |  |
| 1dtk | lmds | A | 57 | `data/raw/decoy_cache/decoys_r_us/lmds/dd/multiple/lmds/doc/pdb_orig/1dtk.pdb` |  |  |  |  |  |
| 1shf-A | lmds | A | 59 | `data/raw/decoy_cache/decoys_r_us/lmds/dd/multiple/lmds/1shf-A/1shf-A.pdb` |  |  |  |  |  |

## Inputs

- FASTA (copy/paste into servers): `docs/runs/quark_itasser_homology_ablation/targets_6.fasta`

## Submission protocol reminder

- QUARK-A: default fragments
- QUARK-B: exclude fragments from proteins with >30% sequence identity
- I-TASSER-A: default
- I-TASSER-B: exclude homologous templates (benchmark option)

## Download checklist (per job)

### QUARK (QAxxxx)
- Download `result.tar.bz2`
- Download `model1.pdb` ... `model5.pdb`
- Copy/paste into notes:
  - Predicted Secondary Structure line (H/S/C string)
  - Predicted Solvent Accessibility line (0-9 digits)

### I-TASSER (Sxxxxx)
- Download `Sxxxxx_results.tar.bz2` (tarball)
- Download `model1.pdb` ... `model5.pdb`
- Record from result page:
  - C-score for models 1–5
  - Estimated TM-score + RMSD (usually for model1)
  - Top 10 threading templates table (copy key rows or screenshot)

## Local folder layout (not committed)

We do **not** commit large server bundles. Suggested local cache layout:

- `data/raw/quark_cache/<target_id>/QAxxxx/`
  - `result.tar.bz2`
  - `models/model1.pdb` ... `models/model5.pdb`
- `data/raw/itasser_cache/<target_id>/Sxxxxx/`
  - `Sxxxxx_results.tar.bz2`
  - `models/model1.pdb` ... `models/model5.pdb`


