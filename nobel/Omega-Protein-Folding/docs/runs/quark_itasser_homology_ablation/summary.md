# QUARK vs I-TASSER homology/fragment ablation benchmark

This run folder is the audit entrypoint for the final reviewer-facing benchmark:

- **Method control**: QUARK vs I‑TASSER on the same sequence (do they converge to the same topology?)
- **Homology dependence controls**:
  - QUARK: default vs *exclude fragments from proteins sharing >30% identity*
  - I‑TASSER: default vs *exclude homologous templates*

## Artifacts

- `targets_6.fasta`: sequences submitted to servers (derived from native PDB cache)
- `targets.md`: target metadata + job ID tracker (QAxxxx, Sxxxxx)
- `tmalign_summary.csv`: TM-align/US-align comparisons (filled after results are downloaded)

## How to reproduce (local prep)

Generate FASTA and target metadata from local Decoys 'R' Us cache:

```bash
python scripts/benchmarks/extract_fasta_from_pdb_cache.py
```

Server submissions and result downloads are manual; record job IDs in `targets.md`.

