# Full1k Hybrid Auric Run (PCA ρA + parity ρB)

Tag: `n1000_seed0_triple232_hybrid_pcaA_parityB_shuffle10_blk8x5_r1`

## What was run

- IDs: `data/processed/pdb_bench/pdb_ids_n1000_seed0.txt` (1000 chains, N in [60,350])
- Alphabet: `triple232`
- Nulls:
  - random reps = 3
  - shuffle reps = 10
  - blockshuffle reps = 5 with k=8
- Auric protocol (hybrid):
  - ρA: `--auric-u-mode pca` + `--auric-rhoA-threshold median`
  - ρB: `--auric-rhoB-mode parity`
  - m: `6,8,10`

Command:

```bash
python scripts/pdb_bench/phason_stats.py \
  --ids data/processed/pdb_bench/pdb_ids_n1000_seed0.txt \
  --tag n1000_seed0_triple232_hybrid_pcaA_parityB_shuffle10_blk8x5_r1 \
  --alphabet triple232 --seed 0 --min-len 60 --max-len 350 \
  --random-reps 3 --perturb-reps 0 \
  --shuffle-reps 10 \
  --blockshuffle-reps 5 --blockshuffle-k 8 \
  --checkpoint-every 10 --resume \
  --auric --auric-m 6,8,10 --auric-readouts all \
  --auric-u-mode pca --auric-rhoA-threshold median \
  --auric-rhoB-mode parity
```

## Outputs

- **Download manifest (audit that mmCIF cache exists for all IDs)**:
  - `pdb_download_summary_n1000_seed0.md`
  - `pdb_download_manifest_n1000_seed0.csv`
- **Paper-facing summary slice (recommended entry point)**:
  - `auric_summary_slice_n1000_seed0_triple232_hybrid_pcaA_parityB_shuffle10_blk8x5_r1.md`
- **Figures**:
  - `figures/` (boxplots, ECDFs, multi-scale curves)

Note: the large CSVs and auto-generated report markdowns are written under
`data/processed/pdb_bench/` and `artifacts/reports/` respectively (those folders are git-ignored).

