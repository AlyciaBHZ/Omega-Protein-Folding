# Auric summary slice (n1000_seed0_triple232_hybrid_pcaA_parityB_shuffle10_blk8x5_r1)

- Source: `data/processed/pdb_bench/auric_stats_summary_n1000_seed0_triple232_hybrid_pcaA_parityB_shuffle10_blk8x5_r1.csv`
- Alphabet: `triple232`
- n proteins (unique pdb_id): 1000

## Metric: `type_entropy`

| m | readout | mean δ(real,shuffle) | median δ(real,shuffle) | mean pct(real,shuffle) | median pct(real,shuffle) | mean δ(real,blk8) | median δ(real,blk8) | mean pct(real,blk8) | median pct(real,blk8) |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 6 | A | -0.281 | -0.400 | 0.361 | 0.300 | -0.247 | -0.200 | 0.379 | 0.400 |
| 6 | B | -0.360 | -0.600 | 0.320 | 0.200 | -0.136 | -0.200 | 0.432 | 0.400 |
| 8 | A | -0.274 | -0.400 | 0.365 | 0.300 | -0.242 | -0.200 | 0.380 | 0.400 |
| 8 | B | -0.366 | -0.600 | 0.317 | 0.200 | -0.119 | -0.200 | 0.441 | 0.400 |
| 10 | A | -0.278 | -0.400 | 0.363 | 0.300 | -0.243 | -0.200 | 0.381 | 0.400 |
| 10 | B | -0.335 | -0.400 | 0.333 | 0.300 | -0.089 | -0.200 | 0.456 | 0.400 |

### Cross-readout agreement (shuffle z-scores)

| m | corr(z_A, z_B) | n pairs |
| --- | --- | --- |
| 6 | -0.004 | 1000 |
| 8 | -0.017 | 1000 |
| 10 | -0.028 | 1000 |

## Metric: `smb_rate_hat`

| m | readout | mean δ(real,shuffle) | median δ(real,shuffle) | mean pct(real,shuffle) | median pct(real,shuffle) | mean δ(real,blk8) | median δ(real,blk8) | mean pct(real,blk8) | median pct(real,blk8) |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 6 | A | -0.298 | -0.400 | 0.352 | 0.300 | -0.278 | -0.200 | 0.362 | 0.400 |
| 6 | B | -0.494 | -0.800 | 0.253 | 0.100 | -0.326 | -0.600 | 0.337 | 0.200 |
| 8 | A | -0.298 | -0.400 | 0.352 | 0.300 | -0.278 | -0.200 | 0.362 | 0.400 |
| 8 | B | -0.494 | -0.800 | 0.253 | 0.100 | -0.326 | -0.600 | 0.337 | 0.200 |
| 10 | A | -0.298 | -0.400 | 0.352 | 0.300 | -0.278 | -0.200 | 0.362 | 0.400 |
| 10 | B | -0.494 | -0.800 | 0.253 | 0.100 | -0.326 | -0.600 | 0.337 | 0.200 |

### Cross-readout agreement (shuffle z-scores)

| m | corr(z_A, z_B) | n pairs |
| --- | --- | --- |
| 6 | -0.083 | 1000 |
| 8 | -0.083 | 1000 |
| 10 | -0.083 | 1000 |

