# Auric decoy benchmark (decoys_4state_reduced_triple232_hybrid)

- Manifest: `docs/runs/decoys_4state_auric/decoy_manifest.csv`
- Alphabet: `triple232`
- m list: 6, 8, 10
- Targets used: 7
- Samples CSV (gitignored): `data/processed/decoy_bench/auric_decoys_samples_decoys_4state_reduced_triple232_hybrid.csv`
- Summary CSV (gitignored): `data/processed/decoy_bench/auric_decoys_summary_decoys_4state_reduced_triple232_hybrid.csv`
- Plots: `docs/runs/decoys_4state_auric/plots/`

## Native vs Decoy separation (summary)

Using m=8. δ<0 means native has lower metric (more ordered) than decoys.

| target_id | n_decoys_used | pct_native_type_entropy_rhoA_m8 | delta_native_vs_decoys_type_entropy_rhoA_m8 | pct_native_smb_rate_hat_rhoA_m8 | delta_native_vs_decoys_smb_rate_hat_rhoA_m8 |
| --- | --- | --- | --- | --- | --- |
| 1ctf | 631 | 0.544 | 0.086 | 0.456 | -0.089 |
| 1r69 | 676 | 0.030 | -0.970 | 0.037 | -0.963 |
| 1sn3 | 660 | 0.285 | -0.432 | 0.302 | -0.398 |
| 2cro | 674 | 0.016 | -0.984 | 0.036 | -0.964 |
| 3icb | 654 | 0.170 | -0.661 | 0.318 | -0.364 |
| 4pti | 687 | 0.159 | -0.689 | 0.226 | -0.557 |
| 4rxn | 677 | 0.171 | -0.663 | 0.258 | -0.487 |

- Fraction targets with native in lowest 10% of decoys (rhoA type_entropy, m=8): 0.286

