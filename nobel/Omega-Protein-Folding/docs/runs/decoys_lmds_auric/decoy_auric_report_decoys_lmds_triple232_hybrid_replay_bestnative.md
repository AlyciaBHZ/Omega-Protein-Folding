# Auric decoy benchmark (decoys_lmds_triple232_hybrid_replay_bestnative)

- Manifest: `docs/runs/decoys_lmds_auric/decoy_manifest.csv`
- Alphabet: `triple232`
- m list: 6, 8, 10
- Axis mode: `best_native_axis`
- FailedGeometry (phason): ph_rms <= 100.0
- Targets used: 9
- Samples CSV (gitignored): `data/processed/decoy_bench/auric_decoys_samples_decoys_lmds_triple232_hybrid_replay_bestnative.csv`
- Summary CSV (gitignored): `data/processed/decoy_bench/auric_decoys_summary_decoys_lmds_triple232_hybrid_replay_bestnative.csv`
- Plots: `docs/runs/decoys_lmds_auric/plots/`
- Geometry QC: max CA-CA <= 10.0 Å and outlier_frac <= 0.20

- Targets with native FailedGeometry under this QC: 0/9

## Native vs Decoy separation (summary)

Using m=8. δ<0 means native has lower metric (more ordered) than decoys.

| target_id | n_decoys_used | pct_native_type_entropy_rhoA_m8 | delta_native_vs_decoys_type_entropy_rhoA_m8 | pct_native_smb_rate_hat_rhoA_m8 | delta_native_vs_decoys_smb_rate_hat_rhoA_m8 | pct_native_type_entropy_rhoB_m8 | pct_native_type_entropy_rhoBvel_m8 | pct_combined_min(A,B)_type_entropy_m8 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 1ctf | 496 (failQC=2) | 0.385 | -0.232 | 0.488 | -0.026 | 0.665 | 0.740 | 0.385 |
| 1dtk | 216 (failQC=0) | 0.051 | -0.917 | 0.282 | -0.449 | 0.954 | 0.472 | 0.051 |
| 1fc2 | 501 (failQC=0) | 0.042 | -0.916 | 0.000 | -1.000 | 0.862 | 0.976 | 0.042 |
| 1igd | 501 (failQC=0) | 0.250 | -0.503 | 0.218 | -0.567 | 0.415 | 0.882 | 0.250 |
| 1shf-A | 436 (failQC=1) | 0.236 | -0.528 | 0.181 | -0.638 | 0.995 | 0.489 | 0.236 |
| 2cro | 501 (failQC=0) | 0.018 | -0.982 | 0.048 | -0.952 | 0.555 | 0.024 | 0.018 |
| 2ovo | 348 (failQC=0) | 0.009 | -0.986 | 0.009 | -0.986 | 0.305 | 0.931 | 0.009 |
| 4pti | 344 (failQC=0) | 0.038 | -0.962 | 0.067 | -0.933 | 0.334 | 0.561 | 0.038 |
| unk | 499 (failQC=0) | 0.144 | -0.711 | 0.289 | -0.423 | 0.972 | 0.760 | 0.144 |

- Fraction targets with native in lowest 10% of decoys (rhoA type_entropy, m=8): 0.556
- Fraction targets with native in lowest 10% of decoys (rhoB=parity type_entropy, m=8): 0.000
- Fraction targets with native in lowest 10% of decoys (rhoB=vel type_entropy, m=8): 0.111
- Fraction targets with native in lowest 10% by min(A,Bparity) (type_entropy, m=8): 0.556

