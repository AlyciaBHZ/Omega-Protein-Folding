# Auric decoy benchmark (decoys_4state_reduced_triple232_hybrid_replay_sharedpca)

- Manifest: `docs/runs/decoys_4state_auric/decoy_manifest.csv`
- Alphabet: `triple232`
- m list: 6, 8, 10
- Axis mode: `shared_pca`
- FailedGeometry (phason): ph_rms <= 100.0
- Targets used: 7
- Samples CSV (gitignored): `data/processed/decoy_bench/auric_decoys_samples_decoys_4state_reduced_triple232_hybrid_replay_sharedpca.csv`
- Summary CSV (gitignored): `data/processed/decoy_bench/auric_decoys_summary_decoys_4state_reduced_triple232_hybrid_replay_sharedpca.csv`
- Plots: `docs/runs/decoys_4state_auric/plots/`
- Geometry QC: max CA-CA <= 10.0 Å and outlier_frac <= 0.20

- Targets with native FailedGeometry under this QC: 1/7

## Native vs Decoy separation (summary)

Using m=8. δ<0 means native has lower metric (more ordered) than decoys.

| target_id | n_decoys_used | pct_native_type_entropy_rhoA_m8 | delta_native_vs_decoys_type_entropy_rhoA_m8 | pct_native_smb_rate_hat_rhoA_m8 | delta_native_vs_decoys_smb_rate_hat_rhoA_m8 | pct_native_type_entropy_rhoB_m8 | pct_native_type_entropy_rhoBvel_m8 | pct_combined_min(A,B)_type_entropy_m8 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 1ctf | 630 (failQC=1) | 0.543 | 0.084 | 0.456 | -0.090 | 0.710 | 0.073 | 0.543 |
| 1r69 | 676 (failQC=0) | 0.030 | -0.970 | 0.037 | -0.963 | 0.562 | 0.482 | 0.030 |
| 1sn3 | 660 (failQC=0) | 0.285 | -0.432 | 0.302 | -0.398 | 0.382 | 0.642 | 0.285 |
| 2cro | 673 (failQC=1) | 0.016 | -0.984 | 0.036 | -0.964 | 0.502 | 0.016 | 0.016 |
| 3icb | 654 (failQC=0) | 0.170 | -0.661 | 0.318 | -0.364 | 0.865 | 0.213 | 0.170 |
| 4pti | 686 (failQC=1) | 0.159 | -0.688 | 0.226 | -0.557 | 0.351 | 0.227 | 0.159 |
| 4rxn | 677 (failQC=0) | 0.171 | -0.663 | 0.258 | -0.487 | 0.619 | 0.109 | 0.171 |

- Fraction targets with native in lowest 10% of decoys (rhoA type_entropy, m=8): 0.286
- Fraction targets with native in lowest 10% of decoys (rhoB=parity type_entropy, m=8): 0.000
- Fraction targets with native in lowest 10% of decoys (rhoB=vel type_entropy, m=8): 0.286
- Fraction targets with native in lowest 10% by min(A,Bparity) (type_entropy, m=8): 0.286

