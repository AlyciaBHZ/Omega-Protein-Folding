# Auric decoy benchmark (toy_v1)

- Manifest: `docs/runs/toy_decoy_auric_smoke/decoy_manifest.csv`
- Alphabet: `triple232`
- m list: 6, 8
- Axis mode: `best_native_axis`
- FailedGeometry (phason): ph_rms <= 100.0
- Targets used: 1
- Samples CSV (gitignored): `data/processed/decoy_bench/auric_decoys_samples_toy_v1.csv`
- Summary CSV (gitignored): `data/processed/decoy_bench/auric_decoys_summary_toy_v1.csv`
- Plots: `docs/runs/toy_decoy_auric_smoke/plots/`
- Geometry QC: max CA-CA <= 10.0 Å and outlier_frac <= 0.20

- Targets with native FailedGeometry under this QC: 0/1

## Native vs Decoy separation (summary)

Using m=8. δ<0 means native has lower metric (more ordered) than decoys.

| target_id | n_decoys_used | pct_native_type_entropy_rhoA_m8 | delta_native_vs_decoys_type_entropy_rhoA_m8 | pct_native_smb_rate_hat_rhoA_m8 | delta_native_vs_decoys_smb_rate_hat_rhoA_m8 | pct_native_type_entropy_rhoB_m8 | pct_native_type_entropy_rhoBvel_m8 | pct_combined_min(A,B)_type_entropy_m8 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| TOY1 | 9 (failQC=0) | 0.333 | -0.333 | 0.222 | -0.556 | 0.333 | 0.778 | 0.333 |

- Fraction targets with native in lowest 10% of decoys (rhoA type_entropy, m=8): 0.000
- Fraction targets with native in lowest 10% of decoys (rhoB=parity type_entropy, m=8): 0.000
- Fraction targets with native in lowest 10% of decoys (rhoB=vel type_entropy, m=8): 0.000
- Fraction targets with native in lowest 10% by min(A,Bparity) (type_entropy, m=8): 0.000

