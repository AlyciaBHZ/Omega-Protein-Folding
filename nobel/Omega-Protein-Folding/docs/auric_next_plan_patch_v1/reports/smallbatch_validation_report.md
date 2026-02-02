# Small-batch validation (auric shared-PCA + rhoB=vel)

This report compares two configurations on the same small PDB set (see `ids_small.txt`).

Controls: shuffle (5 reps) + block-shuffle k=8 (3 reps) + random (3) + perturbed (2).


## Configurations

- **old_fixed+parity**: `rho_A` fixed u=(1,0,0), `rho_B` step-parity; (matches current integration default)

- **new_sharedpca_vel_median**: per-protein u from PCA(real y_perp), shared across all controls; `rho_B` uses perp-velocity with median threshold.


## Key metric: A/B consistency (z-score vs shuffle)

| metric       |   m |   corr_old |   corr_new |
|:-------------|----:|-----------:|-----------:|
| type_entropy |   6 |  -0.47312  |   0.598517 |
| type_entropy |   8 |  -0.203393 |   0.539331 |
| type_entropy |  10 |  -0.124472 |   0.564538 |
| smb_rate_hat |   6 |  -0.143518 |   0.239334 |
| smb_rate_hat |   8 |  -0.143518 |   0.239334 |
| smb_rate_hat |  10 |  -0.143518 |   0.239334 |


## Effect magnitude (mean |z| vs shuffle)

| metric       | readout   |   m |      new |      old |
|:-------------|:----------|----:|---------:|---------:|
| smb_rate_hat | A         |   6 | 0.882989 | 1.36673  |
| smb_rate_hat | A         |   8 | 0.882989 | 1.36673  |
| smb_rate_hat | A         |  10 | 0.882989 | 1.36673  |
| smb_rate_hat | B         |   6 | 1.67204  | 1.89525  |
| smb_rate_hat | B         |   8 | 1.67204  | 1.89525  |
| smb_rate_hat | B         |  10 | 1.67204  | 1.89525  |
| type_entropy | A         |   6 | 1.00622  | 1.80334  |
| type_entropy | A         |   8 | 0.982895 | 1.86483  |
| type_entropy | A         |  10 | 1.01619  | 1.766    |
| type_entropy | B         |   6 | 1.26546  | 0.801475 |
| type_entropy | B         |   8 | 1.04747  | 0.807742 |
| type_entropy | B         |  10 | 0.926315 | 0.984483 |


## Takeaway

- The **main failure mode** is not a numerical bug, but a *protocol issue*: the previous `rho_B` (parity) is not geometrically aligned with `rho_A`, so A/B can disagree.

- With **shared PCA axis** and **geometric rho_B (velocity)**, A/B agreement improves substantially (correlations turn positive for both metrics).

- We recommend using `smb_rate_hat` as the **primary** auric summary statistic and keeping `type_entropy` as a secondary check.

