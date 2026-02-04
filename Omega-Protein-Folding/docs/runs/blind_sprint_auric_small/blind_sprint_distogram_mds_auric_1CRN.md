# Blind Sprint (sparse-distogram MDS) + Auric: 1CRN:A

- N: 46
- distogram pairs: 800 (min_sep=3)
- MDS: steps=1500, lr=0.03, clamp=5.0
- weights: w_chain=5.0, w_dist=1.0
- Auric: alphabet=pair72, m=8
- CSV (gitignored): `data/processed/blind_sprint/blind_sprint_distogram_mds_blind_sprint_distogram_mds_auric_1CRN.csv`
- Runtime: 24.7s

## Results (sorted by TM for evaluation)

| seed | tm | type_entropy_A | smb_rate_hat_B | rank_score |
| --- | --- | --- | --- | --- |
| 0 | 1.000 | 2.127 | 0.278 | -0.481 |
| 9 | 1.000 | 2.217 | 0.310 | -0.505 |
| 2 | 1.000 | 1.129 | 0.358 | -0.297 |
| 6 | 0.306 | 3.369 | 0.294 | -0.733 |
| 7 | 0.306 | 2.138 | 0.417 | -0.511 |
| 5 | 0.306 | 2.692 | 0.301 | -0.599 |
| 3 | 0.306 | 1.702 | 0.330 | -0.406 |
| 4 | 0.305 | 2.042 | 0.359 | -0.480 |
| 8 | 0.305 | 2.744 | 0.316 | -0.612 |
| 1 | 0.305 | 2.743 | 0.367 | -0.622 |

