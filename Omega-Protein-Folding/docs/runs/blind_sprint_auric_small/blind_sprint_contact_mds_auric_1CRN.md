# Blind Sprint (contact-MDS) + Auric rerank: 1CRN:A

- N: 46
- seeds: 0,1,2,3,4,5,6,7,8,9
- contact cutoff: 8.0 Å (min_sep=3)
- MDS: steps=1200, lr=0.03, clamp=5.0
- weights: w_chain=5.0, w_contact=1.0, contact_dist=6.0
- Auric: alphabet=pair72, m=8
- CSV (gitignored): `data/processed/blind_sprint/blind_sprint_contact_mds_blind_sprint_contact_mds_auric_1CRN.csv`
- Runtime: 5.3s

## Results (sorted by rank_score)

| seed | tm | contact_f1 | type_entropy_A | smb_rate_hat_B | rank_score |
| --- | --- | --- | --- | --- | --- |
| 3 | 0.086 | 0.187 | 1.216 | 0.340 | -0.125 |
| 9 | 0.073 | 0.187 | 2.003 | 0.351 | -0.284 |
| 0 | 0.096 | 0.192 | 2.138 | 0.371 | -0.310 |
| 7 | 0.080 | 0.188 | 2.127 | 0.365 | -0.310 |
| 4 | 0.084 | 0.187 | 2.189 | 0.370 | -0.324 |
| 8 | 0.073 | 0.187 | 2.254 | 0.331 | -0.330 |
| 2 | 0.075 | 0.187 | 2.286 | 0.403 | -0.350 |
| 5 | 0.071 | 0.185 | 2.514 | 0.319 | -0.382 |
| 1 | 0.119 | 0.185 | 2.875 | 0.354 | -0.461 |
| 6 | 0.089 | 0.192 | 3.343 | 0.322 | -0.541 |

