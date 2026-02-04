# Blind Sprint suite (sparse-distogram MDS + Auric) — shortN10_v1

- IDs: `scripts/blind_sprint/ids_shortN_10.txt`
- seeds: 0,1,2,3,4
- n_pairs: 300 (capped per target)
- MDS: steps=1500, lr=0.03, clamp=5.0
- Auric: alphabet=pair72, m=8
- CSV (gitignored): `data/processed/blind_sprint/blind_sprint_suite_shortN10_v1.csv`
- Runtime: 78.3s

## Best-of-seeds per target (TM is evaluation-only)

| pdb_id | chain | N | seed | tm | n_pairs | type_entropy_A | smb_rate_hat_B |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 1R69 | A | 63 | 3 | 1.000 | 300 | 1.850 | 0.459 |
| 2CRO | A | 65 | 3 | 1.000 | 300 | 1.519 | 0.449 |
| 1CRN | A | 46 | 0 | 1.000 | 300 | 2.720 | 0.364 |
| 4PTI | A | 58 | 0 | 1.000 | 300 | 1.109 | 0.361 |
| 5PTI | A | 58 | 0 | 1.000 | 300 | 1.109 | 0.409 |
| 1CTF | A | 68 | 4 | 0.998 | 300 | 1.266 | 0.519 |
| 4RXN | A | 54 | 0 | 0.995 | 300 | 1.998 | 0.414 |
| 1UBQ | A | 76 | 2 | 0.993 | 300 | 3.713 | 0.514 |
| 2CI2 | I | 65 | 0 | 0.982 | 300 | 2.430 | 0.418 |
| 1SN3 | ? | 65 | 4 | 0.530 | 300 | 2.325 | 0.428 |

