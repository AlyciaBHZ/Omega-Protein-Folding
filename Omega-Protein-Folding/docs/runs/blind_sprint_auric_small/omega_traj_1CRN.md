# Blind Sprint + Auric (target=1CRN:A)

- N: 46
- codec: `pair72`
- beam: 64
- K per state: 32
- contact cutoff: 8.0 Å (min_sep=3)
- auric_m: 8 (every 5 steps)
- weights: w_contact=1.0, wA=0.2, wB=0.2

- CSV (gitignored): `data/processed/blind_sprint/blind_sprint_auric_runs_omega_traj_1CRN.csv`
- Runtime: 1714.0s

## Results

| seed | codec | beam | K | cutoff | min_sep | auric_m | auric_every | w_contact | w_rmse | wA | wB | tm | contact_f1 | dist_rmse | auric_penalty | bond_len | uvec_x | uvec_y | uvec_z | pdb_id | chain | N |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 0 | pair72 | 64 | 32 | 8.0000 | 3 | 8 | 5 | 1.0000 | 1.0000 | 0.2000 | 0.2000 | 0.2871 | 0.8508 | 1.0658 | 0.4892 | 3.8185 | 0.2838 | -0.3598 | -0.8888 | 1CRN | A | 46 |

