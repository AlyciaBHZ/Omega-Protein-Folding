# Blind Sprint ablations (tag=blind_sprint_ablations_4state_v3)

## Protocol

- targets: 1CTF, 1R69, 1SN3, 2CRO, 3ICB, 4PTI, 4RXN
- seeds: 0, 1, 2, 3, 4
- codec: `pair72`
- beam: 128
- K: 64
- workers: 8

## Best-of-seeds TM per target

(TM is evaluation-only; not used for acceptance.)

| pdb_id | A | A+Bparity | A+Bvel | Bparity | Bvel | baseline |
| --- | --- | --- | --- | --- | --- | --- |
| 1CTF | 0.890 | 0.893 | 0.891 | 0.835 | 0.841 | 0.855 |
| 1R69 | 0.880 | 0.887 | 0.887 | 0.881 | 0.875 | 0.884 |
| 1SN3 | 0.882 | 0.882 | 0.885 | 0.869 | 0.839 | 0.313 |
| 2CRO | 0.868 | 0.865 | 0.863 | 0.862 | 0.870 | 0.856 |
| 3ICB | 0.219 | 0.216 | 0.216 | 0.786 | 0.785 | 0.752 |
| 4PTI | 0.781 | 0.816 | 0.816 | 0.815 | 0.837 | 0.758 |
| 4RXN | 0.818 | 0.831 | 0.803 | 0.839 | 0.843 | 0.800 |

## Acceptance proxy (best-of-seeds TM>0.6 fraction)

| group | frac_best_tm_gt_0p6 | n_targets |
| --- | --- | --- |
| A | 0.857 | 7 |
| A+Bparity | 0.857 | 7 |
| A+Bvel | 0.857 | 7 |
| Bparity | 1.000 | 7 |
| Bvel | 1.000 | 7 |
| baseline | 0.857 | 7 |

## Artifacts

- CSV (gitignored): `data/processed/blind_sprint/blind_sprint_ablations_blind_sprint_ablations_4state_v3.csv`
- Plots: `docs/runs/blind_sprint_ablations_4state_v3/plots/`

## Meta

```json
{
  "ids": [
    "1CTF",
    "1R69",
    "1SN3",
    "2CRO",
    "3ICB",
    "4PTI",
    "4RXN"
  ],
  "seeds": [
    0,
    1,
    2,
    3,
    4
  ],
  "codec": "pair72",
  "beam": 128,
  "K": 64,
  "workers": 8,
  "groups": [
    {
      "name": "baseline",
      "wA": 0.0,
      "wB": 0.0,
      "rhoB_mode": "parity",
      "rhoB_threshold": "median"
    },
    {
      "name": "A",
      "wA": 0.2,
      "wB": 0.0,
      "rhoB_mode": "parity",
      "rhoB_threshold": "median"
    },
    {
      "name": "Bparity",
      "wA": 0.0,
      "wB": 0.2,
      "rhoB_mode": "parity",
      "rhoB_threshold": "median"
    },
    {
      "name": "A+Bparity",
      "wA": 0.2,
      "wB": 0.2,
      "rhoB_mode": "parity",
      "rhoB_threshold": "median"
    },
    {
      "name": "Bvel",
      "wA": 0.0,
      "wB": 0.2,
      "rhoB_mode": "vel",
      "rhoB_threshold": "median"
    },
    {
      "name": "A+Bvel",
      "wA": 0.2,
      "wB": 0.2,
      "rhoB_mode": "vel",
      "rhoB_threshold": "median"
    }
  ]
}
```

