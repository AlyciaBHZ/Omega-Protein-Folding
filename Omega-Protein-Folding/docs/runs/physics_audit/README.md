# Physics audit layer (CA-level) — run artifacts

This folder contains **audit-only, lightweight** artifacts for a coarse-grained
physical consistency check. It does **not** modify any search algorithm.

## Artifacts

- `physics_vs_auric_toy_decoy.csv`
  - Per-PDB physics proxies + Auric geometry certificate scores
  - Data source: `docs/runs/toy_decoy_auric_smoke/pdbs/*.pdb`

## How to reproduce

From repo root:

```bash
python scripts/audit/audit_physics_vs_auric.py
python scripts/audit/make_physics_audit_figure.py
```

The paper-facing figure is written to:

- `paper/assets/figures/figS_physics_vs_auric.png`

# Physics audit layer (CA-level) — run artifacts

This folder contains **audit-only, lightweight** artifacts for a coarse-grained
physical consistency check. It does **not** modify any search algorithm.

## Artifacts

- `physics_vs_auric_toy_decoy.csv`
  - Per-PDB physics proxies + Auric geometry certificate scores
  - Data source: `docs/runs/toy_decoy_auric_smoke/pdbs/*.pdb`

## How to reproduce

From repo root:

```bash
python scripts/audit/audit_physics_vs_auric.py
python scripts/audit/make_physics_audit_figure.py
```

The paper-facing figure is written to:

- `paper/assets/figures/figS_physics_vs_auric.png`

