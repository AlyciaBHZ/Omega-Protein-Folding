# Omega engine interface spec (v1 evidence pack)

This evidence repository focuses on *audit artifacts* (CSVs/plots/logs) and the analysis scripts that reproduce the paper figures.
The full Omega engine implementation (search + operator execution) may live in a separate private repository.

This file provides a minimal **interface contract** so the engine can be reattached to this evidence pack for reproduction.

## Core objects

### Target
- `target_id`: string (e.g., "1TIM", "2O5P_tail250")
- `coords_native`: (N,3) float array (Cα coordinates)
- optional: `mask` / segment selection

### Codec (alphabet)
- `name`: one of `axis12`, `pair72`, `triple232`
- `steps`: list of integer 6D vectors (k,6), each with entries in {-1,0,1}
- `project_step_length(s, B_parallel)`: returns 3D step length distribution

### Phason module
- `B_parallel`, `B_perp`: (3,6) and (3,6) matrices
- `w0`: (3,) perpendicular-space slice offset
- `theta`: scalar threshold
- `phason_violation(n_path) -> ph_max, ph_series`
- `projector_repair(state, theta, budget) -> state_repaired, delta_w0, ph_after`

### Operator bank
Operators act on a `State` containing at minimum:
- `n_path`: (N,6) integer path (or incremental steps)
- `codec`: alphabet
- `audit`: append-only trace

Operators (suggested):
- `LocalStep(i, step_id)`
- `Wormhole_GlobalPatch(region, step_ids)`
- `Wormhole_HotspotPatch(region, step_ids)`
- `PhasonRelax(budget)`
- `PiecewiseProjector(segment_id, budget)`

### Impedance-guided policy (Step 1)
Inputs:
- residual geometry D: derived from current readout vs target (contacts, distogram RMSE, etc.)
Outputs:
- a distribution over operator families (Local / Expander_Global / Expander_Hot / Random)
- and a ranked list of concrete operator candidates

## Runner loop (simplified)

1. Initialize state (codec + random steps or pair72 baseline).
2. For t in 1..T:
   a. compute readout and residual D
   b. select operator family via Imp(p;D)
   c. propose K candidates; score
   d. apply phason gate:
      - filter: reject if ph_max>theta
      - projector: attempt `projector_repair` then accept if feasible
      - piecewise: projector per segment
   e. accept best candidate under objective
   f. emit audit row

3. Return best state + full audit.

## Repro hooks

To reproduce the paper *exactly*, the engine should emit at least these CSV columns (see `data/sprint_to_0p9_AB_summary.csv` etc):
- `protein`, `seed`, `mode`, `tm`, `rmse`, `f1`
- `theta_star`, `ph_max`, `ph_ratio`
- `rejects_total`, `projector_calls`, `w0_delta_norm`
- optional: membrane belt metrics (`belt_corr`, `belt_gap`)

