# Auric integration: next-step patch (shared-PCA + geometric rho_B)

This patch is a *surgical* update to the existing PDB bench so that the
Auric protocol becomes **(i) reproducible**, **(ii) readout-consistent**, and
**(iii) harder to game** by the nulls.

It is designed to be dropped into the existing repo layout:

- `scripts/pdb_bench/phason_stats.py`
- `scripts/auric/readout.py`

The patch does **not** change your Omega architecture, alphabet, operator, or
lift; it only refines the **observation protocol** (how we binarize and audit
the phason/step stream) and keeps the protocol fixed across controls.

---

## What problem this solves

In the current Auric integration, you can observe **readout sensitivity**:

- `rho_A` (from `y_perp`) and `rho_B` (from `n_path` steps) can lead to
  opposite conclusions depending on how the 0/1 readout is defined.
- If you use any adaptive direction choice (e.g., PCA), and you allow that
  direction to be re-estimated *per control replicate*, the nulls are given an
  unfair advantage: each null replicate gets its "best" projection axis,
  washing out the separation.

This patch enforces: **one protein → one observation axis u → all controls share it.**

---

## Key idea

### 1) Protocol-fixed axis `u` per protein (shared-PCA)

Add `--auric-u-mode pca`:

- compute a deterministic PCA axis **once** from the **real** chain's `y_perp`
- then reuse that *same* `u` for:
  - `rho_A` on the real chain
  - `rho_A` on random / perturbed / shuffle / block-shuffle controls
  - `rho_B` (see below) on all groups

If you set `--auric-u-mode fixed`, we fall back to the canonical axis `u=(1,0,0)`.

### 2) Replace/extend `rho_B`

Add `--auric-rhoB-mode`:

- `parity` (legacy): step-parity majority vote from 6D step deltas
- `vel` (new, recommended): **perp-space velocity**
  - compute `Δy_perp[t] = Δn[t] · B_perp`
  - binarize by `sign(u · Δy_perp[t])` with a threshold (see next section)
- `pos` (optional): use `n_path` → `y_perp` and then apply the same `rho_A` rule

### 3) Thresholding that improves stability

Add `--auric-rhoB-threshold`:

- for `vel`, we recommend `median` (balances bits; improves cross-readout agreement)
- `zero` is also supported (velocity sign relative to 0)

---

## What I tested in this sandbox

I ran a small-batch test on real PDB chains (downloaded mmCIF) with
shuffle + block-shuffle controls, using the existing Omega bench pipeline.

**PDB IDs used:** see `reports/ids_small.txt`.

### Compared configurations

1. **Old (baseline in your current integration):**
   - `rho_A`: fixed axis, median threshold
   - `rho_B`: `parity`

2. **New (this patch, recommended):**
   - `--auric-u-mode pca` (shared across all controls)
   - `--auric-rhoB-mode vel`
   - `--auric-rhoB-threshold median`

### Result highlight (consistency gain)

Using **z-scores vs shuffle** (computed from the replicate distribution), the
agreement between readouts improves substantially.

Example (small batch, `m ∈ {6,8,10}`):

- `type_entropy`: corr(z_A, z_B) changes from negative (~ -0.47) to positive (~ +0.55 to +0.60)
- `smb_rate_hat`: corr(z_A, z_B) changes from negative (~ -0.14) to positive (~ +0.24)

Full numbers are saved in:

- `reports/ab_z_compare_smallbatch.csv`
- `reports/auric_stats_samples_*.csv`
- `reports/auric_stats_summary_*.csv`

This directly addresses the current failure mode where A/B readouts can disagree.

---

## How to run in your repo

### Smoke (n≈20)

```bash
python scripts/pdb_bench/phason_stats.py \
  --ids data/processed/pdb_bench/pdb_ids_smoke20.txt \
  --tag auric_sharedpca_vel_smoke \
  --alphabet triple232 \
  --seed 0 \
  --random-reps 5 --perturb-reps 3 --shuffle-reps 10 --blockshuffle-reps 5 --blockshuffle-k 8 \
  --auric --auric-m 6,8,10 \
  --auric-readouts all \
  --auric-u-mode pca \
  --auric-rhoB-mode vel \
  --auric-rhoB-threshold median
```

### Dev (n≈200)

Same command, but point `--ids` to your dev list and reduce reps if needed.

### Full (n=1000)

Same command; recommend keeping `m` small (6,8,10) for runtime.

---

## How to interpret the output

Prefer **shuffle** and **block-shuffle** as the primary nulls.

For each `(pdb_id, readout, m)` look at:

- `pct_*_real_vs_shuffle` (percentile of the real vs shuffle distribution)
- `delta_*_real_vs_shuffle` (Cliff's delta)
- `smb_rate_hat_*` (more stable) and `type_entropy_*` (secondary)

Additionally, compute **cross-readout agreement**:

- corr(z_A, z_B) across proteins should be **positive** if both readouts are
  capturing the same underlying regularity.

---

## Files in this patch

- `scripts/pdb_bench/phason_stats.py` (adds shared-PCA axis + rho_B modes)
- `scripts/auric/readout.py` (adds rho_B=vel/pos + PCA utility)
- `reports/*` (small-batch validation artifacts)

