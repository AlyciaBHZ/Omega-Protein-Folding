# Omega v1 — Supplementary Notes (v1.0)

**Project Omega Consortium**  
**Supplementary Notes — 2026-02-01**

This document provides detailed experiment-by-experiment context, additional tables, and a complete index of artifacts included in this repository.

> **Reading guide.** The main manuscript highlights the v1 core claims: (i) **phason feasibility requires a projector**, (ii) **alphabet expressivity** unlocks the 0.9 regime, and (iii) **Blind Sprint** performance without oracle losses.  
> This Supplementary records the full method evolution from early contact-map proxies to long-chain and membrane case studies.

---

## Supplementary Note 1 | Experiment timeline and method evolution

Omega evolved through the following stages (each stage is backed by audit logs and CSV outputs):

1. **Discrete folding proxy (Contact-map / Z128 low-bit embedding)**  
   Goal: verify that wormhole macro-operators and expander/Ramanujan shaping measurably change search trajectories under a low-resolution readout.

2. **6D cut-and-project feasibility diagnostics (phason proxy tests)**  
   Goal: test whether 6D lifting produces a discriminative signal between native proteins and controls; identify why naïve lifting underperforms without global constraints.

3. **Impedance-guided wormhole selection (\(\mathrm{Imp}(p;D)\))**  
   Goal: turn Omega from a static compiler into an adaptive agent that chooses macro-operators based on residual geometry, with audit-grade traces.

4. **Phason gating and the “hallucination gap”**  
   Goal: show that contact-map success can be hallucinated; introduce semantic residuals and gate-induced separation.

5. **From phason filter to phason projector**  
   Goal: resolve “hard reject freezes the search” by adding an inner-loop projector that adjusts \(w_0\) to restore feasibility.

6. **Scaling: long chains and a membrane case**  
   Goal: stress-test topology and packing; introduce piecewise \(w_0\) drift for long-range strain accumulation.

7. **Sprint experiments**  
   - *Oracle Sprint (upper bound)*: demonstrate an attainable ceiling under native-derived guidance.  
   - *Blind Sprint*: remove oracle TM and native membrane terms; replace with semantic residual proxies + sequence-only belt priors.

---

## Supplementary Note 2 | Early contact-map proxy benchmarks (wormholes / expander shaping / zeta shaping)

**Purpose.** Establish that (a) wormhole macro-operators change hit-rate and trajectory quality under a budget, and (b) spectral shaping influences the search distribution.

**Targets.** 1UBQ(A) and 2CI2(I), with \(N=16\) (contact bits = 120) and a contact definition Cα–Cα < 8Å.

**Key readout.** Contact map F1 between the discrete proxy and the target contact map.

**Summary table.**

| Variant             |   Mean best F1 |   s.d. |   Exact-hit rate |   n |
|:--------------------|---------------:|-------:|-----------------:|----:|
| wormhole_random     |         0.972  | 0.017  |              0.1 |  10 |
| omega_full_random   |         0.9655 | 0.0158 |              0   |  10 |
| baseline            |         0.9646 | 0.0203 |              0   |  10 |
| wormhole_expander   |         0.9627 | 0.0251 |              0.1 |  10 |
| zeta_only           |         0.959  | 0.025  |              0   |  10 |
| omega_full_expander |         0.9515 | 0.0184 |              0   |  10 |

**Interpretation.** Wormholes increase exact-hit probability under small budgets, even when the mean F1 uplift is modest. Zeta shaping stabilizes the trajectory distribution but can become overly conservative under this low-resolution readout.

Artifacts:
- `../../artifacts/reports/small_sample_hpa_omega_protein_experiments.md`
- `../../data/processed/small_sample_contactmap_proxy_results.csv`

---

## Supplementary Note 3 | Icosa-direction codec: 12-direction quantization advantage

**Purpose.** Test whether icosahedral vertex directions form an effective discrete direction codec for protein backbone steps (Cα–Cα vectors).

**Setup.** For each protein, quantize 15 backbone step directions using either:
- 12 icosahedral vertex directions; or
- 200 random 12-direction dictionaries (baseline).

**Result.**

| PDB   |   icosa_mean_deg |   rand12_mean_deg_mean |   rand12_mean_deg_std |   icosa_better_than_frac |
|:------|-----------------:|-----------------------:|----------------------:|-------------------------:|
| 1UBQ  |            24.29 |                  29.86 |                  5.38 |                     87   |
| 1CRN  |            24.12 |                  29.62 |                  4.65 |                     88.5 |
| 2CI2  |            23.87 |                  29.01 |                  4.36 |                     91   |
| 5PTI  |            21.04 |                  28.81 |                  5.09 |                     95.5 |

**Interpretation.** This does **not** prove “proteins are icosahedral”, but it does show that if a finite direction alphabet is desired, an icosahedral dictionary yields systematically smaller quantization error than random dictionaries of the same size—supporting its use as an **inner (non-orthogonal) codec layer**.

Artifacts:
- `../../data/processed/icosa_quantization_results.csv`

---

## Supplementary Note X | Rosetta constrained-relax “physical audit” of Omega backbones

**Purpose.** Use Rosetta/PyRosetta as an external **physical auditor** to verify that Omega’s geometrically generated Cα backbones lie within attractive basins of a standard all-atom energy function (i.e., are physically realizable under constrained relaxation).

**Protocol summary.** For each target, we generated a full-atom pose from FASTA and applied harmonic Cα coordinate restraints to the Omega Cα trace (coord_sd = 1.0 Å, coord_weight = 1.0). We then ran FastRelax for \(nstruct=20\) independent trajectories and recorded best/median total score (REU) and the best Kabsch-aligned Cα RMSD to the original trace. Heavy PDB outputs are cached locally (gitignored); the summary table below is committed.

**Supplementary Table (Rosetta constrained-relax audit, Omega-only; \(nstruct=20\)).**

| Target | Best total score (REU) | Median total score (REU) | Median–best gap (REU) | Best Cα RMSD to Omega trace (Å) |
|---|---:|---:|---:|---:|
| 1R69 | 398.6311 | 551.7821 | 153.1510 | 1.7814 |
| 2CRO | 579.8710 | 675.7528 | 95.8818 | 2.1368 |
| 4PTI | 501.0837 | 805.3716 | 304.2879 | 2.0392 |
| 1CTF | 564.3416 | 603.6578 | 39.3162 | 2.2019 |
| 1DTK | 553.4316 | 735.3389 | 181.9073 | 2.0952 |
| 1SHF-A | 400.2520 | 475.0400 | 74.7880 | 1.8886 |

Artifacts:
- Summary CSV: `../runs/quark_itasser_homology_ablation/rosetta_relax_summary.csv`
- Plan / protocol notes: `../runs/quark_itasser_homology_ablation/rosetta_relax_plan.md`

---

## Supplementary Note Y | Three-way physical audit on 1R69 (Omega vs QUARK vs I‑TASSER)

**Purpose.** Provide a reviewer-facing, target-level comparison showing that Omega’s geometric solver can yield backbones that are not only topologically plausible but also exhibit strong **compatibility with atomic relaxation**.

**Models.** We used committable copies of external server models (model1):
- QUARK job `QA16125`: `../runs/quark_itasser_homology_ablation/server_models/1r69/quark_QA16125_model1.pdb`
- I‑TASSER job `S820065`: `../runs/quark_itasser_homology_ablation/server_models/1r69/itasser_S820065_model1.pdb`
- Omega Cα trace: `../runs/quark_itasser_homology_ablation/omega_models/omega_pred_1r69_seed0.pdb`

**Topological accuracy (TM-score to native).** TM-scores are recorded in `../runs/quark_itasser_homology_ablation/server_summaries.csv` using the repository’s **Kabsch single-pass TM-score** (`scripts/pdb_bench/bench_utils.py`, `tm_score`), computed under 1:1 residue correspondence against RCSB native mmCIF (downloaded into `data/raw/pdb_cache/`, gitignored). This is not TM-align/US-align (no alignment search), but provides a consistent readout for same-length, index-aligned comparisons.

**Physical stability under constrained FastRelax (nstruct=20).** Figure 5B summarizes the drift distribution (Cα RMSD to each method’s input trace) across 20 independent relaxation trajectories. Omega shows a lower drift distribution and a shorter tail than QUARK/I‑TASSER under the same restraint strength (coord_sd=1.0 Å, coord_weight=1.0).

Artifacts:
- Figure 5 panels: `../runs/quark_itasser_homology_ablation/figures/fig5_panelAB_1r69.png`, `../runs/quark_itasser_homology_ablation/figures/fig5_panelC_overlay_1r69.png`
- nstruct summaries: `../runs/quark_itasser_homology_ablation/rosetta_relax_server_models_1r69_n5.csv`, `../runs/quark_itasser_homology_ablation/rosetta_relax_server_models_1r69_n20.csv`
- per-trajectory drift/energy (committable):  
  `../runs/quark_itasser_homology_ablation/rosetta_relax_omega_1r69_n20_per_struct.csv`,  
  `../runs/quark_itasser_homology_ablation/rosetta_relax_quark_1r69_n20_per_struct.csv`,  
  `../runs/quark_itasser_homology_ablation/rosetta_relax_itasser_1r69_n20_per_struct.csv`

---

## Supplementary Note Z | Additional three-way physical audit on 2CRO

To test whether the 1R69 three-way result generalizes beyond a single target, we repeated the same **three-way constrained FastRelax audit** on target **2CRO** (Omega Cα trace vs QUARK model1 vs I‑TASSER model1; \(nstruct=20\)). As in 1R69, Omega exhibits a **left-shifted drift distribution** (lower median drift and a lower best drift) compared to QUARK and I‑TASSER under the same restraint strength, while also reaching a lower Rosetta energy floor.

**Key numbers (2CRO, nstruct=20):**
- **Omega**: best drift **2.01 Å**, median drift **2.27 Å**; best score **560.94 REU**
- **QUARK**: best drift **2.21 Å**, median drift **2.48 Å**; best score **665.05 REU**
- **I‑TASSER**: best drift **2.48 Å**, median drift **2.75 Å**; best score **928.21 REU**

Artifacts:
- Three-way summary: `../runs/quark_itasser_homology_ablation/rosetta_relax_threeway_2cro_n20.csv`
- Per-trajectory drift/energy (committable):  
  `../runs/quark_itasser_homology_ablation/rosetta_relax_omega_2cro_n20_per_struct.csv`,  
  `../runs/quark_itasser_homology_ablation/rosetta_relax_quark_2cro_n20_per_struct.csv`,  
  `../runs/quark_itasser_homology_ablation/rosetta_relax_itasser_2cro_n20_per_struct.csv`
- Figure: `../runs/quark_itasser_homology_ablation/figures/threeway_2cro_panelAB.png`

**Upgraded Figure 5 (multi-target).** A combined multi-target view (1R69 + 2CRO) is provided at:
- `../runs/quark_itasser_homology_ablation/figures/fig5_multi_target.png`

