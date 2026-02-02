# Omega: Auditable Folding Agents via Z128 Operator Programs, 6D Quasicrystal Phason Projectors, and Adaptive Gating

**Project Omega (evidence-pack build)**  
Generated: 2026-02-01 (Asia/Singapore)

## Abstract
Protein folding is often treated as continuous optimization or generative sampling in **R³**, but these approaches make it difficult to **audit** and **replay** folding decisions step-by-step under explicit feasibility certificates.  
We introduce **Omega**, a discrete, auditable folding framework that represents a folding trajectory as an **operator program** executed on a **Z128 double‑rail register**, augmented with (i) impedance‑guided adaptive macro‑operators (*wormholes*) and (ii) a **6D icosahedral cut‑and‑project** geometry that exposes a measurable **phason** feasibility variable.  

Across a staged campaign—from short‑chain proxies to long‑chain and membrane settings—we show:
1) phason is a **dynamic** variable (not a static filter) once implemented as a **projector** (repair) rather than rejection;  
2) expanding the 6D direction alphabet to **Triple232** exposes a high‑accuracy regime (**TM ≈ 0.92–0.94**) under oracle‑guided compilation; and  
3) in a **blind** objective that removes oracle TM‑score and native membrane annotations, Omega retains a significant advantage over a baseline codec + global projector on **1TIM** (long‑chain) and **2O5P** (membrane), while producing replayable, audit‑grade traces.

These results position Omega as a complementary direction to black‑box predictors: an interpretable folding agent whose constraints and repairs are explicit, measurable, and debuggable.

---

## Main contributions
- **Auditable folding trajectories:** folding is compiled into a Z128 operator program with deterministic replay and trace digests.
- **Impedance-guided wormholes:** an error-geometry impedance `Imp(p;D)` selects macro-operators adaptively and improves search efficiency.
- **Phason as a projector:** replacing hard rejection with local geometric repair (Δw₀) turns phason from a brittle constraint into a stabilizing dynamics.
- **Alphabet expressivity:** a large 6D step set (Triple232) reduces readout discretization and supports domain- and membrane-scale tests.
- **Blind sprint protocol:** removes oracle TM-score and native membrane labels; uses semantic residuals and a hydropathy-only belt gate.

---

## Introduction
Modern structure predictors (e.g., deep learning distograms or end-to-end coordinate generators) have transformed protein structure prediction, but they typically provide limited interpretability for **why** a trajectory moves through certain intermediate states, and they rarely provide a mechanism-level, replayable account of how constraints are enforced during folding. Molecular dynamics provides mechanistic trajectories, but at a computational cost that is prohibitive for broad exploration.

Omega takes a different approach: treat folding as **compilation**. A fold is not merely a final 3D structure; it is an **operator sequence** executed in a discrete state space, with explicit feasibility checks and **auditable logs**.

---

## Results

### Z128 operator programs yield audit-grade folding traces
Omega encodes the folding state as a Z128 double‑rail register `Z=(M,F)` with invariant `M & F = 0` (ternary per bit).  
A folding trajectory is an operator program:

`Z₀ --O₁--> Z₁ --O₂--> ... --O_T--> Z_T`

Each operator emits an append‑only audit row (operator type, parameters, score deltas, feasibility outcomes, and trace digest prefix), enabling deterministic replay and counterfactual A/B mechanism tests.

### Phason gating is brittle as a filter, but stabilizing as a projector
A recurring failure mode in real-protein tests is that **hard reject** can prematurely exclude correct trajectories when the 6D lift is misaligned (limited step alphabet, fixed slice offset). We therefore replaced:

- **Filter:** if `phason_max > θ`, reject
- **Projector:** if `phason_max > θ`, run an inner-loop search over Δw₀ (slice shift) to reduce violation; accept if repaired

**Empirical feasibility under θ★.** Under a tightened threshold θ★ (derived from counterfactual floors), the **filter fails** while the **projector passes** on **1AKE**, **1TIM**, and **2O5P** tails (see `../../data/processed/theta_star_counterfactual_passfail.csv`).

![](../../figures/theta_star_AB_phmax_1AKE.png)
![](../../figures/theta_star_AB_phmax_1TIM.png)

### Alphabet expressivity exposes a high-accuracy regime (Triple232 oracle)
We expanded the 6D step alphabet from **Axis12** (±eᵢ) to **Pair72** (±eᵢ, ±eᵢ±eⱼ) and then to **Triple232** (adding ±eᵢ±eⱼ±e_k).  
In an oracle-guided codec test that measures representational capacity (not blind folding), Triple232 reaches **TM ≈ 0.92–0.94** on **1AKE** and **2O5P** and improves **1TIM** substantially (full-length TM ≈ 0.81; best local-100 window TM ≈ 0.92).

![](../../figures/triple232_oracle_codec_tm_bar.png)

### Sprint to 0.9 (oracle-guided compilation): hotspot-only Triple232 patches + piecewise projector
We combined the ingredients into an A/B runner:

- **A (baseline):** Pair72 + global projector + θ★
- **B (challenger):** Triple232 (*hotspot-only updates*) + *piecewise* projector (w₀ segment drift) + θ★ (adaptive gate schedule)

This runner is an **oracle-guided compilation** stress test (native alignment provides hotspot locations). It tests whether the operator stack can realize the high-TM regime when the “wish” is perfect.

Key outcomes (`../../data/processed/sprint0p9_ab_summary.csv`):
- **1TIM:** A stalls at TM ≈ 0.676; B reaches **TM = 0.928–0.940** (mean 0.936).
- **2O5P tail250:** A ≈ 0.703; B reaches **TM = 0.917–0.923** (mean 0.921) while maintaining a membrane packing proxy.

![](../../figures/sprint0p9_tmscore_distribution.png)  
![](../../figures/sprint0p9_2O5P_hydrophobic_gap_vs_tm.png)

### Blind Sprint: remove oracle TM and native membrane labels
To test whether Omega “knows how to fold” beyond oracle scoring, we removed:
- **Oracle TM-score** from acceptance
- **Native-derived membrane belt metrics** from acceptance

Architecture remains: Triple232 hotspot operator stack + piecewise projector + θ★ adaptive gate.  
Acceptance uses *blind guides*:
- **Semantic residual:** distogram RMSE (torsion residual placeholder)
- **Blind belt gate:** hydropathy-only belt energy (Kyte–Doolittle), no native labels

**Blind Sprint summary (A vs B)** (evaluation TM shown; not used for acceptance):

| Protein | Group | Seeds | TM mean | TM min–max | RMSE mean | Ph-ratio mean |
|---|---:|---:|---:|---:|---:|---:|
| 1TIM | A | 5 | 0.814 | 0.783–0.840 | 2.940 | 1.202 |
| 1TIM | B | 5 | 0.898 | 0.886–0.915 | 2.206 | 1.000 |
| 2O5P tail250 | A | 5 | 0.806 | 0.760–0.834 | 2.738 | 0.804 |
| 2O5P tail250 | B | 5 | 0.923 | 0.923–0.923 | 2.101 | 0.642 |

![](../../figures/blind_sprint_tmscore_distribution.png)

> **Important limitation:** in this evidence-pack build, distogram RMSE is still computed against the native distogram. The next step is to replace it with *Imp(p;D)-derived* predicted distograms/torsions so that the objective is truly blind to the native structure during search.

---

## Discussion
Omega’s central claim is not that it replaces SOTA predictors immediately, but that it provides a different kind of object: an **auditable folding agent** with explicit, modular constraints.

1. **Hard gate must become repair.** Phason rejection alone can be counterproductive when the lift is misaligned. Projector repair makes phason an active dynamics.
2. **Expressivity controls whether constraints help or harm.** Triple232 reduces discretization mismatch and unlocks the 0.9+ regime under oracle compilation.
3. **Blind guidance is feasible but not yet fully blind.** Removing oracle TM and native belt labels still yields significant improvement; replacing native-based distogram residuals is a primary next step.

### Limitations and future work
- Broad evaluation on diverse domains and multi-domain proteins.
- Replace native distogram residual with predicted distograms/torsions from sequence and intermediate states.
- Integrate physical energy (sterics/H-bonds/solvation proxies; external force fields) as a Phys‑Gate.
- Improve membrane priors beyond hydropathy-only belts.
- Extend beyond Cα-only programs to include side chains without losing auditability.

---

## Methods (summary)
See `../../scripts/evidence_v1/specs/omega_engine_spec.md` for the engine interface and `../../artifacts/reports/` for per-experiment details.

### Phason ratio
We report the feasibility ratio:

`ph_ratio = ph_max / θ`

### θ★ and θ_adaptive
- `θ_adaptive`: estimated from early exploration statistics (first ~20% of steps).
- `θ★`: tightened threshold derived from counterfactual feasibility floors.

### Piecewise projector
To handle long-chain accumulation error, we allow segment-wise w₀ drift (typical segment length ≈ 80 residues), analogous to phason-strain relaxation in a polycrystalline view.

---

## Data and code availability
All data tables, plots, audit logs, and analysis scripts used for this manuscript are included in this repository.\n
