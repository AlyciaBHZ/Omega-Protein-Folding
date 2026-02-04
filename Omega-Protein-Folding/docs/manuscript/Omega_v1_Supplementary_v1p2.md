# Supplementary Information — Omega v1 Evidence Pack (All Tests)

This document is a *structured index + synthesis* of all experiments executed in the Omega v0.2 → v1 campaign.  
Raw CSVs, plots, and audit logs are stored in `../../data/processed/`, `../../figures/`, and `../../artifacts/logs/` with detailed per-family narratives in `../../artifacts/reports/`.

> **Scope note:** The evidence pack focuses on *audit artifacts* and *analysis scripts*. A full standalone Omega engine implementation is specified in `../../scripts/evidence_v1/specs/omega_engine_spec.md`.

---

## S0. Artifact map (what to open first)

### Core manuscripts
- `Omega_v1_NatureStyle_EN_v1p2.md` (main text; this folder)
- `../../paper/source/Omega_v1_NatureStyle_EN_v1p2.tex` (LaTeX source)

### Key A/B campaigns
- **Blind Sprint (v1):** `../../artifacts/reports/blind_sprint_AB_report_v1.md`
  - summary: `../../artifacts/reports/blind_sprint_ab_summary_v1.csv`
  - effects: `../../artifacts/reports/blind_sprint_ab_effects_v1.csv`
  - plots: `../../figures/blind_sprint_*.png`

- **Sprint to 0.9 (oracle-guided compilation):** `../../artifacts/reports/sprint0p9_AB_runner_report.md`
  - runs: `../../data/processed/sprint0p9_ab_runs.csv`
  - summary: `../../data/processed/sprint0p9_ab_summary.csv`
  - audit: `../../artifacts/logs/sprint0p9_*_audit.csv`
  - plots: `../../figures/sprint0p9_*.png`

- **θ★ feasibility (filter vs projector):** `../../artifacts/reports/step1_theta_star_step2_triple232_evidence_pack.md`
  - success table: `../../data/processed/theta_star_AB_success_1AKE_1TIM.csv`
  - counterfactual pass/fail: `../../data/processed/theta_star_counterfactual_passfail.csv`
  - plots: `../../figures/theta_star_AB_phmax_*.png`

### Scaling and stress
- **Long-chain L1→L2→L3:** `../../artifacts/reports/omega_long_L1_L2_L3_report.md`
  - runs: `../../data/processed/omega_long_L1_runs.csv`, `../../data/processed/omega_long_L2_runs.csv`, `../../data/processed/omega_long_L3_runs.csv`
  - summaries: `../../data/processed/omega_long_L3_final_summary.csv`, `../../data/processed/omega_long_L3_pairwise_AminusB.csv`

- **Projector vs filter on GFP (1EMA):** `../../artifacts/reports/omega_phason_projector_L3_1EMA_report.md`
  - summary: `../../data/processed/omega_projector_L3_1EMA_summary.csv`
  - per-seed audits: `../../data/processed/omega_projector_L3_1EMA_seed0_*_audit.csv`
  - trajectories: `../../figures/omega_projector_1EMA_seed0_*.png`

- **Projector extension (1TIM/1AKE/2O5P):** `../../artifacts/reports/omega_projector_extension_1TIM_1AKE_2O5P_report.md`

---

## S1. Early proxy experiments (short-chain contact-map proxy)

### S1.1 Contact-map proxy A/B (wormholes + ζ shaping)
- Narrative: `../../artifacts/reports/small_sample_hpa_omega_protein_experiments.md`
- Raw: `../../data/processed/small_sample_contactmap_proxy_results.csv` *(added in v1p2 repo)*
- Purpose: establish that **wormhole macro-operators** measurably change search behavior under strict budgets.

### S1.2 Contact-only Phason A/B (low-resolution readout)
- Narrative: `../../artifacts/reports/contact_only_phason_ab_report.md`
- Raw: `../../data/processed/contact_only_phason_ab_runs.csv`, `../../data/processed/contact_only_phason_ab_summary.csv`, `../../data/processed/contact_only_phason_ab_final_states.csv`
- Plots: `../../figures/contact_only_ab_hallgap_*.png`, `../../figures/contact_only_ab_phratio_*.png`
- Takeaway: in a low-res contact objective, phason mostly acts as a **regularizer**; hallucinations exist but can be hidden by coarse readouts.

---

## S2. Geometry priors: Icosa-direction codec test

- Raw: `../../data/processed/icosa_quantization_results.csv`
- Question: is an icosahedral direction dictionary systematically better than a random 12-direction dictionary for backbone step quantization?
- Takeaway (small sample): yes; icosahedral dictionary reduces angular quantization error across multiple PDB examples, supporting its role as an *inner non-orthogonal codec*.

---

## S3. Phason tests: from “proxy” to “necessity” to “projector”

### S3.1 Early phason proxy on PDB vs controls
- Narrative: `../../artifacts/reports/exp_phason_small_sample_report.md`
- Raw: `../../data/processed/phason_test_small_sample.csv`
- Takeaway: weak but nonzero separation; highlighted the need for *global consistency constraints* in lifting.

### S3.2 Synthetic necessity stress test
- Narrative: `../../artifacts/reports/exp_synthetic_phason_necessity_ab_report.md`
- Raw: `../../data/processed/exp_synthetic_phason_necessity_ab_runs.csv`, `../../data/processed/exp_synthetic_phason_necessity_ab_summary.csv`
- Plot: `../../figures/exp_synthetic_phason_necessity_scatter.png`
- Takeaway: in a deliberately adversarial synthetic setting, phason gating can be **decisive**.

### S3.3 Semantic residual + hard phason gate: hallucination-gap forcing
- Narrative: `../../artifacts/reports/semantic_residual_hard_phason_ab_report.md`
- Raw: `../../data/processed/semantic_residual_hard_phason_ab_summary.csv`, `../../data/processed/semantic_residual_hard_phason_ab_paired.csv`
- Plot: `../../figures/semantic_residual_scatter_all.png`
- Takeaway: phason can reject contact-high / semantic-low hallucinations, but early hard-gate suffered misalignment on some real cases → motivates projector.

### S3.4 Projector (repair) vs filter
- Key reports: `../../artifacts/reports/omega_phason_projector_L3_1EMA_report.md`, `../../artifacts/reports/projector_extension_theta_feasibility_addendum.md`
- Key takeaway: **phason becomes practically useful only when it is an active repair operator**, not a brittle filter.

---

## S4. Impedance-guided operator selection (Imp agent)

- Narrative (progress): `../../artifacts/reports/imp_adaptive_wormhole_selector_report.md`, `../../artifacts/reports/imp_adaptive_v2_benchmark_report.md`, `../../artifacts/reports/imp_agent_v5_v7_report.md`
- Raw: `../../data/processed/imp_adaptive_v2_benchmark_runs.csv`, `../../data/processed/imp_adaptive_v2_benchmark_summary.csv`, `../../data/processed/imp_agent_v5_v7_runs.csv`, `../../data/processed/imp_agent_v5_v7_summary.csv`, etc.
- Purpose: show Omega can evolve from a static compiler to a **dynamic, explainable agent** by routing operator families based on error geometry.

---

## S5. Stage I/II integration (θ_adaptive + alphabet extension + TM-score)

- Narrative: `../../artifacts/reports/phason_energy_correlation_analysis_q08.md`
- Key outputs:
  - TM benchmarks: `../../data/processed/omega_tmscore_benchmark_runs_q08.csv`, `../../data/processed/omega_tmscore_benchmark_summary_q08.csv`
  - Hallucination stats: `../../data/processed/hallucination_stats_q08.csv`
  - Plots: `../../figures/ab_*_q08.png`, `../../figures/chain_lift_phason_vs_random_ext.png`

---

## S6. Scaling: long-chain feasibility (N≈100) and L1→L2→L3

- Long-chain L1→L2→L3:
  - Report: `../../artifacts/reports/omega_long_L1_L2_L3_report.md`
  - Data: `../../data/processed/omega_long_L1_runs.csv`, `../../data/processed/omega_long_L2_runs.csv`, `../../data/processed/omega_long_L3_runs.csv`
  - Summary: `../../data/processed/omega_long_L3_final_summary.csv`, `../../data/processed/omega_long_L3_pairwise_AminusB.csv`

---

## S7. Final capability & intelligence tests

### S7.1 Sprint-to-0.9 (oracle-guided compilation)
- Report: `../../artifacts/reports/sprint0p9_AB_runner_report.md`
- Data: `../../data/processed/sprint0p9_ab_runs.csv`, `../../data/processed/sprint0p9_ab_summary.csv`
- Audits: `../../artifacts/logs/sprint0p9_*_audit.csv`
- Takeaway: operator stack can realize TM>0.9 when “wish” is perfect.

### S7.2 Blind Sprint (remove oracle TM and native belt)
- Report: `../../artifacts/reports/blind_sprint_AB_report_v1.md`
- Summary: `../../artifacts/reports/blind_sprint_ab_summary_v1.csv`
- Effects: `../../artifacts/reports/blind_sprint_ab_effects_v1.csv`
- Plots: `../../figures/blind_sprint_*.png`
- Takeaway: Omega retains advantage without optimizing TM-score; indicates **algorithmic intelligence** rather than pure metric hacking.
- Remaining limitation: distogram residual still computed against native distogram in this build → needs predicted residuals.

---

## Reproduction (analysis-only)
From repository root:

```bash
python -m venv .venv
source .venv/bin/activate
pip install -r scripts/evidence_v1/requirements.txt

python scripts/evidence_v1/make_figures.py --outdir figures_repro
python scripts/evidence_v1/make_tables.py  --outdir tables_repro
```

This regenerates key figures/tables from saved CSVs (it does **not** rerun Omega search).

