# Auric 整合实验总结（2026-02-02）

本文档总结本次 Auric multi-scale certificate metrics（Fold\_m）在 `Omega-Protein-Folding` 中的整合、运行与对比结论，并将结果明确对齐你提出的痛点：在大规模 PDB（例如 1000k）上，单一标量（`ph_rms/ph_max`）信号弱且对定义/协议细节敏感，导致难以区分“理论不成立”还是“观测协议不对”。

## 0. 核心目标（与痛点对齐）

- 把 bench 从“单标量显著性”升级为 **可审计的协议对比系统**：同一 lifted 对象，多 readout（ρA/ρB）、多尺度 m、多种 null（random / perturbed / shuffle / blockshuffle），并报告哪里结论稳定、哪里脆弱。
- 引入 **shuffle / block-shuffle** 作为更强、更接近“保持边际统计但破坏长程相关”的主 null，从而更快暴露协议敏感性。

## 1. 我做了什么（实现内容）

### 1.1 新增/升级的核心模块

- `scripts/auric/readout.py`
  - **ρA thresholding**：默认使用 `median` 阈值（可选 `mean` / `zero` 作为 ablation），以降低“硬 sign 分割”的定义敏感性。
  - ρB 保持确定性规则（来自 6D step stream 的多数符号 + tie-break）。

- `scripts/pdb_bench/bench_utils.py`
  - 新增 bond-direction **shuffle** 与 **blockshuffle(block_len)**：保持每步长度与方向边际分布，但打乱顺序/相关长度。

### 1.2 bench 集成（让 shuffle/null 与 auric 变成一等公民）

- `scripts/pdb_bench/phason_stats.py`
  - 新增 CLI：
    - `--shuffle-reps`
    - `--blockshuffle-reps`
    - `--blockshuffle-k 4,8,16`
    - `--auric-rhoA-threshold {median,mean,zero}`
  - `samples` 中新增 group：`shuffle`、`blockshuffle_k{K}`
  - `auric_stats_summary_*.csv` 中新增：
    - 对 `type_entropy` 与 `smb_rate_hat`：δ(real vs random/pert/shuffle/blockshuffle) 以及 **native percentile among null reps**（更易解释、对 reps 更稳定）

- `scripts/pdb_bench/beam_lift_bench.py`
  - 同样新增 shuffle/blockshuffle 控制（**shuffle-before-lift**），并在 `beam_lift_auric_summary_*.csv` 里输出对应 δ 与 percentile（含 blockshuffle 家族）。

### 1.3 可视化与论文说明

- `scripts/pdb_bench/auric_plots.py`
  - 新增对 **shuffle/blockshuffle** 的 box/ECDF/multi-scale 曲线绘图（并兼容旧版本 matplotlib：`boxplot(labels=...)`）。
- `paper/sections/90_supplementary.tex`
  - 更新：ρA 默认 median 阈值与 CLI；新增 shuffle/blockshuffle null 的描述。

### 1.4 版本可追溯与 push blocker 处理

- 识别并从历史中移除 `legacy-docs/**/data/_release_cache/*.tar.gz` 的 >100MB 对象（否则 GitHub push 会被拒）。
- 增补 `.gitignore`，防止该类数据再次进入版本史。
- 已 push：
  - `origin/dev`
  - `origin/feature-auric`

## 2. 我实际跑了哪些实验（含产物路径）

### 2.1 Baseline（未整合/原本 dev 文档里的 phason-only）

对应你当前打开的报告/CSV（oracle triple232 + pmg random + pw64 + lin）：

- Report：
  - `artifacts/reports/pdb_phason_stats_n200_seed0_triple232_pmg_pw64_lin.md`
- Summary CSV：
  - `data/processed/pdb_bench/phason_stats_summary_n200_seed0_triple232_pmg_pw64_lin.csv`

### 2.2 dev200（之前的 auric 整合跑）

- Reports：
  - `artifacts/reports/pdb_phason_stats_dev200_auric.md`
  - `artifacts/reports/pdb_auric_stats_dev200_auric.md`
- CSV：
  - `data/processed/pdb_bench/phason_stats_summary_dev200_auric.csv`
  - `data/processed/pdb_bench/auric_stats_summary_dev200_auric.csv`
  - `data/processed/pdb_bench/beam_lift_bench_summary_dev200_auric_pair72_beam8.csv`
  - `data/processed/pdb_bench/beam_lift_auric_summary_dev200_auric_pair72_beam8.csv`

### 2.3 smoke20\_shuffle（我刚刚直接执行的“auric + shuffle + blockshuffle”）

IDs 自动从本机 `data/raw/pdb_cache/*.cif` 取前 20：

- IDs：
  - `data/processed/pdb_bench/pdb_ids_smoke20_from_cache.txt`
- Reports：
  - `artifacts/reports/pdb_phason_stats_smoke20_shuffle.md`
  - `artifacts/reports/pdb_auric_stats_smoke20_shuffle.md`
  - `artifacts/reports/pdb_beam_lift_bench_smoke20_shuffle.md`
- CSV：
  - `data/processed/pdb_bench/phason_stats_summary_smoke20_shuffle.csv`
  - `data/processed/pdb_bench/auric_stats_summary_smoke20_shuffle.csv`
  - `data/processed/pdb_bench/beam_lift_bench_summary_smoke20_shuffle.csv`
  - `data/processed/pdb_bench/beam_lift_auric_summary_smoke20_shuffle.csv`
- Figures：
  - `figures/pdb_bench/auric_*_smoke20_shuffle.png`

## 3. 结果对比（未整合 baseline vs auric 整合后）

下面所有 δ 都是“每个蛋白先算 δ，再对蛋白取均值”的口径。

### 3.1 Baseline（n200，phason-only）

来自 `phason_stats_summary_n200_seed0_triple232_pmg_pw64_lin.csv` 的均值：

- `cliffs_delta_real_vs_random`（global ph_rms）：**-0.070**
- `cliffs_delta_real_vs_perturbed`（global ph_rms）：**-0.090**
- `cliffs_delta_lin_vs_random`：**+0.035**
- `cliffs_delta_piecewise_vs_perturbed`：**+0.150**

**要点**：同一套数据里，global/piecewise/lin 的方向会变化，说明单标量对“观测协议定义”敏感。

### 3.2 dev200（auric 整合后的对照结果）

phason（来自 `pdb_phason_stats_dev200_auric.md`）：

- `cliffs_delta_real_vs_random`（global ph_rms）：**-0.145**
- `cliffs_delta_real_vs_perturbed`（global ph_rms）：**-0.118**

auric（来自 `pdb_auric_stats_dev200_auric.md`，指标为 `type_entropy`）：

- ρA：m=6/8/10 的 mean δ(real,random) 约 **0 附近（±0.013）**
- ρB：m=6/8/10 的 mean δ(real,random) 约 **-0.12 ~ -0.16**

**要点**：

- 在同一 lift+alphabet 下，ρA 与 ρB 的分离强度/方向不同 → 不确定性主要在 **readout/protocol 层**。
- dev200 的 `random-reps=1` 使 δ 更“量化”（台阶效应），因此后续需要把 shuffle/null + percentiles 做成主口径。

### 3.3 smoke20\_shuffle（加入 shuffle/blockshuffle 后的“协议敏感性显性化”）

phason（来自 `pdb_phason_stats_smoke20_shuffle.md`）：

- `cliffs_delta_real_vs_random`（global ph_rms）：**-0.222**
- 同时输出了 `delta_phmax_real_vs_shuffle` 与 `delta_phmax_real_vs_blockshuffle_k*`（用于观测 ph\_max 的协议敏感性）

auric（来自 `pdb_auric_stats_smoke20_shuffle.md`，指标为 `type_entropy`，并且可用 percentile 解读）：

- ρA：
  - m=6：mean δ(real,random)=**+0.611**；mean percentile(real vs random)≈**0.806**
  - m=8：mean δ(real,random)=**+0.667**；mean percentile(real vs random)≈**0.833**
- ρB：
  - m=6：mean δ(real,random)=**-0.500**；mean percentile(real vs random)≈**0.250**
  - m=8：mean δ(real,random)=**-0.444**；mean percentile(real vs random)≈**0.278**

**要点**：

- 加入更强的 null（shuffle/blockshuffle）后，ρA 与 ρB 在同一数据上出现 **方向相反** 的结论，且 percentile 提供了更直观解释：
  - ρA：real 通常处于 null 的高分位（更“结构化/有序”）
  - ρB：real 通常处于 null 的低分位（方向相反）
- 这正是你要解决的“分不清理论 vs 协议”的问题：现在可以明确定位到 **readout 选择会主导结论**，因此下一步应把 readout 与 null 家族作为审计轴。

## 4. 对比结论（回答“未整合 dev vs auric 整合后”）

1. **未整合（phason-only）**：信号整体偏弱，且对 piecewise/lin 等定义敏感；在大规模上会放大“难以解释”的问题。
2. **auric 整合后**：不是简单“更显著”，而是把问题拆解为可审计维度：
   - readout（ρA/ρB）、尺度（m）、null（random/pert/shuffle/blockshuffle）会系统性改变结论；
   - 当结论不一致时，能明确归因到哪个协议组件，而不是归咎于“样本不够/跑得不够多”。
3. **shuffle/blockshuffle 的价值**：比 random/pert 更能作为“主 null”揭示长程相关的敏感性；配合 percentile 能减少 δ 在小 reps 时的量化问题。

## 5. 下一步建议（直接面向 200/1000 规模）

- 用同一 ID 列表（例如 `pdb_ids_n200_seed0.txt` 或 `pdb_ids_n1000_seed0.txt`）重新跑：
  - `random-reps` 提升到 3+
  - `shuffle-reps` 8–16（主 null）
  - `blockshuffle-reps` 8（k=4/8/16）
  - m 取 `[6,8,10]`
- 报告中把 `type_entropy` 与 `smb_rate_hat` 的 **percentile(real among null)** 作为主读数，δ 作为补充。

