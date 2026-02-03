# Rosetta / PyRosetta 最小验证计划（Omega vs QUARK/I‑TASSER）

## 分支（请在服务器上切到这里查看/运行）

- **开发分支**：`feature-rosetta-relax-validation`
- 目的：把 Rosetta 作为**物理/几何验证器**接入，不做 Rosetta 从零预测（先做最小可复现实验）。

## 为什么需要 PyRosetta（关键约束）

当前 Omega 在仓库内输出的预测结构是 **CA-only**（只有 Cα，且常为 poly-ALA），Rosetta 的 `relax` **不能直接接受 CA-only PDB**。\n因此采用 PyRosetta：用 FASTA 建全原子 pose，再用 **Cα 坐标约束**把 pose “贴合”Omega 的 Cα trace，然后做 `FastRelax`。\n\n## 目标集合（本轮最小实验）

与 QUARK/I‑TASSER 同源排除 benchmark 保持一致的 6 个 target：

- `1R69`, `2CRO`, `4PTI`, `1CTF`, `1DTK`, `1SHF-A`

## 现有输入（已经在 repo 里准备好）

### 1) 序列（服务器提交/本地建模）
- `docs/runs/quark_itasser_homology_ablation/targets_6.fasta`

### 2) Omega 的 CA-only 输出（用于 Rosetta 约束 relax）
- 生成脚本：`scripts/benchmarks/generate_omega_models_for_rosetta.py`
- 输出目录：`docs/runs/quark_itasser_homology_ablation/omega_models/`
  - `omega_pred_<target_id>_seed0.pdb`（CA-only）
  - `omega_audit_<target_id>_seed0.csv`（每步 audit）
- 清单（给 relax 脚本用）：`docs/runs/quark_itasser_homology_ablation/omega_models_manifest.csv`

## Rosetta / PyRosetta 运行方式（两步走）

### Step A：先 dry-run（不需要 PyRosetta）
用途：检查路径、长度是否匹配，生成/覆盖 summary CSV（标记 `dry_run`）。\n
```bash
python scripts/benchmarks/pyrosetta_relax_validate.py --dry-run
```

输出（可提交、轻量）：
- `docs/runs/quark_itasser_homology_ablation/rosetta_relax_summary.csv`

### Step B：Omega-only relax（需要 PyRosetta）
你已选择先做 **smoke test**（nstruct=5），通过后再扩到 nstruct=20。\n
smoke test：

```bash
python scripts/benchmarks/pyrosetta_relax_validate.py --nstruct 5 --coord-sd 1.0 --coord-weight 1.0
```

正式版（按计划）：

```bash
python scripts/benchmarks/pyrosetta_relax_validate.py --nstruct 20 --coord-sd 1.0 --coord-weight 1.0
```

输出（可提交、轻量）：
- `docs/runs/quark_itasser_homology_ablation/rosetta_relax_summary.csv`

输出（不提交、大文件；gitignored）：
- `data/raw/rosetta_cache/quark_itasser_homology_ablation/omega/<target_id>/relax_<target_id>_best.pdb`

## 产物解释（最小可回答的问题）

对每个 target（Omega-only）至少记录：

- **能量**：best / median `total_score`（来自 PyRosetta scorefxn）\n- **稳定性**：`best_rmsd_ca`（relax 后的 pose 的 CA 与输入 CA trace 的 Kabsch RMSD）\n\n判读规则（先验、最小版）：\n- **稳定**：RMSD 小（对 N≈60 蛋白通常 <2–3Å 量级）且能量分布合理\n- **不稳定**：RMSD 大或能量异常（说明输入 CA trace 可能不在 Rosetta basin 附近）\n\n## 后续扩展（等 QUARK / I‑TASSER 回来再做）

当 QUARK / I‑TASSER 结果下载到本地 cache 后（不进 repo）：\n- `data/raw/quark_cache/<target_id>/QAxxxx/models/model1.pdb`\n- `data/raw/itasser_cache/<target_id>/Sxxxxx/models/model1.pdb`\n\n可在 `pyrosetta_relax_validate.py` 里追加“同协议 relax”以公平比较：\n- Omega vs QUARK vs I‑TASSER 的能量/漂移分布\n- 再结合 `docs/runs/quark_itasser_homology_ablation/tmalign_summary.csv` 的 TM-align 共识\n\n## 你需要申请/准备的 PyRosetta（当前阻塞点）

目前本机 Python 环境 **未安装 `pyrosetta`**。\n你需要申请/拿到 PyRosetta 后提供其中之一：\n- PyRosetta 的 whl 安装包路径，或\n- 一个已安装 PyRosetta 的 conda/env（告诉我 env 名称/启动方式）\n\n## 相关脚本索引

- 生成 Omega 输入：`scripts/benchmarks/generate_omega_models_for_rosetta.py`\n- PyRosetta relax 验证：`scripts/benchmarks/pyrosetta_relax_validate.py`\n- 报告模板：`docs/runs/quark_itasser_homology_ablation/rosetta_relax_report.md`\n
