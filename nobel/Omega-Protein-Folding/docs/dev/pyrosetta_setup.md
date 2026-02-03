# 本机安装/配置 PyRosetta（用于 `pyrosetta_relax_validate.py`）

本仓库的 `scripts/benchmarks/pyrosetta_relax_validate.py` 会在 **PyRosetta 已安装**时运行 coordinate‑constrained `FastRelax`；未安装时可用 `--dry-run` 做输入校验并生成 `rosetta_relax_summary.csv`。

> 备注：你当前系统 `python` 是 3.12（`python -V`），通常**不匹配** PyRosetta 的 wheel 版本；建议用 conda 单独建一个与 wheel 匹配的 Python 版本环境（常见是 3.9/3.10）。

## 1) 建 conda 环境（示例：Python 3.10）

```bash
cd nobel/Omega-Protein-Folding
conda create -n omega-rosetta python=3.10 -y
conda activate omega-rosetta
python -m pip install -U pip
python -m pip install numpy pandas biopython
```

如你的 PyRosetta wheel 是 `cp39`，把 `python=3.10` 改为 `python=3.9`（以 wheel 文件名为准）。

## 2) 安装 PyRosetta（从 wheel 本地文件）

拿到 licensed 的 PyRosetta `.whl` 后（不要提交进 repo），在该环境里安装：

```bash
python -m pip install /path/to/PyRosetta-*.whl
```

快速检查：

```bash
python -c "import pyrosetta; pyrosetta.init('-mute all'); print('pyrosetta_ok')"
```

## 3) 跑计划里的 dry‑run（不需要 PyRosetta）

```bash
python scripts/benchmarks/pyrosetta_relax_validate.py --dry-run
```

会写出（可提交、轻量）：

- `docs/runs/quark_itasser_homology_ablation/rosetta_relax_summary.csv`

## 4) PyRosetta smoke test（建议先只跑 1 个 target）

脚本支持 `--targets` 过滤（逗号分隔，大小写不敏感）。例如先跑 `1r69`：

```bash
python scripts/benchmarks/pyrosetta_relax_validate.py --targets 1r69 --nstruct 5 --coord-sd 1.0 --coord-weight 1.0
```

跑通后再去掉 `--targets` 跑全量 6 个 target。

## 5) 大文件输出位置（已 gitignore）

best relaxed PDB 会写到（不提交）：

- `data/raw/rosetta_cache/quark_itasser_homology_ablation/omega/<target_id>/relax_<target_id>_best.pdb`

后续如果你把 QUARK / I‑TASSER 结果下载到本机，也建议放到：

- `data/raw/quark_cache/...`
- `data/raw/itasser_cache/...`

这两类路径也已在该子仓 `.gitignore` 中忽略。

