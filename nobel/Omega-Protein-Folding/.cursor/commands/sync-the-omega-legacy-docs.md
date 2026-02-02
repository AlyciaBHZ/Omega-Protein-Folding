请执行 **the-omega 文档同步** 到本仓库的 `legacy-docs/`：

## 目标

- 从 `C:\Users\zwl62\OneDrive - National University of Singapore\Desktop\nobel\the-omega\docs` 拉取最新内容（`git pull`）
- 将 `docs/` 下的 **文档源码与脚本** 同步到本仓库的 `legacy-docs/`：
  - 文档：`*.tex`、`*.md`、`*.bib` 等
  - 脚本：`*.py`、`*.sh`、`*.ps1` 等
  - **默认只同步论文相关内容**（排除 `books/`、`assets/`）
  - **默认不引入**：图片/epub/pdf/zip 等二进制资产，以及 `data/` 等数据目录

## 执行方式

在仓库根目录运行下面命令（PowerShell）：

```powershell
powershell -ExecutionPolicy Bypass -File .\scripts\sync_the_omega_legacy_docs.ps1
```

默认会跳过 `*.pdf` 这类大文件（避免 OneDrive “云文件超时”），如果你确实需要把 PDF 也同步进来：

```powershell
powershell -ExecutionPolicy Bypass -File .\scripts\sync_the_omega_legacy_docs.ps1 -IncludePdf
```

如果你的 `the-omega/docs` 路径不同，用参数覆盖：

```powershell
powershell -ExecutionPolicy Bypass -File .\scripts\sync_the_omega_legacy_docs.ps1 -TheOmegaDocsPath "D:\path\to\the-omega\docs"
```

如果你想关闭“只同步源码”过滤（不推荐，容易把大量资产同步进来）：

```powershell
powershell -ExecutionPolicy Bypass -File .\scripts\sync_the_omega_legacy_docs.ps1 -SourcesOnly:$false
```

## 验证

- 运行结束后检查 `legacy-docs/` 下对应目录和文件已更新
- 如需更“严格”的镜像同步（会删除目标端多余文件），加 `-Mirror`：

```powershell
powershell -ExecutionPolicy Bypass -File .\scripts\sync_the_omega_legacy_docs.ps1 -Mirror
```
