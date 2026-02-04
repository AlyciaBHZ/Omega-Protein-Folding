# PhasonFold — Paper (LaTeX)

This folder contains a **sectioned LaTeX submission scaffold** for the PhasonFold manuscript.

## Build

- Main entry: `paper/main.tex`
- Compile (recommended): XeLaTeX

## Assets layout

- Figures: `paper/assets/figures/` (e.g. `fig1.pdf`, `fig2.png`, ...)
- Tables: `paper/assets/tables/` (e.g. `table1.tex`, `table2.tex`, ...)

The LaTeX uses helper macros that will **auto-include** `figN.(pdf|png|jpg)` if present; otherwise it shows a placeholder box so the manuscript still compiles.

## Evidence integration status

The repo is organized so that the *paper-ready* artifacts live in:

- `paper/assets/`: figures + tables for LaTeX
- `artifacts/reports/`: experiment reports referenced by the manuscript
- `data/processed/`: CSVs used by the manuscript/tables
- `scripts/evidence_v1/`: analysis scripts + specs for reproducibility

Raw/extracted source packs can be removed after integration to keep the repository lightweight.

## Regenerate (optional, requires external sources)

If you have a DOCX source available locally, you can regenerate the LaTeX scaffold from repo root:

```powershell
python tools/docx_to_paper_latex.py --docx path\to\manuscript.docx --out paper
```

## Import ZIP evidence assets (optional)

When `Omega_v1_evidence_repository.zip` and/or `sprint_to_0p9_outputs.zip` are available on disk:

```powershell
python tools/import_zip_assets.py --zip path\to\Omega_v1_evidence_repository.zip --paper paper
python tools/import_zip_assets.py --zip path\to\sprint_to_0p9_outputs.zip --paper paper
```

