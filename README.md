## PhasonFold evidence repository

This repository is the **open-source evidence pack** for the PhasonFold manuscript (code, scripts, and run artifacts).

- **Project root**: `Omega-Protein-Folding/`
- **Paper (LaTeX)**: `Omega-Protein-Folding/paper/`

### Quickstart

```bash
cd Omega-Protein-Folding

# build the paper
xelatex -interaction=nonstopmode -halt-on-error -output-directory paper paper/main.tex
xelatex -interaction=nonstopmode -halt-on-error -output-directory paper paper/main.tex
```

### Code & data availability

The full source tree and manuscript artifacts are available on GitHub:

- `https://github.com/AlyciaBHZ/Omega-Protein-Folding`
- `https://github.com/AlyciaBHZ/Omega-Protein-Folding.git`

### Repo hygiene

We intentionally do **not** version-control local editor state (e.g. `.cursor/`) or local environments (`.venv/`).
