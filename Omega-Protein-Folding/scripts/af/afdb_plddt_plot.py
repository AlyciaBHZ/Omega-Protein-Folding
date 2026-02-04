from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402


def ensure_dir(p: Path) -> None:
    p.mkdir(parents=True, exist_ok=True)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--plddt-csv", required=True, help="CSV from afdb_download.py")
    ap.add_argument("--out-png", required=True, help="Output PNG path")
    ap.add_argument("--title", default="AlphaFold pLDDT (B-factor)")
    args = ap.parse_args()

    root = Path(__file__).resolve().parents[2]
    plddt_csv = Path(args.plddt_csv)
    if not plddt_csv.is_absolute():
        plddt_csv = (root / plddt_csv).resolve()
    out_png = Path(args.out_png)
    if not out_png.is_absolute():
        out_png = (root / out_png).resolve()
    ensure_dir(out_png.parent)

    df = pd.read_csv(plddt_csv)
    if "pLDDT" not in df.columns:
        raise ValueError("Expected column pLDDT in CSV")
    p = df["pLDDT"].to_numpy(dtype=np.float64)
    x = np.arange(len(p), dtype=np.float64)

    fig, ax = plt.subplots(1, 1, figsize=(8.0, 2.6), constrained_layout=True)
    ax.plot(x, p, lw=1.8, color="#1f77b4")
    ax.axhline(50.0, lw=1.0, color="#d62728", alpha=0.7)
    ax.axhline(70.0, lw=1.0, color="#ff7f0e", alpha=0.7)
    ax.set_ylim(0.0, 100.0)
    ax.set_xlim(0, max(1, len(p) - 1))
    ax.set_xlabel("residue index")
    ax.set_ylabel("pLDDT")
    ax.set_title(str(args.title))
    ax.grid(True, alpha=0.25)
    fig.savefig(out_png, dpi=220)
    print(f"[ok] wrote: {out_png}")


if __name__ == "__main__":
    main()

