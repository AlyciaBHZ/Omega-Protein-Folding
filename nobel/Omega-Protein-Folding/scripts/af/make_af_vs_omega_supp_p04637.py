from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.image as mpimg  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402


def ensure_dir(p: Path) -> None:
    p.mkdir(parents=True, exist_ok=True)


def _read_image(path: Path) -> np.ndarray:
    img = mpimg.imread(str(path))
    if img.ndim == 2:
        img = np.stack([img, img, img], axis=-1)
    return img


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--debugger-png", default="docs/runs/af_vs_omega/debugger_afdb_P04637.png")
    ap.add_argument("--shadow-png", default="docs/runs/af_vs_omega/shadow_afdb_P04637.png")
    ap.add_argument("--out-png", default="paper/assets/figures/fig_af_vs_omega_supp_p04637.png")
    ap.add_argument("--title", default="Supplement: second AFDB example (P04637; mixed confidence)")
    args = ap.parse_args()

    root = Path(__file__).resolve().parents[2]
    dbg = (root / args.debugger_png).resolve()
    sh = (root / args.shadow_png).resolve()
    out = (root / args.out_png).resolve()
    ensure_dir(out.parent)

    fig = plt.figure(figsize=(10.5, 4.2), constrained_layout=True)
    fig.suptitle(str(args.title), fontsize=12)
    gs = fig.add_gridspec(1, 2, width_ratios=[1.2, 1.0])

    ax0 = fig.add_subplot(gs[0, 0])
    ax0.set_axis_off()
    ax0.set_title("AFDB P04637 pLDDT profile", fontsize=10)
    ax0.imshow(_read_image(dbg))

    ax1 = fig.add_subplot(gs[0, 1])
    ax1.set_axis_off()
    ax1.set_title("y_perp shadow projection (triple232 lift)", fontsize=10)
    ax1.imshow(_read_image(sh))

    fig.savefig(out, dpi=220)
    print(f"[ok] wrote: {out}")


if __name__ == "__main__":
    main()

