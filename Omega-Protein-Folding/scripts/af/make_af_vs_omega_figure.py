from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.image as mpimg  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402


def ensure_dir(p: Path) -> None:
    p.mkdir(parents=True, exist_ok=True)


def _read_image(path: Path) -> np.ndarray:
    img = mpimg.imread(str(path))
    if img.ndim == 2:
        img = np.stack([img, img, img], axis=-1)
    return img


def draw_af_blackbox(ax: plt.Axes, *, plddt_csv: Path, title: str) -> None:
    """
    Draw a compact AF schematic and a pLDDT "confidence ribbon".
    """
    df = pd.read_csv(plddt_csv)
    p = df["pLDDT"].to_numpy(dtype=np.float64)
    p = np.clip(p, 0.0, 100.0) / 100.0

    ax.set_axis_off()
    ax.set_title(title, fontsize=10)

    # Schematic: Sequence -> [AF black box] -> Structure
    ax.text(0.05, 0.75, "sequence", fontsize=10, va="center", ha="left")
    ax.annotate("", xy=(0.32, 0.75), xytext=(0.18, 0.75), arrowprops=dict(arrowstyle="->", lw=1.6))
    rect = plt.Rectangle((0.32, 0.63), 0.26, 0.24, facecolor="black", alpha=0.85)
    ax.add_patch(rect)
    ax.text(0.45, 0.75, "AF", fontsize=12, color="white", va="center", ha="center")
    ax.annotate("", xy=(0.80, 0.75), xytext=(0.58, 0.75), arrowprops=dict(arrowstyle="->", lw=1.6))
    ax.text(0.82, 0.75, "structure", fontsize=10, va="center", ha="left")

    # Confidence ribbon (pLDDT as a colormap strip)
    strip = p[None, :]
    ax.imshow(strip, extent=(0.05, 0.95, 0.25, 0.38), aspect="auto", cmap="viridis", vmin=0.0, vmax=1.0)
    ax.text(0.05, 0.18, "pLDDT ribbon (B-factor)", fontsize=9, ha="left", va="center")
    ax.text(0.95, 0.18, f"median={np.median(p)*100.0:.1f}", fontsize=9, ha="right", va="center")


def draw_image(ax: plt.Axes, img_path: Path, title: str) -> None:
    ax.set_axis_off()
    ax.set_title(title, fontsize=10)
    img = _read_image(img_path)
    ax.imshow(img)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--af-plddt-csv", required=True)
    ap.add_argument("--omega-trajectory-png", required=True)
    ap.add_argument("--shadow-af-png", required=True)
    ap.add_argument("--shadow-omega-png", required=True)
    ap.add_argument("--debugger-af-png", required=True)
    ap.add_argument("--debugger-omega-png", required=True)
    ap.add_argument("--out-png", required=True)
    ap.add_argument("--title", default="AlphaFold vs Omega: black-box vs auditable folding program")
    args = ap.parse_args()

    root = Path(__file__).resolve().parents[2]
    def R(p: str) -> Path:
        pp = Path(p)
        return pp if pp.is_absolute() else (root / pp).resolve()

    af_plddt_csv = R(args.af_plddt_csv)
    omega_traj_png = R(args.omega_trajectory_png)
    shadow_af_png = R(args.shadow_af_png)
    shadow_omega_png = R(args.shadow_omega_png)
    dbg_af_png = R(args.debugger_af_png)
    dbg_omega_png = R(args.debugger_omega_png)
    out_png = R(args.out_png)
    ensure_dir(out_png.parent)

    fig = plt.figure(figsize=(12.5, 11.0), constrained_layout=True)
    fig.suptitle(str(args.title), fontsize=13)
    gs = fig.add_gridspec(3, 2, height_ratios=[1.1, 1.0, 0.75])

    # Row A: Process (black-box vs trajectory)
    axA0 = fig.add_subplot(gs[0, 0])
    draw_af_blackbox(axA0, plddt_csv=af_plddt_csv, title="(A) Process: AlphaFold (static output + confidence)")
    axA1 = fig.add_subplot(gs[0, 1])
    draw_image(axA1, omega_traj_png, "(A) Process: Omega (replayable dynamics + audit)")

    # Row B: Geometry (shadow projection)
    axB0 = fig.add_subplot(gs[1, 0])
    draw_image(axB0, shadow_af_png, "(B) Geometry: AFDB shadow projection (y_perp PCA)")
    axB1 = fig.add_subplot(gs[1, 1])
    draw_image(axB1, shadow_omega_png, "(B) Geometry: Omega shadow projection (y_perp PCA)")

    # Row C: Audit/debugger views
    axC0 = fig.add_subplot(gs[2, 0])
    draw_image(axC0, dbg_af_png, "(C) Audit: AF confidence profile (pLDDT)")
    axC1 = fig.add_subplot(gs[2, 1])
    draw_image(axC1, dbg_omega_png, "(C) Audit: Omega trajectory debugger")

    fig.savefig(out_png, dpi=220)
    print(f"[ok] wrote: {out_png}")


if __name__ == "__main__":
    main()

