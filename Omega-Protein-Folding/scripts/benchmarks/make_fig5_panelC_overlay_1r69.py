#!/usr/bin/env python3
from __future__ import annotations

"""
Figure 5 panel C (draft): CA-trace overlay on 1R69.

Produces a simple 3D overlay plot (native vs Omega vs QUARK vs I-TASSER)
aligned to native via Kabsch on Cα coordinates.
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import requests


def repo_root() -> Path:
    return Path(__file__).resolve().parents[2]


def _kabsch_align(P: np.ndarray, Q: np.ndarray) -> np.ndarray:
    """Align P onto Q and return aligned P."""
    root = repo_root()
    import sys

    sys.path.insert(0, str((root / "scripts" / "pdb_bench").resolve()))
    from bench_utils import kabsch  # type: ignore

    P_aligned, _, _ = kabsch(P, Q)
    return np.asarray(P_aligned, dtype=np.float64)


def _load_ca(path: Path) -> np.ndarray:
    root = repo_root()
    import sys

    sys.path.insert(0, str((root / "scripts" / "pdb_bench").resolve()))
    from coords_io import load_ca_coords  # type: ignore

    _, ca = load_ca_coords(path)
    return np.asarray(ca, dtype=np.float64)


def main() -> None:
    root = repo_root()
    run_dir = root / "docs" / "runs" / "quark_itasser_homology_ablation"

    # Native: if local decoy cache is absent, download mmCIF from RCSB into pdb_cache (gitignored).
    native = root / "data" / "raw" / "pdb_cache" / "1R69.cif"
    if not native.exists():
        native.parent.mkdir(parents=True, exist_ok=True)
        url = "https://files.rcsb.org/download/1R69.cif"
        r = requests.get(url, timeout=60)
        r.raise_for_status()
        native.write_bytes(r.content)

    omega = run_dir / "omega_models" / "omega_pred_1r69_seed0.pdb"
    quark = run_dir / "server_models" / "1r69" / "quark_QA16125_model1.pdb"
    itasser = run_dir / "server_models" / "1r69" / "itasser_S820065_model1.pdb"

    ca_native = _load_ca(native)
    ca_omega = _kabsch_align(_load_ca(omega), ca_native)
    ca_quark = _kabsch_align(_load_ca(quark), ca_native)
    ca_it = _kabsch_align(_load_ca(itasser), ca_native)

    fig = plt.figure(figsize=(7.6, 6.4), dpi=200)
    ax = fig.add_subplot(111, projection="3d")

    def plot_line(ca: np.ndarray, *, color: str, label: str, lw: float = 2.0, alpha: float = 1.0):
        ax.plot(ca[:, 0], ca[:, 1], ca[:, 2], color=color, lw=lw, alpha=alpha, label=label)

    plot_line(ca_native, color="#7f7f7f", label="native", lw=2.2, alpha=0.9)
    plot_line(ca_omega, color="#F28E2B", label="omega (CA trace)", lw=2.0, alpha=0.95)
    plot_line(ca_quark, color="#4E79A7", label="quark model1", lw=1.8, alpha=0.85)
    plot_line(ca_it, color="#59A14F", label="i-tasser model1", lw=1.8, alpha=0.85)

    ax.set_title("Figure 5C (draft) | 1R69 Cα overlay (aligned to native)")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_zlabel("z")
    ax.legend(loc="upper left", fontsize=9)
    ax.grid(False)

    # Make equal aspect for 3D
    pts = np.vstack([ca_native, ca_omega, ca_quark, ca_it])
    mn = pts.min(axis=0)
    mx = pts.max(axis=0)
    ctr = (mn + mx) / 2.0
    span = float(np.max(mx - mn))
    for setter, c in zip([ax.set_xlim, ax.set_ylim, ax.set_zlim], ctr):
        setter(c - span / 2.0, c + span / 2.0)

    out_dir = (run_dir / "figures").resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    out_png = out_dir / "fig5_panelC_overlay_1r69.png"
    fig.savefig(out_png, bbox_inches="tight")
    plt.close(fig)
    print(f"[ok] wrote: {out_png}")


if __name__ == "__main__":
    main()

