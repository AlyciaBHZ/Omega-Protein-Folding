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
    ap.add_argument("--audit-csv", required=True, help="CSV from blind_sprint_auric.py --audit-out")
    ap.add_argument("--out-png", required=True, help="Output PNG path")
    ap.add_argument("--title", default="Omega audit debugger view")
    args = ap.parse_args()

    root = Path(__file__).resolve().parents[2]
    audit_csv = Path(args.audit_csv)
    if not audit_csv.is_absolute():
        audit_csv = (root / audit_csv).resolve()
    out_png = Path(args.out_png)
    if not out_png.is_absolute():
        out_png = (root / out_png).resolve()
    ensure_dir(out_png.parent)

    df = pd.read_csv(audit_csv)
    if len(df) == 0:
        raise ValueError(f"Empty audit CSV: {audit_csv}")
    t = df["t"].to_numpy(dtype=np.int32)

    auric_m = int(df["auric_m"].iloc[0]) if "auric_m" in df.columns else 0
    auric_every = int(df["auric_every"].iloc[0]) if "auric_every" in df.columns else 0
    checkpoints = []
    if auric_every > 0:
        for tt in t.tolist():
            if tt >= auric_m and (tt % auric_every == 0):
                checkpoints.append(tt)

    fig, ax = plt.subplots(1, 1, figsize=(8.0, 2.6), constrained_layout=True)
    ax.plot(t, df["ph_max"], lw=1.8, color="#d62728", label="ph_max")
    ax.plot(t, df["auric_penalty"], lw=1.8, color="#2ca02c", label="Auric penalty")
    for tt in checkpoints:
        ax.axvline(tt, lw=0.8, color="k", alpha=0.12)
    ax.set_xlabel("step t")
    ax.set_ylabel("value")
    ax.set_title(str(args.title))
    ax.grid(True, alpha=0.25)
    ax.legend(loc="upper right", frameon=False, ncol=2)
    fig.savefig(out_png, dpi=220)
    print(f"[ok] wrote: {out_png}")


if __name__ == "__main__":
    main()

