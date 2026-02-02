from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import pandas as pd  # noqa: E402


def ensure_dir(p: Path) -> None:
    p.mkdir(parents=True, exist_ok=True)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--audit-csv", required=True, help="CSV from blind_sprint_auric.py --audit-out")
    ap.add_argument("--out-png", required=True, help="Output PNG path")
    ap.add_argument("--title", default="Omega (white-box) trajectory audit")
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

    t = df["t"].to_numpy()

    fig, axs = plt.subplots(3, 1, figsize=(8.5, 8.0), sharex=True, constrained_layout=True)
    fig.suptitle(str(args.title))

    # Panel 1: contact / distogram objectives
    ax = axs[0]
    ax.plot(t, df["contact_f1"], lw=2.0, label="contact F1 (higher=better)", color="#1f77b4")
    ax.set_ylabel("F1")
    ax.grid(True, alpha=0.3)
    ax2 = ax.twinx()
    ax2.plot(t, df["dist_rmse"], lw=2.0, label="dist RMSE (lower=better)", color="#ff7f0e")
    ax2.set_ylabel("RMSE")

    # Make a combined legend
    lines = ax.get_lines() + ax2.get_lines()
    labels = [l.get_label() for l in lines]
    ax.legend(lines, labels, loc="upper right", frameon=False)

    # Panel 2: certificate / phason proxies
    ax = axs[1]
    ax.plot(t, df["auric_penalty"], lw=2.0, label="Auric penalty (lower=better)", color="#2ca02c")
    ax.set_ylabel("Auric")
    ax.grid(True, alpha=0.3)
    ax2 = ax.twinx()
    ax2.plot(t, df["ph_max"], lw=2.0, label="ph_max (lower=better)", color="#d62728")
    ax2.set_ylabel("ph_max")
    lines = ax.get_lines() + ax2.get_lines()
    labels = [l.get_label() for l in lines]
    ax.legend(lines, labels, loc="upper right", frameon=False)

    # Panel 3: total score
    ax = axs[2]
    ax.plot(t, df["score"], lw=2.0, label="total score (higher=better)", color="#9467bd")
    ax.set_xlabel("step t")
    ax.set_ylabel("score")
    ax.grid(True, alpha=0.3)
    ax.legend(loc="lower right", frameon=False)

    fig.savefig(out_png, dpi=220)
    print(f"[ok] wrote: {out_png}")


if __name__ == "__main__":
    main()

