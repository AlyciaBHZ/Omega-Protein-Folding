from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import List, Tuple

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402


def ensure_dir(p: Path) -> None:
    p.mkdir(parents=True, exist_ok=True)


def parse_pdb_ca_coords(pdb_path: Path) -> np.ndarray:
    """
    Parse CA coordinates from a PDB (single model).
    Works for both AFDB PDBs and our minimal CA-only predicted PDB exports.
    """
    coords: List[List[float]] = []
    res_ids: List[Tuple[str, int, str]] = []
    seen = set()
    for line in pdb_path.read_text(encoding="utf-8", errors="replace").splitlines():
        if not line.startswith("ATOM") or len(line) < 54:
            continue
        atom = line[12:16].strip()
        if atom != "CA":
            continue
        altloc = line[16:17]
        if altloc not in {" ", "A"}:
            continue
        chain = line[21:22].strip() or "?"
        resseq_str = line[22:26].strip()
        icode = line[26:27].strip()
        try:
            resseq = int(resseq_str)
        except ValueError:
            continue
        rid = (chain, resseq, icode)
        if rid in seen:
            continue
        seen.add(rid)
        try:
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
        except ValueError:
            continue
        coords.append([x, y, z])
        res_ids.append(rid)
    if len(coords) < 2:
        raise ValueError(f"No CA coords parsed from {pdb_path}")
    return np.asarray(coords, dtype=np.float64)


def pca_2d(Y: np.ndarray) -> np.ndarray:
    Y = np.asarray(Y, dtype=np.float64)
    Yc = Y - Y.mean(axis=0, keepdims=True)
    if np.allclose(Yc, 0.0):
        return np.zeros((Y.shape[0], 2), dtype=np.float64)
    _, _, Vt = np.linalg.svd(Yc, full_matrices=False)
    return Yc @ Vt[:2].T


def auric_type_entropy_from_yperp(y_perp: np.ndarray, *, m: int = 8) -> float:
    """
    Compute a single certificate scalar for annotation: type_entropy of rho_A (PCA u, median threshold).
    """
    root = Path(__file__).resolve().parents[2]
    sys.path.insert(0, str((root / "scripts" / "auric").resolve()))
    from fold_m import fold_sliding_windows_bits  # noqa: E402
    from metrics import metrics_for_stream  # noqa: E402
    from readout import rho_A_from_yperp  # noqa: E402

    bits_A = rho_A_from_yperp(y_perp, u_mode="pca", threshold="median")
    folded = fold_sliding_windows_bits(bits_A, m=int(m))
    met = metrics_for_stream(bits_A, folded)
    return float(met.get("type_entropy", float("nan")))


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--pdb", required=True, help="Input PDB path")
    ap.add_argument("--out-png", required=True, help="Output PNG path")
    ap.add_argument("--label", default="", help="Label for plot title")
    ap.add_argument("--alphabet", default="triple232", choices=["axis12", "pair72", "triple232"])
    ap.add_argument("--m", type=int, default=8, help="Auric window length for annotation")
    args = ap.parse_args()

    root = Path(__file__).resolve().parents[2]
    pdb_path = Path(args.pdb)
    if not pdb_path.is_absolute():
        pdb_path = (root / pdb_path).resolve()
    out_png = Path(args.out_png)
    if not out_png.is_absolute():
        out_png = (root / out_png).resolve()
    ensure_dir(out_png.parent)

    coords = parse_pdb_ca_coords(pdb_path)

    # Use the same reconstruction used elsewhere in this repo to obtain y_perp.
    sys.path.insert(0, str((root / "scripts" / "pdb_bench").resolve()))
    from bench_utils import oracle_direction_reconstruct  # noqa: E402

    _, _, y_perp, ph_rms, ph_max = oracle_direction_reconstruct(coords, alphabet=str(args.alphabet))
    Z = pca_2d(y_perp)
    H = auric_type_entropy_from_yperp(y_perp, m=int(args.m))

    t = np.arange(Z.shape[0], dtype=np.float64)
    fig, ax = plt.subplots(1, 1, figsize=(6.0, 5.2), constrained_layout=True)
    sc = ax.scatter(Z[:, 0], Z[:, 1], c=t, s=12, cmap="viridis", alpha=0.85, linewidths=0.0)
    ax.set_xlabel("PC1(y_perp)")
    ax.set_ylabel("PC2(y_perp)")
    ttl = str(args.label).strip() or pdb_path.name
    ax.set_title(f"{ttl}\nH(type)@m={int(args.m)}={H:.3f}  ph_rms={ph_rms:.3f}  ph_max={ph_max:.3f}")
    ax.grid(True, alpha=0.25)
    cb = fig.colorbar(sc, ax=ax, fraction=0.046, pad=0.04)
    cb.set_label("residue index")
    fig.savefig(out_png, dpi=220)
    print(f"[ok] wrote: {out_png}")


if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        print(f"[error] {e}", file=sys.stderr)
        raise

