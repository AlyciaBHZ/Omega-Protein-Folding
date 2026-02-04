from __future__ import annotations

from dataclasses import dataclass
from typing import Optional, Tuple

import numpy as np


@dataclass
class PseudoMoltenParams:
    """
    A simple "hard" control: WLC-like persistence + SAW-like excluded volume + compactness confinement.

    This is intentionally minimal and fast to generate, not a full physical model.
    """

    persistence: float = 0.7  # 0..1 correlation with previous direction
    min_dist: float = 3.5  # Å excluded-volume for non-neighbor CA
    min_sep: int = 3  # ignore clashes for |i-j| < min_sep
    conf_k: float = 1.0  # strength of soft confinement outside radius
    conf_radius_scale: float = (5.0 / 3.0) ** 0.5  # sphere radius from Rg: R ≈ Rg*sqrt(5/3)
    candidate_dirs: int = 24  # proposals per step (choose best)
    max_tries_per_step: int = 200
    max_restarts: int = 50


def _sample_persistent_dir(rng: np.random.Generator, prev: Optional[np.ndarray], persistence: float) -> np.ndarray:
    v = rng.normal(size=3)
    v /= max(1e-12, float(np.linalg.norm(v)))
    if prev is None:
        return v
    p = float(np.clip(persistence, 0.0, 0.9999))
    u = p * prev + (1.0 - p) * v
    u /= max(1e-12, float(np.linalg.norm(u)))
    return u


def _clash_ok(
    coords: np.ndarray,
    i_next: int,
    x_new: np.ndarray,
    *,
    min_dist2: float,
    min_sep: int,
) -> bool:
    # Only check against positions with index <= i_next - min_sep
    jmax = i_next - int(min_sep)
    if jmax <= 0:
        return True
    d = coords[:jmax] - x_new[None, :]
    ds2 = np.sum(d**2, axis=1)
    return bool(np.all(ds2 > min_dist2))


def generate_pseudo_molten_globule(
    rng: np.random.Generator,
    bond_lengths: np.ndarray,
    *,
    target_rg: float,
    params: PseudoMoltenParams,
) -> np.ndarray:
    """
    Generate a self-avoiding, locally persistent, compact chain with given bond lengths.

    Uses greedy best-of-K proposals under a soft confinement penalty.
    """
    L = np.asarray(bond_lengths, dtype=np.float64)
    N = int(len(L) + 1)
    if N <= 1:
        return np.zeros((N, 3), dtype=np.float64)

    R_conf = float(target_rg) * float(params.conf_radius_scale)
    min_dist2 = float(params.min_dist) ** 2

    for _ in range(int(params.max_restarts)):
        coords = np.zeros((N, 3), dtype=np.float64)
        prev_dir = None
        ok = True
        for t in range(N - 1):
            best = None
            best_e = None

            # try multiple proposals; pick best energy among those that satisfy SAW
            for _try in range(int(params.max_tries_per_step)):
                u = _sample_persistent_dir(rng, prev_dir, params.persistence)
                x_new = coords[t] + u * float(L[t])

                if not _clash_ok(coords, t + 1, x_new, min_dist2=min_dist2, min_sep=params.min_sep):
                    continue

                r = float(np.linalg.norm(x_new))
                # only penalize outside confinement radius (pseudo-molten globule)
                dr = max(0.0, r - R_conf)
                e = float(params.conf_k) * (dr * dr)
                if best_e is None or e < best_e:
                    best_e = e
                    best = (x_new, u)
                    # early stop if perfect (inside conf radius)
                    if e <= 1e-12:
                        break

            if best is None:
                ok = False
                break

            coords[t + 1] = best[0]
            prev_dir = best[1]

        if ok:
            return coords

    # If generation fails repeatedly, fall back to a simple random walk (still returns a chain)
    coords = np.zeros((N, 3), dtype=np.float64)
    for t in range(N - 1):
        v = rng.normal(size=3)
        v /= max(1e-12, float(np.linalg.norm(v)))
        coords[t + 1] = coords[t] + v * float(L[t])
    return coords

