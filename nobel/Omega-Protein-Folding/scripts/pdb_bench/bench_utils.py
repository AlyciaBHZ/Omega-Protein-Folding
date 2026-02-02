from __future__ import annotations

import math
from pathlib import Path
from typing import List, Tuple

import numpy as np
from Bio.PDB.MMCIFParser import MMCIFParser

PHI = (1.0 + 5.0**0.5) / 2.0


def ensure_dir(p: Path) -> None:
    p.mkdir(parents=True, exist_ok=True)


def icosa_B(phi: float) -> np.ndarray:
    """
    Build a 3x6 icosahedral projection basis B (columns correspond to +e_i).
    """
    cols = [
        np.array([0.0, 1.0, phi]),
        np.array([0.0, 1.0, -phi]),
        np.array([1.0, phi, 0.0]),
        np.array([1.0, -phi, 0.0]),
        np.array([phi, 0.0, 1.0]),
        np.array([phi, 0.0, -1.0]),
    ]
    cols = [c / np.linalg.norm(c) for c in cols]
    return np.stack(cols, axis=1)  # (3,6)


def alphabet_steps(name: str) -> np.ndarray:
    """
    Return (K,6) int steps for axis12/pair72/triple232.
    """
    if name not in {"axis12", "pair72", "triple232"}:
        raise ValueError(f"Unknown alphabet: {name}")

    steps: List[np.ndarray] = []
    # axis steps
    for i in range(6):
        e = np.zeros(6, dtype=np.int32)
        e[i] = 1
        steps.append(e.copy())
        steps.append((-e).copy())

    if name in {"pair72", "triple232"}:
        for i in range(6):
            for j in range(i + 1, 6):
                for si in (-1, 1):
                    for sj in (-1, 1):
                        v = np.zeros(6, dtype=np.int32)
                        v[i] = si
                        v[j] = sj
                        steps.append(v)

    if name == "triple232":
        for i in range(6):
            for j in range(i + 1, 6):
                for k in range(j + 1, 6):
                    for si in (-1, 1):
                        for sj in (-1, 1):
                            for sk in (-1, 1):
                                v = np.zeros(6, dtype=np.int32)
                                v[i] = si
                                v[j] = sj
                                v[k] = sk
                                steps.append(v)

    out = np.stack(steps, axis=0)
    expected = {"axis12": 12, "pair72": 72, "triple232": 232}[name]
    if out.shape[0] != expected:
        raise RuntimeError(f"Alphabet {name} size mismatch: got {out.shape[0]}, expected {expected}")
    return out


def load_ca_coords_longest_chain(cif_path: Path) -> Tuple[str, np.ndarray]:
    """
    Return (chain_id, coords) where coords is (N,3) CA trace for the longest chain.
    """
    parser = MMCIFParser(QUIET=True)
    structure_id = cif_path.stem
    s = parser.get_structure(structure_id, str(cif_path))
    model = next(iter(s.get_models()))

    best_chain = None
    best_coords = None
    for chain in model:
        coords = []
        for res in chain:
            if "CA" in res:
                coords.append(res["CA"].get_coord())
        if len(coords) >= 2:
            coords_arr = np.asarray(coords, dtype=np.float64)
            if best_coords is None or len(coords_arr) > len(best_coords):
                best_coords = coords_arr
                best_chain = chain.id

    if best_coords is None or best_chain is None:
        raise ValueError("No CA trace found")
    return str(best_chain), best_coords


def kabsch(P: np.ndarray, Q: np.ndarray) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Return aligned P_aligned, rotation R, translation t such that:
      P_aligned ≈ P @ R + t  matches Q (least-squares).
    """
    assert P.shape == Q.shape and P.shape[1] == 3
    Pc = P - P.mean(axis=0, keepdims=True)
    Qc = Q - Q.mean(axis=0, keepdims=True)
    C = Pc.T @ Qc
    V, _, Wt = np.linalg.svd(C)
    d = np.sign(np.linalg.det(V @ Wt))
    D = np.diag([1.0, 1.0, d])
    R = V @ D @ Wt
    t = Q.mean(axis=0) - P.mean(axis=0) @ R
    P_aligned = P @ R + t
    return P_aligned, R, t


def tm_d0(L: int) -> float:
    # Standard TM-score d0 definition (Zhang & Skolnick).
    if L <= 0:
        return 1.0
    d0 = 1.24 * (L - 15) ** (1.0 / 3.0) - 1.8
    return float(max(d0, 0.5))


def tm_score(P: np.ndarray, Q: np.ndarray) -> float:
    """
    TM-score after optimal Kabsch alignment (single-pass).
    """
    P_aligned, _, _ = kabsch(P, Q)
    d = np.linalg.norm(P_aligned - Q, axis=1)
    d0 = tm_d0(len(d))
    return float(np.mean(1.0 / (1.0 + (d / d0) ** 2)))


def phason_stats_from_yperp(y_perp: np.ndarray) -> Tuple[float, float]:
    """
    Given y_perp (N,3), compute (ph_rms, ph_max) using w0 = mean(y_perp).
    """
    w0 = y_perp.mean(axis=0, keepdims=True)
    r = np.linalg.norm(y_perp - w0, axis=1)
    ph_rms = float(np.sqrt(np.mean(r**2)))
    ph_max = float(np.max(r))
    return ph_rms, ph_max

def phason_stats_piecewise(y_perp: np.ndarray, *, segment_len: int) -> Tuple[float, float]:
    """
    Piecewise proxy: choose w0 per segment as mean(y_perp in segment),
    then compute global ph_rms/ph_max over centered residuals.
    """
    N = y_perp.shape[0]
    seg = int(max(1, segment_len))
    rs: List[np.ndarray] = []
    for s in range(0, N, seg):
        y = y_perp[s : s + seg]
        w0 = y.mean(axis=0, keepdims=True)
        rs.append(np.linalg.norm(y - w0, axis=1))
    r = np.concatenate(rs, axis=0) if rs else np.zeros((0,), dtype=np.float64)
    ph_rms = float(np.sqrt(np.mean(r**2))) if len(r) else float("nan")
    ph_max = float(np.max(r)) if len(r) else float("nan")
    return ph_rms, ph_max

def w0_piecewise_jumps(y_perp: np.ndarray, *, segment_len: int) -> Tuple[float, float]:
    """
    Piecewise w0 trajectory: w0_s = mean(y_perp in segment s).
    Return (mean_jump, max_jump) over ||w0_{s+1}-w0_s||.
    """
    y = np.asarray(y_perp, dtype=np.float64)
    N = y.shape[0]
    seg = int(max(1, segment_len))
    w0s = []
    for s in range(0, N, seg):
        w0s.append(y[s : s + seg].mean(axis=0))
    if len(w0s) <= 1:
        return 0.0, 0.0
    w0s = np.stack(w0s, axis=0)
    d = np.linalg.norm(w0s[1:] - w0s[:-1], axis=1)
    return float(np.mean(d)), float(np.max(d))

def phason_stats_linear_drift(y_perp: np.ndarray) -> Tuple[float, float, float]:
    """
    Fit a global linear drift w0(i) = b + a*i in perp-space (least squares per component),
    then compute residual ph_rms/ph_max. Also return drift_speed = ||a||_2.
    """
    y = np.asarray(y_perp, dtype=np.float64)
    N = y.shape[0]
    if N == 0:
        return float("nan"), float("nan"), float("nan")

    i = np.arange(N, dtype=np.float64)
    X = np.stack([np.ones_like(i), i], axis=1)  # (N,2)
    # Solve for each dimension independently: beta=(b,a)
    # beta shape (2,3)
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    b = beta[0]  # (3,)
    a = beta[1]  # (3,)
    yhat = b[None, :] + i[:, None] * a[None, :]
    r = np.linalg.norm(y - yhat, axis=1)
    ph_rms = float(np.sqrt(np.mean(r**2)))
    ph_max = float(np.max(r))
    drift_speed = float(np.linalg.norm(a))
    return ph_rms, ph_max, drift_speed


def oracle_direction_reconstruct(
    coords_native: np.ndarray,
    *,
    alphabet: str,
    phi: float = PHI,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, float, float]:
    """
    Oracle direction quantization (per-bond). Returns:
      - coords_pred (N,3) before Kabsch
      - n_path (N,6) int walk
      - y_perp (N,3)
      - ph_rms, ph_max (computed from y_perp)

    Important: we scale chosen step direction to the observed bond length, so
    multi-axis steps provide angular expressivity without forcing longer bonds.
    """
    N = coords_native.shape[0]
    d = coords_native[1:] - coords_native[:-1]  # (N-1,3)
    bond_lengths = np.linalg.norm(d, axis=1)
    median_len = float(np.median(bond_lengths))

    Bpar = icosa_B(phi) * median_len
    Bperp = icosa_B(-1.0 / phi) * median_len

    steps = alphabet_steps(alphabet)  # (K,6)
    proj = (Bpar @ steps.T).T  # (K,3)
    proj_norm = np.linalg.norm(proj, axis=1)
    proj_dir = proj / np.clip(proj_norm[:, None], 1e-12, None)

    coords_pred = np.zeros((N, 3), dtype=np.float64)
    n_path = np.zeros((N, 6), dtype=np.int32)

    for t in range(N - 1):
        dt = d[t]
        dt_len = float(np.linalg.norm(dt))
        if dt_len < 1e-8:
            k = 0
            step_vec = proj_dir[k] * 0.0
        else:
            dt_dir = dt / dt_len
            cos = proj_dir @ dt_dir
            k = int(np.argmax(cos))
            step_vec = proj_dir[k] * dt_len
        step = steps[k]
        coords_pred[t + 1] = coords_pred[t] + step_vec
        n_path[t + 1] = n_path[t] + step

    y_perp = (Bperp @ n_path.T).T  # (N,3)
    ph_rms, ph_max = phason_stats_from_yperp(y_perp)
    return coords_pred, n_path, y_perp, ph_rms, ph_max


def random_chain_from_lengths(rng: np.random.Generator, bond_lengths: np.ndarray) -> np.ndarray:
    """
    Generate a random 3D chain with given bond lengths.
    """
    N = len(bond_lengths) + 1
    coords = np.zeros((N, 3), dtype=np.float64)
    for t, L in enumerate(bond_lengths):
        v = rng.normal(size=3)
        v /= max(1e-12, float(np.linalg.norm(v)))
        coords[t + 1] = coords[t] + v * float(L)
    return coords


def perturb_chain_directions(
    rng: np.random.Generator,
    coords_native: np.ndarray,
    *,
    noise: float = 0.15,
) -> np.ndarray:
    """
    Perturb bond directions while preserving per-bond lengths.
    noise: std of additive Gaussian noise in direction space.
    """
    d = coords_native[1:] - coords_native[:-1]
    L = np.linalg.norm(d, axis=1)
    dirs = d / np.clip(L[:, None], 1e-12, None)
    dirs2 = dirs + rng.normal(scale=noise, size=dirs.shape)
    dirs2 = dirs2 / np.clip(np.linalg.norm(dirs2, axis=1)[:, None], 1e-12, None)
    d2 = dirs2 * L[:, None]

    coords = np.zeros_like(coords_native)
    for t in range(len(d2)):
        coords[t + 1] = coords[t] + d2[t]
    return coords


def shuffle_bond_directions(
    rng: np.random.Generator,
    coords_native: np.ndarray,
) -> np.ndarray:
    """
    Shuffle the sequence of 3D bond directions while preserving per-bond lengths.

    This is a correlation-destroying null that preserves the marginal distributions
    of step directions and the original length schedule.
    """
    x = np.asarray(coords_native, dtype=np.float64)
    d = x[1:] - x[:-1]
    L = np.linalg.norm(d, axis=1)
    dirs = d / np.clip(L[:, None], 1e-12, None)
    perm = rng.permutation(dirs.shape[0])
    dirs2 = dirs[perm]
    d2 = dirs2 * L[:, None]
    out = np.zeros_like(x)
    for t in range(d2.shape[0]):
        out[t + 1] = out[t] + d2[t]
    return out


def blockshuffle_bond_directions(
    rng: np.random.Generator,
    coords_native: np.ndarray,
    *,
    block_len: int,
) -> np.ndarray:
    """
    Block-shuffle the bond-direction stream to probe correlation length.

    Split the direction sequence into contiguous blocks of length block_len,
    shuffle the blocks, keep intra-block order; preserve the original length schedule.
    """
    x = np.asarray(coords_native, dtype=np.float64)
    d = x[1:] - x[:-1]
    L = np.linalg.norm(d, axis=1)
    dirs = d / np.clip(L[:, None], 1e-12, None)
    n = dirs.shape[0]
    k = int(max(1, block_len))
    blocks = [dirs[i : i + k] for i in range(0, n, k)]
    perm = rng.permutation(len(blocks))
    dirs2 = np.concatenate([blocks[i] for i in perm], axis=0)
    d2 = dirs2 * L[:, None]
    out = np.zeros_like(x)
    for t in range(d2.shape[0]):
        out[t + 1] = out[t] + d2[t]
    return out


def cliffs_delta_one_vs_many(x: float, ys: np.ndarray) -> float:
    """
    Cliff's delta comparing a single x to a sample ys.
    """
    ys = np.asarray(ys, dtype=np.float64)
    gt = float(np.sum(x > ys))
    lt = float(np.sum(x < ys))
    n = float(len(ys)) if len(ys) else 1.0
    return (gt - lt) / n


def radius_of_gyration(coords: np.ndarray) -> float:
    """
    Rg for a CA trace.
    """
    x = np.asarray(coords, dtype=np.float64)
    c = x.mean(axis=0, keepdims=True)
    return float(np.sqrt(np.mean(np.sum((x - c) ** 2, axis=1))))


def contact_density(coords: np.ndarray, *, cutoff: float = 8.0, min_sep: int = 3) -> float:
    """
    Simple CA contact density: (# pairs within cutoff, excluding |i-j|<min_sep) / N.
    O(N^2) but fine for N<=350.
    """
    x = np.asarray(coords, dtype=np.float64)
    N = x.shape[0]
    if N <= 1:
        return float("nan")
    c2 = float(cutoff) ** 2
    cnt = 0
    for i in range(N):
        xi = x[i]
        j0 = i + int(min_sep)
        if j0 >= N:
            continue
        d = x[j0:] - xi[None, :]
        ds2 = np.sum(d**2, axis=1)
        cnt += int(np.sum(ds2 <= c2))
    return float(cnt) / float(N)

