from __future__ import annotations

from typing import Iterable, Optional, Tuple

import numpy as np


def _pca_uvec(Y: np.ndarray) -> np.ndarray:
    """
    Deterministic first principal component for a (N,3) cloud.

    Uses SVD on centered data. Fixes sign ambiguity by forcing the first nonzero
    component to be positive, making outputs reproducible across runs.
    """
    Y = np.asarray(Y, dtype=np.float64)
    if Y.ndim != 2 or Y.shape[1] != 3:
        raise ValueError(f"Expected (N,3) array, got {Y.shape}")
    if Y.shape[0] == 0:
        return np.array([1.0, 0.0, 0.0], dtype=np.float64)

    Yc = Y - Y.mean(axis=0, keepdims=True)
    if np.allclose(Yc, 0.0):
        return np.array([1.0, 0.0, 0.0], dtype=np.float64)

    # Vt rows are principal directions
    _, _, Vt = np.linalg.svd(Yc, full_matrices=False)
    u = np.asarray(Vt[0], dtype=np.float64)
    n = float(np.linalg.norm(u))
    if n == 0.0:
        return np.array([1.0, 0.0, 0.0], dtype=np.float64)
    u = u / n

    # deterministic sign
    for k in range(3):
        if abs(u[k]) > 1e-12:
            if u[k] < 0:
                u = -u
            break
    return u


def _threshold_bits(s: np.ndarray, threshold: str) -> np.ndarray:
    s = np.asarray(s, dtype=np.float64)
    thr = str(threshold).lower().strip()
    if thr == "median":
        t = float(np.median(s))
    elif thr == "mean":
        t = float(np.mean(s))
    elif thr in {"zero", "0", "sign"}:
        t = 0.0
    else:
        raise ValueError(f"Unknown threshold={threshold!r}; expected one of: median, mean, zero")
    return (s >= t).astype(np.uint8)


def _as_Bperp_6x3(Bperp: np.ndarray) -> np.ndarray:
    """
    Normalize Bperp to shape (6,3) for mapping Δn (6D) -> Δy_perp (3D).
    Accepts either (6,3) or (3,6).
    """
    B = np.asarray(Bperp, dtype=np.float64)
    if B.shape == (6, 3):
        return B
    if B.shape == (3, 6):
        return B.T
    raise ValueError(f"Bperp must be (6,3) or (3,6), got {B.shape}")


def rho_A_from_yperp(
    y_perp: np.ndarray,
    *,
    u: Tuple[float, float, float] = (1.0, 0.0, 0.0),
    threshold: str = "median",
    u_mode: str = "fixed",
) -> np.ndarray:
    """
    ρ_A(t) from perp-space.

    Reuses Omega's convention w0 = mean(y_perp), then binarizes along a direction u:
      s_t = u · (y_perp[t] - w0)

    If u_mode="pca", u is chosen deterministically from the input y_perp cloud.

    Thresholding (to reduce definition sensitivity):
      - threshold="median" (default): b_t = 1{s_t >= median(s)}
      - threshold="mean":             b_t = 1{s_t >= mean(s)}
      - threshold="zero":             b_t = 1{s_t >= 0}   (legacy ablation)

    Returns uint8 array of shape (N,).
    """
    y = np.asarray(y_perp, dtype=np.float64)
    if y.ndim != 2 or y.shape[1] != 3:
        raise ValueError("y_perp must have shape (N,3)")
    if y.shape[0] == 0:
        return np.zeros((0,), dtype=np.uint8)
    w0 = y.mean(axis=0, keepdims=True)
    yc = y - w0

    um = str(u_mode).lower().strip()
    if um == "pca":
        uvec = _pca_uvec(yc)
    elif um == "fixed":
        uvec = np.asarray(u, dtype=np.float64).reshape(3)
        nu = float(np.linalg.norm(uvec))
        uvec = (uvec / nu) if nu != 0.0 else np.array([1.0, 0.0, 0.0], dtype=np.float64)
    else:
        raise ValueError(f"Unknown u_mode={u_mode!r}; expected one of: fixed, pca")

    s = yc @ uvec
    return _threshold_bits(s, threshold)


def rho_B_from_npath(
    n_path: np.ndarray,
    *,
    mode: str = "parity",
    Bperp: Optional[np.ndarray] = None,
    u: Tuple[float, float, float] = (1.0, 0.0, 0.0),
    u_mode: str = "fixed",
    threshold: str = "zero",
) -> np.ndarray:
    """
    ρ_B(t) from 6D step stream.

    Given a 6D integer walk n_path[t] (shape (N,6)), define Δn_t = n_path[t+1]-n_path[t]
    for t=0..N-2, and binarize in an alphabet-agnostic, deterministic way:

      pos = count(Δn_t > 0), neg = count(Δn_t < 0)
      if pos != neg: b_t = 1{pos > neg}
      else: tie-break by sign of first nonzero component (index order); if all zero -> 0

    Returns uint8 array of shape (N-1,).
    """
    n = np.asarray(n_path)
    if n.ndim != 2 or n.shape[1] != 6:
        raise ValueError("n_path must have shape (N,6)")
    if n.shape[0] < 2:
        return np.zeros((0,), dtype=np.uint8)
    m = str(mode).lower().strip()

    # Legacy n_path majority-sign / tie-break (kept for backward-compat).
    if m in {"parity", "majority", "npath_majority"}:
        d = (n[1:] - n[:-1]).astype(np.int32)  # (N-1,6)
        out = np.zeros((d.shape[0],), dtype=np.uint8)
        for i in range(d.shape[0]):
            step = d[i]
            pos = int(np.sum(step > 0))
            neg = int(np.sum(step < 0))
            if pos != neg:
                out[i] = 1 if pos > neg else 0
                continue
            # tie-breaker: first nonzero component sign
            b = 0
            for v in step:
                if v > 0:
                    b = 1
                    break
                if v < 0:
                    b = 0
                    break
            out[i] = b
        return out

    if Bperp is None:
        raise ValueError(f"Bperp is required for mode={mode!r}")
    B = _as_Bperp_6x3(Bperp)

    if m in {"vel", "yperp_vel"}:
        dn = (n[1:] - n[:-1]).astype(np.int64)  # (N-1,6)
        if dn.shape[0] == 0:
            return np.zeros((0,), dtype=np.uint8)
        dy = dn @ B  # (N-1,3)

        um = str(u_mode).lower().strip()
        if um == "pca":
            uvec = _pca_uvec(dy)
        elif um == "fixed":
            uvec = np.asarray(u, dtype=np.float64).reshape(3)
            nu = float(np.linalg.norm(uvec))
            uvec = (uvec / nu) if nu != 0.0 else np.array([1.0, 0.0, 0.0], dtype=np.float64)
        else:
            raise ValueError(f"Unknown u_mode={u_mode!r}; expected one of: fixed, pca")

        s = dy @ uvec
        return _threshold_bits(s, threshold)

    if m in {"pos", "yperp_pos"}:
        # Reconstruct y_perp positions from n_path and apply rho_A rule.
        y = n @ B  # (N,3)
        return rho_A_from_yperp(y, u=u, threshold=threshold, u_mode=u_mode)

    raise ValueError(f"Unknown mode={mode!r}; expected one of: parity, vel, pos")


def rho_B_from_yperp_velocity(
    y_perp: np.ndarray,
    *,
    u: Tuple[float, float, float] = (1.0, 0.0, 0.0),
    u_mode: str = "fixed",
    threshold: str = "median",
) -> np.ndarray:
    """
    Geometrically bound rho_B: perp-space velocity bitstream.

    Define Δy_perp[t] = y_perp[t+1] - y_perp[t], and binarize:
      s_t = u · Δy_perp[t]
      b_t = 1{s_t >= threshold(s)}

    Returns uint8 array of shape (N-1,).
    """
    y = np.asarray(y_perp, dtype=np.float64)
    if y.ndim != 2 or y.shape[1] != 3:
        raise ValueError("y_perp must have shape (N,3)")
    if y.shape[0] < 2:
        return np.zeros((0,), dtype=np.uint8)
    dy = y[1:] - y[:-1]

    um = str(u_mode).lower().strip()
    if um == "pca":
        uvec = _pca_uvec(dy)
    elif um == "fixed":
        uvec = np.asarray(u, dtype=np.float64).reshape(3)
        nu = float(np.linalg.norm(uvec))
        uvec = (uvec / nu) if nu != 0.0 else np.array([1.0, 0.0, 0.0], dtype=np.float64)
    else:
        raise ValueError(f"Unknown u_mode={u_mode!r}; expected one of: fixed, pca")

    s = dy @ uvec
    return _threshold_bits(s, threshold)


def rho_B_from_yperp_position(
    y_perp: np.ndarray,
    *,
    u: Tuple[float, float, float] = (1.0, 0.0, 0.0),
    u_mode: str = "fixed",
    threshold: str = "median",
) -> np.ndarray:
    """
    Geometrically bound rho_B: perp-space position bitstream (rho_A-equivalent).

    This is intentionally identical to rho_A on the same y_perp/u/threshold.
    """
    return rho_A_from_yperp(y_perp, u=u, u_mode=u_mode, threshold=threshold)


def to_uint8_bits(xs: Iterable[int]) -> np.ndarray:
    """Normalize an iterable of bits into a uint8 numpy array (values 0/1)."""
    arr = np.fromiter((int(x) & 1 for x in xs), dtype=np.uint8)
    return arr

