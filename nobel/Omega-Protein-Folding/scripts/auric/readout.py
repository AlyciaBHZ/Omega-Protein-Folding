from __future__ import annotations

from typing import Iterable, Tuple

import numpy as np


def rho_A_from_yperp(y_perp: np.ndarray, *, u: Tuple[float, float, float] = (1.0, 0.0, 0.0)) -> np.ndarray:
    """
    ρ_A(t) from perp-space.

    Reuses Omega's convention w0 = mean(y_perp), then binarizes along a fixed direction u:
      b_t = 1{ u · (y_perp[t] - w0) >= 0 }.

    Returns uint8 array of shape (N,).
    """
    y = np.asarray(y_perp, dtype=np.float64)
    if y.ndim != 2 or y.shape[1] != 3:
        raise ValueError("y_perp must have shape (N,3)")
    w0 = y.mean(axis=0, keepdims=True)
    uvec = np.asarray(u, dtype=np.float64).reshape(3)
    dots = (y - w0) @ uvec
    return (dots >= 0.0).astype(np.uint8)


def rho_B_from_npath(n_path: np.ndarray) -> np.ndarray:
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


def to_uint8_bits(xs: Iterable[int]) -> np.ndarray:
    """Normalize an iterable of bits into a uint8 numpy array (values 0/1)."""
    arr = np.fromiter((int(x) & 1 for x in xs), dtype=np.uint8)
    return arr

