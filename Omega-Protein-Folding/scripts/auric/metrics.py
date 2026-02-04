from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Iterable, List, Tuple

import numpy as np


def shannon_entropy_from_counts(counts: np.ndarray) -> float:
    """Shannon entropy in bits (log2)."""
    c = np.asarray(counts, dtype=np.float64)
    total = float(np.sum(c))
    if total <= 0.0:
        return float("nan")
    p = c / total
    p = p[p > 0.0]
    return float(-np.sum(p * np.log2(p)))


def type_entropy(types: np.ndarray) -> float:
    """Entropy (bits) of a discrete type stream."""
    x = np.asarray(types).reshape(-1)
    if x.size == 0:
        return float("nan")
    _, counts = np.unique(x, return_counts=True)
    return shannon_entropy_from_counts(counts)


def type_support(types: np.ndarray) -> int:
    """Number of observed distinct types."""
    x = np.asarray(types).reshape(-1)
    if x.size == 0:
        return 0
    return int(np.unique(x).size)


@dataclass(frozen=True)
class RunStats:
    run1_mean: float
    run1_max: int
    run1_count: int


def runlength_1(bits: np.ndarray) -> RunStats:
    """Run-length statistics for consecutive 1s in a binary stream."""
    b = np.asarray(bits, dtype=np.uint8).reshape(-1)
    n = int(b.size)
    runs: List[int] = []
    cur = 0
    for i in range(n):
        if int(b[i]) == 1:
            cur += 1
        else:
            if cur > 0:
                runs.append(cur)
                cur = 0
    if cur > 0:
        runs.append(cur)
    if not runs:
        return RunStats(run1_mean=0.0, run1_max=0, run1_count=0)
    return RunStats(run1_mean=float(np.mean(runs)), run1_max=int(np.max(runs)), run1_count=len(runs))


def _pack_blocks_as_uint64(bits: np.ndarray, k: int) -> np.ndarray:
    """
    Pack all length-k blocks into uint64 codes using low-to-high within each block.
    Returns array length (n-k+1).
    """
    b = np.asarray(bits, dtype=np.uint8).reshape(-1)
    n = int(b.size)
    if k <= 0 or k > 63:
        raise ValueError("k must be in [1,63]")
    if n < k:
        return np.zeros((0,), dtype=np.uint64)
    code = np.uint64(0)
    for i in range(k):
        code |= np.uint64(int(b[i]) & 1) << np.uint64(i)
    out = np.empty((n - k + 1,), dtype=np.uint64)
    out[0] = code
    mask = (np.uint64(1) << np.uint64(k)) - np.uint64(1)
    for t in range(1, n - k + 1):
        new_bit = np.uint64(int(b[t + k - 1]) & 1)
        code = ((code >> np.uint64(1)) | (new_bit << np.uint64(k - 1))) & mask
        out[t] = code
    return out


def block_entropy(bits: np.ndarray, k: int) -> float:
    """Block entropy H_k (bits) for length-k blocks."""
    codes = _pack_blocks_as_uint64(bits, k)
    if codes.size == 0:
        return float("nan")
    _, counts = np.unique(codes, return_counts=True)
    return shannon_entropy_from_counts(counts)


def entropy_rate_block_slope(bits: np.ndarray, *, k_min: int = 4, k_max: int = 8) -> float:
    """
    SMB-style practical entropy-rate proxy via block entropy slope.

    Compute H_k for k in [k_min..k_max] (bits), then fit a line H_k ≈ a + h*k.
    Return h as smb_rate_hat.
    """
    b = np.asarray(bits, dtype=np.uint8).reshape(-1)
    if b.size < (k_max + 1):
        return float("nan")
    ks = []
    Hs = []
    for k in range(int(k_min), int(k_max) + 1):
        hk = block_entropy(b, k)
        if not np.isfinite(hk):
            continue
        ks.append(float(k))
        Hs.append(float(hk))
    if len(ks) < 2:
        return float("nan")
    x = np.asarray(ks, dtype=np.float64)
    y = np.asarray(Hs, dtype=np.float64)
    A = np.stack([np.ones_like(x), x], axis=1)
    beta, *_ = np.linalg.lstsq(A, y, rcond=None)
    h = float(beta[1])
    return h


def metrics_for_stream(bits: np.ndarray, folded_types: np.ndarray) -> Dict[str, float]:
    """
    Convenience bundle:
    - run-length stats on raw bits
    - entropy/support on folded type stream
    - entropy-rate proxy on raw bits
    """
    rs = runlength_1(bits)
    return {
        "type_entropy": type_entropy(folded_types),
        "type_support": float(type_support(folded_types)),
        "run1_mean": float(rs.run1_mean),
        "run1_max": float(rs.run1_max),
        "run1_count": float(rs.run1_count),
        "smb_rate_hat": float(entropy_rate_block_slope(bits)),
    }

