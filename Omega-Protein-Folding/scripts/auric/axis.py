from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable, List, Tuple

import numpy as np

from fold_m import fold_sliding_windows_bits
from metrics import metrics_for_stream
from readout import rho_A_from_yperp


def _normalize_u(u: np.ndarray) -> np.ndarray:
    u = np.asarray(u, dtype=np.float64).reshape(3)
    nu = float(np.linalg.norm(u))
    if nu == 0.0 or not np.isfinite(nu):
        u = np.array([1.0, 0.0, 0.0], dtype=np.float64)
    else:
        u = u / nu

    # Deterministic sign: first non-trivial component is positive.
    for k in range(3):
        if abs(u[k]) > 1e-12:
            if u[k] < 0:
                u = -u
            break
    return u


def shared_pca_axis(y_perp: np.ndarray) -> np.ndarray:
    """
    Shared-PCA axis selection used across native + decoys.

    Protocol rule: the PCA axis is computed from the NATIVE y_perp only,
    and then reused for all decoys (prevents decoy leakage).
    """
    Y = np.asarray(y_perp, dtype=np.float64)
    if Y.ndim != 2 or Y.shape[1] != 3 or Y.shape[0] == 0:
        return np.array([1.0, 0.0, 0.0], dtype=np.float64)

    Yc = Y - Y.mean(axis=0, keepdims=True)
    if np.allclose(Yc, 0.0):
        return np.array([1.0, 0.0, 0.0], dtype=np.float64)
    _, _, Vt = np.linalg.svd(Yc, full_matrices=False)
    return _normalize_u(np.asarray(Vt[0], dtype=np.float64))


@dataclass(frozen=True)
class AxisCandidate:
    name: str
    u: np.ndarray


def _candidate_axes(y_perp: np.ndarray) -> List[AxisCandidate]:
    """
    Candidate u vectors for native-only axis selection.

    This is intentionally small and deterministic to keep the ablation protocol stable.
    """
    u_pca = shared_pca_axis(y_perp)
    u_mean_step = np.diff(np.asarray(y_perp, dtype=np.float64), axis=0).mean(axis=0) if len(y_perp) >= 2 else u_pca
    cands = [
        AxisCandidate("pca", _normalize_u(u_pca)),
        AxisCandidate("mean_step", _normalize_u(u_mean_step)),
        AxisCandidate("ex", np.array([1.0, 0.0, 0.0], dtype=np.float64)),
        AxisCandidate("ey", np.array([0.0, 1.0, 0.0], dtype=np.float64)),
        AxisCandidate("ez", np.array([0.0, 0.0, 1.0], dtype=np.float64)),
    ]
    return cands


def best_native_axis(
    y_perp_native: np.ndarray,
    *,
    m: int = 8,
    threshold: str = "median",
    metric: str = "type_entropy",
) -> np.ndarray:
    """
    Strict ablation: choose u using ONLY the native structure.

    We pick the candidate axis that makes the native rhoA stream maximally "ordered"
    under the chosen metric (default: minimize type entropy at window length m).

    This is protocol-safe (no decoy leakage), and gives a reviewer-facing robustness check:
    even with native-only axis selection, native-vs-decoy separation should persist.
    """
    if int(m) <= 0:
        raise ValueError("m must be positive")
    metric = str(metric).strip()
    if metric not in {"type_entropy", "smb_rate_hat"}:
        raise ValueError("metric must be one of: type_entropy, smb_rate_hat")

    Y = np.asarray(y_perp_native, dtype=np.float64)
    if Y.ndim != 2 or Y.shape[1] != 3 or Y.shape[0] == 0:
        return np.array([1.0, 0.0, 0.0], dtype=np.float64)

    best_u = np.array([1.0, 0.0, 0.0], dtype=np.float64)
    best_val = float("inf")

    for cand in _candidate_axes(Y):
        bits_A = rho_A_from_yperp(Y, u=tuple(cand.u.tolist()), u_mode="fixed", threshold=str(threshold))
        folded = fold_sliding_windows_bits(bits_A, m=int(m))
        met = metrics_for_stream(bits_A, folded)
        v = float(met.get(metric, float("nan")))
        if not np.isfinite(v):
            continue
        if v < best_val - 1e-12:
            best_val = v
            best_u = cand.u

    return _normalize_u(best_u)

