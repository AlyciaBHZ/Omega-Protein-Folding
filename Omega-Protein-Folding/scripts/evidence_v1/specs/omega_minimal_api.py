"""
A minimal Python API skeleton for reattaching the Omega engine to this evidence pack.

This is *not* the full engine. It defines data classes and expected entry points.
"""
from dataclasses import dataclass
from typing import Optional, Dict, Any, Tuple
import numpy as np

@dataclass
class Target:
    target_id: str
    coords_native: np.ndarray  # (N,3)
    meta: Dict[str, Any] = None

@dataclass
class RunConfig:
    seed: int = 0
    steps: int = 200
    K: int = 32
    codec: str = "pair72"
    mode: str = "baseline"
    theta: Optional[float] = None
    projector: str = "global"  # "none" | "global" | "piecewise"
    hotspot_only: bool = False

def run_omega(target: Target, cfg: RunConfig) -> Dict[str, Any]:
    """
    Run Omega on a given target. Implementations should return:
      - metrics: dict with keys tm, rmse, f1, ph_max, theta, ph_ratio, etc.
      - audit: optional per-step audit table
      - state: optional final state object
    """
    raise NotImplementedError("Attach your Omega engine implementation here.")

def compute_tm_score(coords_pred: np.ndarray, coords_native: np.ndarray) -> float:
    """Placeholder. Use a standard TM-score implementation."""
    raise NotImplementedError

