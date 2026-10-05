"""Small numerical helpers shared by modal FBIS modules"""
from __future__ import annotations
import numpy as np

_EPS = np.finfo(float).eps

def _positive_int(name: str, value: int, minimum: int = 1) -> int:
    """Return an integer that meets the requested inclusive minimum"""
    ivalue = int(value)
    if ivalue < minimum:
        raise ValueError(f"{name} must be >= {minimum}")
    
    return ivalue

def _safe_interp_rows(x_old: np.ndarray, y_rows: np.ndarray, x_new: np.ndarray, *, left: float = 0.0, right: float = 0.0) -> np.ndarray:
    """Interpolate every row from sorted source coordinates with explicit endpoint fill values"""
    result = np.empty((y_rows.shape[0], x_new.size), dtype=float)
    order = np.argsort(x_old)
    xs = x_old[order]
    for i in range(y_rows.shape[0]):
        ys = y_rows[i, order]
        result[i] = np.interp(x_new, xs, ys, left=left, right=right)

    return result