"""Cached numerical quadrature rules"""
from __future__ import annotations
from functools import lru_cache
import numpy as np

@lru_cache(maxsize=32)
def gauss_legendre_rule(order: int) -> tuple[np.ndarray, np.ndarray]:
    """Return one immutable Gauss Legendre rule"""
    count = int(order)
    if count < 1:
        raise ValueError("Gauss Legendre quadrature order must be positive")
    nodes, weights = np.polynomial.legendre.leggauss(count)
    nodes.setflags(write=False)
    weights.setflags(write=False)
    
    return nodes, weights