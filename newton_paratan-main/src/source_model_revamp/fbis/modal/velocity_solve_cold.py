"""Cold ion Egedal modal velocity solve, Eq 14/15 reference path"""
from __future__ import annotations
import numpy as np
from numpy.typing import ArrayLike
from source_model_revamp.fbis.velocity_space_grid import SpeedGrid

def cold_ion_modal_distribution_eq14(*, speed_grid: SpeedGrid, source_coefficients_j: ArrayLike, source_speed_m_s: float, eigenvalues: ArrayLike, spitzer_slowing_down_time_s: float, critical_velocity_m_s: float, beta_m: float) -> np.ndarray:
    """Egedal Eq 14/15 cold ion modal velocity solution

    source_coefficients_j and eigenvalues have shape (n_mode,)
    The returned distribution has shape (n_mode, n_speed) and vanishes above the monoenergetic source speed
    """
    v = np.asarray(speed_grid.centers_m_s, dtype=float)
    S = np.asarray(source_coefficients_j, dtype=float)
    lambdas = np.asarray(eigenvalues, dtype=float)
    if S.shape != lambdas.shape:
        raise ValueError("source_coefficients_j and eigenvalues must have matching shapes")
    v0 = float(source_speed_m_s)
    tau_s = float(spitzer_slowing_down_time_s)
    vc = float(critical_velocity_m_s)
    beta = float(beta_m)
    if not all(np.isfinite(x) and x > 0.0 for x in (v0, tau_s, vc, beta)):
        raise ValueError("source speed, tau_s, critical velocity, and beta_m must be positive finite values")
    out = np.zeros((S.size, v.size), dtype=float)
    active = v <= v0
    if not np.any(active):
        return out
    va = v[active]
    common = ((v0**3 + vc**3) / (va**3 + vc**3)) * (va**3 / v0**3)
    common = np.maximum(common, np.finfo(float).tiny)
    for j in range(S.size):
        exponent = beta * float(lambdas[j]) / 3.0
        u = common**exponent
        out[j, active] = tau_s * float(S[j]) * u / (va**3 + vc**3)

    return out

def _source_cell_weights(speed_grid: SpeedGrid, source_speed_m_s: float) -> np.ndarray:
    """Return a one cell representation of a monoenergetic source on the speed grid"""
    weights = np.zeros(speed_grid.centers_m_s.size, dtype=float)
    v0 = float(source_speed_m_s)
    faces = np.asarray(speed_grid.faces_m_s, dtype=float)
    if v0 <= faces[0]:
        weights[0] = 1.0
    elif v0 >= faces[-1]:
        weights[-1] = 1.0
    else:
        idx = int(np.searchsorted(faces, v0, side="right") - 1)
        idx = min(max(idx, 0), weights.size - 1)
        weights[idx] = 1.0
        
    return weights