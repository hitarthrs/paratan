"""
Remapping between mirror invariants and local pitch

File contains the production remap between the mirror invariant distribution F(v, Lambda)
and the local gyrotropic distribution f(v, xi, B_tilde)
"""
from __future__ import annotations
import numpy as np
from numpy.typing import ArrayLike
from source_model_revamp.fbis.mirror_invariant_grid import LambdaGrid
from source_model_revamp.fbis.velocity_space_grid import PitchGrid

def _require_positive_B_tilde(B_tilde: float) -> float:
    """Return a validated positive normalized magnetic field"""
    B = float(B_tilde)
    if not np.isfinite(B) or B <= 0.0:
        raise ValueError("B_tilde must be positive and finite")
    return B

def _overlap_width(a_lower: float, a_upper: float, b_lower: np.ndarray, b_upper: np.ndarray) -> np.ndarray:
    """Overlap width between one interval and many intervals"""
    return np.maximum(0.0, np.minimum(float(a_upper), b_upper) - np.maximum(float(a_lower), b_lower))

def conservative_local_pitch_distribution_from_lambda_distribution(lambda_grid: LambdaGrid, distribution_v_lambda: ArrayLike, pitch_grid: PitchGrid, B_tilde: float, fill_value: float = 0.0) -> np.ndarray:
    """
    project cell averaged F(v, Lambda) to pitch cells

    Result is a local pitch cell average on (v, xi)
    Lambda cells whose interval is inaccessible at the supplied B_tilde contribute no local pitch measure
    """
    B = _require_positive_B_tilde(B_tilde)
    distribution = np.asarray(distribution_v_lambda, dtype=float)
    if distribution.ndim != 2 or distribution.shape[1] != lambda_grid.centers.size:
        raise ValueError("distribution_v_lambda must have shape (n_speed, n_lambda)")
    if not np.all(np.isfinite(distribution)):
        raise ValueError("distribution_v_lambda must be finite")
    lambda_faces = np.asarray(lambda_grid.faces, dtype=float)
    pitch_faces = np.asarray(pitch_grid.faces, dtype=float)
    result = np.full((distribution.shape[0], pitch_grid.centers.size), float(fill_value), dtype=float)
    lambda_lower = np.maximum(lambda_faces[:-1], 0.0)
    lambda_upper = np.minimum(lambda_faces[1:], 1.0 / B)
    valid_lambda = lambda_upper > lambda_lower
    xi_inner = np.sqrt(np.maximum(1.0 - lambda_upper * B, 0.0))
    xi_outer = np.sqrt(np.maximum(1.0 - lambda_lower * B, 0.0))
    pitch_lower = pitch_faces[:-1]
    pitch_upper = pitch_faces[1:]
    pitch_width = pitch_upper - pitch_lower
    pitch_integral = np.zeros_like(result)
    covered_width = np.zeros(pitch_grid.centers.size, dtype=float)
    for j_lambda in range(lambda_grid.centers.size):
        if not valid_lambda[j_lambda]:
            continue
        lo = float(xi_inner[j_lambda])
        hi = float(xi_outer[j_lambda])
        if hi <= lo:
            continue
        plus_overlap = _overlap_width(lo, hi, pitch_lower, pitch_upper)
        minus_overlap = _overlap_width(-hi, -lo, pitch_lower, pitch_upper)
        overlap = plus_overlap + minus_overlap
        if np.any(overlap > 0.0):
            pitch_integral += distribution[:, j_lambda][:, None] * overlap[None, :]
            covered_width += overlap
    for j_pitch in range(pitch_grid.centers.size):
        width = float(pitch_width[j_pitch])
        if width <= 0.0:
            raise ValueError("pitch widths must be positive")
        uncovered = max(width - float(covered_width[j_pitch]), 0.0)
        result[:, j_pitch] = (pitch_integral[:, j_pitch] + float(fill_value) * uncovered) / width

    return result

def conservative_lambda_distribution_from_local_pitch_distribution(pitch_grid: PitchGrid, local_distribution_v_xi: ArrayLike, lambda_grid: LambdaGrid, B_tilde: float, fill_value: float = 0.0) -> np.ndarray:
    """project local pitch cell averages back to F(v,Lambda)"""
    B = _require_positive_B_tilde(B_tilde)
    local = np.asarray(local_distribution_v_xi, dtype=float)
    if local.ndim != 2 or local.shape[1] != pitch_grid.centers.size:
        raise ValueError("local_distribution_v_xi must have shape (n_speed, n_pitch)")
    if not np.all(np.isfinite(local)):
        raise ValueError("local_distribution_v_xi must be finite")
    lambda_faces = np.asarray(lambda_grid.faces, dtype=float)
    pitch_faces = np.asarray(pitch_grid.faces, dtype=float)
    result = np.full((local.shape[0], lambda_grid.centers.size), float(fill_value), dtype=float)
    pitch_lower = pitch_faces[:-1]
    pitch_upper = pitch_faces[1:]
    for j_lambda in range(lambda_grid.centers.size):
        lambda_lower = max(float(lambda_faces[j_lambda]), 0.0)
        lambda_upper = min(float(lambda_faces[j_lambda + 1]), 1.0 / B)
        if lambda_upper <= lambda_lower:
            continue
        xi_inner = float(np.sqrt(max(1.0 - lambda_upper * B, 0.0)))
        xi_outer = float(np.sqrt(max(1.0 - lambda_lower * B, 0.0)))
        pitch_measure = 2.0 * (xi_outer - xi_inner)
        if pitch_measure <= 0.0:
            continue
        plus_overlap = _overlap_width(xi_inner, xi_outer, pitch_lower, pitch_upper)
        minus_overlap = _overlap_width(-xi_outer, -xi_inner, pitch_lower, pitch_upper)
        overlap = plus_overlap + minus_overlap
        result[:, j_lambda] = np.sum(local * overlap[None, :], axis=1) / pitch_measure

    return result