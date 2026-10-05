"""Axisymmetric physical position sampling for correlated neutron events"""
from __future__ import annotations
import numpy as np
from source_model_revamp.radial_profiles import KOTELNIKOV_PARABOLIC_FLUX_K1, canonical_radial_profile_model, radial_probability_cdf, radial_probability_inverse_cdf

def sample_axisymmetric_cell_positions(axial_cell_indices: np.ndarray, z_edges_m: np.ndarray, radial_inner_radius_m_by_z: np.ndarray, radial_outer_radius_m_by_z: np.ndarray, rng: np.random.Generator, *, radial_outer_radius_m_by_edge: np.ndarray | None = None, radial_profile_model: str = KOTELNIKOV_PARABOLIC_FLUX_K1, source_connected_radial_rho_limit_edges: np.ndarray | None = None) -> np.ndarray:
    """
    Sample event positions inside selected axial cells from the configured radial source profile
    
    Axial position is uniform inside each cell
    The radial CDF is truncated between the inner radius and the source connected `rho` limit and azimuth is uniform on `[0, 2 pi)`
    """
    axial = np.asarray(axial_cell_indices, dtype=np.int64)
    z_edges = np.asarray(z_edges_m, dtype=float)
    inner = np.asarray(radial_inner_radius_m_by_z, dtype=float)
    outer = np.asarray(radial_outer_radius_m_by_z, dtype=float)
    if axial.ndim != 1:
        raise ValueError("axial_cell_indices must be a 1D array")
    if z_edges.ndim != 1 or z_edges.size < 2 or np.any(np.diff(z_edges) <= 0.0):
        raise ValueError("z_edges_m must be strictly increasing")
    axial_count = z_edges.size - 1
    if inner.shape != (axial_count,) or outer.shape != (axial_count,):
        raise ValueError("radial profiles must have one value per axial cell")
    if np.any(~np.isfinite(inner)) or np.any(~np.isfinite(outer)):
        raise ValueError("radial profiles must be finite")
    if np.any(inner < 0.0) or np.any(outer <= inner):
        raise ValueError("radial bounds must satisfy zero <= inner < outer")
    if radial_outer_radius_m_by_edge is None:
        outer_edges = None
    else:
        outer_edges = np.asarray(radial_outer_radius_m_by_edge, dtype=float)
        if outer_edges.shape != z_edges.shape or np.any(~np.isfinite(outer_edges)) or np.any(outer_edges <= 0.0):
            raise ValueError("radial outer edge profile must match the axial edges and be positive")
    if np.any(axial < 0) or np.any(axial >= axial_count):
        raise IndexError("axial cell index is outside the source grid")
    model = canonical_radial_profile_model(radial_profile_model)
    if source_connected_radial_rho_limit_edges is None:
        rho_limit_edges = np.ones(z_edges.size, dtype=float)
    else:
        rho_limit_edges = np.asarray(source_connected_radial_rho_limit_edges, dtype=float)
    if rho_limit_edges.shape != z_edges.shape:
        raise ValueError("source connected radial rho support must match the axial edges")
    if np.any(~np.isfinite(rho_limit_edges)) or np.any(rho_limit_edges < 0.0) or np.any(rho_limit_edges > 1.0):
        raise ValueError("source connected radial rho support must lie inside zero to one")

    axial_fraction = rng.random(axial.size)
    z = z_edges[axial] + axial_fraction * (z_edges[axial + 1] - z_edges[axial])
    rho_max = rho_limit_edges[axial] + axial_fraction * (rho_limit_edges[axial + 1] - rho_limit_edges[axial])
    sampled_outer = outer[axial] if outer_edges is None else outer_edges[axial] + axial_fraction * (outer_edges[axial + 1] - outer_edges[axial])
    rho_inner = inner[axial] / sampled_outer
    if np.any(rho_inner < 0.0) or np.any(rho_inner >= rho_max):
        raise ValueError("radial source support is empty inside at least one selected axial cell")
    cdf_inner = radial_probability_cdf(rho_inner, model)
    cdf_outer = radial_probability_cdf(rho_max, model)
    # Truncating in CDF space preserves the configured radial law on the connected support
    radial_probability = cdf_inner + rng.random(axial.size) * (cdf_outer - cdf_inner)
    rho = radial_probability_inverse_cdf(radial_probability, model)
    radius = sampled_outer * rho
    azimuth = rng.uniform(0.0, 2.0 * np.pi, axial.size)

    return np.column_stack((radius * np.cos(azimuth), radius * np.sin(azimuth), z))
