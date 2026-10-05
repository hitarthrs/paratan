"""Linear interpolation helpers for evaluated fusion angular data"""
from __future__ import annotations
import numpy as np
from numpy.typing import ArrayLike

def validate_linear_energy_interpolation(breakpoints: ArrayLike, laws: ArrayLike, point_count: int) -> None:
    """Require one dimensional ENDF interpolation regions using linear linear law `INT = 2`"""
    nbt = np.asarray(breakpoints, dtype=np.int64)
    interpolation_laws = np.asarray(laws, dtype=np.int64)
    if nbt.ndim != 1 or interpolation_laws.ndim != 1:
        raise ValueError("interpolation metadata must be one dimensional")
    if nbt.size != interpolation_laws.size or nbt.size == 0:
        raise ValueError("interpolation metadata must contain equal nonempty arrays")
    if np.any(np.diff(nbt) <= 0) or int(nbt[-1]) != int(point_count):
        raise ValueError("interpolation breakpoints must increase through all points")
    if np.any(interpolation_laws != 2):
        raise ValueError("only ENDF linear linear incident energy interpolation is supported")

def bracket_linear_energy(incident_energy_eV: ArrayLike, energy_grid_eV: ArrayLike) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return lower index, upper index, and linear fraction for each incident energy

    Output arrays preserve the query shape and exact grid knots use zero fraction
    """
    energy = np.asarray(incident_energy_eV, dtype=float)
    grid = np.asarray(energy_grid_eV, dtype=float)
    if grid.ndim != 1 or grid.size < 2 or np.any(np.diff(grid) <= 0.0):
        raise ValueError("energy_grid_eV must be a strictly increasing 1D grid")
    if np.any(~np.isfinite(energy)):
        raise ValueError("incident_energy_eV must be finite")
    if np.any(energy < grid[0]) or np.any(energy > grid[-1]):
        raise ValueError(f"incident energy must lie inside [{grid[0]}, {grid[-1]}] eV")
    upper = np.searchsorted(grid, energy, side="right")
    upper = np.clip(upper, 1, grid.size - 1)
    lower = upper - 1
    exact_upper = energy == grid[upper]
    lower = np.where(exact_upper, upper, lower)
    denominator = grid[upper] - grid[lower]
    fraction = np.divide(energy - grid[lower], denominator, out=np.zeros_like(energy, dtype=float), where=denominator > 0.0)

    return lower.astype(np.int64), upper.astype(np.int64), fraction

def linear_blend(lower_values: ArrayLike, upper_values: ArrayLike, fraction: ArrayLike) -> np.ndarray:
    """Evaluate `lower + fraction * (upper − lower)` with NumPy broadcasting"""
    lower = np.asarray(lower_values, dtype=float)
    upper = np.asarray(upper_values, dtype=float)
    weight = np.asarray(fraction, dtype=float)

    return lower + weight * (upper - lower)

def linear_probability_density(mu: ArrayLike, knot_mu: ArrayLike, knot_probability_density: ArrayLike) -> np.ndarray:
    """Interpolate a tabulated `p(mu)` law on a strictly increasing `−1` to `1` cosine grid"""
    query = np.asarray(mu, dtype=float)
    cosine = np.asarray(knot_mu, dtype=float)
    probability = np.asarray(knot_probability_density, dtype=float)
    if cosine.ndim != 1 or probability.ndim != 1 or cosine.size != probability.size:
        raise ValueError("tabulated angular arrays must be equal length 1D arrays")
    if cosine.size < 2 or cosine[0] != -1.0 or cosine[-1] != 1.0:
        raise ValueError("tabulated cosine grid must span minus one to one")
    if np.any(np.diff(cosine) <= 0.0):
        raise ValueError("tabulated cosine grid must be strictly increasing")
    if np.any(query < -1.0) or np.any(query > 1.0):
        raise ValueError("mu must lie inside [-1, 1]")
    
    return np.interp(query, cosine, probability)