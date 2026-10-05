"""Small numerical helpers shared by integration stage modules"""
from __future__ import annotations
import numpy as np
from numpy.typing import ArrayLike
from source_model_revamp.fbis.species import IonSpecies
from source_model_revamp.fbis.beam_source_definition import beam_speed_from_energy_m_s
from source_model_revamp.fbis.velocity_space_grid import SpeedGrid, _distribution_nonnegativity_metrics, source_aligned_stretched_speed_grid
from source_model_revamp.integration.config.root import SourceModelRunConfig

def _require_positive_profile(values: ArrayLike, name: str) -> np.ndarray:
    """Return one finite positive one dimensional profile"""
    arr = np.asarray(values, dtype=float)
    if arr.ndim != 1:
        raise ValueError(f"{name} must be 1D")
    if not np.all(np.isfinite(arr)) or np.any(arr <= 0.0):
        raise ValueError(f"{name} must contain only positive finite values")
    
    return arr

def _volume_average_spatial_lambda_source(source_z_v_lambda_per_s: ArrayLike, cell_volumes_m3: ArrayLike) -> np.ndarray:
    """
    Return the physical volume average of a source with shape `(n_z, n_v, n_lambda)`
    
    The axial source is weighted by cell volume before summing over z
    """
    source = np.asarray(source_z_v_lambda_per_s, dtype=float)
    volumes = _require_positive_profile(cell_volumes_m3, "cell_volumes_m3")
    if source.ndim != 3 or source.shape[0] != volumes.size:
        raise ValueError("source_z_v_lambda_per_s must have shape (n_z, n_v, n_lambda)")
    total_volume = float(np.sum(volumes))

    return np.sum(source * volumes[:, None, None], axis=0) / total_volume

def distribution_roundoff(distribution: ArrayLike, *, relative_tolerance: float, name: str = "distribution") -> tuple[np.ndarray, int, float]:
    """
    Set negative values within the distribution roundoff tolerance to zero and reject materially negative values
    
    The returned tuple contains the corrected array, corrected cell count, and original minimum value
    """
    arr = np.asarray(distribution, dtype=float).copy()
    _, tolerance, min_value, _, significant = _distribution_nonnegativity_metrics(arr, relative_tolerance=relative_tolerance, scale_floor=1.0, name=name)
    if significant:
        raise ValueError(f"{name} contains negative values below the roundoff tolerance: " f"min={min_value:.6e}, tolerance={tolerance:.6e}")
    mask = arr < 0.0
    count = int(np.count_nonzero(mask))
    if count:
        arr[mask] = 0.0

    return arr, count, min_value

def _speed_grid_for_species(config: SourceModelRunConfig, species: IonSpecies, component_energies_J: ArrayLike) -> SpeedGrid:
    """
    Build a source aligned stretched speed grid for one fast ion species
    
    Every distinct enabled beam component speed is aligned to the speed grid and the configured maximum speed is `speed_max_factor` times the largest component speed
    """
    energies = np.asarray(component_energies_J, dtype=float)
    if energies.ndim != 1 or energies.size < 1 or np.any(~np.isfinite(energies)) or np.any(energies <= 0.0):
        raise ValueError("component_energies_J must contain positive finite values")
    component_speeds = np.asarray([beam_speed_from_energy_m_s(float(energy_J), species.mass_kg) for energy_J in energies], dtype=float)
    maximum_source_speed = float(np.max(component_speeds))
    kinetic = config.kinetic_electrostatic
    sorted_speeds = np.sort(component_speeds)
    tolerance = 64.0 * np.finfo(float).eps * max(maximum_source_speed, 1.0)
    distinct_speed_values: list[float] = []
    for value in sorted_speeds:
        if not distinct_speed_values or value - distinct_speed_values[-1] > tolerance:
            distinct_speed_values.append(float(value))
    # Each distinct source speed must fit as an aligned location inside the requested speed grid
    distinct_speeds = np.asarray(distinct_speed_values, dtype=float)
    if distinct_speeds.size >= int(kinetic.speed_bins):
        raise ValueError("kinetic_electrostatic.speed_bins must exceed the number of distinct enabled beam component speeds")
    
    return source_aligned_stretched_speed_grid(max_speed_m_s=float(kinetic.speed_max_factor) * maximum_source_speed, num_cells=int(kinetic.speed_bins), source_speeds_m_s=distinct_speeds, core_cell_fraction=float(kinetic.modal_speed_core_cell_fraction), tail_stretch_power=float(kinetic.modal_speed_tail_stretch_power))

