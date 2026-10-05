"""Eq 70 electron density closure and quasineutrality residual diagnostics"""
from __future__ import annotations
import numpy as np
from scipy.special import gammainc
from source_model_revamp.constants import ELECTRON_CHARGE_C

def _energy_roundoff_tolerance(*energy_scales_J: float) -> float:
    """Return a floating point energy tolerance scaled to finite inputs in J"""
    finite_scales = [abs(float(value)) for value in energy_scales_J if np.isfinite(value)]
    scale = max(finite_scales, default=np.finfo(float).tiny)

    return max(128.0 * np.finfo(float).eps * scale, np.finfo(float).tiny)

def _clip_potential_roundoff_only(potential_drop_magnitude_J: np.ndarray | float, wall_barrier_energy_J: float, *, context: str) -> tuple[np.ndarray, int]:
    """Clip only roundoff scale excursions outside the physical interval [0, wall] in J"""
    values = np.asarray(potential_drop_magnitude_J, dtype=float)
    wall = float(wall_barrier_energy_J)
    if not np.isfinite(wall) or wall < 0.0:
        raise ValueError("wall_barrier_energy_J must be finite and nonnegative")
    if np.any(~np.isfinite(values)):
        raise ValueError(f"{context} must contain only finite values")
    tolerance = _energy_roundoff_tolerance(wall, float(np.max(np.abs(values))) if values.size else 0.0)
    significant = (values < -tolerance) | (values > wall + tolerance)
    if np.any(significant):
        minimum = float(np.min(values))
        maximum = float(np.max(values))
        raise ValueError(f"{context} lies outside [0, wall_barrier_energy_J]: " f"min={minimum:.17g} J, max={maximum:.17g} J, wall={wall:.17g} J")
    roundoff = (values < 0.0) | (values > wall)

    return np.clip(values, 0.0, wall), int(np.count_nonzero(roundoff))

def _electron_density_fraction_eq70(phi_energy_J: float, wall_barrier_energy_J: float, electron_temperature_J: float) -> float:
    """
    Return the confined Maxwellian electron fraction from Egedal Eq 70

    phi_energy_J and wall_barrier_energy_J are positive energy magnitudes in J
    The original sign convention is -e * ϕ(z) and -e * ϕ_w with ϕ(0) = 0
    The returned fraction is n_e / n0 = exp(−x_ϕ²) P(3/2, x_Δ²)
    """
    Te = float(electron_temperature_J)
    wall = float(wall_barrier_energy_J)
    phi = float(phi_energy_J)
    if not np.isfinite(Te) or Te <= 0.0:
        raise ValueError("electron_temperature_J must be positive and finite")
    if not np.isfinite(wall) or wall < 0.0:
        raise ValueError("wall_barrier_energy_J must be finite and nonnegative")
    tolerance = _energy_roundoff_tolerance(wall, Te)
    if phi < -tolerance or phi > wall + tolerance or not np.isfinite(phi):
        raise ValueError("phi_energy_J must lie between the midplane and wall barriers")
    phi = float(_clip_potential_roundoff_only(phi, wall, context="phi_energy_J")[0])
    x_phi_squared = phi / Te
    x_delta = np.sqrt(max(wall - phi, 0.0) / Te)
    bracket = gammainc(1.5, x_delta**2)

    return float(max(np.exp(-x_phi_squared) * bracket, 0.0))

def _quasineutrality_residual_metrics(*, electron_density_m3: np.ndarray, ion_density_m3: np.ndarray, zeta: np.ndarray, cell_volumes_m3: np.ndarray, reference_density_m3: float, support_relative_floor: float = 1.0e-12, central_region_abs_zeta_max: float = 0.5, throat_region_abs_zeta_min: float = 0.9) -> dict[str, object]:
    """Return supported, absolute, and volume integrated Eq 70 residuals

    Relative residuals use max(n_e, n_i) where the local density exceeds the support floor
    Absolute and volume integrated measures retain low density cells in the closure audit
    """
    ne = np.asarray(electron_density_m3, dtype=float)
    ni = np.asarray(ion_density_m3, dtype=float)
    z = np.asarray(zeta, dtype=float)
    volumes = np.asarray(cell_volumes_m3, dtype=float)
    if ne.ndim != 1 or ni.shape != ne.shape or z.shape != ne.shape or volumes.shape != ne.shape:
        raise ValueError("density, zeta, and cell volume profiles must be matching vectors")
    if np.any(~np.isfinite(ne)) or np.any(ne < 0.0) or np.any(~np.isfinite(ni)) or np.any(ni < 0.0):
        raise ValueError("electron and ion densities must be finite and nonnegative")
    if np.any(~np.isfinite(z)) or np.any(~np.isfinite(volumes)) or np.any(volumes <= 0.0):
        raise ValueError("zeta and cell volumes must be finite and cell volumes positive")
    reference = float(reference_density_m3)
    relative_floor = float(support_relative_floor)
    if not np.isfinite(reference) or reference < 0.0:
        raise ValueError("reference_density_m3 must be finite and nonnegative")
    if not np.isfinite(relative_floor) or relative_floor <= 0.0 or relative_floor >= 1.0:
        raise ValueError("support_relative_floor must lie inside (0, 1)")
    if reference == 0.0:
        reference = max(float(np.max(ne)), float(np.max(ni)), np.finfo(float).tiny)
    density_floor = max(reference * relative_floor, np.finfo(float).tiny)
    absolute = ne - ni
    density_scale = np.maximum(ne, ni)
    support = density_scale >= density_floor
    supported_relative = np.zeros_like(absolute)
    np.divide(absolute, density_scale, out=supported_relative, where=support & (density_scale > 0.0))
    normalized_absolute = absolute / reference
    central = support & (np.abs(z) <= float(central_region_abs_zeta_max))
    throat = support & (np.abs(z) >= float(throat_region_abs_zeta_min))
    signed_particle_mismatch = float(np.sum(absolute * volumes))
    absolute_particle_mismatch = float(np.sum(np.abs(absolute) * volumes))
    total_volume = float(np.sum(volumes))
    masked_volume = float(np.sum(volumes[~support]))
    normalization_particles = max(reference * total_volume, np.finfo(float).tiny)
    maximum_absolute_normalized = (float(np.max(np.abs(normalized_absolute))) if normalized_absolute.size else 0.0)

    return {
        "signed_absolute_density_residual_m3": absolute,
        "absolute_residual_normalized_to_reference": normalized_absolute,
        "supported_relative_residual": supported_relative,
        "density_support_mask": support,
        "reference_density_m3": reference,
        "density_support_relative_floor": relative_floor,
        "density_support_floor_m3": density_floor,
        "density_support_floor_derivation": "actual_midplane_density_times_reported_relative_floor",
        "masked_node_count": int(np.count_nonzero(~support)),
        "masked_volume_m3": masked_volume,
        "masked_volume_fraction": masked_volume / total_volume if total_volume > 0.0 else 0.0,
        "volume_integrated_electron_minus_ion_particles": signed_particle_mismatch,
        "volume_integrated_absolute_particle_mismatch": absolute_particle_mismatch,
        "volume_integrated_signed_particle_mismatch_fraction": (signed_particle_mismatch / normalization_particles),
        "volume_integrated_absolute_particle_mismatch_fraction": (absolute_particle_mismatch / normalization_particles),
        "volume_integrated_signed_charge_mismatch_C": -ELECTRON_CHARGE_C * signed_particle_mismatch,
        "volume_integrated_absolute_charge_mismatch_C": ELECTRON_CHARGE_C * absolute_particle_mismatch,
        "maximum_supported_relative_residual": (float(np.max(np.abs(supported_relative[support]))) if np.any(support) else 0.0),
        "maximum_absolute_density_residual_m3": float(np.max(np.abs(absolute))) if absolute.size else 0.0,
        "maximum_absolute_residual_normalized_to_reference": maximum_absolute_normalized,
        "central_region_maximum_supported_relative_residual": (float(np.max(np.abs(supported_relative[central]))) if np.any(central) else 0.0),
        "throat_region_maximum_supported_relative_residual": (float(np.max(np.abs(supported_relative[throat]))) if np.any(throat) else 0.0),
        "central_region_abs_zeta_max": float(central_region_abs_zeta_max),
        "throat_region_abs_zeta_min": float(throat_region_abs_zeta_min),
    }
