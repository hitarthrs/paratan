"""Apply the Egedal low energy electrostatic eigenvalue approximation

The helper preserves the magnetic well shape while constructing the reference mirror ratio basis used for the first mode substitution
"""
from __future__ import annotations
from dataclasses import dataclass
from typing import Callable
import numpy as np
from numpy.typing import ArrayLike

EGEDAL_2022_PUBLISHED_LOW_ENERGY_LAMBDA1_HEURISTIC = "egedal_2022_published_low_energy_lambda1_heuristic"
EGEDAL_2022_REFERENCE_MIRROR_RATIO = 1.5

@dataclass(frozen=True)
class PublishedLambda1Heuristic:
    """Published low energy first eigenvalue substitution state
    
    eigenvalues_by_speed has shape (n_mode, n_speed) and energy_threshold_J is the invariant energy threshold for substitution
    """
    model: str
    reference: str
    is_exact: bool
    eigenvalues_by_speed: np.ndarray
    active_speed_mask: np.ndarray
    energy_threshold_J: float
    magnetic_lambda1: float
    substituted_lambda1: float
    active_population_fraction: float
    population_fraction_assessed: bool
    geometry_extrapolation: bool

def rescaled_magnetic_profile_for_reference_ratio(B_tilde_function: Callable[[float], float], active_mirror_ratio: float, reference_mirror_ratio: float = EGEDAL_2022_REFERENCE_MIRROR_RATIO) -> Callable[[ArrayLike], float | np.ndarray]:
    """Return a magnetic profile with unchanged normalized well shape and a new throat contrast"""
    active = float(active_mirror_ratio)
    reference = float(reference_mirror_ratio)
    if not np.isfinite(active) or active <= 1.0:
        raise ValueError("active_mirror_ratio must be finite and greater than one")
    if not np.isfinite(reference) or reference <= 1.0:
        raise ValueError("reference_mirror_ratio must be finite and greater than one")
    contrast_scale = (reference - 1.0) / (active - 1.0)

    def reference_profile(zeta: ArrayLike) -> float | np.ndarray:
        """Evaluate the rescaled dimensionless magnetic field profile"""
        active_value = np.asarray(B_tilde_function(zeta), dtype=float)
        if np.any(~np.isfinite(active_value)) or np.any(active_value <= 0.0):
            raise ValueError("B_tilde_function returned a nonpositive or nonfinite value")
        result = 1.0 + contrast_scale * (active_value - 1.0)
        if np.any(~np.isfinite(result)) or np.any(result <= 0.0):
            raise ValueError("rescaled reference magnetic profile is nonpositive or nonfinite")
        return float(result) if np.ndim(zeta) == 0 else result

    return reference_profile

def egedal_2022_published_low_energy_lambda1_heuristic(*, speed_centers_m_s: ArrayLike, particle_mass_kg: float, mirror_ratio: float, throat_potential_drop_magnitude_J: float, magnetic_eigenvalues: ArrayLike, reference_lambda1_at_mirror_ratio_1p5: float, invariant_energy_population_weights: ArrayLike | None = None, geometry_extrapolation: bool = True) -> PublishedLambda1Heuristic:
    """Apply the published low energy substitution to the first physical eigenvalue
    
    The substitution is active where invariant kinetic energy is below the throat potential drop divided by mirror ratio
    Higher retained eigenvalues remain unchanged
    """
    speed = np.asarray(speed_centers_m_s, dtype=float)
    eigenvalues = np.asarray(magnetic_eigenvalues, dtype=float)
    mass = float(particle_mass_kg)
    ratio = float(mirror_ratio)
    throat_drop = float(throat_potential_drop_magnitude_J)
    reference_lambda1 = float(reference_lambda1_at_mirror_ratio_1p5)
    if speed.ndim != 1 or speed.size == 0:
        raise ValueError("speed_centers_m_s must be a nonempty vector")
    if np.any(~np.isfinite(speed)) or np.any(speed < 0.0):
        raise ValueError("speed_centers_m_s must be finite and nonnegative")
    if eigenvalues.ndim != 1 or eigenvalues.size == 0:
        raise ValueError("magnetic_eigenvalues must be a nonempty vector")
    if np.any(~np.isfinite(eigenvalues)) or np.any(eigenvalues <= 0.0):
        raise ValueError("magnetic_eigenvalues must be positive and finite")
    if not np.isfinite(mass) or mass <= 0.0:
        raise ValueError("particle_mass_kg must be positive and finite")
    if not np.isfinite(ratio) or ratio <= 1.0:
        raise ValueError("mirror_ratio must be finite and greater than one")
    if not np.isfinite(throat_drop) or throat_drop < 0.0:
        raise ValueError("throat_potential_drop_magnitude_J must be finite and nonnegative")
    if not np.isfinite(reference_lambda1) or reference_lambda1 <= 0.0:
        raise ValueError("reference_lambda1_at_mirror_ratio_1p5 must be positive and finite")

    energy = 0.5 * mass * speed**2
    threshold = throat_drop / ratio
    active = energy < threshold
    substituted_lambda1 = max(float(eigenvalues[0]), reference_lambda1)
    by_speed = np.broadcast_to(eigenvalues[:, None], (eigenvalues.size, speed.size)).copy()
    by_speed[0, active] = substituted_lambda1
    fraction = 0.0
    assessed = invariant_energy_population_weights is not None
    if invariant_energy_population_weights is not None:
        weights = np.asarray(invariant_energy_population_weights, dtype=float)
        if weights.shape != speed.shape:
            raise ValueError("invariant_energy_population_weights must match the speed grid")
        if np.any(~np.isfinite(weights)) or np.any(weights < 0.0):
            raise ValueError("invariant_energy_population_weights must be finite and nonnegative")
        total = float(np.sum(weights))
        fraction = float(np.sum(weights[active]) / total) if total > 0.0 else 0.0

    return PublishedLambda1Heuristic(
        model=EGEDAL_2022_PUBLISHED_LOW_ENERGY_LAMBDA1_HEURISTIC,
        reference="Egedal_et_al_Nuclear_Fusion_62_126053_2022_p15",
        is_exact=False,
        eigenvalues_by_speed=by_speed,
        active_speed_mask=active,
        energy_threshold_J=float(threshold),
        magnetic_lambda1=float(eigenvalues[0]),
        substituted_lambda1=float(substituted_lambda1),
        active_population_fraction=fraction,
        population_fraction_assessed=assessed,
        geometry_extrapolation=bool(geometry_extrapolation),
    )
