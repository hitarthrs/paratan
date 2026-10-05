"""Reduced cross species fast ion Rosenbluth states for the Eq 59 collision operator"""
from __future__ import annotations
from collections.abc import Mapping
from dataclasses import dataclass
import numpy as np
from scipy.integrate import cumulative_trapezoid
from source_model_revamp.fbis.collision_parameters import FBISCollisionParameterState
from source_model_revamp.fbis.modal.species_state import FastIonSystemState
from source_model_revamp.fbis.modal.types import ModalRosenbluthCoefficients
from source_model_revamp.fbis.species import IonSpecies
from source_model_revamp.fbis.velocity_space_grid import SpeedGrid

def _validated_distribution(speed_grid: SpeedGrid, values: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Validate one reduced fast field speed distribution and clip negative roundoff"""
    speed = np.asarray(speed_grid.centers_m_s, dtype=float)
    distribution = np.asarray(values, dtype=float)
    if speed.ndim != 1 or speed.size < 2 or distribution.shape != speed.shape:
        raise ValueError("external fast field distribution must match its speed grid")
    if np.any(~np.isfinite(speed)) or np.any(speed <= 0.0) or np.any(np.diff(speed) <= 0.0):
        raise ValueError("external fast field speed grid must be positive and increasing")
    if np.any(~np.isfinite(distribution)):
        raise ValueError("external fast field distribution must be finite")
    scale = max(float(np.max(np.abs(distribution))), 1.0)
    tolerance = 256.0 * np.finfo(float).eps * scale
    if float(np.min(distribution)) < -tolerance:
        raise ValueError("external fast field first mode contains material negative values")
   
    return speed, np.where(distribution < 0.0, 0.0, distribution)

def _lower_moment_at(test_speed: np.ndarray, field_speed: np.ndarray, field_distribution: np.ndarray, power: int) -> np.ndarray:
    """Evaluate M_p(v) = ∫₀ᵛ uᵖ f(u) du on the requested test speeds"""
    cumulative = cumulative_trapezoid(field_speed**power * field_distribution, field_speed, initial=0.0)
  
    return np.interp(test_speed, field_speed, cumulative, left=0.0, right=float(cumulative[-1]))

@dataclass(frozen=True)
class ExternalFastFieldCollisionState:
    """Store one ordered reduced fast test species and fast field species state
    
    field_first_mode_distribution has shape (n_speed_field,) and supplies the field speed dependence
    """
    test_species: IonSpecies
    field_species: IonSpecies
    field_speed_grid: SpeedGrid
    field_first_mode_distribution: np.ndarray
    field_physical_density_m3: float
    coulomb_log_value: float
    model: str = "reduced_isotropic_first_mode_Rosenbluth_cross_species"
    qualified: bool = False
    limitation: str = "The cross species operator uses the isotropic first pitch mode and the existing representative ion Coulomb logarithm rather than a full anisotropic multispecies Landau operator"

    def __post_init__(self) -> None:
        """Validate distinct species the field distribution density Coulomb log and model labels"""
        if self.test_species.species_id == self.field_species.species_id:
            raise ValueError("external fast field state requires distinct test and field species")
        _validated_distribution(self.field_speed_grid, self.field_first_mode_distribution)
        if not np.isfinite(self.field_physical_density_m3) or self.field_physical_density_m3 <= 0.0:
            raise ValueError("external fast field physical density must be positive and finite")
        if not np.isfinite(self.coulomb_log_value) or self.coulomb_log_value <= 0.0:
            raise ValueError("external fast field Coulomb logarithm must be positive and finite")
        if not self.model or not self.limitation:
            raise ValueError("external fast field model and limitation must be nonempty")

    @property
    def pair_id(self) -> str:
        """Return the ordered fast test species and fast field species identifier"""
        return f"fast_{self.test_species.symbol}<-fast_{self.field_species.symbol}"

def external_fast_field_rosenbluth_coefficients(state: ExternalFastFieldCollisionState, test_speed_grid: SpeedGrid) -> ModalRosenbluthCoefficients:
    """Evaluate density normalized isotropic field potentials on the test speed grid
    
    The normalization is n_f = 4π ∫ v² f(v) dv and the lower moments M1 M2 and M4 are evaluated from the field grid
    """
    field_speed, field_distribution = _validated_distribution(state.field_speed_grid, state.field_first_mode_distribution)
    test_speed = np.asarray(test_speed_grid.centers_m_s, dtype=float)
    if test_speed.ndim != 1 or test_speed.size < 2 or np.any(~np.isfinite(test_speed)) or np.any(test_speed <= 0.0) or np.any(np.diff(test_speed) <= 0.0):
        raise ValueError("test speed grid must be positive and increasing")
    density = float(4.0 * np.pi * np.trapezoid(field_speed**2 * field_distribution, field_speed))
    if not np.isfinite(density) or density <= 0.0:
        raise ValueError("external fast field first mode must have positive Rosenbluth normalization")
    M2 = _lower_moment_at(test_speed, field_speed, field_distribution, 2)
    M4 = _lower_moment_at(test_speed, field_speed, field_distribution, 4)
    L1 = _lower_moment_at(test_speed, field_speed, field_distribution, 1)
    total_L1 = float(np.trapezoid(field_speed * field_distribution, field_speed))
    N1 = np.maximum(total_L1 - L1, 0.0)
    h_tilde = 4.0 * np.pi * M2 / density
    g_tilde_1 = 4.0 * np.pi * (M2 - M4 / (3.0 * test_speed**2) + 2.0 * test_speed * N1 / 3.0) / density
    H = 4.0 * np.pi * (M2 / test_speed + N1)
    dG = g_tilde_1 * density
    g_tilde_2 = test_speed**2 * (2.0 * H - 2.0 * dG / test_speed) / density
    arrays = []
    for name, values in (("h_tilde", h_tilde), ("g_tilde_1", g_tilde_1), ("g_tilde_2", g_tilde_2)):
        values = np.asarray(values, dtype=float)
        scale = max(float(np.max(np.abs(values))), 1.0)
        tolerance = 512.0 * np.finfo(float).eps * scale
        if np.any(~np.isfinite(values)) or float(np.min(values)) < -tolerance:
            raise ValueError(f"external fast field {name} is invalid")
        arrays.append(np.where(values < 0.0, 0.0, values))

    return ModalRosenbluthCoefficients(h_tilde=arrays[0], g_tilde_1=arrays[1], g_tilde_2=arrays[2], density_normalization_m3=density)

def build_external_fast_field_collision_states(*, system: FastIonSystemState, collision_states_by_species: Mapping[str, FBISCollisionParameterState]) -> dict[str, tuple[ExternalFastFieldCollisionState, ...]]:
    """Build ordered reduced fast cross species field states from one solved multi species system"""
    active = {species_id: state for species_id, state in system.species_states.items() if state.modal_result is not None}
    if len(active) < 2:
        return {}
    result: dict[str, tuple[ExternalFastFieldCollisionState, ...]] = {}
    for test_species_id in sorted(active):
        test_collision = collision_states_by_species.get(test_species_id)
        if test_collision is None:
            raise ValueError("cross species field construction requires each test collision state")
        fields: list[ExternalFastFieldCollisionState] = []
        for field_species_id in sorted(active):
            if field_species_id == test_species_id:
                continue
            field_state = active[field_species_id]
            modal = field_state.modal_result
            if modal is None or modal.rosenbluth_coefficients is None:
                raise ValueError("cross species field construction requires a solved hot ion Rosenbluth state")
            fields.append(ExternalFastFieldCollisionState(
                test_species=active[test_species_id].species,
                field_species=field_state.species,
                field_speed_grid=field_state.speed_grid,
                field_first_mode_distribution=np.asarray(modal.modal_distribution_f_j_v[0], dtype=float).copy(),
                field_physical_density_m3=float(modal.density_m3),
                coulomb_log_value=float(test_collision.coulomb_logs.ion_fast_ion),
            ))
        if fields:
            result[test_species_id] = tuple(fields)
   
    return result

__all__ = [
    "ExternalFastFieldCollisionState",
    "external_fast_field_rosenbluth_coefficients",
    "build_external_fast_field_collision_states",
]
