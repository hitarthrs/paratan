"""Evaluate pair resolved energy moments and the discrete Eq 59 energy identity"""
from __future__ import annotations
from dataclasses import dataclass
import numpy as np
from source_model_revamp.fbis.modal.eq59.collision_operator_state import Eq59CollisionOperatorState, Eq59PairOperatorContribution
from source_model_revamp.fbis.modal.eq59.operator import _face_average, _scharfetter_gummel_face_coefficients
from source_model_revamp.fbis.modal.velocity_solve_cold import _source_cell_weights
from source_model_revamp.fbis.modal.types import ModalFBISResult
from source_model_revamp.fbis.velocity_space_grid import SpeedGrid

@dataclass(frozen=True)
class Eq59PairEnergyMoment:
    """Store speed energy exchange and pitch loss for one ordered field pair
    
    A positive field_population_energy_gain_W is the equal and opposite energy lost by the fast test species
    """
    pair_id: str
    field_population_id: str
    field_species_id: str
    test_species_speed_energy_change_W: float
    field_population_energy_gain_W: float | None
    discrete_pitch_loss_power_W: float
    active: bool
    operator_qualified: bool
    limitation: str | None

    def __post_init__(self) -> None:
        """Validate pair identifiers finite moments nonnegative pitch loss and antisymmetric exchange"""
        if not self.pair_id or not self.field_population_id or not self.field_species_id:
            raise ValueError("Eq 59 pair energy moment identifiers must be nonempty")
        values = (self.test_species_speed_energy_change_W, self.discrete_pitch_loss_power_W)
        if any(not np.isfinite(float(value)) for value in values):
            raise ValueError("Eq 59 pair energy moment values must be finite")
        if self.field_population_energy_gain_W is not None and not np.isfinite(float(self.field_population_energy_gain_W)):
            raise ValueError("Eq 59 field population energy gain must be finite when available")
        if self.discrete_pitch_loss_power_W < 0.0:
            raise ValueError("Eq 59 pair pitch loss power must be nonnegative")
        if self.field_population_energy_gain_W is not None:
            scale = max(abs(float(self.test_species_speed_energy_change_W)), abs(float(self.field_population_energy_gain_W)), 1.0)
            if abs(float(self.test_species_speed_energy_change_W) + float(self.field_population_energy_gain_W)) > 256.0 * np.finfo(float).eps * scale:
                raise ValueError("Eq 59 pair energy exchange must be antisymmetric")

@dataclass(frozen=True)
class Eq59EnergyMomentState:
    """Store pair resolved Eq 59 energy moments and discrete balance diagnostics
    
    The discrete identity is P_source + ΔP_speed − P_pitch − P_CX = residual
    """
    fast_species_id: str
    electron_pair_id: str
    electron_heating_power_W: float
    fast_ion_energy_change_from_electrons_W: float
    pair_energy_moments: tuple[Eq59PairEnergyMoment, ...]
    total_speed_operator_energy_change_W: float
    pair_sum_speed_operator_energy_change_W: float
    pair_energy_identity_relative_error: float
    authoritative_confined_source_power_W: float
    projected_exact_source_power_W: float
    discrete_source_power_W: float
    source_projection_relative_error: float
    source_discretization_relative_error: float
    discrete_pitch_loss_power_W: float
    charge_exchange_sink_power_W: float
    charge_exchange_sink_active: bool
    energy_identity_includes_charge_exchange_sink: bool
    discrete_energy_residual_W: float
    discrete_energy_residual_relative: float
    pair_flux_maximum_absolute_error: float
    pair_flux_maximum_norm_scaled_error: float
    pair_flux_l2_relative_error: float
    pair_flux_identity_relative_tolerance: float
    energy_identity_relative_tolerance: float
    source_power_relative_tolerance: float
    pair_flux_identity_passed: bool
    pair_energy_identity_passed: bool
    discrete_energy_identity_passed: bool
    source_projection_passed: bool
    source_discretization_passed: bool
    electron_pair_active: bool
    electron_pair_operator_qualified: bool
    electron_pair_numerically_valid: bool
    electron_pair_qualified: bool
    full_collision_model_qualified: bool
    numerically_valid: bool
    qualified: bool
    failure_reason: str | None
    model: str = "direct_pair_resolved_conservative_eq59_energy_moment"

    def __post_init__(self) -> None:
        """Validate energy moment values tolerances qualification logic and electron exchange antisymmetry"""
        if not self.fast_species_id or not self.electron_pair_id or not self.model:
            raise ValueError("Eq 59 energy moment identifiers must be nonempty")
        scalar_values = (
            self.electron_heating_power_W,
            self.fast_ion_energy_change_from_electrons_W,
            self.total_speed_operator_energy_change_W,
            self.pair_sum_speed_operator_energy_change_W,
            self.pair_energy_identity_relative_error,
            self.authoritative_confined_source_power_W,
            self.projected_exact_source_power_W,
            self.discrete_source_power_W,
            self.source_projection_relative_error,
            self.source_discretization_relative_error,
            self.discrete_pitch_loss_power_W,
            self.charge_exchange_sink_power_W,
            self.discrete_energy_residual_W,
            self.discrete_energy_residual_relative,
            self.pair_flux_maximum_absolute_error,
            self.pair_flux_maximum_norm_scaled_error,
            self.pair_flux_l2_relative_error,
            self.pair_flux_identity_relative_tolerance,
            self.energy_identity_relative_tolerance,
            self.source_power_relative_tolerance,
        )
        if any(not np.isfinite(float(value)) for value in scalar_values):
            raise ValueError("Eq 59 energy moment values must be finite")
        nonnegative_values = (
            self.pair_energy_identity_relative_error,
            self.authoritative_confined_source_power_W,
            self.projected_exact_source_power_W,
            self.discrete_source_power_W,
            self.source_projection_relative_error,
            self.source_discretization_relative_error,
            self.discrete_pitch_loss_power_W,
            self.charge_exchange_sink_power_W,
            self.discrete_energy_residual_relative,
            self.pair_flux_maximum_absolute_error,
            self.pair_flux_maximum_norm_scaled_error,
            self.pair_flux_l2_relative_error,
        )
        if any(float(value) < 0.0 for value in nonnegative_values):
            raise ValueError("Eq 59 energy moment losses and errors must be nonnegative")
        if self.pair_flux_identity_relative_tolerance <= 0.0 or self.energy_identity_relative_tolerance <= 0.0 or self.source_power_relative_tolerance <= 0.0:
            raise ValueError("Eq 59 energy moment tolerances must be positive")
        scale = max(abs(float(self.electron_heating_power_W)), abs(float(self.fast_ion_energy_change_from_electrons_W)), 1.0)
        if abs(float(self.electron_heating_power_W) + float(self.fast_ion_energy_change_from_electrons_W)) > 256.0 * np.finfo(float).eps * scale:
            raise ValueError("fast ion electron heating and test species energy change must be antisymmetric")
        if self.qualified and not self.numerically_valid:
            raise ValueError("qualified Eq 59 energy moments must be numerically valid")
        if self.qualified and not self.electron_pair_qualified:
            raise ValueError("qualified Eq 59 energy moments require a qualified electron pair")

def _positive_finite_scalar(value: float, name: str) -> float:
    """Validate one strictly positive finite scalar"""
    scalar = float(value)
    if not np.isfinite(scalar) or scalar <= 0.0:
        raise ValueError(f"{name} must be positive and finite")
   
    return scalar

def _validated_inputs(*, speed_grid: SpeedGrid, modal_distribution_f_j_v: np.ndarray, source_coefficients_by_component_j: np.ndarray, source_speeds_m_s: np.ndarray, active_eigenvalues_by_speed: np.ndarray, collision_operator_state: Eq59CollisionOperatorState, mode_eta_integrals: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Validate modal source eigenvalue and η integral shapes used by the energy moment calculation"""
    speed = np.asarray(speed_grid.centers_m_s, dtype=float)
    modal = np.asarray(modal_distribution_f_j_v, dtype=float)
    source_coefficients = np.asarray(source_coefficients_by_component_j, dtype=float)
    source_speeds = np.asarray(source_speeds_m_s, dtype=float)
    eigenvalues = np.asarray(active_eigenvalues_by_speed, dtype=float)
    eta_integrals = np.asarray(mode_eta_integrals, dtype=float)
    if modal.ndim != 2 or modal.shape[1] != speed.size or np.any(~np.isfinite(modal)):
        raise ValueError("modal_distribution_f_j_v must be finite with shape n_mode by n_speed")
    if source_coefficients.ndim != 2 or source_coefficients.shape[1] != modal.shape[0] or np.any(~np.isfinite(source_coefficients)):
        raise ValueError("source_coefficients_by_component_j must be finite with shape n_component by n_mode")
    if source_speeds.shape != (source_coefficients.shape[0],) or np.any(~np.isfinite(source_speeds)) or np.any(source_speeds <= 0.0):
        raise ValueError("source_speeds_m_s must contain one positive finite speed per source component")
    if eigenvalues.ndim == 1:
        eigenvalues = np.broadcast_to(eigenvalues[:, None], modal.shape)
    if eigenvalues.shape != modal.shape or np.any(~np.isfinite(eigenvalues)) or np.any(eigenvalues <= 0.0):
        raise ValueError("active_eigenvalues_by_speed must be positive and match the modal distribution")
    if eta_integrals.shape != (modal.shape[0],) or np.any(~np.isfinite(eta_integrals)):
        raise ValueError("mode_eta_integrals must contain one finite value per mode")
    collision_speed = np.asarray(collision_operator_state.speed_grid.centers_m_s, dtype=float)
    collision_faces = np.asarray(collision_operator_state.speed_grid.faces_m_s, dtype=float)
    if not np.array_equal(collision_speed, speed) or not np.array_equal(collision_faces, np.asarray(speed_grid.faces_m_s, dtype=float)):
        raise ValueError("Eq 59 collision operator grid must match the energy moment grid")
  
    return modal, source_coefficients, source_speeds, eigenvalues, eta_integrals

def _face_flux_and_fitted_state(*, speed_grid: SpeedGrid, modal_distribution: np.ndarray, drift_center: np.ndarray, diffusion_center: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Reconstruct the Scharfetter Gummel face flux distribution and gradient
    
    The fitted face state is the common state used to decompose the total speed flux into ordered pair fluxes
    """
    speed = np.asarray(speed_grid.centers_m_s, dtype=float)
    faces = np.asarray(speed_grid.faces_m_s, dtype=float)
    mode_count, cell_count = modal_distribution.shape
    drift = np.asarray(drift_center, dtype=float)
    diffusion = np.asarray(diffusion_center, dtype=float)
    if cell_count != speed.size or drift.shape != speed.shape or diffusion.shape != speed.shape:
        raise ValueError("Eq 59 fitted face state grid mismatch")
    if np.any(~np.isfinite(drift)) or np.any(~np.isfinite(diffusion)) or np.any(diffusion < 0.0):
        raise ValueError("Eq 59 fitted face coefficients are invalid")
    drift_face = _face_average(drift)
    diffusion_face = _face_average(diffusion)
    flux = np.zeros((mode_count, cell_count + 1), dtype=float)
    fitted_distribution = np.zeros_like(flux)
    fitted_gradient = np.zeros_like(flux)
    for face_index in range(1, cell_count):
        left_index = face_index - 1
        right_index = face_index
        delta_v = float(speed[right_index] - speed[left_index])
        left_to_face = float(faces[face_index] - speed[left_index])
        b = float(drift_face[face_index])
        D = float(diffusion_face[face_index])
        left_coefficient, right_coefficient = _scharfetter_gummel_face_coefficients(b, D, delta_v)
        face_flux = left_coefficient * modal_distribution[:, left_index] + right_coefficient * modal_distribution[:, right_index]
        flux[:, face_index] = face_flux
        if D == 0.0:
            if b > 0.0:
                fitted_distribution[:, face_index] = modal_distribution[:, right_index]
            elif b < 0.0:
                fitted_distribution[:, face_index] = modal_distribution[:, left_index]
            else:
                fitted_distribution[:, face_index] = 0.5 * (modal_distribution[:, left_index] + modal_distribution[:, right_index])
            continue
        peclet = b * delta_v / D
        if abs(peclet) < 1.0e-7 or abs(b) <= np.finfo(float).tiny:
            fraction = left_to_face / delta_v
            face_value = modal_distribution[:, left_index] + fraction * (modal_distribution[:, right_index] - modal_distribution[:, left_index])
        else:
            equilibrium = face_flux / b
            exponent = max(min(-b * left_to_face / D, 0.0), -745.0)
            face_value = equilibrium + (modal_distribution[:, left_index] - equilibrium) * np.exp(exponent)
        face_gradient = (face_flux - b * face_value) / D
        if np.any(~np.isfinite(face_value)) or np.any(~np.isfinite(face_gradient)):
            raise RuntimeError("Eq 59 face state reconstruction produced nonfinite values")
        fitted_distribution[:, face_index] = face_value
        fitted_gradient[:, face_index] = face_gradient
   
    return flux, fitted_distribution, fitted_gradient

def _pair_drift_diffusion_center(*, contribution: Eq59PairOperatorContribution, speed: np.ndarray, electron_pair_id: str) -> tuple[np.ndarray, np.ndarray]:
    """Return center drift and diffusion arrays for one ordered pair contribution"""
    if contribution.pair_id == electron_pair_id:
        return speed**3, np.zeros_like(speed)
  
    return np.asarray(contribution.drag_velocity_cubed_m3_s3, dtype=float), np.asarray(contribution.energy_diffusion_velocity_fourth_m4_s4, dtype=float)

def _pair_fluxes(*, collision: Eq59CollisionOperatorState, speed: np.ndarray, fitted_distribution: np.ndarray, fitted_gradient: np.ndarray) -> dict[str, np.ndarray]:
    """Reconstruct every ordered pair speed flux from the common fitted face state"""
    pair_fluxes: dict[str, np.ndarray] = {}
    cell_count = speed.size
    for pair_id, contribution in collision.pair_contributions_by_id.items():
        drift_center, diffusion_center = _pair_drift_diffusion_center(contribution=contribution, speed=speed, electron_pair_id=collision.electron_pair_id)
        drift_face = _face_average(drift_center)
        diffusion_face = _face_average(diffusion_center)
        pair_flux = np.zeros_like(fitted_distribution)
        pair_flux[:, 1:-1] = drift_face[None, 1:-1] * fitted_distribution[:, 1:-1] + diffusion_face[None, 1:-1] * fitted_gradient[:, 1:-1]
        pair_fluxes[pair_id] = pair_flux
  
    return pair_fluxes

def _speed_flux_energy_change_W(*, face_flux: np.ndarray, speed: np.ndarray, mode_eta_integrals: np.ndarray, mass_kg: float, volume_m3: float, slowing_time_s: float) -> float:
    """Integrate the speed flux divergence against E = m v² / 2 and the modal η measures"""
    energy = 0.5 * float(mass_kg) * speed**2
    divergence = face_flux[:, 1:] - face_flux[:, :-1]
    weighted = mode_eta_integrals[:, None] * energy[None, :] * divergence
  
    return float(4.0 * np.pi * float(volume_m3) * np.sum(weighted) / float(slowing_time_s))

def _discrete_pitch_loss_power_W(*, speed_grid: SpeedGrid, modal_distribution: np.ndarray, active_eigenvalues_by_speed: np.ndarray, mode_eta_integrals: np.ndarray, pitch_scattering_velocity_cubed_m3_s3: np.ndarray, mass_kg: float, volume_m3: float, slowing_time_s: float) -> float:
    """Integrate the Eq 59 modal pitch scattering sink against kinetic energy"""
    speed = np.asarray(speed_grid.centers_m_s, dtype=float)
    widths = np.asarray(speed_grid.widths_m_s, dtype=float)
    pitch = np.asarray(pitch_scattering_velocity_cubed_m3_s3, dtype=float)
    energy = 0.5 * float(mass_kg) * speed**2
    weighted = modal_distribution * active_eigenvalues_by_speed * mode_eta_integrals[:, None] * pitch[None, :] * widths[None, :] * energy[None, :] / speed[None, :]
    value = float(4.0 * np.pi * float(volume_m3) * np.sum(weighted) / float(slowing_time_s))
    scale = max(abs(value), 1.0)
    if value < -256.0 * np.finfo(float).eps * scale:
        raise RuntimeError("Eq 59 discrete pitch loss power is materially negative")
  
    return max(value, 0.0)

def _charge_exchange_sink_power_W(*, speed_grid: SpeedGrid, modal_distribution: np.ndarray, mode_eta_integrals: np.ndarray, loss_frequency_s: np.ndarray | None, mass_kg: float, volume_m3: float) -> tuple[float, bool]:
    """Return the discrete Eq 59 charge exchange energy removal"""
    if loss_frequency_s is None:
        return 0.0, False
    frequency = np.asarray(loss_frequency_s, dtype=float)
    speed = np.asarray(speed_grid.centers_m_s, dtype=float)
    widths = np.asarray(speed_grid.widths_m_s, dtype=float)
    if frequency.shape != speed.shape or np.any(~np.isfinite(frequency)) or np.any(frequency < 0.0):
        raise ValueError("charge exchange loss frequency must be finite and nonnegative on the Eq 59 speed grid")
    active = bool(np.any(frequency > 0.0))
    if not active:
        return 0.0, False
    energy = 0.5 * float(mass_kg) * speed**2
    integrated_distribution = np.sum(modal_distribution * mode_eta_integrals[:, None], axis=0)
    value = float(4.0 * np.pi * float(volume_m3) * np.sum(integrated_distribution * frequency * speed**2 * widths * energy))
    scale = max(abs(value), 1.0)
    if value < -256.0 * np.finfo(float).eps * scale:
        raise RuntimeError("Eq 59 discrete charge exchange sink power is materially negative")
 
    return max(value, 0.0), True

def _source_powers_W(*, speed_grid: SpeedGrid, source_coefficients: np.ndarray, source_speeds: np.ndarray, mode_eta_integrals: np.ndarray, mass_kg: float, volume_m3: float) -> tuple[float, float]:
    """Return source power from exact component speeds and from their represented speed cells"""
    speed = np.asarray(speed_grid.centers_m_s, dtype=float)
    center_energy = 0.5 * float(mass_kg) * speed**2
    source_mode_weight = source_coefficients * mode_eta_integrals[None, :]
    projected_exact = 0.0
    discrete = 0.0
    for component_index, source_speed in enumerate(source_speeds):
        mode_sum = float(np.sum(source_mode_weight[component_index]))
        exact_energy = 0.5 * float(mass_kg) * float(source_speed) ** 2
        represented_energy = float(np.sum(_source_cell_weights(speed_grid, float(source_speed)) * center_energy))
        projected_exact += 4.0 * np.pi * float(volume_m3) * mode_sum * exact_energy
        discrete += 4.0 * np.pi * float(volume_m3) * mode_sum * represented_energy
   
    return float(projected_exact), float(discrete)

def build_eq59_energy_moment_state(*, speed_grid: SpeedGrid, modal_distribution_f_j_v: np.ndarray, source_coefficients_by_component_j: np.ndarray, source_speeds_m_s: np.ndarray, active_eigenvalues_by_speed: np.ndarray, collision_operator_state: Eq59CollisionOperatorState, mode_eta_integrals: np.ndarray, particle_mass_kg: float, volume_m3: float, exact_confined_source_power_W: float, charge_exchange_loss_frequency_s: np.ndarray | None = None, pair_flux_relative_tolerance: float = 1.0e-10, energy_residual_relative_tolerance: float = 1.0e-4, source_power_relative_tolerance: float = 1.0e-2) -> Eq59EnergyMomentState:
    """Build pair energy exchange and the full discrete Eq 59 energy identity
    
    The calculation verifies pair flux reconstruction pair energy reconstruction source representation pitch loss and optional charge exchange removal
    """
    mass = _positive_finite_scalar(particle_mass_kg, "particle_mass_kg")
    volume = _positive_finite_scalar(volume_m3, "volume_m3")
    exact_source = float(exact_confined_source_power_W)
    if not np.isfinite(exact_source) or exact_source < 0.0:
        raise ValueError("exact_confined_source_power_W must be finite and nonnegative")
    flux_tolerance = _positive_finite_scalar(pair_flux_relative_tolerance, "pair_flux_relative_tolerance")
    energy_tolerance = _positive_finite_scalar(energy_residual_relative_tolerance, "energy_residual_relative_tolerance")
    source_tolerance = _positive_finite_scalar(source_power_relative_tolerance, "source_power_relative_tolerance")
    modal, source_coefficients, source_speeds, eigenvalues, eta_integrals = _validated_inputs(speed_grid=speed_grid, modal_distribution_f_j_v=modal_distribution_f_j_v, source_coefficients_by_component_j=source_coefficients_by_component_j, source_speeds_m_s=source_speeds_m_s, active_eigenvalues_by_speed=active_eigenvalues_by_speed, collision_operator_state=collision_operator_state, mode_eta_integrals=mode_eta_integrals)
    speed = np.asarray(speed_grid.centers_m_s, dtype=float)
    slowing_time = _positive_finite_scalar(collision_operator_state.spitzer_slowing_down_time_s, "spitzer_slowing_down_time_s")
    total_drift = speed**3 + np.asarray(collision_operator_state.ion_drag_velocity_cubed_m3_s3, dtype=float)
    total_diffusion = np.asarray(collision_operator_state.ion_energy_diffusion_velocity_fourth_m4_s4, dtype=float)
    total_flux, fitted_distribution, fitted_gradient = _face_flux_and_fitted_state(speed_grid=speed_grid, modal_distribution=modal, drift_center=total_drift, diffusion_center=total_diffusion)
    # Pair fluxes use the same fitted face state as the total flux identity
    pair_fluxes = _pair_fluxes(collision=collision_operator_state, speed=speed, fitted_distribution=fitted_distribution, fitted_gradient=fitted_gradient)
    pair_sum_flux = np.sum(np.stack(tuple(pair_fluxes.values()), axis=0), axis=0)
    flux_difference = pair_sum_flux - total_flux
    internal_total = total_flux[:, 1:-1]
    internal_pair = pair_sum_flux[:, 1:-1]
    internal_difference = flux_difference[:, 1:-1]
    maximum_absolute_error = float(np.max(np.abs(internal_difference))) if internal_difference.size else 0.0
    maximum_scale = max(float(np.max(np.abs(internal_total))) if internal_total.size else 0.0, float(np.max(np.abs(internal_pair))) if internal_pair.size else 0.0, np.finfo(float).tiny)
    maximum_norm_scaled_error = maximum_absolute_error / maximum_scale
    l2_scale = max(float(np.linalg.norm(internal_total.ravel())), float(np.linalg.norm(internal_pair.ravel())), np.finfo(float).tiny)
    l2_relative_error = float(np.linalg.norm(internal_difference.ravel())) / l2_scale
    pair_flux_identity_passed = bool(maximum_norm_scaled_error <= flux_tolerance and l2_relative_error <= flux_tolerance)
    total_energy_change = _speed_flux_energy_change_W(face_flux=total_flux, speed=speed, mode_eta_integrals=eta_integrals, mass_kg=mass, volume_m3=volume, slowing_time_s=slowing_time)
    pair_moments: list[Eq59PairEnergyMoment] = []
    pair_sum_energy_change = 0.0
    for pair_id, contribution in collision_operator_state.pair_contributions_by_id.items():
        test_change = _speed_flux_energy_change_W(face_flux=pair_fluxes[pair_id], speed=speed, mode_eta_integrals=eta_integrals, mass_kg=mass, volume_m3=volume, slowing_time_s=slowing_time)
        pitch_loss = _discrete_pitch_loss_power_W(speed_grid=speed_grid, modal_distribution=modal, active_eigenvalues_by_speed=eigenvalues, mode_eta_integrals=eta_integrals, pitch_scattering_velocity_cubed_m3_s3=np.asarray(contribution.pitch_scattering_velocity_cubed_m3_s3, dtype=float), mass_kg=mass, volume_m3=volume, slowing_time_s=slowing_time)
        pair_sum_energy_change += test_change
        field_population_energy_gain = -test_change if pair_id == collision_operator_state.electron_pair_id else None
        pair_moments.append(Eq59PairEnergyMoment(pair_id=pair_id, field_population_id=contribution.field_population_id, field_species_id=contribution.field_species_id, test_species_speed_energy_change_W=test_change, field_population_energy_gain_W=field_population_energy_gain, discrete_pitch_loss_power_W=pitch_loss, active=bool(contribution.active), operator_qualified=bool(contribution.qualified), limitation=contribution.limitation))
    pair_energy_scale = max(abs(total_energy_change), abs(pair_sum_energy_change), 1.0)
    pair_energy_error = abs(pair_sum_energy_change - total_energy_change) / pair_energy_scale
    pair_energy_identity_passed = bool(pair_energy_error <= flux_tolerance)
    projected_exact_source, discrete_source = _source_powers_W(speed_grid=speed_grid, source_coefficients=source_coefficients, source_speeds=source_speeds, mode_eta_integrals=eta_integrals, mass_kg=mass, volume_m3=volume)
    source_projection_error = abs(projected_exact_source - exact_source) / max(abs(projected_exact_source), abs(exact_source), 1.0)
    source_discretization_error = abs(discrete_source - projected_exact_source) / max(abs(discrete_source), abs(projected_exact_source), 1.0)
    source_projection_passed = bool(source_projection_error <= source_tolerance)
    source_discretization_passed = bool(source_discretization_error <= source_tolerance)
    pitch_loss_power = float(sum(moment.discrete_pitch_loss_power_W for moment in pair_moments))
    charge_exchange_sink_power, charge_exchange_sink_active = _charge_exchange_sink_power_W(speed_grid=speed_grid, modal_distribution=modal, mode_eta_integrals=eta_integrals, loss_frequency_s=charge_exchange_loss_frequency_s, mass_kg=mass, volume_m3=volume)
    energy_residual = discrete_source + total_energy_change - pitch_loss_power - charge_exchange_sink_power
    energy_scale = max(abs(discrete_source), abs(total_energy_change) + abs(pitch_loss_power) + abs(charge_exchange_sink_power), 1.0)
    energy_residual_relative = abs(energy_residual) / energy_scale
    discrete_energy_identity_passed = bool(energy_residual_relative <= energy_tolerance)
    electron_pair = collision_operator_state.contribution(collision_operator_state.electron_pair_id)
    electron_moment = next(moment for moment in pair_moments if moment.pair_id == collision_operator_state.electron_pair_id)
    if electron_moment.field_population_energy_gain_W is None:
        raise RuntimeError("Eq 59 electron pair field energy gain is unavailable")
    electron_heating = float(electron_moment.field_population_energy_gain_W)
    electron_pair_active = bool(electron_pair.active)
    electron_pair_operator_qualified = bool(electron_pair.qualified)
    electron_pair_numerically_valid = bool(pair_flux_identity_passed and pair_energy_identity_passed and discrete_energy_identity_passed and source_projection_passed and source_discretization_passed and np.isfinite(electron_heating))
    electron_pair_qualified = bool(electron_pair_numerically_valid and electron_pair_active and electron_pair_operator_qualified)
    full_collision_model_qualified = bool(collision_operator_state.qualified)
    numerically_valid = bool(pair_flux_identity_passed and pair_energy_identity_passed and discrete_energy_identity_passed and source_projection_passed and source_discretization_passed)
    qualified = bool(numerically_valid and electron_pair_qualified)
    failures: list[str] = []
    if not pair_flux_identity_passed:
        failures.append("pair speed fluxes do not reconstruct the active Eq 59 flux")
    if not pair_energy_identity_passed:
        failures.append("pair speed energy moments do not reconstruct the active Eq 59 speed energy moment")
    if not discrete_energy_identity_passed:
        failures.append("the discrete Eq 59 source collision and pitch loss energy identity failed")
    if not source_projection_passed:
        failures.append("the projected modal source power does not match the confined source power")
    if not source_discretization_passed:
        failures.append("the speed cell source energy does not represent the exact source energy")
    if not electron_pair_active:
        failures.append("electron collision pair is inactive")
    if not electron_pair_operator_qualified:
        failures.append(electron_pair.limitation or "electron collision pair is not qualified")
  
    return Eq59EnergyMomentState(
        fast_species_id=collision_operator_state.test_species.species_id,
        electron_pair_id=collision_operator_state.electron_pair_id,
        electron_heating_power_W=electron_heating,
        fast_ion_energy_change_from_electrons_W=float(electron_moment.test_species_speed_energy_change_W),
        pair_energy_moments=tuple(pair_moments),
        total_speed_operator_energy_change_W=total_energy_change,
        pair_sum_speed_operator_energy_change_W=pair_sum_energy_change,
        pair_energy_identity_relative_error=float(pair_energy_error),
        authoritative_confined_source_power_W=exact_source,
        projected_exact_source_power_W=projected_exact_source,
        discrete_source_power_W=discrete_source,
        source_projection_relative_error=float(source_projection_error),
        source_discretization_relative_error=float(source_discretization_error),
        discrete_pitch_loss_power_W=pitch_loss_power,
        charge_exchange_sink_power_W=charge_exchange_sink_power,
        charge_exchange_sink_active=charge_exchange_sink_active,
        energy_identity_includes_charge_exchange_sink=True,
        discrete_energy_residual_W=float(energy_residual),
        discrete_energy_residual_relative=float(energy_residual_relative),
        pair_flux_maximum_absolute_error=maximum_absolute_error,
        pair_flux_maximum_norm_scaled_error=float(maximum_norm_scaled_error),
        pair_flux_l2_relative_error=float(l2_relative_error),
        pair_flux_identity_relative_tolerance=flux_tolerance,
        energy_identity_relative_tolerance=energy_tolerance,
        source_power_relative_tolerance=source_tolerance,
        pair_flux_identity_passed=pair_flux_identity_passed,
        pair_energy_identity_passed=pair_energy_identity_passed,
        discrete_energy_identity_passed=discrete_energy_identity_passed,
        source_projection_passed=source_projection_passed,
        source_discretization_passed=source_discretization_passed,
        electron_pair_active=electron_pair_active,
        electron_pair_operator_qualified=electron_pair_operator_qualified,
        electron_pair_numerically_valid=electron_pair_numerically_valid,
        electron_pair_qualified=electron_pair_qualified,
        full_collision_model_qualified=full_collision_model_qualified,
        numerically_valid=numerically_valid,
        qualified=qualified,
        failure_reason=None if not failures else "; ".join(failures),
    )

def build_fast_ion_electron_heating_state(*, modal_result: ModalFBISResult, volume_m3: float, pair_flux_relative_tolerance: float = 1.0e-10, energy_residual_relative_tolerance: float = 1.0e-4, source_power_relative_tolerance: float = 1.0e-2) -> Eq59EnergyMomentState:
    """Build the Eq 59 energy moment state from a completed modal species result"""
    warm = modal_result.eq59_warm_start_state
    collision = modal_result.eq59_collision_operator_state
    if warm is None or collision is None or warm.mode_eta_integrals is None:
        raise ValueError("completed Eq 59 energy moments require the warm start and collision operator states")
    exact_source_power = float(sum(component.confined_birth_power_W for component in modal_result.component_results))
    charge_exchange_frequency = None if warm.charge_exchange_loss_frequency_s is None else np.asarray(warm.charge_exchange_loss_frequency_s, dtype=float)
  
    return build_eq59_energy_moment_state(speed_grid=modal_result.speed_grid, modal_distribution_f_j_v=modal_result.modal_distribution_f_j_v, source_coefficients_by_component_j=np.asarray(warm.source_coefficients_by_component_j, dtype=float), source_speeds_m_s=np.asarray(warm.source_speeds_m_s, dtype=float), active_eigenvalues_by_speed=np.asarray(warm.active_eigenvalues_by_speed, dtype=float), collision_operator_state=collision, mode_eta_integrals=np.asarray(warm.mode_eta_integrals, dtype=float), particle_mass_kg=modal_result.species.mass_kg, volume_m3=volume_m3, exact_confined_source_power_W=exact_source_power, charge_exchange_loss_frequency_s=charge_exchange_frequency, pair_flux_relative_tolerance=pair_flux_relative_tolerance, energy_residual_relative_tolerance=energy_residual_relative_tolerance, source_power_relative_tolerance=source_power_relative_tolerance)

def eq59_energy_moment_metadata(state: Eq59EnergyMomentState, prefix: str = "eq59") -> dict[str, object]:
    """Return flattened metadata for one fast ion Eq 59 energy state"""
    key = str(prefix).strip()
    pair_gain = {moment.pair_id: moment.field_population_energy_gain_W for moment in state.pair_energy_moments}
    pair_pitch = {moment.pair_id: moment.discrete_pitch_loss_power_W for moment in state.pair_energy_moments}
    return {
        f"{key}_fast_species_id": state.fast_species_id,
        f"{key}_electron_heating_power_W": state.electron_heating_power_W,
        f"{key}_electron_pair_id": state.electron_pair_id,
        f"{key}_electron_pair_test_species_energy_change_W": state.fast_ion_energy_change_from_electrons_W,
        f"{key}_field_population_energy_gain_W_by_pair": pair_gain,
        f"{key}_discrete_pitch_loss_power_W_by_pair": pair_pitch,
        f"{key}_total_speed_operator_energy_change_W": state.total_speed_operator_energy_change_W,
        f"{key}_pair_sum_speed_operator_energy_change_W": state.pair_sum_speed_operator_energy_change_W,
        f"{key}_pair_energy_identity_relative_error": state.pair_energy_identity_relative_error,
        f"{key}_authoritative_confined_source_power_W": state.authoritative_confined_source_power_W,
        f"{key}_projected_exact_source_power_W": state.projected_exact_source_power_W,
        f"{key}_discrete_source_power_W": state.discrete_source_power_W,
        f"{key}_source_projection_relative_error": state.source_projection_relative_error,
        f"{key}_source_discretization_relative_error": state.source_discretization_relative_error,
        f"{key}_discrete_pitch_loss_power_W": state.discrete_pitch_loss_power_W,
        f"{key}_charge_exchange_sink_power_W": state.charge_exchange_sink_power_W,
        f"{key}_charge_exchange_sink_active": state.charge_exchange_sink_active,
        f"{key}_energy_identity_includes_charge_exchange_sink": state.energy_identity_includes_charge_exchange_sink,
        f"{key}_discrete_energy_residual_W": state.discrete_energy_residual_W,
        f"{key}_discrete_energy_residual_relative": state.discrete_energy_residual_relative,
        f"{key}_pair_flux_maximum_absolute_error": state.pair_flux_maximum_absolute_error,
        f"{key}_pair_flux_maximum_norm_scaled_error": state.pair_flux_maximum_norm_scaled_error,
        f"{key}_pair_flux_l2_relative_error": state.pair_flux_l2_relative_error,
        f"{key}_pair_flux_identity_relative_tolerance": state.pair_flux_identity_relative_tolerance,
        f"{key}_energy_identity_relative_tolerance": state.energy_identity_relative_tolerance,
        f"{key}_source_power_relative_tolerance": state.source_power_relative_tolerance,
        f"{key}_pair_flux_identity_passed": state.pair_flux_identity_passed,
        f"{key}_pair_energy_identity_passed": state.pair_energy_identity_passed,
        f"{key}_discrete_energy_identity_passed": state.discrete_energy_identity_passed,
        f"{key}_source_projection_passed": state.source_projection_passed,
        f"{key}_source_discretization_passed": state.source_discretization_passed,
        f"{key}_electron_pair_active": state.electron_pair_active,
        f"{key}_electron_pair_operator_qualified": state.electron_pair_operator_qualified,
        f"{key}_electron_pair_numerically_valid": state.electron_pair_numerically_valid,
        f"{key}_electron_pair_qualified": state.electron_pair_qualified,
        f"{key}_full_collision_model_qualified": state.full_collision_model_qualified,
        f"{key}_energy_moment_qualification_scope": "active_electron_pair_and_discrete_Eq59_numerical_identities",
        f"{key}_full_collision_model_qualification_role": "separate_model_applicability_diagnostic",
        f"{key}_numerically_valid": state.numerically_valid,
        f"{key}_qualified": state.qualified,
        f"{key}_failure_reason": state.failure_reason,
        f"{key}_model": state.model,
    }

fast_ion_electron_heating_metadata = eq59_energy_moment_metadata
FastIonElectronHeatingState = Eq59EnergyMomentState

__all__ = [
    "Eq59PairEnergyMoment",
    "Eq59EnergyMomentState",
    "FastIonElectronHeatingState",
    "build_eq59_energy_moment_state",
    "build_fast_ion_electron_heating_state",
    "eq59_energy_moment_metadata",
    "fast_ion_electron_heating_metadata",
]
