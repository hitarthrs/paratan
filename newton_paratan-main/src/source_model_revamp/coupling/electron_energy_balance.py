"""Assemble represented electron heating end loss and ambipolar power terms"""
from __future__ import annotations
from collections.abc import Mapping
from dataclasses import dataclass
from types import MappingProxyType
import numpy as np
from source_model_revamp.constants import ELECTRON_CHARGE_C
from source_model_revamp.fbis.modal.species_state import FastIonSystemState

@dataclass(frozen=True)
class ElectronEnergyBalanceState:
    """Validated electron energy balance at one solved operating point"""
    electron_temperature_J: float
    fast_ion_heating_power_W_by_species: Mapping[str, float]
    external_electron_heating_W: float
    total_electron_heating_W: float
    electron_wall_kinetic_power_W: float
    ambipolar_ion_acceleration_power_W: float
    independent_ion_acceleration_power_W: float
    ambipolar_power_identity_relative_error: float
    total_electron_loss_W: float
    residual_W: float
    relative_residual: float
    current_balance_relative_residual: float
    total_current_numerical_fixed_point_converged: bool
    total_current_end_model_applicable: bool
    fast_ion_transfer_available: bool
    fast_ion_transfer_numerically_valid: bool
    fast_ion_transfer_qualified: bool
    ambipolar_power_identity_passed: bool
    represented_terms_available: bool
    represented_terms_evaluable: bool
    represented_terms_qualified: bool
    missing_physics_terms: tuple[str, ...]
    failure_reason: str | None
    model: str = "classical_represented_confined_plasma_electron_energy_balance"

    def __post_init__(self) -> None:
        """Validate component powers identities and qualification consistency"""
        fast = dict(self.fast_ion_heating_power_W_by_species)
        if not self.model:
            raise ValueError("electron energy balance model must be nonempty")
        if any(not key for key in fast):
            raise ValueError("electron energy balance species identifiers must be nonempty")
        scalar_values = (self.electron_temperature_J, self.external_electron_heating_W, self.total_electron_heating_W, self.electron_wall_kinetic_power_W, self.ambipolar_ion_acceleration_power_W, self.independent_ion_acceleration_power_W, self.ambipolar_power_identity_relative_error, self.total_electron_loss_W, self.residual_W, self.relative_residual, self.current_balance_relative_residual, *fast.values())
        if any(not np.isfinite(float(value)) for value in scalar_values):
            raise ValueError("electron energy balance values must be finite")
        if self.electron_temperature_J <= 0.0:
            raise ValueError("electron temperature must be positive")
        if self.external_electron_heating_W < 0.0 or self.electron_wall_kinetic_power_W < 0.0 or self.ambipolar_ion_acceleration_power_W < 0.0 or self.independent_ion_acceleration_power_W < 0.0 or self.total_electron_loss_W < 0.0:
            raise ValueError("electron heating input and loss powers must be nonnegative")
        if self.ambipolar_power_identity_relative_error < 0.0 or self.current_balance_relative_residual < 0.0:
            raise ValueError("electron energy balance errors must be nonnegative")
        heating_sum = float(sum(fast.values()) + self.external_electron_heating_W)
        loss_sum = float(self.electron_wall_kinetic_power_W + self.ambipolar_ion_acceleration_power_W)
        expected_residual = heating_sum - loss_sum
        expected_relative_residual = expected_residual / max(abs(heating_sum), abs(loss_sum), 1.0)
        consistency_scale = max(abs(heating_sum), abs(loss_sum), abs(self.total_electron_heating_W), abs(self.total_electron_loss_W), 1.0)
        consistency_tolerance = 4096.0 * np.finfo(float).eps * consistency_scale
        if abs(self.total_electron_heating_W - heating_sum) > consistency_tolerance:
            raise ValueError("total electron heating is inconsistent with component powers")
        if abs(self.total_electron_loss_W - loss_sum) > consistency_tolerance:
            raise ValueError("total electron loss is inconsistent with component powers")
        if abs(self.residual_W - expected_residual) > consistency_tolerance:
            raise ValueError("electron energy residual is inconsistent with heating and loss powers")
        relative_tolerance = 4096.0 * np.finfo(float).eps * max(abs(expected_relative_residual), abs(self.relative_residual), 1.0)
        if abs(self.relative_residual - expected_relative_residual) > relative_tolerance:
            raise ValueError("relative electron energy residual is inconsistent with heating and loss powers")
        if self.represented_terms_qualified and not self.represented_terms_evaluable:
            raise ValueError("qualified electron energy terms must be evaluable")
        object.__setattr__(self, "fast_ion_heating_power_W_by_species", MappingProxyType(fast))

def _positive_finite_scalar(value: float, name: str) -> float:
    """Return one validated positive finite scalar"""
    scalar = float(value)
    if not np.isfinite(scalar) or scalar <= 0.0:
        raise ValueError(f"{name} must be positive and finite")
   
    return scalar

def _nonnegative_finite_scalar(value: float, name: str) -> float:
    """Return one validated nonnegative finite scalar"""
    scalar = float(value)
    if not np.isfinite(scalar) or scalar < 0.0:
        raise ValueError(f"{name} must be finite and nonnegative")
 
    return scalar

def _fast_ion_transfer_state(system: FastIonSystemState) -> tuple[dict[str, float], bool, bool, bool, list[str]]:
    """Collect direct Eq 59 electron energy transfer powers and qualification flags"""
    powers: dict[str, float] = {}
    available = True
    numerically_valid = True
    qualified = True
    failures: list[str] = []
    for species_id, species_state in system.species_states.items():
        if not species_state.active:
            powers[species_id] = 0.0
            continue
        modal = species_state.modal_result
        energy_state = None if modal is None else modal.fast_ion_electron_heating_state
        if energy_state is None:
            available = False
            numerically_valid = False
            qualified = False
            failures.append(f"fast {species_id} electron transfer state is unavailable")
            continue
        powers[species_id] = float(energy_state.electron_heating_power_W)
        if not energy_state.electron_pair_numerically_valid:
            numerically_valid = False
            failures.append(f"fast {species_id} electron transfer is not numerically valid")
        modal_metadata = getattr(modal, "metadata", {})
        speed_resolution_assessed = modal_metadata.get("modal_rosenbluth_speed_electron_heating_resolution_assessed") is True
        speed_resolution_converged = modal_metadata.get("modal_rosenbluth_speed_electron_heating_resolution_converged") is True
        speed_domain_assessed = modal_metadata.get("modal_rosenbluth_speed_electron_heating_domain_assessed") is True
        speed_domain_converged = modal_metadata.get("modal_rosenbluth_speed_electron_heating_domain_converged") is True
        if not energy_state.qualified:
            qualified = False
            failures.append(f"fast {species_id} electron transfer state is not fully qualified")
        if not speed_resolution_assessed or not speed_resolution_converged:
            qualified = False
            failures.append(f"fast {species_id} electron transfer speed resolution is not qualified")
        if not speed_domain_assessed or not speed_domain_converged:
            qualified = False
            failures.append(f"fast {species_id} electron transfer speed domain is not qualified")

    return powers, available, numerically_valid, qualified, failures

def _independent_ion_acceleration_power_W(fast_system: FastIonSystemState) -> tuple[float, bool, list[str]]:
    """Reconstruct ion electrostatic acceleration power from midplane and wall losses"""
    total = 0.0
    available = True
    failures: list[str] = []
    for species_id, state in fast_system.species_states.items():
        if state.modal_result is not None:
            wall = state.modal_result.ion_wall_power_loss_W
            if wall is None:
                available = False
                failures.append(f"fast {species_id} ion wall power is unavailable")
            else:
                # Wall minus midplane kinetic power isolates electrostatic acceleration
                total += float(wall) - float(state.modal_result.ion_midplane_kinetic_power_loss_W)
        elif state.prompt_only_particle_loss_rate_s > 0.0:
            if state.prompt_only_ion_wall_power_loss_W is None:
                available = False
                failures.append(f"prompt {species_id} ion wall power is unavailable")
            else:
                total += float(state.prompt_only_ion_wall_power_loss_W) - float(state.prompt_only_midplane_kinetic_power_loss_W)
    if total < 0.0:
        scale = max(abs(total), 1.0)
        if total < -256.0 * np.finfo(float).eps * scale:
            available = False
            failures.append("independent ion acceleration power is materially negative")
        total = max(total, 0.0)

    return float(total), available, failures

def build_electron_energy_balance_state(*, electron_temperature_J: float, fast_ion_system_state: FastIonSystemState, external_electron_heating_W: float=0.0, ambipolar_power_relative_tolerance: float=0.001, current_balance_relative_tolerance: float=0.001, terminal_ion_current_A: float | None=None, terminal_independent_ion_acceleration_power_W: float | None=None) -> ElectronEnergyBalanceState:
    """Assemble R = P_heat − P_loss with P_amb = (I_i / e) E_barrier"""
    temperature = _positive_finite_scalar(electron_temperature_J, "electron_temperature_J")
    external = _nonnegative_finite_scalar(external_electron_heating_W, "external_electron_heating_W")
    ambipolar_tolerance = _positive_finite_scalar(ambipolar_power_relative_tolerance, "ambipolar_power_relative_tolerance")
    current_tolerance = _positive_finite_scalar(current_balance_relative_tolerance, "current_balance_relative_tolerance")
    current_balance = fast_ion_system_state.shared_current_balance
    if current_balance is None:
        raise ValueError("electron energy balance requires the shared total current balance")
    electron_balance = current_balance.electron_balance
    barrier = float(electron_balance.barrier_energy_J)
    ion_current = float(current_balance.ion_current_loss.total_ion_current_A) if terminal_ion_current_A is None else _nonnegative_finite_scalar(terminal_ion_current_A, "terminal_ion_current_A")
    electron_wall_power = float(electron_balance.electron_wall_power_loss_W)
    if np.isnan(barrier) or barrier < 0.0 or not np.isfinite(ion_current) or ion_current < 0.0 or not np.isfinite(electron_wall_power) or electron_wall_power < 0.0:
        raise ValueError("electron energy balance current and wall terms must be nonnegative and may use an infinite zero current barrier")
    # Convert matched ion current to particle rate before applying the barrier energy
    if ion_current == 0.0:
        ambipolar_power = 0.0
    elif np.isfinite(barrier):
        ambipolar_power = float((ion_current / ELECTRON_CHARGE_C) * barrier)
    else:
        raise ValueError("a nonzero ion current requires a finite ambipolar barrier")
    if terminal_independent_ion_acceleration_power_W is None:
        independent_power, independent_available, independent_failures = _independent_ion_acceleration_power_W(fast_ion_system_state)
    else:
        independent_power = _nonnegative_finite_scalar(terminal_independent_ion_acceleration_power_W, "terminal_independent_ion_acceleration_power_W")
        independent_available = True
        independent_failures = []
    ambipolar_identity_error = abs(independent_power - ambipolar_power) / max(abs(independent_power), abs(ambipolar_power), 1.0)
    ambipolar_identity_passed = bool(independent_available and ambipolar_identity_error <= ambipolar_tolerance)
    target_current = ion_current
    achieved_current = float(electron_balance.electron_current_A)
    current_error = abs(achieved_current - target_current) / max(abs(achieved_current), abs(target_current), 1.0e-300)
    total_current_numerical = bool(fast_ion_system_state.metadata.get("total_current_numerical_fixed_point_converged", fast_ion_system_state.metadata.get("total_current_balance_numerical_fixed_point_converged", False)))
    total_current_end_model_applicable = bool(fast_ion_system_state.metadata.get("equivalent_end_model_applicable", fast_ion_system_state.metadata.get("total_current_balance_end_model_applicable", False)))
    fast_powers, fast_available, fast_numerical, fast_qualified, fast_failures = _fast_ion_transfer_state(fast_ion_system_state)
    total_heating = float(sum(fast_powers.values()) + external)
    total_loss = float(electron_wall_power + ambipolar_power)
    residual = total_heating - total_loss
    relative_residual = residual / max(abs(total_heating), abs(total_loss), 1.0)
    represented_terms_available = bool(fast_available and independent_available)
    represented_terms_evaluable = bool(represented_terms_available and fast_numerical and total_current_numerical and current_error <= current_tolerance and ambipolar_identity_passed)
    represented_terms_qualified = bool(represented_terms_evaluable and fast_qualified and total_current_end_model_applicable)
    failures = [*fast_failures, *independent_failures]
    if current_error > current_tolerance:
        failures.append("total electron and ion current residual exceeds tolerance")
    if not total_current_numerical:
        failures.append("total ion and electron current numerical fixed point is not converged")
    if not total_current_end_model_applicable:
        failures.append("equivalent two end electron current model is not applicable")
    if not ambipolar_identity_passed:
        failures.append("ambipolar ion acceleration power identity failed")
    if not represented_terms_available:
        failures.append("one or more represented electron energy terms are unavailable")
    missing_terms = ("fusion_alpha_heating", "DD_charged_product_heating", "bremsstrahlung_radiation", "synchrotron_radiation", "electron_radial_or_turbulent_transport", "beam_ionization_electron_energy_cost", "neutral_recycling", "expander_population_maintenance", "external_bias_circuit_power")

    return ElectronEnergyBalanceState(electron_temperature_J=temperature, fast_ion_heating_power_W_by_species=fast_powers, external_electron_heating_W=external, total_electron_heating_W=total_heating, electron_wall_kinetic_power_W=electron_wall_power, ambipolar_ion_acceleration_power_W=ambipolar_power, independent_ion_acceleration_power_W=independent_power, ambipolar_power_identity_relative_error=float(ambipolar_identity_error), total_electron_loss_W=total_loss, residual_W=float(residual), relative_residual=float(relative_residual), current_balance_relative_residual=float(current_error), total_current_numerical_fixed_point_converged=total_current_numerical, total_current_end_model_applicable=total_current_end_model_applicable, fast_ion_transfer_available=fast_available, fast_ion_transfer_numerically_valid=fast_numerical, fast_ion_transfer_qualified=fast_qualified, ambipolar_power_identity_passed=ambipolar_identity_passed, represented_terms_available=represented_terms_available, represented_terms_evaluable=represented_terms_evaluable, represented_terms_qualified=represented_terms_qualified, missing_physics_terms=missing_terms, failure_reason=None if not failures else '; '.join(dict.fromkeys(failures)))

def electron_energy_balance_metadata(state: ElectronEnergyBalanceState) -> dict[str, object]:
    """Return stable metadata for the represented electron energy balance"""
    return {'physical_electron_energy_equation_active': True, 'electron_energy_balance_model': state.model, 'electron_energy_balance_scope': 'classical_represented_confined_plasma', 'electron_energy_balance_fast_ion_heating_power_W_by_species': dict(state.fast_ion_heating_power_W_by_species), 'electron_energy_balance_external_electron_heating_W': state.external_electron_heating_W, 'electron_energy_balance_total_heating_W': state.total_electron_heating_W, 'electron_energy_balance_electron_wall_kinetic_power_W': state.electron_wall_kinetic_power_W, 'electron_energy_balance_ambipolar_ion_acceleration_power_W': state.ambipolar_ion_acceleration_power_W, 'electron_energy_balance_independent_ion_acceleration_power_W': state.independent_ion_acceleration_power_W, 'electron_energy_balance_ambipolar_power_identity_relative_error': state.ambipolar_power_identity_relative_error, 'electron_energy_balance_ambipolar_power_identity_passed': state.ambipolar_power_identity_passed, 'electron_energy_balance_total_loss_W': state.total_electron_loss_W, 'electron_energy_balance_residual_W': state.residual_W, 'electron_energy_balance_relative_residual': state.relative_residual, 'electron_energy_balance_current_relative_residual': state.current_balance_relative_residual, 'electron_energy_balance_total_current_numerical_fixed_point_converged': state.total_current_numerical_fixed_point_converged, 'electron_energy_balance_total_current_end_model_applicable': state.total_current_end_model_applicable, 'fast_ion_electron_transfer_available': state.fast_ion_transfer_available, 'fast_ion_electron_transfer_numerically_valid': state.fast_ion_transfer_numerically_valid, 'fast_ion_electron_transfer_qualified': state.fast_ion_transfer_qualified, 'represented_electron_energy_terms_available': state.represented_terms_available, 'represented_electron_energy_terms_evaluable': state.represented_terms_evaluable, 'represented_electron_energy_terms_qualified': state.represented_terms_qualified, 'electron_energy_balance_missing_physics_terms': list(state.missing_physics_terms), 'electron_energy_balance_failure_reason': state.failure_reason, 'fusion_alpha_heating_included': False, 'radiative_electron_losses_included': False, 'beam_ionization_electron_energy_loss_included': False, 'radial_electron_transport_included': False, 'global_plasma_power_balance_claimed': False, 'global_device_power_balance_claimed': False}

__all__ = [
    "ElectronEnergyBalanceState",
    "build_electron_energy_balance_state",
    "electron_energy_balance_metadata",
]
