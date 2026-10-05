"""Build pair resolved speed dependent collision coefficients for Egedal Eq 59"""
from __future__ import annotations
from collections.abc import Mapping
from dataclasses import dataclass
from types import MappingProxyType
import numpy as np
from source_model_revamp.constants import ELECTRON_CHARGE_C, VACUUM_PERMITTIVITY_F_PER_M
from source_model_revamp.fbis.pairwise_collisions import PairwiseCollisionRecord, PairwiseCollisionState
from source_model_revamp.fbis.species import IonSpecies
from source_model_revamp.fbis.modal.types import ModalRosenbluthCoefficients
from source_model_revamp.fbis.modal.eq59.maxwellian_pairs import maxwellian_rosenbluth_coefficients
from source_model_revamp.fbis.modal.eq59.cross_species import ExternalFastFieldCollisionState, external_fast_field_rosenbluth_coefficients
from source_model_revamp.fbis.velocity_space_grid import SpeedGrid

def _validated_nonnegative_array(values: np.ndarray, shape: tuple[int, ...], name: str) -> np.ndarray:
    """Validate one speed array and clip negative roundoff to zero"""
    array = np.asarray(values, dtype=float)
    if array.shape != shape or np.any(~np.isfinite(array)):
        raise ValueError(f"{name} must contain one finite value per speed cell")
    scale = max(float(np.max(np.abs(array))) if array.size else 0.0, 1.0)
    tolerance = 256.0 * np.finfo(float).eps * scale
    if float(np.min(array)) < -tolerance:
        raise ValueError(f"{name} contains material negative values")
   
    return np.where(array < 0.0, 0.0, array)

def _positive_finite_scalar(value: float, name: str) -> float:
    """Validate one strictly positive finite scalar"""
    scalar = float(value)
    if not np.isfinite(scalar) or scalar <= 0.0:
        raise ValueError(f"{name} must be positive and finite")
   
    return scalar

def _nonnegative_finite_scalar(value: float, name: str) -> float:
    """Validate one finite scalar that may be zero"""
    scalar = float(value)
    if not np.isfinite(scalar) or scalar < 0.0:
        raise ValueError(f"{name} must be finite and nonnegative")
   
    return scalar

@dataclass(frozen=True)
class Eq59PairOperatorContribution:
    """Store one ordered field population contribution to the Eq 59 speed operator
    
    Coefficient arrays have shape (n_speed,) on the test species speed grid
    Drag and pitch arrays carry velocity cubed units and diffusion carries velocity fourth units
    """
    pair_id: str
    test_species_id: str
    field_population_id: str
    field_species_id: str
    field_distribution_model: str
    field_mass_kg: float
    field_charge_number: float
    field_density_m3: float
    field_temperature_J_or_none: float | None
    coulomb_log_value: float
    rosenbluth_g_weight_m3: float
    rosenbluth_h_weight_m3: float
    g_scale_velocity_cubed_m3_s3: float
    net_drag_scale_velocity_cubed_m3_s3: float
    h_tilde: np.ndarray
    g_tilde_1: np.ndarray
    g_tilde_2: np.ndarray
    drag_velocity_cubed_m3_s3: np.ndarray
    energy_diffusion_velocity_fourth_m4_s4: np.ndarray
    pitch_scattering_velocity_cubed_m3_s3: np.ndarray
    active_operator_terms: tuple[str, ...]
    active: bool
    qualified: bool
    limitation: str | None

    def __post_init__(self) -> None:
        """Validate pair identifiers scalar inputs coefficient shapes and active term rules"""
        text_values = (self.pair_id, self.test_species_id, self.field_population_id, self.field_species_id, self.field_distribution_model)
        if any(not value for value in text_values):
            raise ValueError("Eq 59 pair operator text fields must be nonempty")
        if not self.active_operator_terms or any(not str(value).strip() for value in self.active_operator_terms):
            raise ValueError("Eq 59 pair operator terms must be nonempty")
        if not np.isfinite(self.field_mass_kg) or self.field_mass_kg <= 0.0:
            raise ValueError("Eq 59 field mass must be positive and finite")
        if not np.isfinite(self.field_charge_number) or self.field_charge_number == 0.0:
            raise ValueError("Eq 59 field charge must be finite and nonzero")
        if not np.isfinite(self.field_density_m3) or self.field_density_m3 < 0.0:
            raise ValueError("Eq 59 field density must be finite and nonnegative")
        if self.field_temperature_J_or_none is not None and (not np.isfinite(self.field_temperature_J_or_none) or self.field_temperature_J_or_none <= 0.0):
            raise ValueError("Eq 59 field temperature must be positive and finite when available")
        for name, value in (("coulomb_log_value", self.coulomb_log_value), ("rosenbluth_g_weight_m3", self.rosenbluth_g_weight_m3), ("rosenbluth_h_weight_m3", self.rosenbluth_h_weight_m3), ("g_scale_velocity_cubed_m3_s3", self.g_scale_velocity_cubed_m3_s3), ("net_drag_scale_velocity_cubed_m3_s3", self.net_drag_scale_velocity_cubed_m3_s3)):
            if not np.isfinite(value) or value < 0.0:
                raise ValueError(f"{name} must be finite and nonnegative")
        shape = np.asarray(self.h_tilde, dtype=float).shape
        if len(shape) != 1 or shape[0] < 1:
            raise ValueError("Eq 59 pair operator arrays must be nonempty one dimensional arrays")
        for name, values in (("h_tilde", self.h_tilde), ("g_tilde_1", self.g_tilde_1), ("g_tilde_2", self.g_tilde_2), ("drag_velocity_cubed_m3_s3", self.drag_velocity_cubed_m3_s3), ("energy_diffusion_velocity_fourth_m4_s4", self.energy_diffusion_velocity_fourth_m4_s4), ("pitch_scattering_velocity_cubed_m3_s3", self.pitch_scattering_velocity_cubed_m3_s3)):
            _validated_nonnegative_array(values, shape, name)
        if self.active and self.field_population_id != "electrons" and self.g_scale_velocity_cubed_m3_s3 <= 0.0:
            raise ValueError("active Eq 59 ion pairs require a positive g scale")
        if self.active and self.field_population_id == "electrons" and self.active_operator_terms != ("analytic_electron_drag_normalization",):
            raise ValueError("the active electron pair must use the analytic electron drag normalization")

@dataclass(frozen=True)
class Eq59CollisionOperatorState:
    """Store the summed Eq 59 collision coefficients for one fast test species
    
    The electron pair supplies the analytic τ_s normalization
    Explicit ion pair arrays sum thermal fast self and optional fast cross species contributions
    """
    test_species: IonSpecies
    speed_grid: SpeedGrid
    spitzer_slowing_down_time_s: float
    electron_pair_id: str
    pair_contributions_by_id: Mapping[str, Eq59PairOperatorContribution]
    ion_drag_velocity_cubed_m3_s3: np.ndarray
    ion_energy_diffusion_velocity_fourth_m4_s4: np.ndarray
    ion_pitch_scattering_velocity_cubed_m3_s3: np.ndarray
    fast_self_physical_density_m3: float
    fast_self_rosenbluth_density_normalization_m3: float
    operator_model: str
    density_convention: str
    temperature_convention: str
    coulomb_log_convention: str
    qualified: bool
    limitations: tuple[str, ...]

    def __post_init__(self) -> None:
        """Validate the pair map and require each summed coefficient to equal its pair contributions"""
        contributions = dict(self.pair_contributions_by_id)
        if not contributions:
            raise ValueError("Eq 59 collision operator requires pair contributions")
        if not self.electron_pair_id or self.electron_pair_id not in contributions:
            raise ValueError("Eq 59 collision operator requires the electron pair identity")
        electron_contribution = contributions[self.electron_pair_id]
        if electron_contribution.field_population_id != "electrons" or not electron_contribution.active or electron_contribution.active_operator_terms != ("analytic_electron_drag_normalization",):
            raise ValueError("Eq 59 electron pair must be active through the analytic electron drag normalization")
        if any(not value for value in (self.operator_model, self.density_convention, self.temperature_convention, self.coulomb_log_convention)):
            raise ValueError("Eq 59 collision operator convention fields must be nonempty")
        _positive_finite_scalar(self.spitzer_slowing_down_time_s, "spitzer_slowing_down_time_s")
        _nonnegative_finite_scalar(self.fast_self_physical_density_m3, "fast_self_physical_density_m3")
        _nonnegative_finite_scalar(self.fast_self_rosenbluth_density_normalization_m3, "fast_self_rosenbluth_density_normalization_m3")
        shape = np.asarray(self.speed_grid.centers_m_s, dtype=float).shape
        drag = _validated_nonnegative_array(self.ion_drag_velocity_cubed_m3_s3, shape, "ion_drag_velocity_cubed_m3_s3")
        diffusion = _validated_nonnegative_array(self.ion_energy_diffusion_velocity_fourth_m4_s4, shape, "ion_energy_diffusion_velocity_fourth_m4_s4")
        pitch = _validated_nonnegative_array(self.ion_pitch_scattering_velocity_cubed_m3_s3, shape, "ion_pitch_scattering_velocity_cubed_m3_s3")
        sum_drag = np.zeros(shape, dtype=float)
        sum_diffusion = np.zeros(shape, dtype=float)
        sum_pitch = np.zeros(shape, dtype=float)
        for pair_id, contribution in contributions.items():
            if pair_id != contribution.pair_id:
                raise ValueError("Eq 59 collision operator keys must match pair identifiers")
            if contribution.test_species_id != self.test_species.species_id:
                raise ValueError("Eq 59 pair test species must match the operator state")
            if np.asarray(contribution.h_tilde, dtype=float).shape != shape:
                raise ValueError("Eq 59 pair operator speed arrays must match the operator grid")
            sum_drag += np.asarray(contribution.drag_velocity_cubed_m3_s3, dtype=float)
            sum_diffusion += np.asarray(contribution.energy_diffusion_velocity_fourth_m4_s4, dtype=float)
            sum_pitch += np.asarray(contribution.pitch_scattering_velocity_cubed_m3_s3, dtype=float)
        for name, total, expected in (("ion drag", drag, sum_drag), ("ion energy diffusion", diffusion, sum_diffusion), ("ion pitch scattering", pitch, sum_pitch)):
            if not np.allclose(total, expected, rtol=4.0e-15, atol=0.0):
                raise ValueError(f"Eq 59 total {name} array does not equal the sum of pair contributions")
        object.__setattr__(self, "pair_contributions_by_id", MappingProxyType(contributions))

    def contribution(self, pair_id: str) -> Eq59PairOperatorContribution:
        """Return one ordered pair contribution"""
        try:
            return self.pair_contributions_by_id[str(pair_id)]
        except KeyError as exc:
            raise ValueError(f"Eq 59 pair contribution {pair_id!r} is unavailable") from exc

    def structure_fingerprint(self) -> tuple[object, ...]:
        """Return immutable operator structure identity for warm start checks"""
        pair_values = tuple((contribution.pair_id, contribution.field_population_id, contribution.field_species_id, contribution.field_distribution_model, contribution.field_mass_kg, contribution.field_charge_number, contribution.active_operator_terms, contribution.active) for contribution in sorted(self.pair_contributions_by_id.values(), key=lambda item: item.pair_id))
        return (self.operator_model, self.test_species.species_id, self.test_species.mass_kg, self.test_species.charge_number, pair_values)

    def closure_fingerprint(self) -> tuple[object, ...]:
        """Return closure values while excluding endogenous fast self density and shape"""
        pair_values: list[object] = []
        for contribution in sorted(self.pair_contributions_by_id.values(), key=lambda item: item.pair_id):
            common = (contribution.pair_id, contribution.field_temperature_J_or_none, contribution.coulomb_log_value)
            if contribution.field_population_id == "fast_self":
                pair_values.append(common + ("endogenous_fast_self_density_scale_and_shape",))
            else:
                pair_values.append(common + (contribution.field_density_m3, contribution.rosenbluth_g_weight_m3, contribution.rosenbluth_h_weight_m3, contribution.g_scale_velocity_cubed_m3_s3, contribution.net_drag_scale_velocity_cubed_m3_s3))
      
        return (self.spitzer_slowing_down_time_s, tuple(pair_values))

def _zero_contribution(record: PairwiseCollisionRecord, speed: np.ndarray, *, active: bool = False) -> Eq59PairOperatorContribution:
    """Build an inactive or analytic only pair record with zero explicit speed coefficients"""
    zeros = np.zeros_like(speed)
    
    return Eq59PairOperatorContribution(pair_id=record.pair_id, test_species_id=record.test_species_id, field_population_id=record.field_population_id, field_species_id=record.field_species_id, field_distribution_model=record.field_distribution_model, field_mass_kg=record.field_mass_kg, field_charge_number=record.field_charge_number, field_density_m3=record.field_density_m3, field_temperature_J_or_none=record.field_temperature_J_or_none, coulomb_log_value=float(record.coulomb_log_state.value or 0.0), rosenbluth_g_weight_m3=float(record.rosenbluth_g_weight_m3 or 0.0), rosenbluth_h_weight_m3=float(record.rosenbluth_h_weight_m3 or 0.0), g_scale_velocity_cubed_m3_s3=0.0, net_drag_scale_velocity_cubed_m3_s3=0.0, h_tilde=zeros.copy(), g_tilde_1=zeros.copy(), g_tilde_2=zeros.copy(), drag_velocity_cubed_m3_s3=zeros.copy(), energy_diffusion_velocity_fourth_m4_s4=zeros.copy(), pitch_scattering_velocity_cubed_m3_s3=zeros.copy(), active_operator_terms=record.active_operator_terms, active=bool(active), qualified=record.applicability, limitation=record.limitation)

def _pair_functions(*, record: PairwiseCollisionRecord, speed: np.ndarray, fast_self_coefficients: ModalRosenbluthCoefficients | None) -> ModalRosenbluthCoefficients:
    """Select density normalized Rosenbluth functions for one active ion field population"""
    if record.field_population_id == "fast_self":
        if fast_self_coefficients is None:
            raise ValueError("active fast self pair requires density normalized Rosenbluth functions")
        coefficients = fast_self_coefficients
    elif record.field_population_id.startswith("thermal_"):
        if record.field_temperature_J_or_none is None:
            raise ValueError("active thermal pair requires a field temperature")
        coefficients = maxwellian_rosenbluth_coefficients(speed_m_s=speed, field_temperature_J=record.field_temperature_J_or_none, field_mass_kg=record.field_mass_kg)
    else:
        raise ValueError(f"unsupported active Eq 59 ion field population {record.field_population_id!r}")
    shape = speed.shape
    _validated_nonnegative_array(coefficients.h_tilde, shape, "h_tilde")
    _validated_nonnegative_array(coefficients.g_tilde_1, shape, "g_tilde_1")
    _validated_nonnegative_array(coefficients.g_tilde_2, shape, "g_tilde_2")
  
    return coefficients

def build_eq59_collision_operator_state(*, pairwise_collision_state: PairwiseCollisionState, speed_grid: SpeedGrid, spitzer_slowing_down_time_s: float, fast_self_rosenbluth_coefficients: ModalRosenbluthCoefficients | None, fast_self_physical_density_m3: float, external_fast_field_states: tuple[ExternalFastFieldCollisionState, ...] = ()) -> Eq59CollisionOperatorState:
    """Build the active pair resolved speed operator for one fast species
    
    For each ion field pair the implementation forms
    g_scale = τ_s Γ W_g
    drag_scale = τ_s Γ (W_h − W_g)
    then applies h_tilde g_tilde_1 and g_tilde_2 to drag pitch scattering and energy diffusion
    """
    tau_s = _positive_finite_scalar(spitzer_slowing_down_time_s, "spitzer_slowing_down_time_s")
    fast_density = _nonnegative_finite_scalar(fast_self_physical_density_m3, "fast_self_physical_density_m3")
    speed = np.asarray(speed_grid.centers_m_s, dtype=float)
    if speed.ndim != 1 or speed.size < 2 or np.any(~np.isfinite(speed)) or np.any(speed <= 0.0):
        raise ValueError("Eq 59 speed grid centers must be positive and finite")
    species = pairwise_collision_state.test_species
    gamma = species.charge_number**4 * ELECTRON_CHARGE_C**4 / (4.0 * np.pi * VACUUM_PERMITTIVITY_F_PER_M**2 * species.mass_kg**2)
    if not np.isfinite(gamma) or gamma <= 0.0:
        raise ValueError("Rosenbluth collision coefficient must be positive and finite")
    contributions: dict[str, Eq59PairOperatorContribution] = {}
    total_drag = np.zeros_like(speed)
    total_diffusion = np.zeros_like(speed)
    total_pitch = np.zeros_like(speed)
    electron_pair_id = ""
    fast_self_normalization = 0.0
    # The electron pair fixes τ_s while ion pairs supply the explicit h and g terms
    for record in pairwise_collision_state.records_by_pair_id.values():
        if record.field_population_id == "electrons":
            electron_pair_id = record.pair_id
            contributions[record.pair_id] = _zero_contribution(record, speed, active=record.active_in_current_eq59_operator)
            continue
        if not record.active_in_current_eq59_operator or record.field_density_m3 <= 0.0:
            contributions[record.pair_id] = _zero_contribution(record, speed)
            continue
        if record.coulomb_log_state.value is None or record.rosenbluth_g_weight_m3 is None or record.rosenbluth_h_weight_m3 is None:
            raise ValueError("active Eq 59 ion pairs require Coulomb logs and Rosenbluth weights")
        if record.field_population_id == "fast_self":
            density_tolerance = 256.0 * np.finfo(float).eps * max(record.field_density_m3, fast_density, 1.0)
            if abs(record.field_density_m3 - fast_density) > density_tolerance:
                raise ValueError("fast self pair density must match the physical fast species density")
        coefficients = _pair_functions(record=record, speed=speed, fast_self_coefficients=fast_self_rosenbluth_coefficients)
        g_scale = tau_s * gamma * float(record.rosenbluth_g_weight_m3)
        drag_scale = tau_s * gamma * (float(record.rosenbluth_h_weight_m3) - float(record.rosenbluth_g_weight_m3))
        if g_scale <= 0.0 or drag_scale <= 0.0:
            raise ValueError("active Eq 59 ion pair scales must be positive")
        expected_ratio = species.mass_kg / record.field_mass_kg
        if not np.isclose(drag_scale / g_scale, expected_ratio, rtol=4.0e-15, atol=0.0):
            raise ValueError("Eq 59 pair drag and diffusion scales violate the Killeen mass identity")
        h_tilde = np.asarray(coefficients.h_tilde, dtype=float)
        g1 = np.asarray(coefficients.g_tilde_1, dtype=float)
        g2 = np.asarray(coefficients.g_tilde_2, dtype=float)
        drag = drag_scale * h_tilde
        diffusion = 0.5 * g_scale * g2
        pitch = 0.5 * g_scale * g1
        contribution = Eq59PairOperatorContribution(pair_id=record.pair_id, test_species_id=record.test_species_id, field_population_id=record.field_population_id, field_species_id=record.field_species_id, field_distribution_model=record.field_distribution_model, field_mass_kg=record.field_mass_kg, field_charge_number=record.field_charge_number, field_density_m3=record.field_density_m3, field_temperature_J_or_none=record.field_temperature_J_or_none, coulomb_log_value=float(record.coulomb_log_state.value), rosenbluth_g_weight_m3=float(record.rosenbluth_g_weight_m3), rosenbluth_h_weight_m3=float(record.rosenbluth_h_weight_m3), g_scale_velocity_cubed_m3_s3=float(g_scale), net_drag_scale_velocity_cubed_m3_s3=float(drag_scale), h_tilde=h_tilde.copy(), g_tilde_1=g1.copy(), g_tilde_2=g2.copy(), drag_velocity_cubed_m3_s3=drag.copy(), energy_diffusion_velocity_fourth_m4_s4=diffusion.copy(), pitch_scattering_velocity_cubed_m3_s3=pitch.copy(), active_operator_terms=record.active_operator_terms, active=True, qualified=record.applicability, limitation=record.limitation)
        contributions[record.pair_id] = contribution
        total_drag += drag
        total_diffusion += diffusion
        total_pitch += pitch
        if record.field_population_id == "fast_self":
            fast_self_normalization = float(coefficients.density_normalization_m3)
    active_external_pair_ids: set[str] = set()
    for external in tuple(external_fast_field_states):
        if external.test_species != species:
            raise ValueError("external fast field test species must match the Eq 59 operator")
        if external.pair_id in contributions:
            raise ValueError("external fast field pair identifier duplicates an existing contribution")
        coefficients = external_fast_field_rosenbluth_coefficients(external, speed_grid)
        g_weight = float(external.field_physical_density_m3) * (float(external.field_species.charge_number) / float(species.charge_number)) ** 2 * float(external.coulomb_log_value)
        h_weight = g_weight * (float(species.mass_kg) + float(external.field_species.mass_kg)) / float(external.field_species.mass_kg)
        g_scale = tau_s * gamma * g_weight
        drag_scale = tau_s * gamma * (h_weight - g_weight)
        expected_ratio = species.mass_kg / external.field_species.mass_kg
        if g_scale <= 0.0 or drag_scale <= 0.0 or not np.isclose(drag_scale / g_scale, expected_ratio, rtol=4.0e-15, atol=0.0):
            raise ValueError("external fast field scales violate the Killeen mass identity")
        h_tilde = np.asarray(coefficients.h_tilde, dtype=float)
        g1 = np.asarray(coefficients.g_tilde_1, dtype=float)
        g2 = np.asarray(coefficients.g_tilde_2, dtype=float)
        drag = drag_scale * h_tilde
        diffusion = 0.5 * g_scale * g2
        pitch = 0.5 * g_scale * g1
        contribution = Eq59PairOperatorContribution(
            pair_id=external.pair_id,
            test_species_id=species.species_id,
            field_population_id=f"fast_{external.field_species.species_id}",
            field_species_id=external.field_species.species_id,
            field_distribution_model=external.model,
            field_mass_kg=external.field_species.mass_kg,
            field_charge_number=external.field_species.charge_number,
            field_density_m3=float(external.field_physical_density_m3),
            field_temperature_J_or_none=None,
            coulomb_log_value=float(external.coulomb_log_value),
            rosenbluth_g_weight_m3=g_weight,
            rosenbluth_h_weight_m3=h_weight,
            g_scale_velocity_cubed_m3_s3=float(g_scale),
            net_drag_scale_velocity_cubed_m3_s3=float(drag_scale),
            h_tilde=h_tilde.copy(),
            g_tilde_1=g1.copy(),
            g_tilde_2=g2.copy(),
            drag_velocity_cubed_m3_s3=drag.copy(),
            energy_diffusion_velocity_fourth_m4_s4=diffusion.copy(),
            pitch_scattering_velocity_cubed_m3_s3=pitch.copy(),
            active_operator_terms=("ion_drag", "energy_diffusion", "pitch_scattering"),
            active=True,
            qualified=bool(external.qualified),
            limitation=external.limitation,
        )
        contributions[external.pair_id] = contribution
        total_drag += drag
        total_diffusion += diffusion
        total_pitch += pitch
        active_external_pair_ids.add(external.pair_id)
    limitations = []
    for limitation in pairwise_collision_state.missing_pair_physics:
        pair_id = str(limitation).split(":", 1)[0]
        if pair_id not in active_external_pair_ids:
            limitations.append(limitation)
    for contribution in contributions.values():
        if contribution.active and contribution.limitation:
            limitations.append(f"{contribution.pair_id}:{contribution.limitation}")
    limitations.extend(("anisotropic_Rosenbluth_terms_unavailable", "full_multispecies_Landau_operator_not_claimed", "collision_coefficients_use_confined_volume_average_scalar_densities"))
    if any(contribution.active and contribution.field_population_id.startswith("thermal_") for contribution in contributions.values()):
        limitations.append("reference_Maxwellian_fields_are_prescribed")
    limitations = list(dict.fromkeys(limitations))
    qualified = bool(all((not contribution.active) or contribution.qualified for contribution in contributions.values()) and not limitations)
    operator_model = "egedal_eq59_isotropic_pairwise_with_reduced_fast_cross_species" if active_external_pair_ids else "egedal_eq59_isotropic_pairwise_thermal_D_T_and_fast_self"
    coulomb_log_convention = "thermal_pairs_use_pair_logs_fast_self_and_fast_cross_pairs_use_the_existing_representative_ion_log" if active_external_pair_ids else "thermal_D_and_T_use_pair_logs_fast_self_uses_legacy_shared_ion_log"
  
    return Eq59CollisionOperatorState(test_species=species, speed_grid=speed_grid, spitzer_slowing_down_time_s=tau_s, electron_pair_id=electron_pair_id, pair_contributions_by_id=contributions, ion_drag_velocity_cubed_m3_s3=total_drag, ion_energy_diffusion_velocity_fourth_m4_s4=total_diffusion, ion_pitch_scattering_velocity_cubed_m3_s3=total_pitch, fast_self_physical_density_m3=fast_density, fast_self_rosenbluth_density_normalization_m3=fast_self_normalization, operator_model=operator_model, density_convention=pairwise_collision_state.density_convention, temperature_convention=pairwise_collision_state.temperature_convention, coulomb_log_convention=coulomb_log_convention, qualified=qualified, limitations=tuple(limitations))

def eq59_collision_operator_metadata(state: Eq59CollisionOperatorState) -> dict[str, object]:
    """Return compact JSON friendly metadata for the pair resolved Eq 59 operator"""
    active = [contribution for contribution in state.pair_contributions_by_id.values() if contribution.active]
    pair_terms = {contribution.pair_id: list(contribution.active_operator_terms) for contribution in active}
    pair_g_scales = {contribution.pair_id: contribution.g_scale_velocity_cubed_m3_s3 for contribution in active}
    pair_drag_scales = {contribution.pair_id: contribution.net_drag_scale_velocity_cubed_m3_s3 for contribution in active}
    fast_self = next((contribution for contribution in active if contribution.field_population_id == "fast_self"), None)
    fast_cross = [contribution for contribution in active if contribution.field_population_id.startswith("fast_") and contribution.field_population_id != "fast_self"]
    return {
        "eq59_collision_operator_model": state.operator_model,
        "eq59_electron_collision_role": "analytic_spitzer_fast_ion_drag",
        "eq59_thermal_D_active": any(contribution.field_population_id == "thermal_deuterium" for contribution in active),
        "eq59_thermal_T_active": any(contribution.field_population_id == "thermal_tritium" for contribution in active),
        "eq59_fast_self_active": fast_self is not None,
        "eq59_pair_ids": [contribution.pair_id for contribution in state.pair_contributions_by_id.values()],
        "eq59_active_pair_ids": [contribution.pair_id for contribution in active],
        "eq59_pair_operator_terms": pair_terms,
        "eq59_pair_g_scales_velocity_cubed_m3_s3": pair_g_scales,
        "eq59_pair_net_drag_scales_velocity_cubed_m3_s3": pair_drag_scales,
        "eq59_fast_self_physical_density_m3": state.fast_self_physical_density_m3,
        "eq59_fast_self_shape_density_normalization_m3": state.fast_self_rosenbluth_density_normalization_m3,
        "eq59_fast_self_coulomb_log_model": (None if fast_self is None else "legacy_shared_ion_log_without_non_Maxwellian_relative_energy"),
        "eq59_fast_cross_species_collisions_available": bool(fast_cross),
        "eq59_fast_cross_species_pair_ids": [contribution.pair_id for contribution in fast_cross],
        "eq59_fast_cross_species_model": (None if not fast_cross else "reduced_isotropic_first_mode_Rosenbluth_cross_species"),
        "eq59_anisotropic_Rosenbluth_terms_available": False,
        "eq59_full_multispecies_Landau_operator_claimed": False,
        "eq59_collision_operator_qualified": state.qualified,
        "eq59_collision_operator_limitations": list(state.limitations),
        "eq59_max_ion_drag_velocity_cubed_m3_s3": float(np.max(state.ion_drag_velocity_cubed_m3_s3)),
        "eq59_max_ion_energy_diffusion_velocity_fourth_m4_s4": float(np.max(state.ion_energy_diffusion_velocity_fourth_m4_s4)),
        "eq59_max_ion_pitch_scattering_velocity_cubed_m3_s3": float(np.max(state.ion_pitch_scattering_velocity_cubed_m3_s3)),
    }

__all__ = [
    "Eq59PairOperatorContribution",
    "Eq59CollisionOperatorState",
    "build_eq59_collision_operator_state",
    "eq59_collision_operator_metadata",
]
