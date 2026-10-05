"""Pairwise collision state for the species generic FBIS backend"""
from __future__ import annotations
from collections.abc import Mapping
from dataclasses import dataclass
from types import MappingProxyType
import numpy as np
from source_model_revamp.fbis.species import IonSpecies

@dataclass(frozen=True)
class PairCoulombLogState:
    """One ordered pair Coulomb logarithm state"""
    model: str
    value: float | None
    screening_length_m: float
    reduced_mass_kg: float
    relative_energy_J: float | None
    relative_energy_model: str
    b_90_m: float | None
    b_quantum_m: float | None
    b_min_m: float | None
    reference_identity: str
    active_in_current_eq59_operator: bool
    applicability: bool
    limitation: str | None

    def __post_init__(self) -> None:
        if not self.model or not self.relative_energy_model or not self.reference_identity:
            raise ValueError("pair Coulomb log text fields must be nonempty")
        if not np.isfinite(self.screening_length_m) or self.screening_length_m <= 0.0:
            raise ValueError("pair Coulomb log screening length must be positive and finite")
        if not np.isfinite(self.reduced_mass_kg) or self.reduced_mass_kg <= 0.0:
            raise ValueError("pair Coulomb log reduced mass must be positive and finite")
        if self.value is not None and (not np.isfinite(self.value) or self.value <= 0.0):
            raise ValueError("pair Coulomb log value must be positive and finite when available")
        if self.relative_energy_J is not None and (not np.isfinite(self.relative_energy_J) or self.relative_energy_J <= 0.0):
            raise ValueError("pair relative energy must be positive and finite when available")
        for name, value in (("b_90_m", self.b_90_m), ("b_quantum_m", self.b_quantum_m), ("b_min_m", self.b_min_m)):
            if value is not None and (not np.isfinite(value) or value <= 0.0):
                raise ValueError(f"{name} must be positive and finite when available")
        detailed_cutoffs = (self.relative_energy_J, self.b_90_m, self.b_quantum_m, self.b_min_m)
        if any(value is None for value in detailed_cutoffs) and not all(value is None for value in detailed_cutoffs):
            raise ValueError("pair Coulomb log relative energy and cutoff fields must be available together")

@dataclass(frozen=True)
class PairwiseCollisionRecord:
    """One ordered test species and field population record"""
    pair_id: str
    test_species_id: str
    field_population_id: str
    field_species_id: str
    field_distribution_model: str
    field_mass_kg: float
    field_charge_number: float
    field_density_m3: float
    field_temperature_J_or_none: float | None
    field_distribution_reference_or_none: str | None
    coulomb_log_state: PairCoulombLogState
    rosenbluth_g_weight_m3: float | None
    rosenbluth_h_weight_m3: float | None
    legacy_critical_velocity_cubed_contribution_m3_s3_or_none: float | None
    operator_role: str
    active_in_current_eq59_operator: bool
    active_operator_terms: tuple[str, ...]
    applicability: bool
    limitation: str | None

    def __post_init__(self) -> None:
        text_values = (self.pair_id, self.test_species_id, self.field_population_id, self.field_species_id, self.field_distribution_model, self.operator_role)
        if any(not value for value in text_values):
            raise ValueError("pairwise collision record text fields must be nonempty")
        if not self.active_operator_terms or any(not str(value).strip() for value in self.active_operator_terms):
            raise ValueError("pairwise collision record active operator terms must be nonempty")
        if not np.isfinite(self.field_mass_kg) or self.field_mass_kg <= 0.0:
            raise ValueError("field population mass must be positive and finite")
        if not np.isfinite(self.field_charge_number) or self.field_charge_number == 0.0:
            raise ValueError("field population charge number must be finite and nonzero")
        if not np.isfinite(self.field_density_m3) or self.field_density_m3 < 0.0:
            raise ValueError("field population density must be finite and nonnegative")
        if self.field_temperature_J_or_none is not None and (not np.isfinite(self.field_temperature_J_or_none) or self.field_temperature_J_or_none <= 0.0):
            raise ValueError("field population temperature must be positive and finite when available")
        for name, value in (("rosenbluth_g_weight_m3", self.rosenbluth_g_weight_m3), ("rosenbluth_h_weight_m3", self.rosenbluth_h_weight_m3), ("legacy_critical_velocity_cubed_contribution_m3_s3_or_none", self.legacy_critical_velocity_cubed_contribution_m3_s3_or_none)):
            if value is not None and (not np.isfinite(value) or value < 0.0):
                raise ValueError(f"{name} must be finite and nonnegative when available")
        if self.active_in_current_eq59_operator and self.field_density_m3 <= 0.0:
            raise ValueError("zero density field populations cannot be active in the current Eq 59 operator")
        weights = (self.rosenbluth_g_weight_m3, self.rosenbluth_h_weight_m3)
        if any(value is None for value in weights) and not all(value is None for value in weights):
            raise ValueError("Rosenbluth pair weights must be available together")

@dataclass(frozen=True)
class PairwiseCollisionScreeningState:
    """Shared screening state used by all ordered pair records"""
    model: str
    screening_length_m: float
    electron_density_m3: float
    electron_temperature_J: float
    ion_population_ids: tuple[str, ...]
    ion_densities_m3: tuple[float, ...]
    ion_charge_numbers: tuple[float, ...]
    ion_temperature_J: float
    fast_self_included_in_screening: bool
    reference_identity: str

    def __post_init__(self) -> None:
        if not self.model or not self.reference_identity:
            raise ValueError("pairwise screening text fields must be nonempty")
        if not np.isfinite(self.screening_length_m) or self.screening_length_m <= 0.0:
            raise ValueError("pairwise screening length must be positive and finite")
        if not np.isfinite(self.electron_density_m3) or self.electron_density_m3 <= 0.0:
            raise ValueError("pairwise screening electron density must be positive and finite")
        if not np.isfinite(self.electron_temperature_J) or self.electron_temperature_J <= 0.0:
            raise ValueError("pairwise screening electron temperature must be positive and finite")
        if not np.isfinite(self.ion_temperature_J) or self.ion_temperature_J <= 0.0:
            raise ValueError("pairwise screening ion temperature must be positive and finite")
        if not (len(self.ion_population_ids) == len(self.ion_densities_m3) == len(self.ion_charge_numbers)):
            raise ValueError("pairwise screening ion arrays must have matching lengths")
        if any(not np.isfinite(value) or value < 0.0 for value in self.ion_densities_m3):
            raise ValueError("pairwise screening ion densities must be finite and nonnegative")
        if any(not np.isfinite(value) or value == 0.0 for value in self.ion_charge_numbers):
            raise ValueError("pairwise screening ion charges must be finite and nonzero")

@dataclass(frozen=True)
class ReducedEq59CollisionProjection:
    """Legacy scalar collision projection retained for diagnostics and Eq 14"""
    spitzer_slowing_down_time_s: float
    critical_velocity_m_s: float
    beta_m: float
    electron_coulomb_log: float
    legacy_shared_ion_coulomb_log: float
    projection_model: str
    source_pair_ids: tuple[str, ...]
    representative_field_ion_species_id: str
    mixed_background_projection_used: bool
    qualified: bool
    limitation: str | None

    def __post_init__(self) -> None:
        for name, value in (("spitzer_slowing_down_time_s", self.spitzer_slowing_down_time_s), ("critical_velocity_m_s", self.critical_velocity_m_s), ("beta_m", self.beta_m), ("electron_coulomb_log", self.electron_coulomb_log), ("legacy_shared_ion_coulomb_log", self.legacy_shared_ion_coulomb_log)):
            if not np.isfinite(value) or value <= 0.0:
                raise ValueError(f"{name} must be positive and finite")
        if not self.projection_model or not self.representative_field_ion_species_id:
            raise ValueError("reduced Eq 59 projection text fields must be nonempty")
        if not self.source_pair_ids:
            raise ValueError("reduced Eq 59 projection requires source pair identifiers")

@dataclass(frozen=True)
class PairwiseCollisionState:
    """Authoritative ordered collision identity and provenance state"""
    test_species: IonSpecies
    records_by_pair_id: Mapping[str, PairwiseCollisionRecord]
    screening_state: PairwiseCollisionScreeningState
    density_convention: str
    temperature_convention: str
    reduced_eq59_projection: ReducedEq59CollisionProjection
    missing_pair_physics: tuple[str, ...]
    applicability: str

    def __post_init__(self) -> None:
        records = dict(self.records_by_pair_id)
        if not records:
            raise ValueError("pairwise collision state requires at least one record")
        if len(records) != len(set(records)):
            raise ValueError("pairwise collision record identifiers must be unique")
        for pair_id, record in records.items():
            if pair_id != record.pair_id:
                raise ValueError("pairwise collision state keys must match record identifiers")
            if record.test_species_id != self.test_species.species_id:
                raise ValueError("pairwise collision record test species must match the state")
        if not self.density_convention or not self.temperature_convention or not self.applicability:
            raise ValueError("pairwise collision state convention fields must be nonempty")
        object.__setattr__(self, "records_by_pair_id", MappingProxyType(records))

    def record(self, pair_id: str) -> PairwiseCollisionRecord:
        """Return one ordered pair record"""
        try:
            return self.records_by_pair_id[str(pair_id)]
        except KeyError as exc:
            raise ValueError(f"pairwise collision record {pair_id!r} is unavailable") from exc

def pair_coulomb_log_metadata(state: PairCoulombLogState) -> dict[str, object]:
    """Return JSON friendly pair Coulomb log metadata"""
    return {
        "model": state.model,
        "value": state.value,
        "screening_length_m": state.screening_length_m,
        "reduced_mass_kg": state.reduced_mass_kg,
        "relative_energy_J": state.relative_energy_J,
        "relative_energy_model": state.relative_energy_model,
        "b_90_m": state.b_90_m,
        "b_quantum_m": state.b_quantum_m,
        "b_min_m": state.b_min_m,
        "reference_identity": state.reference_identity,
        "active_in_current_eq59_operator": state.active_in_current_eq59_operator,
        "applicability": state.applicability,
        "limitation": state.limitation,
    }

def pairwise_collision_record_metadata(record: PairwiseCollisionRecord) -> dict[str, object]:
    """Return JSON friendly ordered pair metadata"""
    return {
        "pair_id": record.pair_id,
        "test_species_id": record.test_species_id,
        "field_population_id": record.field_population_id,
        "field_species_id": record.field_species_id,
        "field_distribution_model": record.field_distribution_model,
        "field_mass_kg": record.field_mass_kg,
        "field_charge_number": record.field_charge_number,
        "field_density_m3": record.field_density_m3,
        "field_temperature_J": record.field_temperature_J_or_none,
        "field_distribution_reference": record.field_distribution_reference_or_none,
        "coulomb_log_state": pair_coulomb_log_metadata(record.coulomb_log_state),
        "rosenbluth_g_weight_m3": record.rosenbluth_g_weight_m3,
        "rosenbluth_h_weight_m3": record.rosenbluth_h_weight_m3,
        "legacy_critical_velocity_cubed_contribution_m3_s3": record.legacy_critical_velocity_cubed_contribution_m3_s3_or_none,
        "operator_role": record.operator_role,
        "active_in_current_eq59_operator": record.active_in_current_eq59_operator,
        "active_operator_terms": list(record.active_operator_terms),
        "applicability": record.applicability,
        "limitation": record.limitation,
    }

def pairwise_collision_metadata(state: PairwiseCollisionState) -> dict[str, object]:
    """Return compact canonical pairwise collision metadata"""
    projection = state.reduced_eq59_projection
    screening = state.screening_state
    return {
        "pairwise_collision_state_model": "ordered_test_species_field_population_state",
        "pairwise_collision_reference_model": "Rosenbluth_MacDonald_Judd_1957_Killeen_Mirin_Rensink_1976",
        "pairwise_collision_test_species_id": state.test_species.species_id,
        "pairwise_collision_density_convention": state.density_convention,
        "pairwise_collision_temperature_convention": state.temperature_convention,
        "pairwise_collision_applicability": state.applicability,
        "pairwise_collision_missing_pair_physics": list(state.missing_pair_physics),
        "pairwise_collision_screening_state": {
            "model": screening.model,
            "screening_length_m": screening.screening_length_m,
            "electron_density_m3": screening.electron_density_m3,
            "electron_temperature_J": screening.electron_temperature_J,
            "ion_population_ids": list(screening.ion_population_ids),
            "ion_densities_m3": list(screening.ion_densities_m3),
            "ion_charge_numbers": list(screening.ion_charge_numbers),
            "ion_temperature_J": screening.ion_temperature_J,
            "fast_self_included_in_screening": screening.fast_self_included_in_screening,
            "reference_identity": screening.reference_identity,
        },
        "pairwise_collision_records": [pairwise_collision_record_metadata(record) for record in state.records_by_pair_id.values()],
        "reduced_eq59_collision_projection": {
            "spitzer_slowing_down_time_s": projection.spitzer_slowing_down_time_s,
            "critical_velocity_m_s": projection.critical_velocity_m_s,
            "beta_m": projection.beta_m,
            "electron_coulomb_log": projection.electron_coulomb_log,
            "legacy_shared_ion_coulomb_log": projection.legacy_shared_ion_coulomb_log,
            "projection_model": projection.projection_model,
            "source_pair_ids": list(projection.source_pair_ids),
            "representative_field_ion_species_id": projection.representative_field_ion_species_id,
            "mixed_background_projection_used": projection.mixed_background_projection_used,
            "qualified": projection.qualified,
            "limitation": projection.limitation,
            "hot_eq59_authority": False,
            "cold_eq14_reference_authority": True,
            "diagnostic_comparison": True,
        },
    }

__all__ = [
    "PairCoulombLogState",
    "PairwiseCollisionRecord",
    "PairwiseCollisionScreeningState",
    "ReducedEq59CollisionProjection",
    "PairwiseCollisionState",
    "pair_coulomb_log_metadata",
    "pairwise_collision_record_metadata",
    "pairwise_collision_metadata",
]
