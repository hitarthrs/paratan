"""Shared immutable result and continuation state dataclasses for the integration pipeline"""
from __future__ import annotations
from collections.abc import Callable, Mapping
from dataclasses import dataclass
from typing import Any, TYPE_CHECKING
import numpy as np
from numpy.typing import ArrayLike

if TYPE_CHECKING:
    from source_model_revamp.beam.beam_path_geometry import BeamPathOnAxisymmetricGrid
    from source_model_revamp.beam.source_from_attenuation import AttenuatedMultiEnergyBeamSource
    from source_model_revamp.beam.charge_exchange_redistribution import ChargeExchangeTargetSinkState
    from source_model_revamp.export.openmc_file_source import OpenMCFileSourceBundle
    from source_model_revamp.fbis.mirror_invariant_grid import LambdaGrid
    from source_model_revamp.fbis.modal.species_state import FastIonSystemState
    from source_model_revamp.fbis.velocity_space_grid import PitchGrid, SpeedGrid
    from source_model_revamp.fusion.reactions import FusionReaction
    from source_model_revamp.fusion.populations import FusionReactantPopulation
    from source_model_revamp.fusion.pair_kernel import FusionPairKernelResult
    from source_model_revamp.fusion.source_profiles import FusionSourceProfile, MultiReactionFusionSourceProfile
    from source_model_revamp.integration.config.root import SourceModelRunConfig
    from source_model_revamp.integration.modal_stage.types import OperatingPointDensityState
    from source_model_revamp.neutrons.lambda_spectrum import AxialNeutronSpectrumMatrix
    from source_model_revamp.neutrons.events import CorrelatedNeutronEventBank
    from source_model_revamp.expander.types import ExpanderSystemState
    from source_model_revamp.coupling.electron_energy_balance import ElectronEnergyBalanceState
    from source_model_revamp.geometry.background_profiles import IonSeedReferenceProfiles
    from source_model_revamp.geometry.device_domains import DeviceDomains
    from source_model_revamp.geometry.plasma_boundaries import PlasmaBoundarySelection
    from source_model_revamp.integration.assessment import RunAssessment

@dataclass(frozen=True)
class GeometryStageResult:
    """
    Magnetic and flux tube geometry shared by downstream stages
    
    Confined arrays use the throat to throat grid while optional full device arrays extend to the selected material boundaries
    """
    z_edges_m: np.ndarray
    z_centers_m: np.ndarray
    zeta_edges: np.ndarray
    zeta_centers: np.ndarray
    B_tilde_centers: np.ndarray
    B_tilde_midpoints: np.ndarray
    B0_T: float
    mirror_ratio: float
    half_length_m: float
    plasma_radius_m: float
    radial_outer_radius_m_by_z: np.ndarray
    radial_inner_radius_m_by_z: np.ndarray
    midplane_area_m2: float
    cell_volumes_m3: np.ndarray
    volume_m3: float
    magnetic_field_model: str
    linked_to_paratan_geometry: bool
    metadata: dict[str, Any]
    B_tilde_function: Callable[[ArrayLike], float | np.ndarray]
    full_device_z_edges_m: np.ndarray | None = None
    full_device_z_centers_m: np.ndarray | None = None
    full_device_B_T_centers: np.ndarray | None = None
    full_device_B_tilde_centers: np.ndarray | None = None
    full_device_cell_volumes_m3: np.ndarray | None = None
    domains: DeviceDomains | None = None
    plasma_boundaries: PlasmaBoundarySelection | None = None
    background_profiles: IonSeedReferenceProfiles | None = None
    B_T_function: Callable[[ArrayLike], float | np.ndarray] | None = None

@dataclass(frozen=True)
class BeamDepositionResult:
    """
    Attenuation and mapped fast ion source for one enabled neutral beam
    
    `source_z_v_lambda_per_s` has shape `(n_z, n_speed, n_lambda)` and its volume averaged companion removes the axial dimension
    """
    beam_id: str
    projectile_species: str
    path_geometry: BeamPathOnAxisymmetricGrid
    attenuation: Any
    attenuated_source: AttenuatedMultiEnergyBeamSource
    component_energies_J: np.ndarray
    source_z_v_lambda_per_s: np.ndarray
    volume_averaged_source_v_lambda_per_s: np.ndarray
    total_fast_birth_rate_s: float
    total_deposited_birth_power_W: float
    metadata: dict[str, Any]

@dataclass(frozen=True)
class BeamEnsembleResult:
    """
    Combined enabled beam state with per beam diagnostics and species grouped fast ion sources
    
    Species can use distinct speed grids while sharing the magnetic Λ and local pitch grids
    """
    beam_results_by_id: Mapping[str, BeamDepositionResult]
    beam_order: tuple[str, ...]
    attenuated_source: AttenuatedMultiEnergyBeamSource
    attenuated_source_by_species: Mapping[str, AttenuatedMultiEnergyBeamSource]
    speed_grid: SpeedGrid
    speed_grid_by_species: Mapping[str, SpeedGrid]
    lambda_grid: LambdaGrid
    pitch_grid: PitchGrid
    component_energies_J: np.ndarray
    source_z_v_lambda_per_s: np.ndarray
    volume_averaged_source_v_lambda_per_s: np.ndarray
    combined_source_by_species: Mapping[str, np.ndarray]
    combined_volume_averaged_source_by_species: Mapping[str, np.ndarray]
    total_fast_birth_rate_by_species: Mapping[str, float]
    charge_exchange_sink_by_target_species: Mapping[str, ChargeExchangeTargetSinkState] | None
    total_fast_birth_rate_s: float
    total_deposited_birth_power_W: float
    total_injected_power_W: float
    metadata: dict[str, Any]

@dataclass(frozen=True)
class KineticStageResult:
    """
    Modal FBIS result plus local reconstruction, full device mapping, density state, and continuation data
    
    Legacy deuterium fields are retained together with species resolved mappings for coupled D plus T operation
    """
    speed_grid: SpeedGrid
    lambda_grid: LambdaGrid
    pitch_grid: PitchGrid
    final_distribution_v_lambda: np.ndarray
    metadata: dict[str, Any]
    local_speed_grid: SpeedGrid | None = None
    local_distribution_z_v_lambda: np.ndarray | None = None
    local_distribution_z_v_pitch: np.ndarray | None = None
    local_density_m3: np.ndarray | None = None
    full_device_local_distribution_z_v_pitch: np.ndarray | None = None
    eq59_warm_start_state: object | None = None
    operating_point_density_state: OperatingPointDensityState | None = None
    modal_basis_warm_start: object | None = None
    eq70_potential_energy_warm_start_J: np.ndarray | None = None
    eq42_density_profile_warm_start: object | None = None
    fast_ion_system_state: FastIonSystemState | None = None
    speed_grid_by_species: Mapping[str, SpeedGrid] | None = None
    final_distribution_v_lambda_by_species: Mapping[str, np.ndarray] | None = None
    local_speed_grid_by_species: Mapping[str, SpeedGrid] | None = None
    local_distribution_z_v_lambda_by_species: Mapping[str, np.ndarray] | None = None
    local_distribution_z_v_pitch_by_species: Mapping[str, np.ndarray] | None = None
    local_density_m3_by_species: Mapping[str, np.ndarray] | None = None
    full_device_local_distribution_z_v_pitch_by_species: Mapping[str, np.ndarray] | None = None
    eq59_warm_start_state_by_species: Mapping[str, object] | None = None
    scalar_collision_operator_warm_start_by_species: Mapping[str, tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]] | None = None
    scalar_external_fast_field_warm_start_by_test_species: Mapping[str, tuple[object, ...]] | None = None
    scalar_wall_barrier_energy_warm_start_J: float | None = None

@dataclass(frozen=True)
class BeamDensityCouplingIteration:
    """One fixed electron temperature beam and density fixed point iteration with density, beam source, qualification, and consistency metrics"""
    iteration: int
    density_pointwise_relative_change: float
    density_volume_L2_relative_change: float
    density_absolute_reference_change: float
    species_density_pointwise_relative_change: float | None
    species_density_volume_L2_relative_change: float | None
    species_density_absolute_reference_change: float | None
    species_density_relative_change_by_species: Mapping[str, float] | None
    total_fast_birth_rate_relative_change: float | None
    deposited_beam_power_relative_change: float | None
    birth_profile_pointwise_relative_change: float | None
    birth_profile_volume_L2_relative_change: float | None
    birth_profile_absolute_reference_change: float | None
    component_pointwise_relative_change: float | None
    component_L2_relative_change: float | None
    component_absolute_reference_change: float | None
    state_finite: bool
    state_physically_evaluable: bool
    kinetic_numerical_convergence_passed: bool
    end_loss_power_terms_available: bool
    candidate_fixed_point_passed: bool = False
    final_consistency_recompute_passed: bool = False
    species_density_change_metrics_by_species: Mapping[str, Mapping[str, float]] | None = None

    def as_metadata(self) -> dict[str, object]:
        """Return the iteration convergence and qualification fields as serializable metadata"""
        return {
            "iteration": self.iteration,
            "density_pointwise_relative_change": self.density_pointwise_relative_change,
            "density_volume_L2_relative_change": self.density_volume_L2_relative_change,
            "density_absolute_reference_change": self.density_absolute_reference_change,
            "species_density_pointwise_relative_change": self.species_density_pointwise_relative_change,
            "species_density_volume_L2_relative_change": self.species_density_volume_L2_relative_change,
            "species_density_absolute_reference_change": self.species_density_absolute_reference_change,
            "species_density_relative_change_by_species": None if self.species_density_relative_change_by_species is None else dict(self.species_density_relative_change_by_species),
            "species_density_change_metrics_by_species": None if self.species_density_change_metrics_by_species is None else {species_id: dict(metrics) for species_id, metrics in self.species_density_change_metrics_by_species.items()},
            "total_fast_birth_rate_relative_change": self.total_fast_birth_rate_relative_change,
            "deposited_beam_power_relative_change": self.deposited_beam_power_relative_change,
            "birth_profile_pointwise_relative_change": self.birth_profile_pointwise_relative_change,
            "birth_profile_volume_L2_relative_change": self.birth_profile_volume_L2_relative_change,
            "birth_profile_absolute_reference_change": self.birth_profile_absolute_reference_change,
            "component_pointwise_relative_change": self.component_pointwise_relative_change,
            "component_L2_relative_change": self.component_L2_relative_change,
            "component_absolute_reference_change": self.component_absolute_reference_change,
            "state_finite": self.state_finite,
            "state_physically_evaluable": self.state_physically_evaluable,
            "kinetic_numerical_convergence_passed": self.kinetic_numerical_convergence_passed,
            "end_loss_power_terms_available": self.end_loss_power_terms_available,
            "candidate_fixed_point_passed": self.candidate_fixed_point_passed,
            "final_consistency_recompute_passed": self.final_consistency_recompute_passed,
        }

@dataclass(frozen=True)
class OperatingPointWarmStartState:
    """Same run continuation state for beam target density, modal basis, Eq 59, Eq 70, Eq 42, collision fields, and wall barrier"""
    target_density_profile_m3: np.ndarray
    target_density_defined_mask: np.ndarray
    kinetic_target_state: KineticStageResult | None = None
    eq59_warm_start_state: object | None = None
    density_state: OperatingPointDensityState | None = None
    modal_basis_warm_start: object | None = None
    eq70_potential_energy_warm_start_J: np.ndarray | None = None
    eq42_density_profile_warm_start: object | None = None
    eq59_warm_start_state_by_species: Mapping[str, object] | None = None
    scalar_collision_operator_warm_start_by_species: Mapping[str, tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]] | None = None
    scalar_external_fast_field_warm_start_by_test_species: Mapping[str, tuple[object, ...]] | None = None
    scalar_wall_barrier_energy_warm_start_J: float | None = None

@dataclass(frozen=True)
class OperatingPointStageResult:
    """Closed beam, kinetic, expander, density, and represented electron energy state for one electron temperature"""
    beam: BeamEnsembleResult
    kinetic: KineticStageResult
    density_state: OperatingPointDensityState | None
    beam_density_coupling_history: tuple[BeamDensityCouplingIteration, ...]
    converged: bool
    failure_reason: str | None
    warm_start_state: OperatingPointWarmStartState
    metadata: dict[str, Any]
    expander: ExpanderStageResult | None = None
    electron_energy_balance: ElectronEnergyBalanceState | None = None
    state_finite: bool = False
    state_physically_evaluable: bool = False
    beam_density_fixed_point_converged: bool = False
    kinetic_numerical_convergence_passed: bool = False
    end_loss_power_terms_available: bool = False
    electron_temperature_assignment_valid: bool = False
    electron_temperature_power_balance_converged: bool | None = None

@dataclass(frozen=True)
class ExpanderStageResult:
    """Source connected expander system state and associated metadata"""
    system_state: ExpanderSystemState | None
    metadata: dict[str, Any]

@dataclass(frozen=True)
class FusionComponent:
    """One labeled fusion population pair, reaction branch, axial source profile, and optional pair kernel"""
    label: str
    kind: str
    reaction: FusionReaction
    reactant_a_population_id: str
    reactant_b_population_id: str
    profile: FusionSourceProfile
    pair_kernel: FusionPairKernelResult | None = None

@dataclass(frozen=True)
class FusionStageResult:
    """All active fusion components, combined source profile, reactant population registry, and metadata"""
    components: tuple[FusionComponent, ...]
    combined_profile: MultiReactionFusionSourceProfile
    populations_by_id: Mapping[str, FusionReactantPopulation]
    metadata: dict[str, Any]

@dataclass(frozen=True)
class NeutronStageResult:
    """
    Deterministic neutron spectrum matrices and optional correlated event bank
    
    `source_rate_matrix_s` has shape `(n_z, n_energy)` and stores physical neutron rates in axial and energy bins
    """
    component_matrices: tuple[tuple[str, AxialNeutronSpectrumMatrix], ...]
    source_rate_matrix_s: np.ndarray
    energy_edges_J: np.ndarray
    total_neutron_rate_s: float
    metadata: dict[str, Any]
    correlated_event_bank: CorrelatedNeutronEventBank | None = None
    correlated_event_source_rate_matrix_s: np.ndarray | None = None

@dataclass(frozen=True)
class OpenMCExportStageResult:
    """Optional correlated OpenMC file source bundle and export metadata"""
    file_source_bundle: OpenMCFileSourceBundle | None
    metadata: dict[str, Any]

@dataclass(frozen=True)
class SourceModelPipelineResult:
    """Complete source model result with typed stage outputs, assessment, and canonical metadata"""
    config: SourceModelRunConfig
    geometry: GeometryStageResult | None
    beam: BeamEnsembleResult | None
    kinetic: KineticStageResult | None
    expander: ExpanderStageResult | None
    fusion: FusionStageResult | None
    neutrons: NeutronStageResult | None
    openmc_export: OpenMCExportStageResult | None
    metadata: dict[str, Any]
    operating_point: OperatingPointStageResult | None = None
    assessment: RunAssessment | None = None

    @property
    def calculation_valid(self) -> bool:
        """Return whether a final RunAssessment was constructed after direct validation"""
        return self.assessment is not None

    @property
    def numerically_converged(self) -> bool:
        """Return the aggregate numerical convergence state when an assessment exists"""
        return bool(self.assessment is not None and self.assessment.numerically_converged)

    @property
    def model_applicable(self) -> bool:
        """Return the aggregate model applicability state when an assessment exists"""
        return bool(self.assessment is not None and self.assessment.model_applicable)

    @property
    def result_status(self) -> str | None:
        """Return the combined convergence and applicability result status"""
        return None if self.assessment is None else self.assessment.result_status

    @property
    def failed_convergence_checks(self) -> tuple[str, ...]:
        """Return names of convergence checks that explicitly failed"""
        return () if self.assessment is None else self.assessment.failed_convergence_checks

    @property
    def applicability_limitations(self) -> tuple[str, ...]:
        """Return names of applicability checks that explicitly failed"""
        return () if self.assessment is None else self.assessment.applicability_limitations

    @property
    def source_constructible(self) -> bool:
        """Return the final correlated source constructibility flag from canonical metadata"""
        return self.metadata.get("source_constructible") is True

    @property
    def source_generated(self) -> bool:
        """Return whether an OpenMC source bundle was generated"""
        return self.metadata.get("source_generated") is True

    @property
    def source_file_valid(self) -> bool | None:
        """Return the OpenMC source file validity flag when it is a literal boolean"""
        value = self.metadata.get("source_file_valid")
        return value if isinstance(value, bool) else None
