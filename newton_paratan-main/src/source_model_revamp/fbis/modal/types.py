"""Shared dataclasses and physical validity errors for the Egedal modal FBIS backend"""
from __future__ import annotations
from collections.abc import Mapping
from dataclasses import dataclass
from typing import TYPE_CHECKING
import numpy as np
from source_model_revamp.fbis.mirror_invariant_grid import LambdaGrid
from source_model_revamp.fbis.velocity_space_grid import PitchGrid, SpeedGrid
from source_model_revamp.fbis.species import IonSpecies
from source_model_revamp.orbits.eta_mapping import EtaLambdaMap
from source_model_revamp.scattering.physical_eigenbasis import PhysicalEigenbasis
from source_model_revamp.scattering.square_mirror_basis import SquareMirrorBasis

if TYPE_CHECKING:
    from source_model_revamp.fbis.modal.eq59.collision_operator_state import Eq59CollisionOperatorState
    from source_model_revamp.fbis.modal.eq59.types import Eq59WarmStartState
    from source_model_revamp.fbis.modal.eq59.energy_moments import FastIonElectronHeatingState
    from source_model_revamp.electrostatic.current_balance import AmbipolarCurrentBalance

class ModalPhysicalDistributionError(ValueError):
    """Raised when a reconstructed modal distribution violates a hard physical validity requirement"""

@dataclass(frozen=True)
class ModalFBISNumerics:
    """Numerical controls for basis, Eq 59, electrostatic, local reconstruction, loss, source, and Eq 42 convergence"""
    n_lambda_grid: int = 161
    n_eta_grid: int = 121
    n_square_basis_modes: int = 8
    n_physical_modes: int = 4
    square_scan_points: int = 20000
    speed_core_cell_fraction: float = 0.7
    speed_tail_stretch_power: float = 2.0
    velocity_solution_model: str = "hot_ion_rosenbluth_eq59"
    rosenbluth_max_iterations: int = 4
    rosenbluth_relative_tolerance: float = 1.0e-4
    rosenbluth_absolute_tolerance: float = 0.0
    rosenbluth_relaxation: float = 1.0
    rosenbluth_min_iterations: int = 2
    rosenbluth_high_speed_tail_population_tolerance: float = 1.0e-4
    rosenbluth_high_speed_tail_energy_tolerance: float = 1.0e-4
    rosenbluth_speed_convergence_relative_tolerance: float = 1.0e-2
    retained_mode_relative_tolerance: float = 1.0e-2
    local_speed_refinement_relative_tolerance: float = 1.0e-2
    phi_z_iterations: int = 32
    phi_z_root_scan_points: int = 33
    phi_z_relative_tolerance: float = 1.0e-3
    phi_z_relaxation: float = 0.8
    local_velocity_quadrature_order: int = 3
    phi_z_low_energy_weight_fraction_tolerance: float = 1.0e-3
    electrostatic_feedback_iterations: int = 4
    electrostatic_feedback_relative_tolerance: float = 1.0e-3
    basis_convergence_levels: int = 1
    basis_eigenvalue_relative_tolerance: float = 5.0e-3
    basis_eigenfunction_overlap_tolerance: float = 0.995
    loss_convention_relative_tolerance: float = 1.0e-2
    source_projection_relative_tolerance: float = 1.0e-2
    eq42_density_basis_iterations: int = 4
    eq42_density_shape_relative_tolerance: float = 1.0e-3
    eq42_density_shape_relaxation: float = 0.5
    eq42_basis_eigenvalue_relative_tolerance: float = 5.0e-3
    eq42_basis_eigenfunction_overlap_tolerance: float = 0.995
    eq42_density_symmetry_tolerance: float = 1.0e-2
    negative_distribution_roundoff_tolerance: float = 1.0e-12
    reconstruction_negative_particle_fraction_tolerance: float = 1.0e-3
    reconstruction_energy_correction_relative_tolerance: float = 5.0e-3

@dataclass(frozen=True)
class ModalFBISBasis:
    """Orbit averaged modal basis and loss geometry state
    
    η and Λ boundaries are dimensionless and mode_slopes_dI_dlambda_at_boundary contains one value per retained physical mode
    """
    eta_lambda_map: EtaLambdaMap
    square_basis: SquareMirrorBasis
    physical_basis: PhysicalEigenbasis
    eta_boundary: float
    lambda_boundary: float
    mode_slopes_dI_dlambda_at_boundary: np.ndarray
    geometry_factor_G: float
    volume_geometry_integral: float
    eq61_geometry_density_weighted: bool = False
    eq61_geometry_weighting_model: str = "uniform_density_reference"
    eq61_geometry_density_ratio_minimum: float = 1.0
    eq61_geometry_density_ratio_maximum: float = 1.0
    eq42_density_weighting_model: str = "uniform_density_reference"
    eq42_density_weighting_coupled: bool = False
    eq42_nonuniform_density_weighting_active: bool = False
    eq42_density_profile_source: str = "uniform_density_reference"
    eq42_volume_average_density_m3: float | None = None
    eq42_density_ratio_zeta: np.ndarray | None = None
    eq42_density_ratio_values: np.ndarray | None = None
    eq42_density_profile_identity: str | None = None
    eq42_density_profile_symmetry_relative_error: float = 0.0
    eq42_density_profile_symmetry_tolerance: float = 1.0e-2
    eq42_density_profile_symmetric: bool = True

@dataclass(frozen=True)
class ModalRosenbluthCoefficients:
    """Density normalized Eq 59 Rosenbluth coefficients on the species speed grid"""
    h_tilde: np.ndarray
    g_tilde_1: np.ndarray
    g_tilde_2: np.ndarray
    density_normalization_m3: float

@dataclass(frozen=True)
class ModalElectrostaticProfile:
    """Egedal Eq 70-72 quasineutral axial potential profile

    Primary axial arrays have shape (n_z,) with energies in J, potentials in V, and densities in m⁻³
    Additional fields retain iteration, root, support, inventory, exact node, and species resolved qualification diagnostics
    """
    zeta: np.ndarray
    B_tilde: np.ndarray
    potential_energy_J: np.ndarray
    potential_relative_to_midplane_V: np.ndarray
    electron_density_m3: np.ndarray
    ion_density_m3: np.ndarray
    iterations: int
    max_relative_quasineutrality_error: float
    converged: bool
    electron_density_scale_n0_m3: float
    electron_midplane_density_m3: float
    electron_parent_maxwellian_n0_m3: float
    electron_collision_density_m3: float
    electron_volume_average_density_m3: float
    fast_ion_density_m3: np.ndarray
    background_positive_charge_density_m3: np.ndarray
    electron_density_input_semantics: str
    electron_density_closure_model: str
    relative_residual_profile: np.ndarray
    absolute_residual_profile_m3: np.ndarray
    residual_history: tuple[float, ...]
    electron_density_scale_history_m3: tuple[float, ...]
    stagnated: bool = False
    oscillatory: bool = False
    diverged: bool = False
    failure_reason: str | None = None
    energy_mapping_model: str = "conserved_total_ion_energy"
    invariant_mapping_model: str = "egedal_eq71_eq72"
    density_measure: str = "local_2pi_v_squared_dv_dxi"
    roundoff_clip_count: int = 0
    throat_extrapolation_clip_count: int = 0
    throat_potential_energy_J: float = 0.0
    throat_potential_left_energy_J: float = 0.0
    throat_potential_right_energy_J: float = 0.0
    low_energy_approximation_sample_count: int = 0
    nonpositive_total_energy_sample_count: int = 0
    closed_eq71_boundary_sample_count: int = 0
    low_energy_approximation_max_weight_fraction: float = 0.0
    low_energy_approximation_valid: bool = True
    low_energy_approximation_weight_fraction_assessed: bool = False
    low_energy_approximation_weight_model: str = "not_assessed"
    eq71_closed_interval_intersected_population_cell_count: int = 0
    eq71_closed_interval_intersects_distribution_support: bool = False
    eta_to_local_phase_space_normalization: float = 1.0
    ion_inventory_normalization_history: tuple[float, ...] = ()
    phase_space_measure_conversion_factor: float = 1.0
    phase_space_measure_conversion_history: tuple[float, ...] = ()
    phase_space_measure_conversion_model: str = "not_evaluated"
    fast_ion_inventory_particles: float = 0.0
    fast_ion_volume_average_density_m3: float = 0.0
    target_volume_averaged_fast_ion_density_m3: float | None = None
    fast_ion_volume_average_density_absolute_error_m3: float | None = None
    fast_ion_volume_average_density_relative_error: float | None = None
    symmetric_input: bool = False
    midplane_reference_V: float = 0.0
    throat_extrapolation_valid: bool = True
    throat_potential_definition: str = "direct_quasineutral_throat_node_solve"
    potential_clipping_policy: str = "roundoff_only_otherwise_raise"
    eq71_closed_interval_threshold_J: float = 0.0
    effective_potential_check_assessed: bool = False
    effective_potential_throat_boundary_valid: bool = False
    effective_potential_failure_reason: str | None = None
    effective_potential_active_energy_sample_count: int = 0
    effective_potential_open_energy_sample_count: int = 0
    effective_potential_interior_minimum_energy_sample_count: int = 0
    effective_potential_interior_minimum_weight_fraction: float = 0.0
    effective_potential_max_relative_boundary_shortfall: float = 0.0
    effective_potential_max_relative_throat_asymmetry: float = 0.0
    effective_potential_worst_total_energy_J: float = 0.0
    effective_potential_worst_interior_location_zeta: float | None = None
    effective_potential_relative_tolerance: float = 1.0e-9
    absolute_residual_normalized_to_reference: np.ndarray | None = None
    density_support_mask: np.ndarray | None = None
    density_support_reference_m3: float = 0.0
    density_support_relative_floor: float = 0.0
    density_support_floor_m3: float = 0.0
    density_support_floor_derivation: str = "not_evaluated"
    density_support_masked_node_count: int = 0
    density_support_masked_volume_m3: float = 0.0
    density_support_masked_volume_fraction: float = 0.0
    volume_integrated_electron_minus_ion_particles: float = 0.0
    volume_integrated_absolute_particle_mismatch: float = 0.0
    volume_integrated_signed_charge_mismatch_C: float = 0.0
    volume_integrated_absolute_charge_mismatch_C: float = 0.0
    central_region_maximum_supported_relative_residual: float = 0.0
    throat_region_maximum_supported_relative_residual: float = 0.0
    maximum_absolute_density_residual_m3: float = 0.0
    volume_integrated_signed_particle_mismatch_fraction: float = 0.0
    maximum_absolute_density_residual_normalized_to_reference: float = 0.0
    volume_integrated_absolute_particle_mismatch_fraction: float = 0.0
    electrostatic_node_zeta: np.ndarray | None = None
    electrostatic_node_B_tilde: np.ndarray | None = None
    electrostatic_node_potential_energy_J: np.ndarray | None = None
    electrostatic_node_potential_relative_to_midplane_V: np.ndarray | None = None
    electrostatic_node_electron_density_m3: np.ndarray | None = None
    electrostatic_node_ion_density_m3: np.ndarray | None = None
    electrostatic_node_fast_ion_density_m3: np.ndarray | None = None
    electrostatic_node_background_positive_charge_density_m3: np.ndarray | None = None
    electrostatic_node_absolute_residual_profile_m3: np.ndarray | None = None
    electrostatic_node_relative_residual_profile: np.ndarray | None = None
    electrostatic_left_throat_index: int | None = None
    electrostatic_midplane_index: int | None = None
    electrostatic_right_throat_index: int | None = None
    electrostatic_left_throat_node_exact: bool = False
    electrostatic_midplane_node_exact: bool = False
    electrostatic_right_throat_node_exact: bool = False
    electrostatic_midplane_gauge_exact: bool = False
    electrostatic_node_model: str = "not_evaluated"
    floating_nonnegative_reference_active: bool = False
    floating_zero_reference_kind: str = "not_evaluated"
    floating_zero_reference_zeta: float | None = None
    floating_zero_reference_B_tilde: float | None = None
    exact_midplane_potential_energy_J: float = 0.0
    floating_scalar_closure_converged: bool = False
    floating_scalar_closure_history: tuple[Mapping[str, object], ...] = ()
    direct_exact_midplane_root_valid: bool = False
    local_speed_grid: SpeedGrid | None = None
    exact_anchor_maximum_quasineutrality_residual: float = 0.0
    exact_anchor_quasineutrality_valid: bool = False
    exact_anchor_residual_support_model: str = "not_evaluated"
    exact_midplane_quasineutrality_residual: float = 0.0
    exact_left_throat_quasineutrality_residual: float = 0.0
    exact_right_throat_quasineutrality_residual: float = 0.0
    exact_throat_maximum_quasineutrality_residual: float = 0.0
    exact_midplane_quasineutrality_valid: bool = False
    exact_throat_quasineutrality_valid: bool = False
    direct_left_throat_root_valid: bool = False
    direct_right_throat_root_valid: bool = False
    direct_throat_roots_valid: bool = False
    direct_central_roots_valid: bool = False
    direct_central_no_root_cell_count: int = 0
    direct_central_interior_scan_cell_count: int = 0
    direct_central_multiple_root_cell_count: int = 0
    direct_central_maximum_root_count: int = 0
    direct_central_root_scan_points: int = 0
    fast_ion_density_by_species_m3: Mapping[str, np.ndarray] | None = None
    electrostatic_node_fast_ion_density_by_species_m3: Mapping[str, np.ndarray] | None = None
    local_speed_grid_by_species: Mapping[str, SpeedGrid] | None = None
    fast_ion_inventory_particles_by_species: Mapping[str, float] | None = None
    fast_ion_volume_average_density_m3_by_species: Mapping[str, float] | None = None

@dataclass(frozen=True)
class ModalLocalReconstruction:
    """Local confined and lost ion reconstructions after applying ϕ(z)

    Distribution arrays use axial position first and local speed second with Λ or pitch as the final coordinate
    local_density_m3 has one value per axial cell
    """
    local_distribution_z_v_lambda: np.ndarray
    local_distribution_z_v_pitch: np.ndarray
    lost_ion_distribution_z_v_lambda: np.ndarray | None
    lost_ion_parallel_temperature_J: float | None
    lost_ion_geometry_profile_G_z: np.ndarray | None
    local_density_m3: np.ndarray
    local_speed_grid: SpeedGrid | None = None
    local_speed_grid_diagnostics: dict[str, object] | None = None

@dataclass(frozen=True)
class ModalBeamComponentResult:
    """Modal source projection and velocity solution for one attenuated beam energy component
    
    Rates use particles s⁻¹, powers use W, and modal arrays use retained mode first
    """
    energy_J: float
    speed_m_s: float
    birth_rate_s: float
    deposited_birth_power_W: float
    confined_birth_rate_s: float
    confined_birth_power_W: float
    prompt_loss_birth_rate_s: float
    prompt_loss_birth_power_W: float
    modal_rhs: np.ndarray
    modal_source_coefficients: np.ndarray
    modal_distribution_f_j_v: np.ndarray
    contributing_axial_cells: int
    birth_lambda_weighted_mean: float
    birth_eta_weighted_mean: float
    first_mode_weighted_mean: float

@dataclass(frozen=True)
class ModalFBISResult:
    """Completed species level modal FBIS state consumed by the integration pipeline
    
    The result preserves invariant and local representations, loss rates, electrostatic closure, collision provenance, and convergence metadata
    """
    species: IonSpecies
    speed_grid: SpeedGrid
    lambda_grid: LambdaGrid
    pitch_grid: PitchGrid | None
    basis: ModalFBISBasis
    component_results: tuple[ModalBeamComponentResult, ...]
    modal_distribution_f_j_v: np.ndarray
    distribution_v_eta: np.ndarray
    distribution_v_lambda: np.ndarray
    density_m3: float
    inventory_particles: float
    ion_particle_loss_rate_s: float
    ion_midplane_kinetic_power_loss_W: float
    ion_wall_power_loss_W: float | None
    electron_wall_power_loss_W: float | None
    wall_barrier_energy_J: float
    normalized_wall_barrier: float
    electron_wall_potential_relative_to_midplane_V: float | None
    dF_dlambda_boundary_v: np.ndarray
    boundary_flux_density_v: np.ndarray
    confinement_time_s: float | None
    rosenbluth_coefficients: ModalRosenbluthCoefficients | None
    eq59_collision_operator_state: Eq59CollisionOperatorState | None
    electrostatic_profile: ModalElectrostaticProfile | None
    local_reconstruction: ModalLocalReconstruction | None
    metadata: dict[str, object]
    eq59_warm_start_state: Eq59WarmStartState | None = None
    current_balance: AmbipolarCurrentBalance | None = None
    fast_ion_electron_heating_state: FastIonElectronHeatingState | None = None
