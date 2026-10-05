"""Modal FBIS, electrostatic, loss, and expander numerical configuration"""
from __future__ import annotations
from collections.abc import Mapping
from dataclasses import dataclass, fields
from typing import Any
import numpy as np
from source_model_revamp.fbis.modal.density_weighting import EQ42_SELF_CONSISTENT_QUASINEUTRAL_DENSITY, eq42_density_weighting_model
from source_model_revamp.fbis.modal.models import electrostatic_feedback_model as canonical_electrostatic_feedback_model, velocity_solution_model as canonical_velocity_solution_model
from source_model_revamp.integration.config.common import canonical_model_name, _check_unknown, _number, _integer
from source_model_revamp.integration.config.constants import *
from source_model_revamp.radial_profiles import KOTELNIKOV_PARABOLIC_FLUX_K1, canonical_radial_profile_model

@dataclass(frozen=True)
class KineticElectrostaticConfig:
    """
    Numerical controls shared by the modal FBIS and coupled electrostatic path
    
    The fields cover speed and pitch grids, physical basis construction, Eq 59 iteration, Eq 42 density weighting, Eq 70 potential iteration, fixed magnetic boundary losses, current balance, expander closure, and convergence checks
    """
    speed_bins: int = 8
    lambda_bins: int = 10
    pitch_bins: int = 6
    speed_max_factor: float = 1.4
    modal_speed_core_cell_fraction: float = 0.7
    modal_speed_tail_stretch_power: float = 2.0
    negative_distribution_roundoff_tolerance: float = 1.0e-12
    modal_reconstruction_negative_particle_fraction_tolerance: float = 1.0e-3
    modal_reconstruction_energy_correction_relative_tolerance: float = 5.0e-3
    ion_loss_closure_model: str = "egedal_hot_fixed_magnetic_boundary"
    lost_ion_parallel_temperature_model: str = "baldwin_1972_throat_density_closure"
    modal_n_lambda_grid: int = 61
    modal_n_eta_grid: int = 45
    modal_n_square_basis_modes: int = 6
    modal_n_physical_modes: int = 1
    modal_square_scan_points: int = 4000
    modal_velocity_solution_model: str = "hot_ion_rosenbluth_eq59"
    modal_rosenbluth_max_iterations: int = 30
    modal_rosenbluth_relative_tolerance: float = 1.0e-4
    modal_rosenbluth_absolute_tolerance: float = 0.0
    modal_rosenbluth_relaxation: float = 0.8
    modal_rosenbluth_min_iterations: int = 2
    modal_rosenbluth_high_speed_tail_population_tolerance: float = 1.0e-4
    modal_rosenbluth_high_speed_tail_energy_tolerance: float = 1.0e-4
    modal_rosenbluth_speed_convergence_relative_tolerance: float = 1.0e-2
    modal_retained_mode_relative_tolerance: float = 1.0e-2
    modal_local_speed_refinement_relative_tolerance: float = 1.0e-2
    modal_phi_z_iterations: int = 32
    modal_phi_z_root_scan_points: int = 33
    modal_phi_z_relative_tolerance: float = 1.0e-3
    modal_phi_z_relaxation: float = 0.8
    modal_local_velocity_quadrature_order: int = 5
    modal_phi_z_low_energy_weight_fraction_tolerance: float = 1.0e-3
    modal_basis_convergence_levels: int = 3
    modal_basis_eigenvalue_relative_tolerance: float = 5.0e-3
    modal_basis_eigenfunction_overlap_tolerance: float = 0.995
    modal_loss_convention_relative_tolerance: float = 1.0e-2
    modal_source_projection_relative_tolerance: float = 1.0e-2
    modal_electrostatic_feedback_iterations: int = 4
    modal_electrostatic_feedback_relative_tolerance: float = 1.0e-3
    modal_density_closure_model: str = "nbi_supported_stationary"
    modal_density_closure_iterations: int = 4
    modal_density_closure_relative_tolerance: float = 1.0e-3
    modal_density_iteration_relative_tolerance: float | None = None
    modal_density_closure_relaxation: float = 1.0
    modal_eq42_density_weighting_model: str = EQ42_SELF_CONSISTENT_QUASINEUTRAL_DENSITY
    modal_eq42_density_basis_iterations: int = 4
    modal_eq42_density_shape_relative_tolerance: float = 1.0e-3
    modal_eq42_density_shape_relaxation: float = 0.5
    modal_eq42_basis_eigenvalue_relative_tolerance: float = 5.0e-3
    modal_eq42_basis_eigenfunction_overlap_tolerance: float = 0.995
    modal_eq42_density_symmetry_tolerance: float = 1.0e-2
    electrostatic_feedback_model: str = "egedal_2022_published_electrostatic_fbis_approximation"
    total_current_balance_iterations: int = 4
    total_current_balance_relative_tolerance: float = 1.0e-3
    total_current_balance_relaxation: float = 1.0
    total_current_balance_end_asymmetry_relative_tolerance: float = 5.0e-2
    expander_potential_iterations: int = 16
    expander_potential_relative_tolerance: float = 1.0e-3
    expander_fast_throat_pitch_quantile_nodes: int = 64
    expander_radial_profile_model: str = KOTELNIKOV_PARABOLIC_FLUX_K1
    expander_fast_nonterminal_fraction_tolerance: float = 1.0e-3
    expander_unclosed_population_fraction_tolerance: float = 1.0e-3
    fast_fusion_burnup_relative_tolerance: float = 1.0e-2

    def __post_init__(self) -> None:
        """
        Validate cross parameter numerical constraints and canonicalize the expander radial profile model
        
        This includes basis mode ordering, relaxation ranges, overlap bounds, odd Eq 70 root scan size, and nonnegative convergence tolerances
        """
        if self.modal_n_square_basis_modes < self.modal_n_physical_modes:
            raise ValueError("modal_n_square_basis_modes must be at least modal_n_physical_modes")
        for name, value in (("modal_reconstruction_negative_particle_fraction_tolerance", self.modal_reconstruction_negative_particle_fraction_tolerance), ("modal_reconstruction_energy_correction_relative_tolerance", self.modal_reconstruction_energy_correction_relative_tolerance)):
            if not np.isfinite(value) or value < 0.0:
                raise ValueError(f"{name} must be finite and nonnegative")
        if self.modal_rosenbluth_min_iterations > self.modal_rosenbluth_max_iterations:
            raise ValueError("modal_rosenbluth_min_iterations must not exceed modal_rosenbluth_max_iterations")
        if not 0.0 < self.modal_rosenbluth_relaxation <= 1.0:
            raise ValueError("modal_rosenbluth_relaxation must lie inside (0, 1]")
        if self.modal_basis_convergence_levels != 1 and self.modal_basis_convergence_levels < 3:
            raise ValueError("modal_basis_convergence_levels must be one or at least three")
        if not 0.0 < self.modal_basis_eigenfunction_overlap_tolerance <= 1.0:
            raise ValueError("modal_basis_eigenfunction_overlap_tolerance must lie inside (0, 1]")
        if not 0.0 < self.modal_speed_core_cell_fraction < 1.0:
            raise ValueError("modal_speed_core_cell_fraction must lie inside (0, 1)")
        if self.modal_speed_tail_stretch_power < 1.0:
            raise ValueError("modal_speed_tail_stretch_power must be at least one")
        for name, value in (
            ("modal_rosenbluth_high_speed_tail_energy_tolerance", self.modal_rosenbluth_high_speed_tail_energy_tolerance),
            ("modal_rosenbluth_speed_convergence_relative_tolerance", self.modal_rosenbluth_speed_convergence_relative_tolerance),
            ("modal_retained_mode_relative_tolerance", self.modal_retained_mode_relative_tolerance),
            ("modal_local_speed_refinement_relative_tolerance", self.modal_local_speed_refinement_relative_tolerance)):
            if not np.isfinite(value) or value <= 0.0:
                raise ValueError(f"{name} must be finite and positive")
        if self.modal_density_iteration_relative_tolerance is not None and (not np.isfinite(self.modal_density_iteration_relative_tolerance) or self.modal_density_iteration_relative_tolerance < 0.0):
            raise ValueError("modal_density_iteration_relative_tolerance must be finite and nonnegative when specified")
        eq42_density_weighting_model(self.modal_eq42_density_weighting_model)
        if self.modal_eq42_density_basis_iterations < 1:
            raise ValueError("modal_eq42_density_basis_iterations must be at least one")
        if (not np.isfinite(self.modal_eq42_density_shape_relative_tolerance) or self.modal_eq42_density_shape_relative_tolerance < 0.0):
            raise ValueError("modal_eq42_density_shape_relative_tolerance must be finite and nonnegative")
        if not 0.0 < self.modal_eq42_density_shape_relaxation <= 1.0:
            raise ValueError("modal_eq42_density_shape_relaxation must lie inside (0, 1]")
        if (not np.isfinite(self.modal_eq42_basis_eigenvalue_relative_tolerance) or self.modal_eq42_basis_eigenvalue_relative_tolerance < 0.0):
            raise ValueError("modal_eq42_basis_eigenvalue_relative_tolerance must be finite and nonnegative")
        if not 0.0 < self.modal_eq42_basis_eigenfunction_overlap_tolerance <= 1.0:
            raise ValueError("modal_eq42_basis_eigenfunction_overlap_tolerance must lie inside (0, 1]")
        if (not np.isfinite(self.modal_eq42_density_symmetry_tolerance) or self.modal_eq42_density_symmetry_tolerance < 0.0):
            raise ValueError("modal_eq42_density_symmetry_tolerance must be finite and nonnegative")
        if self.modal_phi_z_root_scan_points < 5 or self.modal_phi_z_root_scan_points % 2 == 0:
            raise ValueError("modal_phi_z_root_scan_points must be an odd integer of at least five")
        if self.total_current_balance_iterations < 1:
            raise ValueError("total_current_balance_iterations must be at least one")
        if not np.isfinite(self.total_current_balance_relative_tolerance) or self.total_current_balance_relative_tolerance < 0.0:
            raise ValueError("total_current_balance_relative_tolerance must be finite and nonnegative")
        if not 0.0 < self.total_current_balance_relaxation <= 1.0:
            raise ValueError("total_current_balance_relaxation must lie inside (0, 1]")
        if not np.isfinite(self.total_current_balance_end_asymmetry_relative_tolerance) or self.total_current_balance_end_asymmetry_relative_tolerance < 0.0:
            raise ValueError("total_current_balance_end_asymmetry_relative_tolerance must be finite and nonnegative")
        if self.expander_potential_iterations < 1:
            raise ValueError("expander_potential_iterations must be at least one")
        if not np.isfinite(self.expander_potential_relative_tolerance) or self.expander_potential_relative_tolerance < 0.0:
            raise ValueError("expander_potential_relative_tolerance must be finite and nonnegative")
        object.__setattr__(self, "expander_radial_profile_model", canonical_radial_profile_model(self.expander_radial_profile_model))
        if self.expander_fast_throat_pitch_quantile_nodes < 2:
            raise ValueError("expander_fast_throat_pitch_quantile_nodes must be at least two")
        if not np.isfinite(self.expander_fast_nonterminal_fraction_tolerance) or self.expander_fast_nonterminal_fraction_tolerance < 0.0:
            raise ValueError("expander_fast_nonterminal_fraction_tolerance must be finite and nonnegative")
        if not np.isfinite(self.expander_unclosed_population_fraction_tolerance) or self.expander_unclosed_population_fraction_tolerance < 0.0:
            raise ValueError("expander_unclosed_population_fraction_tolerance must be finite and nonnegative")
        if not np.isfinite(self.fast_fusion_burnup_relative_tolerance) or self.fast_fusion_burnup_relative_tolerance < 0.0:
            raise ValueError("fast_fusion_burnup_relative_tolerance must be finite and nonnegative")

    @classmethod
    def from_mapping(cls, data: Mapping[str, Any], *, strict: bool = False) -> "KineticElectrostaticConfig":
        """
        Parse the modal kinetic and electrostatic numerical controls
        
        Most numerical controls are required explicitly while selected model names and newer optional controls retain canonical defaults
        """
        allowed = {item.name for item in fields(cls)}
        _check_unknown(data, allowed, "source_model.kinetic_electrostatic", strict)
        optional = {"ion_loss_closure_model", "lost_ion_parallel_temperature_model", "modal_density_closure_model", "modal_density_iteration_relative_tolerance", "modal_eq42_density_weighting_model", "electrostatic_feedback_model", "modal_phi_z_root_scan_points", "expander_fast_throat_pitch_quantile_nodes", "expander_radial_profile_model"}
        missing = sorted((allowed - optional) - set(data))
        if missing:
            raise ValueError("source_model.kinetic_electrostatic is missing required fields: " + ", ".join(missing))
        modal_density_closure_model = canonical_model_name(data.get("modal_density_closure_model"), MODAL_DENSITY_CLOSURE_ALIASES, "nbi_supported_stationary")
        modal_eq42_density_weighting_model = eq42_density_weighting_model(data.get("modal_eq42_density_weighting_model", EQ42_SELF_CONSISTENT_QUASINEUTRAL_DENSITY))
        ion_loss_closure_model = str(data.get("ion_loss_closure_model", "egedal_hot_fixed_magnetic_boundary")).strip().lower()
        if ion_loss_closure_model != "egedal_hot_fixed_magnetic_boundary":
            raise ValueError("kinetic_electrostatic.ion_loss_closure_model must be egedal_hot_fixed_magnetic_boundary")
        lost_ion_parallel_temperature_model = str(data.get("lost_ion_parallel_temperature_model", "baldwin_1972_throat_density_closure")).strip().lower()
        if lost_ion_parallel_temperature_model != "baldwin_1972_throat_density_closure":
            raise ValueError("kinetic_electrostatic.lost_ion_parallel_temperature_model must be baldwin_1972_throat_density_closure")
        electrostatic_feedback_model = canonical_electrostatic_feedback_model(data.get("electrostatic_feedback_model", "egedal_2022_published_electrostatic_fbis_approximation"))

        return cls(speed_bins=_integer(data.get('speed_bins'), 'kinetic.speed_bins', minimum=2), lambda_bins=_integer(data.get('lambda_bins'), 'kinetic.lambda_bins', minimum=2), pitch_bins=_integer(data.get('pitch_bins'), 'kinetic.pitch_bins', minimum=3), speed_max_factor=_number(data.get('speed_max_factor'), 'kinetic.speed_max_factor', positive=True), modal_speed_core_cell_fraction=_number(data.get('modal_speed_core_cell_fraction'), 'kinetic.modal_speed_core_cell_fraction', positive=True), modal_speed_tail_stretch_power=_number(data.get('modal_speed_tail_stretch_power'), 'kinetic.modal_speed_tail_stretch_power', positive=True), negative_distribution_roundoff_tolerance=_number(data.get('negative_distribution_roundoff_tolerance'), 'kinetic.negative_distribution_roundoff_tolerance', nonnegative=True), modal_reconstruction_negative_particle_fraction_tolerance=_number(data.get('modal_reconstruction_negative_particle_fraction_tolerance'), 'kinetic.modal_reconstruction_negative_particle_fraction_tolerance', nonnegative=True), modal_reconstruction_energy_correction_relative_tolerance=_number(data.get('modal_reconstruction_energy_correction_relative_tolerance'), 'kinetic.modal_reconstruction_energy_correction_relative_tolerance', nonnegative=True), ion_loss_closure_model=ion_loss_closure_model, lost_ion_parallel_temperature_model=lost_ion_parallel_temperature_model, modal_n_lambda_grid=_integer(data.get('modal_n_lambda_grid'), 'kinetic.modal_n_lambda_grid', minimum=5), modal_n_eta_grid=_integer(data.get('modal_n_eta_grid'), 'kinetic.modal_n_eta_grid', minimum=5), modal_n_square_basis_modes=_integer(data.get('modal_n_square_basis_modes'), 'kinetic.modal_n_square_basis_modes', minimum=1), modal_n_physical_modes=_integer(data.get('modal_n_physical_modes'), 'kinetic.modal_n_physical_modes', minimum=1), modal_square_scan_points=_integer(data.get('modal_square_scan_points'), 'kinetic.modal_square_scan_points', minimum=100), modal_velocity_solution_model=canonical_velocity_solution_model(data.get('modal_velocity_solution_model')), modal_rosenbluth_max_iterations=_integer(data.get('modal_rosenbluth_max_iterations'), 'kinetic.modal_rosenbluth_max_iterations', minimum=1), modal_rosenbluth_relative_tolerance=_number(data.get('modal_rosenbluth_relative_tolerance'), 'kinetic.modal_rosenbluth_relative_tolerance', positive=True), modal_rosenbluth_absolute_tolerance=_number(data.get('modal_rosenbluth_absolute_tolerance'), 'kinetic.modal_rosenbluth_absolute_tolerance', nonnegative=True), modal_rosenbluth_relaxation=_number(data.get('modal_rosenbluth_relaxation'), 'kinetic.modal_rosenbluth_relaxation', positive=True), modal_rosenbluth_min_iterations=_integer(data.get('modal_rosenbluth_min_iterations'), 'kinetic.modal_rosenbluth_min_iterations', minimum=1), modal_rosenbluth_high_speed_tail_population_tolerance=_number(data.get('modal_rosenbluth_high_speed_tail_population_tolerance'), 'kinetic.modal_rosenbluth_high_speed_tail_population_tolerance', nonnegative=True), modal_rosenbluth_high_speed_tail_energy_tolerance=_number(data.get('modal_rosenbluth_high_speed_tail_energy_tolerance'), 'kinetic.modal_rosenbluth_high_speed_tail_energy_tolerance', nonnegative=True), modal_rosenbluth_speed_convergence_relative_tolerance=_number(data.get('modal_rosenbluth_speed_convergence_relative_tolerance'), 'kinetic.modal_rosenbluth_speed_convergence_relative_tolerance', positive=True), modal_retained_mode_relative_tolerance=_number(data.get('modal_retained_mode_relative_tolerance'), 'kinetic.modal_retained_mode_relative_tolerance', positive=True), modal_local_speed_refinement_relative_tolerance=_number(data.get('modal_local_speed_refinement_relative_tolerance'), 'kinetic.modal_local_speed_refinement_relative_tolerance', positive=True), modal_phi_z_iterations=_integer(data.get('modal_phi_z_iterations'), 'kinetic.modal_phi_z_iterations', minimum=1), modal_phi_z_root_scan_points=_integer(data.get('modal_phi_z_root_scan_points', 33), 'kinetic.modal_phi_z_root_scan_points', minimum=5), modal_phi_z_relative_tolerance=_number(data.get('modal_phi_z_relative_tolerance'), 'kinetic.modal_phi_z_relative_tolerance', nonnegative=True), modal_phi_z_relaxation=_number(data.get('modal_phi_z_relaxation'), 'kinetic.modal_phi_z_relaxation', nonnegative=True), modal_local_velocity_quadrature_order=_integer(data.get('modal_local_velocity_quadrature_order'), 'kinetic.modal_local_velocity_quadrature_order', minimum=2), modal_phi_z_low_energy_weight_fraction_tolerance=_number(data.get('modal_phi_z_low_energy_weight_fraction_tolerance'), 'kinetic.modal_phi_z_low_energy_weight_fraction_tolerance', nonnegative=True), modal_basis_convergence_levels=_integer(data.get('modal_basis_convergence_levels'), 'kinetic.modal_basis_convergence_levels', minimum=1), modal_basis_eigenvalue_relative_tolerance=_number(data.get('modal_basis_eigenvalue_relative_tolerance'), 'kinetic.modal_basis_eigenvalue_relative_tolerance', positive=True), modal_basis_eigenfunction_overlap_tolerance=_number(data.get('modal_basis_eigenfunction_overlap_tolerance'), 'kinetic.modal_basis_eigenfunction_overlap_tolerance', positive=True), modal_loss_convention_relative_tolerance=_number(data.get('modal_loss_convention_relative_tolerance'), 'kinetic.modal_loss_convention_relative_tolerance', nonnegative=True), modal_source_projection_relative_tolerance=_number(data.get('modal_source_projection_relative_tolerance'), 'kinetic.modal_source_projection_relative_tolerance', nonnegative=True), modal_electrostatic_feedback_iterations=_integer(data.get('modal_electrostatic_feedback_iterations'), 'kinetic.modal_electrostatic_feedback_iterations', minimum=1), modal_electrostatic_feedback_relative_tolerance=_number(data.get('modal_electrostatic_feedback_relative_tolerance'), 'kinetic.modal_electrostatic_feedback_relative_tolerance', nonnegative=True), modal_density_closure_model=modal_density_closure_model, modal_density_closure_iterations=_integer(data.get('modal_density_closure_iterations'), 'kinetic.modal_density_closure_iterations', minimum=1), modal_density_closure_relative_tolerance=_number(data.get('modal_density_closure_relative_tolerance'), 'kinetic.modal_density_closure_relative_tolerance', nonnegative=True), modal_density_iteration_relative_tolerance=None if data.get('modal_density_iteration_relative_tolerance') is None else _number(data.get('modal_density_iteration_relative_tolerance'), 'kinetic.modal_density_iteration_relative_tolerance', nonnegative=True), modal_density_closure_relaxation=_number(data.get('modal_density_closure_relaxation'), 'kinetic.modal_density_closure_relaxation', positive=True), modal_eq42_density_weighting_model=modal_eq42_density_weighting_model, modal_eq42_density_basis_iterations=_integer(data.get('modal_eq42_density_basis_iterations'), 'kinetic.modal_eq42_density_basis_iterations', minimum=1), modal_eq42_density_shape_relative_tolerance=_number(data.get('modal_eq42_density_shape_relative_tolerance'), 'kinetic.modal_eq42_density_shape_relative_tolerance', nonnegative=True), modal_eq42_density_shape_relaxation=_number(data.get('modal_eq42_density_shape_relaxation'), 'kinetic.modal_eq42_density_shape_relaxation', positive=True), modal_eq42_basis_eigenvalue_relative_tolerance=_number(data.get('modal_eq42_basis_eigenvalue_relative_tolerance'), 'kinetic.modal_eq42_basis_eigenvalue_relative_tolerance', nonnegative=True), modal_eq42_basis_eigenfunction_overlap_tolerance=_number(data.get('modal_eq42_basis_eigenfunction_overlap_tolerance'), 'kinetic.modal_eq42_basis_eigenfunction_overlap_tolerance', positive=True), modal_eq42_density_symmetry_tolerance=_number(data.get('modal_eq42_density_symmetry_tolerance'), 'kinetic.modal_eq42_density_symmetry_tolerance', nonnegative=True), electrostatic_feedback_model=electrostatic_feedback_model, total_current_balance_iterations=_integer(data.get('total_current_balance_iterations'), 'kinetic.total_current_balance_iterations', minimum=1), total_current_balance_relative_tolerance=_number(data.get('total_current_balance_relative_tolerance'), 'kinetic.total_current_balance_relative_tolerance', nonnegative=True), total_current_balance_relaxation=_number(data.get('total_current_balance_relaxation'), 'kinetic.total_current_balance_relaxation', positive=True), total_current_balance_end_asymmetry_relative_tolerance=_number(data.get('total_current_balance_end_asymmetry_relative_tolerance'), 'kinetic.total_current_balance_end_asymmetry_relative_tolerance', nonnegative=True), expander_potential_iterations=_integer(data.get('expander_potential_iterations'), 'kinetic.expander_potential_iterations', minimum=1), expander_potential_relative_tolerance=_number(data.get('expander_potential_relative_tolerance'), 'kinetic.expander_potential_relative_tolerance', nonnegative=True), expander_fast_throat_pitch_quantile_nodes=_integer(data.get('expander_fast_throat_pitch_quantile_nodes', 64), 'kinetic.expander_fast_throat_pitch_quantile_nodes', minimum=2), expander_radial_profile_model=canonical_radial_profile_model(data.get('expander_radial_profile_model', KOTELNIKOV_PARABOLIC_FLUX_K1)), expander_fast_nonterminal_fraction_tolerance=_number(data.get('expander_fast_nonterminal_fraction_tolerance'), 'kinetic.expander_fast_nonterminal_fraction_tolerance', nonnegative=True), expander_unclosed_population_fraction_tolerance=_number(data.get('expander_unclosed_population_fraction_tolerance'), 'kinetic.expander_unclosed_population_fraction_tolerance', nonnegative=True), fast_fusion_burnup_relative_tolerance=_number(data.get('fast_fusion_burnup_relative_tolerance'), 'kinetic.fast_fusion_burnup_relative_tolerance', nonnegative=True))

    def backend_numerics_metadata(self) -> dict[str, Any]:
        """Return the modal, electrostatic, loss, and expander numerical controls recorded with solver metadata"""
        return {'speed_max_factor': self.speed_max_factor, 'modal_speed_core_cell_fraction': self.modal_speed_core_cell_fraction, 'modal_speed_tail_stretch_power': self.modal_speed_tail_stretch_power, 'negative_distribution_roundoff_tolerance': self.negative_distribution_roundoff_tolerance, 'modal_reconstruction_negative_particle_fraction_tolerance': self.modal_reconstruction_negative_particle_fraction_tolerance, 'modal_reconstruction_energy_correction_relative_tolerance': self.modal_reconstruction_energy_correction_relative_tolerance, 'ion_loss_closure_model': self.ion_loss_closure_model, 'lost_ion_parallel_temperature_model': self.lost_ion_parallel_temperature_model, 'modal_n_lambda_grid': self.modal_n_lambda_grid, 'modal_n_eta_grid': self.modal_n_eta_grid, 'modal_n_square_basis_modes': self.modal_n_square_basis_modes, 'modal_n_physical_modes': self.modal_n_physical_modes, 'modal_square_scan_points': self.modal_square_scan_points, 'modal_velocity_solution_model': self.modal_velocity_solution_model, 'modal_rosenbluth_max_iterations': self.modal_rosenbluth_max_iterations, 'modal_rosenbluth_relative_tolerance': self.modal_rosenbluth_relative_tolerance, 'modal_rosenbluth_absolute_tolerance': self.modal_rosenbluth_absolute_tolerance, 'modal_rosenbluth_relaxation': self.modal_rosenbluth_relaxation, 'modal_rosenbluth_min_iterations': self.modal_rosenbluth_min_iterations, 'modal_rosenbluth_high_speed_tail_population_tolerance': self.modal_rosenbluth_high_speed_tail_population_tolerance, 'modal_rosenbluth_high_speed_tail_energy_tolerance': self.modal_rosenbluth_high_speed_tail_energy_tolerance, 'modal_rosenbluth_speed_convergence_relative_tolerance': self.modal_rosenbluth_speed_convergence_relative_tolerance, 'modal_retained_mode_relative_tolerance': self.modal_retained_mode_relative_tolerance, 'modal_local_speed_refinement_relative_tolerance': self.modal_local_speed_refinement_relative_tolerance, 'modal_phi_z_iterations': self.modal_phi_z_iterations, 'modal_phi_z_root_scan_points': self.modal_phi_z_root_scan_points, 'modal_phi_z_relative_tolerance': self.modal_phi_z_relative_tolerance, 'modal_phi_z_relaxation': self.modal_phi_z_relaxation, 'modal_local_velocity_quadrature_order': self.modal_local_velocity_quadrature_order, 'modal_phi_z_low_energy_weight_fraction_tolerance': self.modal_phi_z_low_energy_weight_fraction_tolerance, 'modal_basis_convergence_levels': self.modal_basis_convergence_levels, 'modal_basis_eigenvalue_relative_tolerance': self.modal_basis_eigenvalue_relative_tolerance, 'modal_basis_eigenfunction_overlap_tolerance': self.modal_basis_eigenfunction_overlap_tolerance, 'modal_loss_convention_relative_tolerance': self.modal_loss_convention_relative_tolerance, 'modal_source_projection_relative_tolerance': self.modal_source_projection_relative_tolerance, 'modal_electrostatic_feedback_iterations': self.modal_electrostatic_feedback_iterations, 'modal_electrostatic_feedback_relative_tolerance': self.modal_electrostatic_feedback_relative_tolerance, 'modal_density_closure_model': self.modal_density_closure_model, 'modal_density_closure_iterations': self.modal_density_closure_iterations, 'modal_density_closure_relative_tolerance': self.modal_density_closure_relative_tolerance, 'modal_density_iteration_relative_tolerance': self.modal_density_iteration_relative_tolerance, 'modal_density_closure_relaxation': self.modal_density_closure_relaxation, 'modal_eq42_density_weighting_model': self.modal_eq42_density_weighting_model, 'modal_eq42_density_basis_iterations': self.modal_eq42_density_basis_iterations, 'modal_eq42_density_shape_relative_tolerance': self.modal_eq42_density_shape_relative_tolerance, 'modal_eq42_density_shape_relaxation': self.modal_eq42_density_shape_relaxation, 'modal_eq42_basis_eigenvalue_relative_tolerance': self.modal_eq42_basis_eigenvalue_relative_tolerance, 'modal_eq42_basis_eigenfunction_overlap_tolerance': self.modal_eq42_basis_eigenfunction_overlap_tolerance, 'modal_eq42_density_symmetry_tolerance': self.modal_eq42_density_symmetry_tolerance, 'electrostatic_feedback_model': self.electrostatic_feedback_model, 'total_current_balance_iterations': self.total_current_balance_iterations, 'total_current_balance_relative_tolerance': self.total_current_balance_relative_tolerance, 'total_current_balance_relaxation': self.total_current_balance_relaxation, 'total_current_balance_end_asymmetry_relative_tolerance': self.total_current_balance_end_asymmetry_relative_tolerance, 'expander_potential_iterations': self.expander_potential_iterations, 'expander_potential_relative_tolerance': self.expander_potential_relative_tolerance, 'expander_fast_throat_pitch_quantile_nodes': self.expander_fast_throat_pitch_quantile_nodes, 'expander_radial_profile_model': self.expander_radial_profile_model, 'expander_fast_nonterminal_fraction_tolerance': self.expander_fast_nonterminal_fraction_tolerance, 'expander_unclosed_population_fraction_tolerance': self.expander_unclosed_population_fraction_tolerance, 'fast_fusion_burnup_relative_tolerance': self.fast_fusion_burnup_relative_tolerance}
