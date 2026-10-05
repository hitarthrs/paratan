"""Eq 70 local reconstruction and fixed magnetic boundary Eq 63 branch assembly"""
from __future__ import annotations
import numpy as np
from source_model_revamp.fbis.modal.electrostatic.iteration import _solve_phi_profile_quasineutrality
from source_model_revamp.fbis.modal.local.speed import assess_local_physical_speed_grid_convergence, derive_local_physical_speed_grid
from source_model_revamp.fbis.modal.lost_ion_distribution import build_fixed_boundary_eq63_branch, calculate_baldwin_throat_parallel_temperature
from source_model_revamp.fbis.modal.types import ModalElectrostaticProfile, ModalLocalReconstruction
from collections.abc import Mapping
from typing import Any

_FIXED_BOUNDARY_RECONSTRUCTION_REQUIRED = "required"
_FIXED_BOUNDARY_RECONSTRUCTION_DEFERRED = "deferred_eq42_iteration"

def solve_modal_fixed_boundary_reconstruction(context: Mapping[str, Any]) -> dict[str, Any]:
    """Solve the local electrostatic state and build left and right fixed boundary Eq 63 branches when active"""
    B_tilde_midpoints = context['B_tilde_midpoints']
    _eq70_initial_potential_energy_J = context['_eq70_initial_potential_energy_J']
    _evaluation_mode = context['_evaluation_mode']
    _run_numerical_convergence_assessments = context['_run_numerical_convergence_assessments']
    background_left_throat_positive_charge_density_m3 = context['background_left_throat_positive_charge_density_m3']
    background_midplane_positive_charge_density_m3 = context['background_midplane_positive_charge_density_m3']
    background_positive_charge_density_m3 = context['background_positive_charge_density_m3']
    background_right_throat_positive_charge_density_m3 = context['background_right_throat_positive_charge_density_m3']
    barrier = context['barrier']
    basis = context['basis']
    boundary_flux_single_density_v = context['boundary_flux_single_density_v']
    cell_volumes_m3 = context['cell_volumes_m3']
    collision_state = context['collision_state']
    dfdlambda = context['dfdlambda']
    eq59_collision_operator_state = context['eq59_collision_operator_state']
    electron_collision_density_m3 = context['electron_collision_density_m3']
    eq70_prescribed_electron_midplane_density_m3 = context['eq70_prescribed_electron_midplane_density_m3']
    electron_temperature_J = context['electron_temperature_J']
    f_lambda = context['f_lambda']
    fixed_boundary_production = context['fixed_boundary_production']
    fast_ion_species = context['fast_ion_species']
    fixed_boundary_reconstruction_policy = context['fixed_boundary_reconstruction_policy']
    invariant_energy_cell_population_weights = context['invariant_energy_cell_population_weights']
    inventory = context['inventory']
    lambda_grid = context['lambda_grid']
    midplane_area_m2 = context['midplane_area_m2']
    mirror_ratio = context['mirror_ratio']
    numerics = context['numerics']
    pitch_grid = context['pitch_grid']
    speed_grid = context['speed_grid']
    volume_m3 = context['volume_m3']
    zeta_faces = context['zeta_faces']
    electrostatic_profile: ModalElectrostaticProfile | None = None
    local_reconstruction: ModalLocalReconstruction | None = None
    directed_eq63_throat_rate_v_lambda_s: np.ndarray | None = None
    left_throat_temperature = None
    right_throat_temperature = None
    left_eq63_branch = None
    right_eq63_branch = None
    fixed_boundary_reconstruction_status = ("not_applicable_non_fixed_boundary_model")
    fixed_boundary_reconstruction_failure_reason: str | None = None
    # Local Eq 70 to Eq 72 reconstruction is available only when pitch and axial volume grids are supplied
    if pitch_grid is not None and cell_volumes_m3 is not None:
        zeta_centers = 0.5 * (np.asarray(zeta_faces, dtype=float)[:-1] + np.asarray(zeta_faces, dtype=float)[1:])
        eta_to_local_normalization = float(basis.volume_geometry_integral / basis.eta_lambda_map.lambda_normalization)
        electrostatic_profile, local_lambda, local_pitch, local_density = _solve_phi_profile_quasineutrality(
            speed_grid=speed_grid,
            lambda_grid=lambda_grid,
            pitch_grid=pitch_grid,
            base_distribution_v_lambda=f_lambda,
            zeta=zeta_centers,
            B_tilde=np.asarray(B_tilde_midpoints, dtype=float),
            cell_volumes_m3=np.asarray(cell_volumes_m3, dtype=float),
            mirror_ratio=mirror_ratio,
            wall_barrier_energy_J=barrier,
            electron_temperature_J=electron_temperature_J,
            particle_mass_kg=fast_ion_species.mass_kg,
            iterations=numerics.phi_z_iterations,
            relaxation=numerics.phi_z_relaxation,
            relative_tolerance=numerics.phi_z_relative_tolerance,
            electron_midplane_density_m3=eq70_prescribed_electron_midplane_density_m3,
            electron_collision_density_m3=electron_collision_density_m3,
            background_positive_charge_density_m3=background_positive_charge_density_m3,
            background_midplane_positive_charge_density_m3=background_midplane_positive_charge_density_m3,
            background_left_throat_positive_charge_density_m3=background_left_throat_positive_charge_density_m3,
            background_right_throat_positive_charge_density_m3=background_right_throat_positive_charge_density_m3,
            target_volume_averaged_ion_density_m3=None,
            eta_to_local_phase_space_normalization=eta_to_local_normalization,
            local_velocity_quadrature_order=numerics.local_velocity_quadrature_order,
            low_energy_weight_fraction_tolerance=numerics.phi_z_low_energy_weight_fraction_tolerance,
            invariant_energy_cell_population_weights=invariant_energy_cell_population_weights,
            initial_potential_energy_J=_eq70_initial_potential_energy_J,
        )
        local_speed_grid = electrostatic_profile.local_speed_grid
        if local_speed_grid is None:
            raise RuntimeError("Eq. 70 local reconstruction did not expose its physical speed grid")
        _, local_grid_diagnostics = derive_local_physical_speed_grid(invariant_speed_grid=speed_grid, maximum_potential_drop_magnitude_J=barrier, particle_mass_kg=fast_ion_species.mass_kg)
        if _run_numerical_convergence_assessments:
            local_refinement = assess_local_physical_speed_grid_convergence(
                invariant_speed_grid=speed_grid,
                local_speed_grid=local_speed_grid,
                lambda_grid=lambda_grid,
                pitch_grid=pitch_grid,
                base_distribution_v_lambda=f_lambda,
                mirror_ratio=mirror_ratio,
                zeta=zeta_centers,
                B_tilde=np.asarray(B_tilde_midpoints, dtype=float),
                potential_drop_magnitude_J=np.asarray(electrostatic_profile.potential_energy_J, dtype=float),
                left_throat_potential_drop_magnitude_J=float(electrostatic_profile.throat_potential_left_energy_J),
                right_throat_potential_drop_magnitude_J=float(electrostatic_profile.throat_potential_right_energy_J),
                eta_to_local_phase_space_normalization=(eta_to_local_normalization),
                cell_volumes_m3=np.asarray(cell_volumes_m3, dtype=float),
                quadrature_order=numerics.local_velocity_quadrature_order,
                relative_tolerance=(numerics.local_speed_refinement_relative_tolerance),
                reference_local_distribution_z_v_pitch=local_pitch,
                particle_mass_kg=fast_ion_species.mass_kg,
            )
        else:
            local_refinement = {"local_grid_refinement_assessed": False, "local_grid_refinement_converged": False, "local_grid_refinement_status": ("deferred_during_intermediate_closure_iteration"), "local_grid_refinement_failure_reason": ("intermediate_closure_state"), "local_grid_refinement_relative_tolerance": float(numerics.local_speed_refinement_relative_tolerance), "local_grid_refinement_history": [], "eq71_72_population_change_included_in_grid_error": False}
        local_grid_diagnostics = {**local_grid_diagnostics, "remap_inventory_relative_error": ((float(electrostatic_profile.fast_ion_inventory_particles) - float(inventory)) / float(inventory) if float(inventory) > 0.0 else 0.0), "remap_energy_relative_error": None, "remap_energy_status": "not_a_conservative_energy_remap_Eq71_72_is_an_energy_mapping", **local_refinement}
        # The fixed magnetic boundary Eq 63 branches are matched independently to the one end Eq 61 spectrum
        if fixed_boundary_production:
            if (fixed_boundary_reconstruction_policy == _FIXED_BOUNDARY_RECONSTRUCTION_DEFERRED):
                fixed_boundary_reconstruction_status = ("deferred_during_eq42_density_basis_iteration")
            else:
                one_end_eq61_rate_v_s = np.asarray(boundary_flux_single_density_v, dtype=float) * np.asarray(speed_grid.shell_volumes_m3_s3, dtype=float) * float(volume_m3)
                common = dict(
                    speed_grid=speed_grid,
                    one_end_eq61_rate_v_s=one_end_eq61_rate_v_s,
                    dfdlambda_boundary_v=np.asarray(dfdlambda, dtype=float),
                    mirror_ratio=float(mirror_ratio),
                    midplane_area_m2=float(midplane_area_m2),
                    half_length_m=float(context['half_length_m']),
                    geometry_factor_G=float(basis.geometry_factor_G),
                    particle_mass_kg=float(fast_ion_species.mass_kg),
                    collision_state=collision_state,
                    collision_operator_state=eq59_collision_operator_state,
                    volume_geometry_integral=float(basis.volume_geometry_integral),
                    volume_m3=float(volume_m3),
                )
                left_throat_temperature = calculate_baldwin_throat_parallel_temperature(side="left", **common)
                right_throat_temperature = calculate_baldwin_throat_parallel_temperature(side="right", **common)
                left_eq63_branch = build_fixed_boundary_eq63_branch(side="left", speed_grid=speed_grid, lambda_grid=lambda_grid, one_end_eq61_rate_v_s=one_end_eq61_rate_v_s, mirror_ratio=mirror_ratio, midplane_area_m2=midplane_area_m2, geometry_factor_G=basis.geometry_factor_G, parallel_temperature_J=left_throat_temperature.eq63_parallel_temperature_J, particle_mass_kg=fast_ion_species.mass_kg)
                right_eq63_branch = build_fixed_boundary_eq63_branch(side="right", speed_grid=speed_grid, lambda_grid=lambda_grid, one_end_eq61_rate_v_s=one_end_eq61_rate_v_s, mirror_ratio=mirror_ratio, midplane_area_m2=midplane_area_m2, geometry_factor_G=basis.geometry_factor_G, parallel_temperature_J=right_throat_temperature.eq63_parallel_temperature_J, particle_mass_kg=fast_ion_species.mass_kg)
                directed_eq63_throat_rate_v_lambda_s = right_eq63_branch.throat_rate_v_lambda_s
                fixed_boundary_reconstruction_status = "represented"
        representative_T_L = (None if left_throat_temperature is None or right_throat_temperature is None else 0.5 * (left_throat_temperature.eq63_parallel_temperature_J + right_throat_temperature.eq63_parallel_temperature_J))
        local_reconstruction = ModalLocalReconstruction(local_distribution_z_v_lambda=local_lambda, local_distribution_z_v_pitch=local_pitch, lost_ion_distribution_z_v_lambda=None, lost_ion_parallel_temperature_J=representative_T_L, lost_ion_geometry_profile_G_z=None, local_density_m3=local_density, local_speed_grid=local_speed_grid, local_speed_grid_diagnostics=local_grid_diagnostics)
    
    return {
        'directed_eq63_throat_rate_v_lambda_s': directed_eq63_throat_rate_v_lambda_s,
        'electrostatic_profile': electrostatic_profile,
        'fixed_boundary_reconstruction_failure_reason': fixed_boundary_reconstruction_failure_reason,
        'fixed_boundary_reconstruction_status': fixed_boundary_reconstruction_status,
        'left_eq63_branch': left_eq63_branch,
        'left_throat_temperature': left_throat_temperature,
        'local_reconstruction': local_reconstruction,
        'right_eq63_branch': right_eq63_branch,
        'right_throat_temperature': right_throat_temperature,
        'operating_point_evaluation_mode': _evaluation_mode.value,
    }
