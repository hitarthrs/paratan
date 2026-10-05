"""
Modal numerical controls, Eq 42 density profiles, basis construction, and collision state helpers

This file connects integration configuration and geometry to the reusable modal basis used by the single species and coupled species stage paths
"""
from __future__ import annotations
from dataclasses import replace
import numpy as np
from source_model_revamp.fbis.collision_parameters import FBISCollisionParameterState, build_fbis_collision_parameter_state
from source_model_revamp.fbis.species import DEUTERON, IonSpecies
from source_model_revamp.fbis.modal.basis import assess_modal_fbis_basis_convergence
from source_model_revamp.fbis.modal.density_weighting import EQ42_PRESCRIBED_AXIAL_DENSITY_PROFILE, EQ42_SELF_CONSISTENT_QUASINEUTRAL_DENSITY, Eq42DensityProfile, build_eq42_density_profile, eq42_density_weighting_model, uniform_eq42_density_profile
from source_model_revamp.fbis.modal.electrostatic_feedback import rescaled_magnetic_profile_for_reference_ratio
from source_model_revamp.fbis.modal import ModalFBISNumerics, build_modal_fbis_basis
from source_model_revamp.integration.pipeline_types import GeometryStageResult
from source_model_revamp.integration.config.root import SourceModelRunConfig
from source_model_revamp.integration.modal_stage.types import ModalBasisCache, _ModalBasisCacheEntry, _ReusableModalBasisBundle

def _modal_numerics(config: SourceModelRunConfig) -> ModalFBISNumerics:
    """
    Build ModalFBISNumerics from the kinetic configuration
    
    The returned controls cover Λ and η grids, Eq 59 iteration, Eq 70 reconstruction, Eq 42 density coupling, basis checks, loss checks, and local reconstruction checks
    """
    k = config.kinetic_electrostatic
    defaults = ModalFBISNumerics()
    return ModalFBISNumerics(
        n_lambda_grid=k.modal_n_lambda_grid,
        n_eta_grid=k.modal_n_eta_grid,
        n_square_basis_modes=k.modal_n_square_basis_modes,
        n_physical_modes=k.modal_n_physical_modes,
        square_scan_points=k.modal_square_scan_points,
        speed_core_cell_fraction=float(getattr(k, "modal_speed_core_cell_fraction", 0.7)),
        speed_tail_stretch_power=float(getattr(k, "modal_speed_tail_stretch_power", 2.0)),
        velocity_solution_model=k.modal_velocity_solution_model,
        rosenbluth_max_iterations=int(k.modal_rosenbluth_max_iterations),
        rosenbluth_relative_tolerance=float(getattr(k, "modal_rosenbluth_relative_tolerance", defaults.rosenbluth_relative_tolerance)),
        rosenbluth_absolute_tolerance=float(getattr(k, "modal_rosenbluth_absolute_tolerance", defaults.rosenbluth_absolute_tolerance)),
        rosenbluth_relaxation=float(getattr(k, "modal_rosenbluth_relaxation", defaults.rosenbluth_relaxation)),
        rosenbluth_min_iterations=int(k.modal_rosenbluth_min_iterations),
        rosenbluth_high_speed_tail_population_tolerance=float(getattr(k, "modal_rosenbluth_high_speed_tail_population_tolerance", defaults.rosenbluth_high_speed_tail_population_tolerance)),
        rosenbluth_high_speed_tail_energy_tolerance=float(getattr(k, "modal_rosenbluth_high_speed_tail_energy_tolerance", defaults.rosenbluth_high_speed_tail_energy_tolerance)),
        rosenbluth_speed_convergence_relative_tolerance=float(getattr(k, "modal_rosenbluth_speed_convergence_relative_tolerance", defaults.rosenbluth_speed_convergence_relative_tolerance)),
        retained_mode_relative_tolerance=float(getattr(k, "modal_retained_mode_relative_tolerance", defaults.retained_mode_relative_tolerance)),
        local_speed_refinement_relative_tolerance=float(getattr(k, "modal_local_speed_refinement_relative_tolerance", defaults.local_speed_refinement_relative_tolerance)),
        phi_z_iterations=k.modal_phi_z_iterations,
        phi_z_root_scan_points=int(getattr(k, "modal_phi_z_root_scan_points", defaults.phi_z_root_scan_points)),
        phi_z_relative_tolerance=k.modal_phi_z_relative_tolerance,
        phi_z_relaxation=k.modal_phi_z_relaxation,
        local_velocity_quadrature_order=int(getattr(k, "modal_local_velocity_quadrature_order", defaults.local_velocity_quadrature_order)),
        phi_z_low_energy_weight_fraction_tolerance=float(getattr(k, "modal_phi_z_low_energy_weight_fraction_tolerance", defaults.phi_z_low_energy_weight_fraction_tolerance)),
        electrostatic_feedback_iterations=int(getattr(k, "modal_electrostatic_feedback_iterations", defaults.electrostatic_feedback_iterations)),
        electrostatic_feedback_relative_tolerance=float(getattr(k, "modal_electrostatic_feedback_relative_tolerance", defaults.electrostatic_feedback_relative_tolerance)),
        basis_convergence_levels=int(getattr(k, "modal_basis_convergence_levels", defaults.basis_convergence_levels)),
        basis_eigenvalue_relative_tolerance=float(getattr(k, "modal_basis_eigenvalue_relative_tolerance", defaults.basis_eigenvalue_relative_tolerance)),
        basis_eigenfunction_overlap_tolerance=float(getattr(k, "modal_basis_eigenfunction_overlap_tolerance", defaults.basis_eigenfunction_overlap_tolerance)),
        loss_convention_relative_tolerance=float(getattr(k, "modal_loss_convention_relative_tolerance", defaults.loss_convention_relative_tolerance)),
        source_projection_relative_tolerance=float(getattr(k, "modal_source_projection_relative_tolerance", defaults.source_projection_relative_tolerance)),
        eq42_density_basis_iterations=int(getattr(k, "modal_eq42_density_basis_iterations", defaults.eq42_density_basis_iterations)),
        eq42_density_shape_relative_tolerance=float(getattr(k, "modal_eq42_density_shape_relative_tolerance", defaults.eq42_density_shape_relative_tolerance)),
        eq42_density_shape_relaxation=float(getattr(k, "modal_eq42_density_shape_relaxation", defaults.eq42_density_shape_relaxation)),
        eq42_basis_eigenvalue_relative_tolerance=float(getattr(k, "modal_eq42_basis_eigenvalue_relative_tolerance", defaults.eq42_basis_eigenvalue_relative_tolerance)),
        eq42_basis_eigenfunction_overlap_tolerance=float(getattr(k, "modal_eq42_basis_eigenfunction_overlap_tolerance", defaults.eq42_basis_eigenfunction_overlap_tolerance)),
        eq42_density_symmetry_tolerance=float(getattr(k, "modal_eq42_density_symmetry_tolerance", defaults.eq42_density_symmetry_tolerance)),
        negative_distribution_roundoff_tolerance=float(k.negative_distribution_roundoff_tolerance),
        reconstruction_negative_particle_fraction_tolerance=float(k.modal_reconstruction_negative_particle_fraction_tolerance),
        reconstruction_energy_correction_relative_tolerance=float(k.modal_reconstruction_energy_correction_relative_tolerance),
    )

def _seed_reference_positive_charge_at(geometry: GeometryStageResult, z_m: np.ndarray) -> np.ndarray:
    """Evaluate configured D plus T positive charge number density at machine coordinates in m⁻³"""
    background = geometry.background_profiles
    if background is None:
        raise ValueError("seed or reference positive charge requires geometry owned ion profiles")
    coordinates = np.asarray(z_m, dtype=float)
    density = np.asarray(background.positive_charge_density_m3_at(coordinates), dtype=float)
    if density.shape != coordinates.shape:
        raise ValueError("seed or reference positive charge evaluation must match its coordinates")
    if np.any(~np.isfinite(density)) or np.any(density < 0.0):
        raise ValueError("seed or reference positive charge must be finite and nonnegative")

    return density

def _exact_background_positive_charge_reference_values(geometry: GeometryStageResult) -> tuple[float, float, float]:
    """Return configured positive charge number density at the magnetic midplane and both confined throat coordinates"""
    background = geometry.background_profiles
    if background is None:
        return 0.0, 0.0, 0.0
    z_edges = np.asarray(geometry.z_edges_m, dtype=float)
    if z_edges.ndim != 1 or z_edges.size < 2:
        raise ValueError("confined cell edges are required for exact background values")
    reference_z = np.asarray([background.magnetic_midplane_z_m, z_edges[0], z_edges[-1]], dtype=float)
    values = _seed_reference_positive_charge_at(geometry, reference_z)

    return float(values[0]), float(values[1]), float(values[2])

def _prescribed_eq42_density_profile(config: SourceModelRunConfig, geometry: GeometryStageResult, numerics: ModalFBISNumerics) -> Eq42DensityProfile:
    """
    Build the explicit Eq 42 density profile from the configured D plus T reference state
    
    Cell and node densities are normalized by physical flux tube cell volume inside build_eq42_density_profile
    """
    background = geometry.background_profiles
    if background is None:
        raise ValueError("prescribed_axial_density_profile requires geometry owned seed or reference ion profiles")
    density_cells = _seed_reference_positive_charge_at(geometry, np.asarray(geometry.z_centers_m, dtype=float))
    zeta_cells = np.asarray(geometry.zeta_centers, dtype=float)
    zeta_edges = np.asarray(geometry.zeta_edges, dtype=float)
    z_edges_m = np.asarray(geometry.z_edges_m, dtype=float)
    if zeta_edges.shape != z_edges_m.shape:
        raise ValueError("confined normalized and machine coordinate edges must match")
    zeta_nodes = np.unique(np.concatenate((zeta_edges, np.asarray([0.0]))))
    z_nodes_m = np.interp(zeta_nodes, zeta_edges, z_edges_m)
    midplane = np.isclose(zeta_nodes, 0.0, rtol=0.0, atol=1.0e-14)
    z_nodes_m[midplane] = float(background.magnetic_midplane_z_m)
    density_nodes = _seed_reference_positive_charge_at(geometry, z_nodes_m)

    return build_eq42_density_profile(
        model=EQ42_PRESCRIBED_AXIAL_DENSITY_PROFILE,
        source=("prescribed_prescribed_reference_D_plus_T_positive_charge_profile_" f"used_as_density_shape_diagnostic_{background.model}"),
        zeta_cells=zeta_cells,
        density_cells_m3=density_cells,
        cell_volumes_m3=geometry.cell_volumes_m3,
        zeta_nodes=zeta_nodes,
        density_nodes_m3=density_nodes,
        symmetry_tolerance=numerics.eq42_density_symmetry_tolerance,
    )

def _initial_eq42_density_profile(config: SourceModelRunConfig, geometry: GeometryStageResult, numerics: ModalFBISNumerics) -> Eq42DensityProfile:
    """
    Return the density profile used for the first Eq 42 basis construction
    
    The self consistent path begins from `n(z) / <n> = 1` while retaining a physical reference density amplitude for metadata
    """
    model = eq42_density_weighting_model(config.kinetic_electrostatic.modal_eq42_density_weighting_model)
    if model == EQ42_PRESCRIBED_AXIAL_DENSITY_PROFILE:
        return _prescribed_eq42_density_profile(config, geometry, numerics)
    reference_density = 1.0
    source_suffix = "arbitrary_unit_amplitude_without_prescribed_reference"
    if geometry.background_profiles is not None:
        density_cells = _seed_reference_positive_charge_at(geometry, np.asarray(geometry.z_centers_m, dtype=float))
        volumes = np.asarray(geometry.cell_volumes_m3, dtype=float)
        if volumes.shape != density_cells.shape or np.any(~np.isfinite(volumes)) or np.any(volumes <= 0.0):
            raise ValueError("physical confined cell volumes must be finite and positive")
        maintained_average = float(np.sum(density_cells * volumes) / np.sum(volumes))
        if np.isfinite(maintained_average) and maintained_average > 0.0:
            reference_density = maintained_average
            source_suffix = "prescribed_reference_D_plus_T_confined_physical_volume_average_amplitude"

    return uniform_eq42_density_profile(
        zeta_cells=geometry.zeta_centers,
        cell_volumes_m3=geometry.cell_volumes_m3,
        reference_density_m3=reference_density,
        model=model,
        symmetry_tolerance=numerics.eq42_density_symmetry_tolerance,
        source=("Egedal_2022_section_4_5_first_iteration_n_over_nbar_equals_one_" f"{source_suffix}" if model == EQ42_SELF_CONSISTENT_QUASINEUTRAL_DENSITY else f"explicit_uniform_density_reference_{source_suffix}"),
    )

def _basis_numerics_without_refinement_assessment(numerics: ModalFBISNumerics) -> ModalFBISNumerics:
    """Return modal numerics with basis refinement assessment disabled during the Eq 42 fixed point iteration"""
    return replace(numerics, basis_convergence_levels=1)

def _modal_basis_cache_key(*, geometry: GeometryStageResult, numerics: ModalFBISNumerics, profile: Eq42DensityProfile, assess_grid_convergence: bool, feedback_model: str) -> tuple[object, ...]:
    """Return a run local cache key containing every magnetic grid, numerical, density profile, and feedback input that defines the basis"""
    zeta_faces = np.asarray(geometry.zeta_edges, dtype=float)
    field_midpoints = np.asarray(geometry.B_tilde_midpoints, dtype=float)
  
    return (id(geometry.B_tilde_function), geometry.magnetic_field_model, float(geometry.mirror_ratio), zeta_faces.shape, zeta_faces.tobytes(), field_midpoints.shape, field_midpoints.tobytes(), numerics, profile.cache_key, bool(assess_grid_convergence), str(feedback_model))

def _basis_with_current_profile_metadata(basis, profile: Eq42DensityProfile):
    """Attach the current Eq 42 density profile identity and normalized density ratio to an existing basis result"""
    return replace(basis, eq42_density_weighting_model=profile.model, eq42_density_weighting_coupled=True, eq42_nonuniform_density_weighting_active=profile.nonuniform_weighting_active, eq42_density_profile_source=profile.source, eq42_volume_average_density_m3=profile.volume_average_density_m3, eq42_density_ratio_zeta=np.asarray(profile.zeta_positive_nodes, dtype=float), eq42_density_ratio_values=np.asarray(profile.density_ratio_positive_nodes, dtype=float), eq42_density_profile_identity=profile.profile_identity, eq42_density_profile_symmetry_relative_error=profile.symmetry_relative_error, eq42_density_profile_symmetry_tolerance=profile.symmetry_tolerance, eq42_density_profile_symmetric=profile.symmetric_within_tolerance)

def _cached_bundle(entry: _ModalBasisCacheEntry, *, profile: Eq42DensityProfile, basis_cache: ModalBasisCache, published_reference_lambda1: float | None) -> _ReusableModalBasisBundle:
    """Return a reusable basis bundle with the current Eq 42 profile metadata attached to a cached basis and convergence history"""
    assessment = entry.convergence_assessment
    if assessment is None:
        basis = _basis_with_current_profile_metadata(entry.basis, profile)
    else:
        history = tuple(_basis_with_current_profile_metadata(item, profile) for item in assessment.basis_history)
        assessment = replace(assessment, basis_history=history)
        basis = assessment.basis_history[0]
    reference_lambda1 = entry.published_reference_lambda1 if published_reference_lambda1 is None else float(published_reference_lambda1)
  
    return _ReusableModalBasisBundle(basis=basis, convergence_assessment=assessment, published_reference_lambda1=reference_lambda1, eq42_density_profile=profile, basis_cache=basis_cache)

def _build_modal_basis_bundle_for_eq42_profile(*, config: SourceModelRunConfig, geometry: GeometryStageResult, profile: Eq42DensityProfile, assess_grid_convergence: bool, published_reference_lambda1: float | None = None, basis_cache: ModalBasisCache | None = None, bypass_cache: bool = False) -> _ReusableModalBasisBundle:
    """
    Build or reuse one density weighted modal basis bundle for an Eq 42 profile
    
    Final qualification can assess basis grid convergence while ordinary Eq 42 iterations use the configured base grid
    """
    numerics = _modal_numerics(config)
    active_cache = ModalBasisCache() if basis_cache is None else basis_cache
    feedback_model = str(config.kinetic_electrostatic.electrostatic_feedback_model).strip().lower()
    cache_key = _modal_basis_cache_key(geometry=geometry, numerics=numerics, profile=profile, assess_grid_convergence=assess_grid_convergence, feedback_model=feedback_model)
    if bypass_cache:
        active_cache.record_bypass()
    else:
        cached = active_cache.lookup(cache_key)
        if cached is not None:
            return _cached_bundle(cached, profile=profile, basis_cache=active_cache, published_reference_lambda1=published_reference_lambda1)
    # Reserve the independent grid refinement assessment for the selected final profile
    if assess_grid_convergence and numerics.basis_convergence_levels >= 3:
        assessment = assess_modal_fbis_basis_convergence(mirror_ratio=geometry.mirror_ratio, B_tilde_function=geometry.B_tilde_function, zeta_faces=geometry.zeta_edges, B_tilde_midpoints=geometry.B_tilde_midpoints, numerics=numerics, eq42_density_profile=profile, eigenvalue_relative_tolerance=numerics.basis_eigenvalue_relative_tolerance, eigenfunction_minimum_overlap=numerics.basis_eigenfunction_overlap_tolerance)
        basis = assessment.basis_history[0]
    else:
        assessment = None
        basis = build_modal_fbis_basis(mirror_ratio=geometry.mirror_ratio, B_tilde_function=geometry.B_tilde_function, zeta_faces=geometry.zeta_edges, B_tilde_midpoints=geometry.B_tilde_midpoints, numerics=_basis_numerics_without_refinement_assessment(numerics), eq42_density_profile=profile)
    reference_lambda1 = published_reference_lambda1
    # The feedback option needs a same shape reference basis evaluated at mirror ratio 1.5
    if (feedback_model not in {"magnetic_only_fbis"} and reference_lambda1 is None):
        reference_profile = rescaled_magnetic_profile_for_reference_ratio(geometry.B_tilde_function, geometry.mirror_ratio)
        active_midpoints = np.asarray(geometry.B_tilde_midpoints, dtype=float)
        reference_midpoints = 1.0 + 0.5 * (active_midpoints - 1.0) / (float(geometry.mirror_ratio) - 1.0)
        reference_basis = build_modal_fbis_basis(mirror_ratio=1.5, B_tilde_function=reference_profile, zeta_faces=geometry.zeta_edges, B_tilde_midpoints=reference_midpoints, numerics=numerics)
        reference_lambda1 = float(reference_basis.physical_basis.eigenvalues[0])
    active_cache.store(cache_key, _ModalBasisCacheEntry(basis=basis, convergence_assessment=assessment, published_reference_lambda1=reference_lambda1))
   
    return _ReusableModalBasisBundle(basis=basis, convergence_assessment=assessment, published_reference_lambda1=reference_lambda1, eq42_density_profile=profile, basis_cache=active_cache)

def build_reusable_modal_basis(config: SourceModelRunConfig, geometry: GeometryStageResult) -> _ReusableModalBasisBundle:
    """Build the initial Eq 42 modal basis and run local basis cache for one source model evaluation"""
    numerics = _modal_numerics(config)
    profile = _initial_eq42_density_profile(config, geometry, numerics)
    model = eq42_density_weighting_model(config.kinetic_electrostatic.modal_eq42_density_weighting_model)
    
    return _build_modal_basis_bundle_for_eq42_profile(config=config, geometry=geometry, profile=profile, assess_grid_convergence=(model != EQ42_SELF_CONSISTENT_QUASINEUTRAL_DENSITY), basis_cache=ModalBasisCache())

def _collision_state_for_density(*, fast_ion_species: IonSpecies, electron_density_m3: float, ion_densities_m3: np.ndarray, ion_charge_numbers: np.ndarray, ion_masses_kg: np.ndarray, electron_temperature_J: float, ion_temperature_J: float, fast_ion_density_m3: float = 0.0, collision_scope: str = "maxwellian_reference_with_fast_ion") -> FBISCollisionParameterState:
    """
    Build the FBIS collision parameter state for supplied electron and ion densities
    
    The representative background ion mass is the unique active ion mass when only one mass is present and otherwise defaults to the deuteron reference used by the collision helper
    """
    densities = np.asarray(ion_densities_m3, dtype=float)
    masses = np.asarray(ion_masses_kg, dtype=float)
    active_masses = np.unique(masses[densities > 0.0])
    representative_mass = (float(active_masses[0]) if active_masses.size == 1 else DEUTERON.mass_kg)
  
    return build_fbis_collision_parameter_state(
        electron_density_m3=float(electron_density_m3),
        electron_temperature_J=float(electron_temperature_J),
        ion_densities_m3=densities,
        ion_charge_numbers=np.asarray(ion_charge_numbers, dtype=float),
        ion_masses_kg=masses,
        ion_temperature_J=float(ion_temperature_J),
        fast_ion_species=fast_ion_species,
        background_ion_mass_kg=representative_mass,
        fast_ion_density_m3=float(fast_ion_density_m3),
        collision_scope=str(collision_scope),
    )

def _background_positive_charge_on_confined_grid(geometry: GeometryStageResult,) -> tuple[np.ndarray | None, float, float, str]:
    """Return the configured D plus T positive charge profile on confined cells plus exact left and right throat values"""
    background = geometry.background_profiles
    if background is None:
        return None, 0.0, 0.0, "not_available_manual_or_legacy_geometry"
    z = np.asarray(geometry.z_centers_m, dtype=float)
    profile = _seed_reference_positive_charge_at(geometry, z)
    _, left_throat, right_throat = _exact_background_positive_charge_reference_values(geometry)

    return profile, left_throat, right_throat, "geometry_owned_prescribed_reference_D_plus_T_profile_exactly_evaluated_on_confined_grid"
