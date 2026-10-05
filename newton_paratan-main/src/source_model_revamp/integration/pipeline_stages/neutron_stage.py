"""
Deterministic neutron spectrum and correlated event stage

The deterministic matrix uses the full device fusion populations and exact isotropic center of momentum energy marginal while the event path supplies evaluated angular correlations
"""
from __future__ import annotations
from time import perf_counter
import numpy as np
from source_model_revamp.constants import EV_TO_J
from source_model_revamp.fusion.reactions import MEV_TO_J
from source_model_revamp.fusion.populations import fusion_population_display_label
from source_model_revamp.integration.full_device_populations import build_full_device_population_grid
from source_model_revamp.integration.validation import NEUTRON_RAW_TO_FUSION_RELATIVE_TOLERANCE
from source_model_revamp.integration.pipeline_types import FusionComponent, FusionStageResult, GeometryStageResult, KineticStageResult, NeutronStageResult
from source_model_revamp.integration.pipeline_stages.neutron_event_stage import build_correlated_neutron_event_stage
from source_model_revamp.integration.config.root import SourceModelRunConfig
from source_model_revamp.neutrons.lambda_spectrum import AxialNeutronSpectrumMatrix, local_gyrotropic_pair_neutron_spectrum_matrix
from source_model_revamp.integration.pipeline_stages.neutron_event_stage import _component_pair_gyroangles
from source_model_revamp.runtime_profile import finalize_runtime_profile, new_runtime_profile, record_runtime

def _neutron_energy_edges(config: SourceModelRunConfig) -> np.ndarray:
    """Return configured neutron energy bin edges converted from MeV to J"""
    return np.linspace(config.neutrons.energy_min_MeV * MEV_TO_J, config.neutrons.energy_max_MeV * MEV_TO_J, config.neutrons.energy_bins + 1)

def _component_spectrum_matrix(*, component: FusionComponent, config: SourceModelRunConfig, coordinate: np.ndarray, volumes: np.ndarray, populations_by_id, energy_edges: np.ndarray, max_pair_states_per_batch: int, exact_isotropic_energy_marginal: bool) -> AxialNeutronSpectrumMatrix:
    """
    Build one axial deterministic neutron spectrum matrix from a fusion component
    
    The pair uses the same full device reactant populations and relative gyro angle count as the fusion kernel and does not apply posterior target rate normalization
    """
    ncfg = config.neutrons
    num_gyroangle_points = _component_pair_gyroangles(component)
   
    try:
        population_a = populations_by_id[component.reactant_a_population_id]
        population_b = populations_by_id[component.reactant_b_population_id]
    except KeyError as exc:
        raise ValueError(f"fusion component {component.label!r} references an unavailable neutron reactant population") from exc
   
    identical_population = component.reactant_a_population_id == component.reactant_b_population_id
  
    return local_gyrotropic_pair_neutron_spectrum_matrix(coordinate=coordinate, speed_grid_a=population_a.speed_grid, pitch_grid_a=population_a.pitch_grid, local_distribution_a_z_v_pitch=population_a.local_distribution_z_v_pitch, speed_grid_b=population_b.speed_grid, pitch_grid_b=population_b.pitch_grid, local_distribution_b_z_v_pitch=population_b.local_distribution_z_v_pitch, reaction=component.reaction, energy_edges_J=energy_edges, cell_volumes_m3=volumes, target_rate_density_m3_s=None, identical_population=identical_population, num_gyroangle_points=num_gyroangle_points, num_emission_directions=ncfg.num_emission_directions, allow_out_of_range=ncfg.allow_out_of_range, max_pair_states_per_batch=max_pair_states_per_batch, exact_isotropic_energy_marginal=exact_isotropic_energy_marginal)

def _component_rate_diagnostics(component: FusionComponent, matrix: AxialNeutronSpectrumMatrix, volumes: np.ndarray) -> dict[str, object]:
    """Compare one raw deterministic spectrum rate with its fusion pair rate and retain represented and out of range accounting"""
    raw_density = np.asarray(matrix.raw_total_rate_density_m3_s, dtype=float)
    out_of_range_density = np.asarray(matrix.out_of_range_rate_density_m3_s, dtype=float,)
    represented_density = np.asarray(matrix.axial_rate_density_m3_s, dtype=float,)
    fusion_density = np.asarray(component.profile.neutron_source_density_m3_s, dtype=float,)
    expected_shape = volumes.shape
   
    for name, values in (("raw_density", raw_density), ("out_of_range_density", out_of_range_density), ("represented_density", represented_density), ("fusion_density", fusion_density),):
        if (values.shape != expected_shape or np.any(~np.isfinite(values)) or np.any(values < 0.0)):
            raise ValueError(f"{component.label} {name} must be finite and nonnegative")
  
    tolerance = (64.0 * np.finfo(float).eps * np.maximum(raw_density, np.finfo(float).tiny))
   
    if np.any(out_of_range_density > raw_density + tolerance):
        raise ValueError(f"{component.label} out of range rate exceeds the raw rate")
   
    raw_axial = raw_density * volumes
    out_of_range_axial = out_of_range_density * volumes
    represented_axial = represented_density * volumes
    fusion_axial = fusion_density * volumes
    raw_total = float(np.sum(raw_axial))
    out_of_range_total = float(np.sum(out_of_range_axial))
    represented_total = float(np.sum(represented_axial))
    fusion_total = float(np.sum(fusion_axial))
    fusion_denominator = max(fusion_total, np.finfo(float).tiny)
    raw_denominator = max(raw_total, np.finfo(float).tiny)
    raw_global_error = (raw_total - fusion_total) / fusion_denominator
    raw_axial_error = float(np.sum(np.abs(raw_axial - fusion_axial)) / fusion_denominator)
    represented_global_error = (represented_total - fusion_total) / fusion_denominator
    represented_axial_error = float(np.sum(np.abs(represented_axial - fusion_axial)) / fusion_denominator)
    represented_accounting_error = float(np.sum(np.abs(represented_axial - (raw_axial - out_of_range_axial))) / raw_denominator)
    kernel_model = (component.pair_kernel.kernel_model if component.pair_kernel is not None else "missing_canonical_pair_kernel")
  
    return {
        "label": component.label,
        "component_kind": component.kind,
        "reaction": component.reaction.key,
        "reactant_a_population_id": component.reactant_a_population_id,
        "reactant_b_population_id": component.reactant_b_population_id,
        "reactant_a_population_label": fusion_population_display_label(component.reactant_a_population_id),
        "reactant_b_population_label": fusion_population_display_label(component.reactant_b_population_id),
        "fusion_rate_kernel_model": kernel_model,
        "raw_neutron_rate_kernel_model": kernel_model,
        "raw_to_fusion_same_rate_kernel": True,
        "fusion_target_total_rate_s": fusion_total,
        "fusion_target_axial_bin_rates_s": fusion_axial,
        "raw_total_rate_s": raw_total,
        "raw_axial_bin_rates_s": raw_axial,
        "represented_total_rate_s": represented_total,
        "represented_axial_bin_rates_s": represented_axial,
        "out_of_range_rate_s": out_of_range_total,
        "out_of_range_axial_bin_rates_s": out_of_range_axial,
        "raw_to_fusion_global_relative_error": raw_global_error,
        "raw_to_fusion_axial_l1_relative_error": raw_axial_error,
        "represented_to_fusion_global_relative_error": represented_global_error,
        "represented_to_fusion_axial_l1_relative_error": represented_axial_error,
        "represented_to_raw_in_range_axial_l1_relative_error": (represented_accounting_error),
        "raw_quadrature_diagnostics_available": True,
        "raw_quadrature_diagnostics_unavailable_reason": None,
        "posterior_rate_normalization_applied": False,
        "num_gyroangle_points": _component_pair_gyroangles(component),
    }

def build_neutron_stage(config: SourceModelRunConfig, geometry: GeometryStageResult, kinetic: KineticStageResult, fusion: FusionStageResult) -> NeutronStageResult:
    """
    Build deterministic D D and D T neutron spectra and the optional correlated event bank
    
    Scalar conservation is checked against the fusion components
    The deterministic energy marginal remains independent of the evaluated angular law, which is applied only to correlated events
    """
    neutron_runtime_started = perf_counter()
    runtime_profile = new_runtime_profile()
    ncfg = config.neutrons
    max_pair_states_per_batch = int(getattr(ncfg, "max_pair_states_per_batch", 131072))
    exact_isotropic_energy_marginal = True
    energy_edges = _neutron_energy_edges(config)
    population_grid_started = perf_counter()
    population_grid = build_full_device_population_grid(geometry, kinetic)
  
    record_runtime(runtime_profile, "neutron_population_grid", perf_counter() - population_grid_started)
  
    coordinate = np.asarray(population_grid.z_centers_m, dtype=float)
    volumes = np.asarray(population_grid.cell_volumes_m3, dtype=float)
    populations_by_id = fusion.populations_by_id
    primary_fast = populations_by_id.get("fast_D") or populations_by_id.get("fast_T")
    local_fast = None if primary_fast is None else primary_fast.local_distribution_z_v_pitch
    use_local_fast = primary_fast is not None
    local_fast_speed_grid = kinetic.speed_grid if primary_fast is None else primary_fast.speed_grid

    matrices: list[tuple[str, AxialNeutronSpectrumMatrix]] = []
    neutron_components: list[FusionComponent] = []
   
    for component in fusion.components:
        if component.reaction.neutron_yield_per_reaction <= 0.0:
            continue
        target_density = np.asarray(component.profile.neutron_source_density_m3_s, dtype=float,)
       
        if np.all(target_density == 0.0):
            continue 
        component_spectrum_started = perf_counter()
        matrix = _component_spectrum_matrix(component=component, config=config, coordinate=coordinate, volumes=volumes, populations_by_id=populations_by_id, energy_edges=energy_edges, max_pair_states_per_batch=max_pair_states_per_batch, exact_isotropic_energy_marginal=exact_isotropic_energy_marginal)
        component_spectrum_runtime_s = perf_counter() - component_spectrum_started
      
        record_runtime(runtime_profile, "neutron_spectrum_matrices", component_spectrum_runtime_s)
        record_runtime(runtime_profile, f"neutron_spectrum_matrix_{component.label}", component_spectrum_runtime_s, nested=True)
       
        if matrix.posterior_rate_normalization_applied:
            raise RuntimeError("source neutron spectra must not be normalized to a fusion target rate")
        matrices.append((component.label, matrix))
        neutron_components.append(component)
    if matrices:
        source_rate_matrix_s = np.sum([matrix.source_rate_matrix_s for _, matrix in matrices], axis=0,)
    else:
        source_rate_matrix_s = np.zeros((coordinate.size, ncfg.energy_bins), dtype=float,)
    correlated_event_started = perf_counter()
   
    # Evaluated CM angular data enter only through the correlated event bank
    event_stage = build_correlated_neutron_event_stage(config=config, neutron_components=tuple(neutron_components), populations_by_id=populations_by_id, population_grid=population_grid, energy_edges_J=energy_edges)
    record_runtime(runtime_profile, "correlated_event_generation", perf_counter() - correlated_event_started)
    neutron_postprocessing_started = perf_counter()
    correlated_event_bank = event_stage.event_bank
    correlated_event_source_rate_matrix_s = event_stage.source_rate_matrix_s
    correlated_event_out_of_range_rate_s = event_stage.out_of_range_rate_s
    component_diagnostics = [ _component_rate_diagnostics(component, matrix, volumes) for component, (_, matrix) in zip(neutron_components, matrices, strict=True)]
  
    if component_diagnostics:
        raw_axial_bin_rates_s = np.sum([item["raw_axial_bin_rates_s"] for item in component_diagnostics], axis=0,)
        out_of_range_axial_bin_rates_s = np.sum([ item["out_of_range_axial_bin_rates_s"] for item in component_diagnostics], axis=0,)
        fusion_profile_axial_bin_rates_s = np.sum([ item["fusion_target_axial_bin_rates_s"] for item in component_diagnostics], axis=0,)
        component_max_global_error = max(abs(float(item["raw_to_fusion_global_relative_error"])) for item in component_diagnostics)
        component_max_axial_error = max(float(item["raw_to_fusion_axial_l1_relative_error"]) for item in component_diagnostics)
        pair_gyroangle_counts = sorted({int(item["num_gyroangle_points"]) for item in component_diagnostics})
        pair_gyroangle_metadata: int | list[int] | None = (pair_gyroangle_counts[0] if len(pair_gyroangle_counts) == 1 else pair_gyroangle_counts)
    else:
        raw_axial_bin_rates_s = np.zeros(coordinate.size, dtype=float)
        out_of_range_axial_bin_rates_s = np.zeros(coordinate.size, dtype=float,)
        fusion_profile_axial_bin_rates_s = np.zeros(coordinate.size, dtype=float,)
        component_max_global_error = 0.0
        component_max_axial_error = 0.0
        pair_gyroangle_metadata = None
    
    represented_axial_bin_rates_s = np.sum(source_rate_matrix_s, axis=1)
    raw_total_rate_s = float(np.sum(raw_axial_bin_rates_s))
    out_of_range_rate_s = float(np.sum(out_of_range_axial_bin_rates_s))
    fusion_profile_total_rate_s = float(np.sum(fusion_profile_axial_bin_rates_s))
    total = float(np.sum(source_rate_matrix_s))
    denominator = max(fusion_profile_total_rate_s, np.finfo(float).tiny)
    raw_to_fusion_relative_error = (raw_total_rate_s - fusion_profile_total_rate_s) / denominator
    raw_to_fusion_axial_l1_relative_error = float(np.sum(np.abs(raw_axial_bin_rates_s - fusion_profile_axial_bin_rates_s)) / denominator)
    spectrum_relative_error = (total - fusion_profile_total_rate_s) / denominator
    spectrum_axial_l1_relative_error = float(np.sum(np.abs(represented_axial_bin_rates_s - fusion_profile_axial_bin_rates_s)) / denominator)
    represented_accounting_error = float(np.sum(np.abs(represented_axial_bin_rates_s  - (raw_axial_bin_rates_s - out_of_range_axial_bin_rates_s))) / max(raw_total_rate_s, np.finfo(float).tiny))
    same_kernel_max_relative_error = max(abs(raw_to_fusion_relative_error), raw_to_fusion_axial_l1_relative_error, component_max_global_error, component_max_axial_error)
    conservative_rate_error = max(same_kernel_max_relative_error, abs(spectrum_relative_error), spectrum_axial_l1_relative_error, represented_accounting_error)
    out_of_range_fraction = out_of_range_rate_s / max(raw_total_rate_s, np.finfo(float).tiny)
    axial_bin_rates_s = np.sum(source_rate_matrix_s, axis=1)
    energy_bin_rates_s = np.sum(source_rate_matrix_s, axis=0)
    dt_source_rate_matrix_s = np.sum([ matrix.source_rate_matrix_s for component, (_, matrix) in zip(neutron_components, matrices, strict=True,) if component.reaction.key == "dt_n"], axis=0,) if any(component.reaction.key == "dt_n" for component in neutron_components) else np.zeros_like(source_rate_matrix_s)
    dt_axial_bin_rates_s = np.sum(dt_source_rate_matrix_s, axis=1)
    dt_total_rate_s = float(np.sum(dt_source_rate_matrix_s))
    fusion_profile_dt_rate_s = float(sum(float(np.sum(component.profile.neutron_source_density_m3_s * volumes) ) for component in neutron_components if component.reaction.key == "dt_n"))
    dt_relative_error = (dt_total_rate_s - fusion_profile_dt_rate_s) / max(fusion_profile_dt_rate_s, np.finfo(float).tiny)
    evaluated_angular_active = correlated_event_bank is not None
    event_rate_s = (None if correlated_event_bank is None else float(correlated_event_bank.physical_total_rate_s))
    event_rate_relative_error = (None if event_rate_s is None else (event_rate_s - fusion_profile_total_rate_s) / denominator)
    event_matrix_rate_s = (None if correlated_event_source_rate_matrix_s is None else float(np.sum(correlated_event_source_rate_matrix_s)))
    event_matrix_relative_error = (None if event_matrix_rate_s is None else (event_matrix_rate_s - fusion_profile_total_rate_s) / denominator)
    event_metadata = (None if correlated_event_bank is None else dict(correlated_event_bank.metadata))
    event_reaction_validation = (None if event_metadata is None else event_metadata.get("endf_b_viii1_reaction_validation") )
    direct_angular_validation_fields = ("reaction_key", "reaction_identity_valid", "reference_frame", "reference_frame_valid", "angular_law_representation_valid", "angular_law_validation_passed", "angular_normalization_measured", "maximum_angular_normalization_error", "angular_normalization_tolerance", "minimum_angular_probability_density_per_mu", "minimum_angular_probability_tolerance")
    direct_angular_validation_measured = bool(isinstance(event_reaction_validation, dict) and event_reaction_validation and all(isinstance(validation, dict) and all(field in validation for field in direct_angular_validation_fields) for validation in event_reaction_validation.values()))
    reaction_identity_passed = bool(direct_angular_validation_measured and all(validation.get("reaction_key") == reaction_key and validation.get("reaction_identity_valid") is True for reaction_key, validation in event_reaction_validation.items()))
    reference_frame_passed = bool(direct_angular_validation_measured and all(validation.get("reference_frame") == "center_of_mass" and validation.get("reference_frame_valid") is True for validation in event_reaction_validation.values()))
    angular_representation_passed = bool(direct_angular_validation_measured and all(validation.get("angular_law_representation_valid") is True for validation in event_reaction_validation.values()))
    angular_normalization_passed = bool(isinstance(event_reaction_validation, dict) and event_reaction_validation and all(isinstance(validation, dict) and validation.get("angular_law_validation_passed") is True and validation.get("angular_normalization_measured") is True and float(validation["maximum_angular_normalization_error"]) <= float(validation["angular_normalization_tolerance"]) for validation in event_reaction_validation.values()))
    angular_nonnegativity_passed = bool(direct_angular_validation_measured and all("minimum_angular_probability_density_per_mu" in validation and "minimum_angular_probability_tolerance" in validation and float(validation["minimum_angular_probability_density_per_mu"]) >= float(validation["minimum_angular_probability_tolerance"]) for validation in event_reaction_validation.values()))
    angular_coverage_passed = bool(event_metadata is not None and event_metadata.get("angular_coverage_measured") is True and event_metadata.get("angular_coverage_passed") is True)
    event_probability_passed = bool(event_metadata is not None and event_metadata.get("event_probability_normalization_measured") is True and event_metadata.get("event_probability_normalization_passed") is True)
    event_rate_conservation_passed = bool(event_rate_relative_error is not None and abs(float(event_rate_relative_error)) <= 5.0e-13)
    relativistic_kinematics_passed = bool(event_metadata is not None and float(event_metadata["maximum_residual_mass_shell_relative_error"]) <= 1.0e-10)
    differential_fidelity_passed = (None if event_metadata is None else bool(reaction_identity_passed and reference_frame_passed and angular_representation_passed and angular_normalization_passed and angular_nonnegativity_passed and angular_coverage_passed and event_rate_conservation_passed and event_probability_passed and relativistic_kinematics_passed and int(event_metadata["positive_weight_event_count"]) > 0))
    metadata = {'neutron_spectrum_model': ncfg.model, 'neutron_angular_model': ncfg.angular_model, 'neutron_radial_source_profile_model': ncfg.radial_source_profile_model, 'neutron_radial_source_profile_role': 'prescribed_position_only_source_shaping_after_one_dimensional_fusion_rate', 'neutron_radial_source_profile_is_radially_self_consistent_FBIS': False, 'neutron_radial_source_profile_reference': 'Kotelnikov_et_al_Journal_of_Plasma_Physics_2025_Eq3_11_k1', 'neutron_radial_source_profile_reference_role': 'published_pressure_profile_used_as_a_prescribed_reduced_neutron_spatial_source_closure_not_a_published_neutron_emissivity_law', 'neutron_radial_coordinate_definition': 'rho_equals_r_over_nominal_a_of_z_and_psi_tilde_equals_rho_squared', 'neutron_radial_source_shape_definition': 'G_of_rho_equals_2_times_one_minus_rho_squared', 'neutron_radial_geometry_conditioned_source_shape_definition': 'G_of_rho_divided_by_F_of_rho_max_inside_source_connected_support', 'neutron_source_connected_radial_rho_limit_edges': np.asarray(population_grid.source_connected_radial_rho_max_edges, dtype=float), 'neutron_radial_geometry_truncation_active': bool(np.any(np.asarray(population_grid.source_connected_radial_rho_max_edges, dtype=float) < 1.0 - 1e-12)), 'neutron_radial_sample_support_check_passed': None if event_metadata is None else event_metadata.get('radial_sample_support_check_passed'), 'neutron_reaction_cross_section_model': 'bosch_hale_table_iv_total_cross_section', 'neutron_rate_kernel_model': 'bosch_hale_table_iv_total_cross_section_quadrature', 'neutron_scalar_and_spectrum_share_rate_kernel': True, 'neutron_posterior_rate_normalization_applied': False, 'neutron_posterior_rate_normalization_absence_measured': True, 'neutron_posterior_rate_normalization_absence_passed': True, 'neutron_energy_marginal_model': 'deterministic_exact_isotropic_cm_rate_benchmark_plus_correlated_endf_event_bank' if evaluated_angular_active else 'exact_uniform_lab_energy_interval_from_isotropic_cm_emission', 'neutron_exact_isotropic_cm_energy_marginal_active': exact_isotropic_energy_marginal, 'neutron_cm_emission_direction_sampling_model': 'endf_b_viii1_mf6_law2_lct2_conditional_mu_plus_uniform_azimuth' if evaluated_angular_active else 'not_used_for_exact_energy_only_isotropic_cm_marginal', 'neutron_differential_cross_section_model': 'bosch_hale_table_iv_total_sigma_times_normalized_endf_b_viii1_conditional_cm_angular_law' if evaluated_angular_active else None, 'neutron_relativistic_two_body_kinematics_active': True, 'neutron_differential_cross_section_fidelity_required': True, 'neutron_differential_cross_section_fidelity_available': bool(evaluated_angular_active), 'neutron_differential_cross_section_fidelity_assessed': bool(evaluated_angular_active), 'neutron_differential_cross_section_fidelity_passed': differential_fidelity_passed, 'neutron_differential_cross_section_unavailability_reason': None if evaluated_angular_active else 'isotropic_cm_reference_does_not_use_evaluated_angular_data', 'neutron_endf_reaction_validation': event_reaction_validation, 'neutron_endf_reaction_identity_measured': direct_angular_validation_measured, 'neutron_endf_reaction_identity_passed': reaction_identity_passed if evaluated_angular_active else None, 'neutron_endf_reference_frame_measured': direct_angular_validation_measured, 'neutron_endf_reference_frame_passed': reference_frame_passed if evaluated_angular_active else None, 'neutron_endf_angular_law_representation_measured': direct_angular_validation_measured, 'neutron_endf_angular_law_representation_passed': angular_representation_passed if evaluated_angular_active else None, 'neutron_endf_angular_normalization_measured': bool(evaluated_angular_active), 'neutron_endf_angular_normalization_passed': angular_normalization_passed if evaluated_angular_active else None, 'neutron_endf_angular_nonnegativity_measured': direct_angular_validation_measured, 'neutron_endf_angular_nonnegativity_passed': angular_nonnegativity_passed if evaluated_angular_active else None, 'neutron_endf_angular_coverage_measured': bool(evaluated_angular_active), 'neutron_endf_angular_coverage_passed': angular_coverage_passed if evaluated_angular_active else None, 'neutron_correlated_event_bank_available': bool(evaluated_angular_active), 'neutron_correlated_event_bank_metadata': event_metadata, 'neutron_correlated_event_count': None if correlated_event_bank is None else correlated_event_bank.event_count, 'neutron_correlated_event_seed': None if event_metadata is None else event_metadata['random_seed'], 'neutron_correlated_event_physical_rate_s': event_rate_s, 'neutron_correlated_event_to_fusion_rate_relative_error': event_rate_relative_error, 'neutron_correlated_event_rate_conservation_measured': bool(evaluated_angular_active), 'neutron_correlated_event_rate_conservation_passed': event_rate_conservation_passed if evaluated_angular_active else None, 'neutron_correlated_event_histogram_rate_s': event_matrix_rate_s, 'neutron_correlated_event_histogram_to_fusion_relative_error': event_matrix_relative_error, 'neutron_correlated_event_out_of_range_rate_s': correlated_event_out_of_range_rate_s, 'neutron_correlated_event_probability_weight_sum': None if correlated_event_bank is None else float(np.sum(correlated_event_bank.normalized_weights)), 'neutron_correlated_event_probability_normalization_measured': bool(evaluated_angular_active), 'neutron_correlated_event_probability_normalization_passed': event_probability_passed if evaluated_angular_active else None, 'neutron_correlated_event_positive_weight_count': None if event_metadata is None else event_metadata['positive_weight_event_count'], 'neutron_correlated_event_zero_weight_count': None if event_metadata is None else event_metadata['zero_weight_event_count'], 'neutron_correlated_event_maximum_mass_shell_relative_error': None if event_metadata is None else event_metadata['maximum_residual_mass_shell_relative_error'], 'neutron_relativistic_event_kinematics_measured': bool(evaluated_angular_active), 'neutron_relativistic_event_kinematics_passed': relativistic_kinematics_passed if evaluated_angular_active else None, 'neutron_correlated_event_angular_coverage_excluded_fraction': None if event_metadata is None else event_metadata['angular_coverage_excluded_importance_fraction'], 'neutron_lab_energy_angle_correlation_preserved_in_event_bank': bool(evaluated_angular_active), 'neutron_dd_identical_population_projectile_label_model': 'random_equal_probability_exchange_symmetrization' if evaluated_angular_active else None, 'neutron_energy_min_MeV': ncfg.energy_min_MeV, 'neutron_energy_max_MeV': ncfg.energy_max_MeV, 'neutron_energy_bins': ncfg.energy_bins, 'neutron_requested_num_gyroangle_points': ncfg.num_gyroangle_points, 'neutron_num_gyroangle_points': pair_gyroangle_metadata, 'neutron_pair_kernel_num_gyroangle_points': pair_gyroangle_metadata, 'neutron_num_emission_directions': ncfg.num_emission_directions, 'neutron_num_emission_directions_used_for_energy_marginal': False, 'neutron_max_pair_states_per_batch': max_pair_states_per_batch, 'neutron_allow_out_of_range': bool(ncfg.allow_out_of_range), 'neutron_component_count': len(matrices), 'neutron_component_population_pairs': {component.label: [component.reactant_a_population_id, component.reactant_b_population_id] for component in neutron_components}, 'neutron_population_display_label_by_id': dict(fusion.metadata.get('fusion_population_display_label_by_id', {})), 'neutron_nbi_supported_kinetic_population_only': bool(fusion.metadata.get('fusion_nbi_supported_kinetic_population_only', False)), 'neutron_startup_seed_authority_fraction': fusion.metadata.get('fusion_startup_seed_authority_fraction'), 'neutron_DT_canonical_deuteron_first_ordering': True, 'neutron_DT_equivalent_incident_energy_model': 'pair_invariant_to_equivalent_deuteron_on_stationary_triton_lab_energy', 'neutron_component_raw_diagnostics': component_diagnostics, 'neutron_raw_same_kernel_to_fusion_component_max_relative_error': same_kernel_max_relative_error, 'neutron_thermal_cross_approximant_comparison_present': bool(fusion.metadata.get('fusion_thermal_table_iv_to_table_vii_validation_present', False)), 'neutron_fast_distribution_source': fusion.metadata.get('fusion_fast_distribution_source', 'modal_phi_z_local_z_v_pitch' if use_local_fast else 'global_v_lambda_remapped_by_B_tilde'), 'neutron_local_phi_z_distribution_used': bool(use_local_fast), 'neutron_invariant_speed_grid_max_m_s': float(kinetic.speed_grid.faces_m_s[-1]), 'neutron_local_speed_grid_max_m_s': None if not use_local_fast else float(local_fast_speed_grid.faces_m_s[-1]), 'neutron_distinct_local_speed_grid_used': bool(use_local_fast and local_fast_speed_grid is not kinetic.speed_grid), 'neutron_domain_model': 'full_device_population_coexistence_grid', 'neutron_electron_density_required': False, 'neutron_electron_density_used': False, 'neutron_solved_electron_profile_scope': population_grid.electron_profile_scope, 'neutron_expander_electron_profile_available': bool(fusion.metadata.get('fusion_expander_electron_profile_available', False)), 'neutron_lost_fast_ion_population_included': bool(population_grid.includes_lost_fast_ions), 'neutron_lost_fast_ion_population_included_by_species': dict(population_grid.includes_lost_fast_ions_by_species), 'neutron_lost_fast_ion_population_status': 'included_Egedal_Eq63_central_and_source_connected_expander_population' if fusion.metadata.get('fusion_expander_population_included') and population_grid.includes_lost_fast_ions else 'included_explicitly_represented_population' if population_grid.includes_lost_fast_ions else 'unavailable_not_substituted', 'neutron_background_profile_source': fusion.metadata.get('fusion_background_profile_source', 'geometry_derived_central_cell_Maxwellian'), 'neutron_expander_population_included': bool(fusion.metadata.get('fusion_expander_population_included', False)), 'neutron_source_disconnected_expander_population_included': False, 'neutron_total_rate_from_spectrum_s': total, 'neutron_scalar_table_iv_fusion_rate_s': fusion_profile_total_rate_s, 'neutron_raw_deterministic_spectrum_rate_s': raw_total_rate_s, 'neutron_raw_deterministic_to_scalar_rate_difference_s': raw_total_rate_s - fusion_profile_total_rate_s, 'neutron_correlated_event_to_scalar_rate_difference_s': None if event_rate_s is None else event_rate_s - fusion_profile_total_rate_s, 'neutron_fusion_profile_total_rate_s': fusion_profile_total_rate_s, 'neutron_spectrum_to_fusion_rate_relative_error': spectrum_relative_error, 'neutron_spectrum_to_fusion_axial_l1_relative_error': spectrum_axial_l1_relative_error, 'neutron_raw_total_rate_s': raw_total_rate_s, 'neutron_raw_axial_bin_rates_s': raw_axial_bin_rates_s, 'neutron_raw_to_fusion_rate_relative_error': raw_to_fusion_relative_error, 'neutron_raw_to_fusion_axial_l1_relative_error': raw_to_fusion_axial_l1_relative_error, 'neutron_raw_to_fusion_component_max_global_relative_error': component_max_global_error, 'neutron_raw_to_fusion_component_max_axial_l1_relative_error': component_max_axial_error, 'neutron_raw_to_fusion_conservative_max_relative_error': same_kernel_max_relative_error, 'neutron_represented_to_raw_in_range_relative_error': represented_accounting_error, 'neutron_out_of_range_rate_s': out_of_range_rate_s, 'neutron_out_of_range_axial_bin_rates_s': out_of_range_axial_bin_rates_s, 'neutron_out_of_range_fraction_of_raw': out_of_range_fraction, 'neutron_raw_quadrature_diagnostics_available': bool(component_diagnostics), 'neutron_raw_quadrature_is_pre_target_normalization': True, 'neutron_dt_rate_from_spectrum_s': dt_total_rate_s, 'neutron_dt_rate_from_fusion_profile_s': fusion_profile_dt_rate_s, 'neutron_dt_spectrum_to_fusion_rate_relative_error': dt_relative_error, 'neutron_fusion_profile_axial_bin_rates_s': fusion_profile_axial_bin_rates_s, 'neutron_fusion_to_spectrum_conservation_tolerance': NEUTRON_RAW_TO_FUSION_RELATIVE_TOLERANCE, 'neutron_fusion_to_spectrum_conservation_max_relative_error': conservative_rate_error, 'neutron_fusion_to_spectrum_conservation_check_passed': bool(conservative_rate_error <= NEUTRON_RAW_TO_FUSION_RELATIVE_TOLERANCE), 'neutron_dt_source_rate_matrix_s': dt_source_rate_matrix_s, 'neutron_dt_axial_bin_rates_s': dt_axial_bin_rates_s, 'source_rate_matrix_s': source_rate_matrix_s, 'energy_edges_J': energy_edges, 'energy_edges_MeV': energy_edges / MEV_TO_J, 'energy_edges_eV': energy_edges / EV_TO_J, 'axial_bin_rates_s': axial_bin_rates_s, 'energy_bin_rates_s': energy_bin_rates_s}
  
    record_runtime(runtime_profile, "neutron_postprocessing", perf_counter() - neutron_postprocessing_started)
  
    metadata = {**metadata, "runtime_neutron_profile": finalize_runtime_profile(runtime_profile, total_s=perf_counter() - neutron_runtime_started)}

    return NeutronStageResult(
        component_matrices=tuple(matrices),
        source_rate_matrix_s=source_rate_matrix_s,
        energy_edges_J=energy_edges,
        total_neutron_rate_s=total,
        metadata=metadata,
        correlated_event_bank=correlated_event_bank,
        correlated_event_source_rate_matrix_s=(correlated_event_source_rate_matrix_s),
    )
