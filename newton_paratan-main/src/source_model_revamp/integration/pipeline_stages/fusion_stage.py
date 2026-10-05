"""
Full device fusion source assembly from solved fast ion populations

The stage forms active D D and D T population pairs on the common full device grid and evaluates the Bosch and Hale Table IV pair kernel
"""
from __future__ import annotations
from time import perf_counter
import numpy as np
from source_model_revamp.fbis.species import DEUTERON, TRITON
from source_model_revamp.fusion.pair_kernel import FusionPairKernelResult, bosch_hale_table_iv_pair_kernel_profiles
from source_model_revamp.fusion.populations import FAST_D, FAST_T, FusionReactantPopulation, active_population_pair_definitions, fusion_population_display_label
from source_model_revamp.fusion.reactions import FusionReaction
from source_model_revamp.fusion.source_profiles import combine_fusion_source_profiles, fusion_source_profile_from_rate_density
from source_model_revamp.integration.full_device_populations import FullDevicePopulationGrid, build_full_device_population_grid
from source_model_revamp.integration.pipeline_types import ExpanderStageResult, FusionComponent, FusionStageResult, GeometryStageResult, KineticStageResult
from source_model_revamp.integration.config.root import SourceModelRunConfig
from source_model_revamp.integration.plasma_closure import uses_startup_seed_only
from source_model_revamp.runtime_profile import finalize_runtime_profile, new_runtime_profile, record_runtime

def _component_from_pair_kernel(*, label: str, kind: str, reaction: FusionReaction, reactant_a_population_id: str, reactant_b_population_id: str, pair_kernel: FusionPairKernelResult, coordinate: np.ndarray, cell_volumes_m3: np.ndarray) -> FusionComponent:
    """Convert one fusion pair kernel into an axial FusionComponent with reaction rate, neutron rate, and fusion power profiles"""
    profile = fusion_source_profile_from_rate_density(coordinate, reaction, pair_kernel.rate_density_m3_s, cell_volumes_m3)
  
    return FusionComponent(label=label, kind=kind, reaction=reaction, reactant_a_population_id=reactant_a_population_id, reactant_b_population_id=reactant_b_population_id, profile=profile, pair_kernel=pair_kernel)

def _pair_kernel_metadata(component: FusionComponent, cell_volumes_m3: np.ndarray) -> dict[str, object]:
    """Return one pair kernel rate, fit interval, identical pair, and reactant kinetic energy removal diagnostic record"""
    pair_kernel = component.pair_kernel
  
    if pair_kernel is None:
        raise ValueError(f"fusion component {component.label!r} has no pair kernel")
   
    total_rate_s = float(np.sum(pair_kernel.rate_density_m3_s * cell_volumes_m3))
    inside_rate_s = float(np.sum(pair_kernel.inside_fit_domain_rate_density_m3_s * cell_volumes_m3))
    outside_rate_s = float(np.sum(pair_kernel.outside_fit_domain_rate_density_m3_s * cell_volumes_m3))
    reactant_a_energy_removal_W = float(np.sum(pair_kernel.reactant_a_kinetic_energy_removal_density_W_m3 * cell_volumes_m3))
    reactant_b_energy_removal_W = float(np.sum(pair_kernel.reactant_b_kinetic_energy_removal_density_W_m3 * cell_volumes_m3))
    outside_fraction = outside_rate_s / max(total_rate_s, np.finfo(float).tiny)
   
    return {'label': component.label, 'kind': component.kind, 'reaction': component.reaction.key, 'reactant_a_population_id': component.reactant_a_population_id, 'reactant_b_population_id': component.reactant_b_population_id, 'kernel_model': pair_kernel.kernel_model, 'identical_population': pair_kernel.identical_population, 'pair_counting_factor': 0.5 if pair_kernel.identical_population else 1.0, 'num_gyroangle_points': pair_kernel.num_gyroangle_points, 'max_pair_states_per_batch': pair_kernel.max_pair_states_per_batch, 'fit_domain_E_cm_keV': list(pair_kernel.fit_domain_E_cm_keV), 'total_rate_s': total_rate_s, 'inside_fit_domain_rate_s': inside_rate_s, 'outside_fit_domain_rate_s': outside_rate_s, 'outside_fit_domain_fraction': outside_fraction, 'reactant_a_kinetic_energy_removal_W': reactant_a_energy_removal_W, 'reactant_b_kinetic_energy_removal_W': reactant_b_energy_removal_W, 'reactant_a_reaction_weighted_mean_kinetic_energy_J': reactant_a_energy_removal_W / max(total_rate_s, np.finfo(float).tiny), 'reactant_b_reaction_weighted_mean_kinetic_energy_J': reactant_b_energy_removal_W / max(total_rate_s, np.finfo(float).tiny)}

def _expander_population_authority_metadata(*, population_grid: FullDevicePopulationGrid, kinetic: KineticStageResult, expander: ExpanderStageResult | None, startup_seed_only: bool) -> dict[str, object]:
    """
    Check that fusion uses the source connected full device population state required by the active closure
    
    The startup seed path requires the mapped fast populations and does not admit an earlier linear expander approximation
    """
    if expander is None or expander.system_state is None:
        return {
            "fusion_expander_population_authority_measured": False,
            "fusion_expander_population_authority_passed": False,
            "fusion_expander_population_authority_model": "none",
            "fusion_legacy_linear_expander_population_used": False,
            "fusion_expander_population_authority_failure_reasons": (),
        }
   
    expected_expander_model = "source_connected_U_mu_conserving_residence_mapping"
    expected_lost_model = "Egedal_Eq63_source_connected_U_mu_mapping"
    reasons: list[str] = []
  
    if kinetic.metadata.get("full_device_expander_population_model") != expected_expander_model:
        reasons.append("kinetic_full_device_expander_model_is_not_source_connected_U_mu_mapping")
    if kinetic.metadata.get("full_device_lost_fast_ion_population_model") != expected_lost_model:
        reasons.append("kinetic_lost_fast_population_model_is_not_Egedal_Eq63_source_connected_U_mu_mapping")
    if kinetic.metadata.get("full_device_linear_expander_potential_is_authoritative") is not False:
        reasons.append("legacy_linear_expander_potential_remains_authoritative")
    if kinetic.metadata.get("expander_population_included_in_fusion_input") is not True:
        reasons.append("Pass13_expander_population_not_marked_for_fusion_input")
   
    directed_by_population = expander.system_state.full_device_distribution_by_population
    source_connected_by_species = kinetic.metadata.get("full_device_fast_ion_distribution_includes_source_connected_expander_by_species", {})
   
    if not isinstance(source_connected_by_species, dict):
        reasons.append("kinetic_source_connected_expander_species_flags_are_unavailable")
        source_connected_by_species = {}
 
    population_checks: dict[str, dict[str, object]] = {}
  
    for species_id, population_id in ((DEUTERON.species_id, FAST_D), (TRITON.species_id, FAST_T)):
        directed = directed_by_population.get(population_id)
        full_distribution = population_grid.fast_distribution_z_v_pitch_by_species.get(species_id)
        directed_present = bool(directed is not None and np.any(np.asarray(directed, dtype=float) > 0.0))
        shape_matches = bool(directed is None or full_distribution is not None and np.asarray(directed).shape == np.asarray(full_distribution).shape)
        source_connected_flag = bool(source_connected_by_species.get(species_id, False))
      
        if not shape_matches:
            reasons.append(f"{population_id}_expander_distribution_shape_mismatch")
        if directed_present and full_distribution is None:
            reasons.append(f"{population_id}_full_device_distribution_unavailable")
        if directed_present and not source_connected_flag:
            reasons.append(f"{population_id}_source_connected_expander_population_not_tagged")
     
        population_checks[population_id] = {"directed_distribution_present": directed_present, "full_device_distribution_available": full_distribution is not None, "shape_matches": shape_matches, "source_connected_expander_population_tagged": source_connected_flag}
   
    passed = not reasons
   
    if not passed:
        raise ValueError("Pass 13 expander population authority failed: " + ";".join(reasons))
    
    return {
        "fusion_expander_population_authority_measured": True,
        "fusion_expander_population_authority_passed": True,
        "fusion_expander_population_authority_model": ("Egedal_Eq63_kinetic_D_T_source_connected_U_mu_mapping_with_Eq70_center_nodal_potential" if startup_seed_only else "Egedal_Eq63_and_Sato_boundary_face_source_connected_U_mu_mapping_with_Eq70_center_nodal_potential"),
        "fusion_legacy_linear_expander_population_used": False,
        "fusion_expander_population_authority_by_population": population_checks,
        "fusion_expander_population_authority_failure_reasons": (),
    }

def _build_population_registry(*, population_grid: FullDevicePopulationGrid, kinetic: KineticStageResult) -> dict[str, FusionReactantPopulation]:
    """Build the fast D and fast T fusion reactant registry from full device local distributions"""
    populations: dict[str, FusionReactantPopulation] = {}
    speed_by_species = kinetic.local_speed_grid_by_species or {}
    fast_distribution_by_species = population_grid.fast_distribution_z_v_pitch_by_species
    fast_definition = ((DEUTERON.species_id, FAST_D, DEUTERON), (TRITON.species_id, FAST_T, TRITON))
   
    for species_id, population_id, species in fast_definition:
        distribution = fast_distribution_by_species.get(species_id)
        speed_grid = speed_by_species.get(species_id)
        if distribution is None or speed_grid is None:
            continue
        populations[population_id] = FusionReactantPopulation(population_id=population_id, species=species, population_kind="fast", speed_grid=speed_grid, pitch_grid=kinetic.pitch_grid, local_distribution_z_v_pitch=np.asarray(distribution, dtype=float), full_device_population=True, includes_directed_lost_population=bool(population_grid.includes_lost_fast_ions_by_species.get(species_id, False)))
  
    system = kinetic.fast_ion_system_state
  
    if system is not None:
        active_ids = set(system.active_species_ids)
        missing = [species_id for species_id, population_id, _ in fast_definition if species_id in active_ids and population_id not in populations]
        if missing:
            raise ValueError("active fast fusion species require full device distributions and species owned local speed grids")
        
    return populations

def _reaction_totals(components: tuple[FusionComponent, ...]) -> tuple[dict[str, float], dict[str, float], dict[str, float]]:
    """Sum fusion reaction rates, neutron rates, and fusion powers by reaction key"""
    reaction_rate_s: dict[str, float] = {}
    neutron_rate_s: dict[str, float] = {}
    power_W: dict[str, float] = {}
  
    for component in components:
        key = component.reaction.key
        reaction_rate_s[key] = reaction_rate_s.get(key, 0.0) + float(component.profile.total_reaction_rate_s or 0.0)
        neutron_rate_s[key] = neutron_rate_s.get(key, 0.0) + float(component.profile.total_neutron_rate_s or 0.0)
        power_W[key] = power_W.get(key, 0.0) + float(component.profile.total_fusion_power_W or 0.0)

    return reaction_rate_s, neutron_rate_s, power_W

def _population_burnup_totals(components: tuple[FusionComponent, ...], cell_volumes_m3: np.ndarray) -> tuple[dict[str, float], dict[str, float]]:
    """Return population resolved fusion consumption rates and reactant kinetic energy removal powers"""
    rates: dict[str, float] = {}
    powers: dict[str, float] = {}
    volumes = np.asarray(cell_volumes_m3, dtype=float)
   
    for component in components:
        pair = component.pair_kernel
        if pair is None:
            raise ValueError(f"fusion component {component.label!r} has no pair kernel")
        reaction_rate = float(np.sum(pair.rate_density_m3_s * volumes))
        energy_a = float(np.sum(pair.reactant_a_kinetic_energy_removal_density_W_m3 * volumes))
        energy_b = float(np.sum(pair.reactant_b_kinetic_energy_removal_density_W_m3 * volumes))
        rates[component.reactant_a_population_id] = rates.get(component.reactant_a_population_id, 0.0) + reaction_rate
        rates[component.reactant_b_population_id] = rates.get(component.reactant_b_population_id, 0.0) + reaction_rate
        powers[component.reactant_a_population_id] = powers.get(component.reactant_a_population_id, 0.0) + energy_a
        powers[component.reactant_b_population_id] = powers.get(component.reactant_b_population_id, 0.0) + energy_b

    return rates, powers

def build_fusion_stage(config: SourceModelRunConfig, geometry: GeometryStageResult, kinetic: KineticStageResult, expander: ExpanderStageResult | None = None) -> FusionStageResult:
    """
    Evaluate every active fusion population pair on the full device grid
    
    The same Table IV total cross section kernel supplies scalar reaction rates and downstream neutron component rates without posterior rate rescaling
    """
    fusion_runtime_started = perf_counter()
    runtime_profile = new_runtime_profile()
    fcfg = config.fusion
    max_pair_states_per_batch = int(getattr(fcfg, "max_pair_states_per_batch", 262144))
    startup_seed_only = uses_startup_seed_only(config.plasma_closure)
    population_grid_started = perf_counter()
    population_grid = build_full_device_population_grid(geometry, kinetic)
    record_runtime(runtime_profile, "fusion_population_grid", perf_counter() - population_grid_started)
   

    population_registry_started = perf_counter()
    expander_authority_metadata = _expander_population_authority_metadata(population_grid=population_grid, kinetic=kinetic, expander=expander, startup_seed_only=startup_seed_only)
    coordinate = np.asarray(population_grid.z_centers_m, dtype=float)
    volumes = np.asarray(population_grid.cell_volumes_m3, dtype=float)
    populations = _build_population_registry(population_grid=population_grid, kinetic=kinetic)
    components: list[FusionComponent] = []
   
    if startup_seed_only:
        unexpected_populations = tuple(sorted(set(populations) - {FAST_D, FAST_T}))
        if unexpected_populations:
            raise ValueError("NBI supported fusion requires kinetic D and T populations only")
        definitions = active_population_pair_definitions(include_fast_fast=True)
    else:
        definitions = active_population_pair_definitions(include_fast_fast=fcfg.include_fast_fast)
  
    record_runtime(runtime_profile, "fusion_population_registry", perf_counter() - population_registry_started)
 
    for definition in definitions:
        population_a = populations.get(definition.population_a_id)
        population_b = populations.get(definition.population_b_id)
      
        if population_a is None or population_b is None:
            continue
        for reaction in definition.reactions:
            identical_population = definition.identical_population
            label = f"{population_a.population_id}__{population_b.population_id}__{reaction.key}"
            pair_kernel_started = perf_counter()
            pair_kernel = bosch_hale_table_iv_pair_kernel_profiles(population_a.speed_grid, population_a.pitch_grid, population_a.local_distribution_z_v_pitch, population_a.species.mass_kg, population_b.speed_grid, population_b.pitch_grid, population_b.local_distribution_z_v_pitch, population_b.species.mass_kg, reaction, identical_population=identical_population, num_gyroangle_points=fcfg.num_gyroangle_points, max_pair_states_per_batch=max_pair_states_per_batch)
            components.append(_component_from_pair_kernel(label=label, kind=definition.component_kind, reaction=reaction, reactant_a_population_id=population_a.population_id, reactant_b_population_id=population_b.population_id, pair_kernel=pair_kernel, coordinate=coordinate, cell_volumes_m3=volumes))
            pair_kernel_runtime_s = perf_counter() - pair_kernel_started
         
            record_runtime(runtime_profile, "fusion_pair_kernels", pair_kernel_runtime_s)
            record_runtime(runtime_profile, f"fusion_pair_kernel_{label}", pair_kernel_runtime_s, nested=True)
 
    fusion_postprocessing_started = perf_counter()
   
    if not components:
        raise ValueError("fusion configuration produced no active components")
    component_tuple = tuple(components)
    combined = combine_fusion_source_profiles(tuple(component.profile for component in component_tuple))
    channel_axial_neutron_rate_s: dict[str, np.ndarray] = {}
    channel_axial_neutron_rate_density_m3_s: dict[str, np.ndarray] = {}
  
    for component in component_tuple:
        channel_axial_neutron_rate_density_m3_s[component.label] = component.profile.neutron_source_density_m3_s
        channel_axial_neutron_rate_s[component.label] = component.profile.neutron_source_density_m3_s * volumes
   
    pair_kernel_diagnostics = [_pair_kernel_metadata(component, volumes) for component in component_tuple]
    total_pair_rate_s = float(sum(float(item["total_rate_s"]) for item in pair_kernel_diagnostics))
    total_outside_rate_s = float(sum(float(item["outside_fit_domain_rate_s"]) for item in pair_kernel_diagnostics))
    outside_fraction = total_outside_rate_s / max(total_pair_rate_s, np.finfo(float).tiny)
    total_axial_neutron_rate_s = combined.total_neutron_source_density_m3_s * volumes
    component_axial_neutron_rate_sum_s = np.sum(tuple(channel_axial_neutron_rate_s.values()), axis=0) if channel_axial_neutron_rate_s else np.zeros_like(total_axial_neutron_rate_s)
    component_axial_neutron_rate_identity_error = float(np.max(np.abs(component_axial_neutron_rate_sum_s - total_axial_neutron_rate_s)) / max(float(np.max(np.abs(total_axial_neutron_rate_s))), np.finfo(float).tiny))
    reaction_rate_s, neutron_rate_s, power_W = _reaction_totals(component_tuple)
    population_burnup_rate_s, population_reactant_energy_removal_W = _population_burnup_totals(component_tuple, volumes)
    particle_balance_by_species = kinetic.metadata.get("nbi_supported_particle_balance_by_species", {})
  
    if not isinstance(particle_balance_by_species, dict):
        particle_balance_by_species = {}
   
    kinetic_system = kinetic.fast_ion_system_state
    charge_exchange_sink_active_by_species = {} if kinetic_system is None else {species_id: bool(state.modal_result is not None and state.modal_result.metadata.get("modal_charge_exchange_target_sink_active") is True) for species_id, state in kinetic_system.species_states.items()}
    burnup_rate_by_species = {DEUTERON.species_id: float(population_burnup_rate_s.get(FAST_D, 0.0)), TRITON.species_id: float(population_burnup_rate_s.get(FAST_T, 0.0))}
    burnup_fraction_by_species: dict[str, float | None] = {}
  
    for species_id, burnup_rate in burnup_rate_by_species.items():
        balance = particle_balance_by_species.get(species_id, {})
        if not isinstance(balance, dict):
            balance = {}
        net_source = float(balance.get("net_nbi_source_rate_s", 0.0) or 0.0)
        terminal_loss = float(balance.get("terminal_loss_rate_s", 0.0) or 0.0)
        if burnup_rate <= 0.0:
            burnup_fraction_by_species[species_id] = 0.0
        elif net_source > 0.0 and terminal_loss > 0.0:
            burnup_fraction_by_species[species_id] = max(burnup_rate / net_source, burnup_rate / terminal_loss)
        else:
            burnup_fraction_by_species[species_id] = None
  
    burnup_present = bool(any(value > 0.0 for value in burnup_rate_by_species.values()))
    active_burnup_species = tuple(species_id for species_id, value in burnup_rate_by_species.items() if value > 0.0)
    burnup_tolerance = float(config.kinetic_electrostatic.fast_fusion_burnup_relative_tolerance)
    burnup_assessed = bool(all(burnup_fraction_by_species[species_id] is not None for species_id in active_burnup_species))
    burnup_passed = bool(burnup_assessed and all(float(burnup_fraction_by_species[species_id]) <= burnup_tolerance for species_id in active_burnup_species))
    burnup_maximum_fraction = max((float(burnup_fraction_by_species[species_id]) for species_id in active_burnup_species if burnup_fraction_by_species[species_id] is not None), default=0.0 if not active_burnup_species else None)
    fast_population_ids = tuple(population_id for population_id in (FAST_D, FAST_T) if population_id in populations)
    primary_fast = populations.get(FAST_D) or populations.get(FAST_T)
    expander_available = expander is not None and expander.system_state is not None
    fast_distribution_source = "not_applicable_no_confined_fast_population" if primary_fast is None else "modal_phi_z_central_confined_plus_central_Eq63_lost_plus_Egedal_Eq63_source_connected_expander" if expander_available else "modal_phi_z_central_confined_plus_Egedal_Eq63_central_lost" if population_grid.includes_lost_fast_ions else "modal_phi_z_local_z_v_pitch"
    metadata = {**expander_authority_metadata, 'fusion_model': 'distribution_integrated_dd_dt', 'fusion_plasma_closure_model': config.plasma_closure.model, 'fusion_startup_seed_population_included': False, 'fusion_startup_seed_authority_fraction': 0.0 if startup_seed_only else None, 'fusion_nbi_supported_kinetic_population_only': bool(startup_seed_only), 'fusion_charge_exchange_target_sink_active_by_species': charge_exchange_sink_active_by_species if startup_seed_only else {}, 'fusion_rate_kernel_model': 'bosch_hale_table_iv_total_cross_section_quadrature', 'fusion_all_components_use_canonical_table_iv': True, 'fusion_table_iv_canonical_kernel_measured': True, 'fusion_table_iv_canonical_kernel_passed': bool(pair_kernel_diagnostics and all((item['kernel_model'] == 'bosch_hale_table_iv_total_cross_section_quadrature' for item in pair_kernel_diagnostics))), 'fusion_posterior_rate_normalization_applied': False, 'fusion_posterior_rate_normalization_absence_measured': True, 'fusion_posterior_rate_normalization_absence_passed': True, 'fusion_num_gyroangle_points': fcfg.num_gyroangle_points, 'fusion_max_pair_states_per_batch': max_pair_states_per_batch, 'fusion_pair_kernel_component_diagnostics': pair_kernel_diagnostics, 'fusion_fast_distribution_source': fast_distribution_source, 'fusion_fast_distribution_source_by_species': {populations[population_id].species.species_id: fast_distribution_source for population_id in fast_population_ids}, 'fusion_local_phi_z_distribution_used': bool(fast_population_ids), 'fusion_local_phi_z_fast_fast_used': bool(fast_population_ids and any((component.kind == 'fast_fast' for component in component_tuple))), 'fusion_local_distribution_z_bins': None if primary_fast is None else int(primary_fast.local_distribution_z_v_pitch.shape[0]), 'fusion_invariant_speed_grid_max_m_s': float(kinetic.speed_grid.faces_m_s[-1]), 'fusion_local_speed_grid_max_m_s': None if primary_fast is None else float(primary_fast.speed_grid.faces_m_s[-1]), 'fusion_distinct_local_speed_grid_used': bool(primary_fast is not None and primary_fast.speed_grid is not kinetic.speed_grid), 'fusion_local_speed_grid_max_m_s_by_population': {population_id: float(populations[population_id].speed_grid.faces_m_s[-1]) for population_id in fast_population_ids}, 'fusion_domain_model': 'full_device_population_coexistence_grid', 'fusion_expander_population_model': 'source_connected_U_mu_conserving_residence_mapping_with_Eq70_center_nodal_potential' if expander_available else 'none', 'fusion_expander_population_included': expander_available, 'fusion_expander_population_qualified': bool(expander_available and expander.system_state.status == 'qualified'), 'fusion_source_disconnected_expander_population_included': False, 'fusion_electron_density_required': False, 'fusion_electron_density_used': False, 'fusion_solved_electron_profile_scope': population_grid.electron_profile_scope, 'fusion_expander_electron_profile_available': bool(expander is not None and expander.system_state is not None and expander.system_state.potential_by_side), 'fusion_lost_fast_ion_population_included': bool(primary_fast is not None and primary_fast.includes_directed_lost_population), 'fusion_lost_fast_ion_population_included_by_population': {population_id: populations[population_id].includes_directed_lost_population for population_id in fast_population_ids}, 'fusion_lost_fast_ion_population_status': 'included_Egedal_Eq63_central_and_source_connected_expander' if expander_available and primary_fast is not None and primary_fast.includes_directed_lost_population else 'included_Egedal_Eq63_central' if primary_fast is not None and primary_fast.includes_directed_lost_population else 'unavailable_not_replaced_by_another_population', 'fusion_cross_section_model': 'bosch_hale_table_iv_total_cross_section', 'fusion_bosch_hale_fit_domain_coverage_required': True, 'fusion_bosch_hale_fit_domain_coverage_measured': True, 'fusion_bosch_hale_fit_domain_coverage_assessed': False, 'fusion_bosch_hale_fit_domain_coverage_passed': None, 'fusion_bosch_hale_fit_domain_coverage_assessment_method': 'reaction_rate_weighted_table_iv_pair_kernel_partition', 'fusion_bosch_hale_rate_weighted_out_of_fit_domain_fraction': outside_fraction, 'fusion_bosch_hale_rate_weighted_out_of_fit_domain_rate_s': total_outside_rate_s, 'fusion_bosch_hale_rate_weighted_total_rate_s': total_pair_rate_s, 'fusion_bosch_hale_rate_weighted_out_of_fit_domain_fraction_tolerance': None, 'fusion_bosch_hale_fit_domain_authoritative_tolerance_available': False, 'fusion_bosch_hale_cross_section_fit_domains_E_cm_keV': {item['reaction']: item['fit_domain_E_cm_keV'] for item in pair_kernel_diagnostics}, 'fusion_bosch_hale_fit_domain_reference': 'Bosch_and_Hale_1992_Nuclear_Fusion_32_611_Tables_IV_and_VII', 'active_fusion_channels': [component.label for component in component_tuple], 'fusion_population_ids': tuple(populations), 'fusion_population_species_by_id': {population_id: population.species.species_id for population_id, population in populations.items()}, 'fusion_population_display_label_by_id': {population_id: fusion_population_display_label(population_id) for population_id in populations}, 'fusion_population_kind_by_id': {population_id: population.population_kind for population_id, population in populations.items()}, 'fusion_population_full_device_available_by_id': {population_id: population.full_device_population for population_id, population in populations.items()}, 'fusion_population_includes_directed_lost_by_id': {population_id: population.includes_directed_lost_population for population_id, population in populations.items()}, 'fusion_component_population_pairs': {component.label: [component.reactant_a_population_id, component.reactant_b_population_id] for component in component_tuple}, 'fusion_component_kind_by_label': {component.label: component.kind for component in component_tuple}, 'fusion_component_reaction_by_label': {component.label: component.reaction.key for component in component_tuple}, 'fusion_component_pair_counting_factors': {component.label: 0.5 if component.pair_kernel.identical_population else 1.0 for component in component_tuple}, 'fusion_component_rates_s': {component.label: float(component.profile.total_reaction_rate_s or 0.0) for component in component_tuple}, 'fusion_component_neutron_rates_s': {component.label: float(component.profile.total_neutron_rate_s or 0.0) for component in component_tuple}, 'fusion_component_reactant_a_kinetic_energy_removal_W': {component.label: float(np.sum(component.pair_kernel.reactant_a_kinetic_energy_removal_density_W_m3 * volumes)) for component in component_tuple}, 'fusion_component_reactant_b_kinetic_energy_removal_W': {component.label: float(np.sum(component.pair_kernel.reactant_b_kinetic_energy_removal_density_W_m3 * volumes)) for component in component_tuple}, 'fusion_population_burnup_rate_s': population_burnup_rate_s, 'fusion_population_reactant_kinetic_energy_removal_W': population_reactant_energy_removal_W, 'fusion_fast_deuterium_burnup_rate_s': float(population_burnup_rate_s.get(FAST_D, 0.0)), 'fusion_fast_tritium_burnup_rate_s': float(population_burnup_rate_s.get(FAST_T, 0.0)), 'fusion_fast_deuterium_reactant_kinetic_energy_removal_W': float(population_reactant_energy_removal_W.get(FAST_D, 0.0)), 'fusion_fast_tritium_reactant_kinetic_energy_removal_W': float(population_reactant_energy_removal_W.get(FAST_T, 0.0)), 'fusion_fast_population_burnup_in_Eq59_sink': False, 'fusion_fast_population_burnup_negligibility_tolerance_available': bool(startup_seed_only), 'fusion_fast_population_burnup_negligibility_tolerance': burnup_tolerance if startup_seed_only else None, 'fusion_fast_population_burnup_negligibility_tolerance_source': 'kinetic.fast_fusion_burnup_relative_tolerance' if startup_seed_only else None, 'fusion_fast_population_burnup_relative_fraction_by_species': burnup_fraction_by_species if startup_seed_only else {}, 'fusion_fast_population_burnup_maximum_relative_fraction': burnup_maximum_fraction if startup_seed_only else None, 'fusion_fast_population_burnup_applicability_assessed': burnup_assessed if startup_seed_only else False, 'fusion_fast_population_burnup_applicability_passed': burnup_passed if startup_seed_only else False, 'fusion_fast_population_burnup_applicability_failure_reason': None if not startup_seed_only or not burnup_present or burnup_passed else 'fast_ion_fusion_burnup_not_negligible' if burnup_assessed else 'fast_ion_fusion_burnup_not_assessable_against_source_and_terminal_loss', 'fusion_reaction_rate_s_by_key': reaction_rate_s, 'fusion_neutron_rate_s_by_reaction_key': neutron_rate_s, 'fusion_power_W_by_reaction_key': power_W, 'fusion_DD_neutron_component_rate_sum_s': neutron_rate_s.get('dd_n', 0.0), 'fusion_DD_proton_component_rate_sum_s': reaction_rate_s.get('dd_p', 0.0), 'fusion_DT_neutron_component_rate_sum_s': neutron_rate_s.get('dt_n', 0.0), 'fusion_fast_T_channels_available': True, 'fusion_TT_channels_available': False, 'fusion_TT_channels_active': False, 'fusion_DT_canonical_deuteron_first_ordering': True, 'fusion_total_neutron_rate_density_m3_s': combined.total_neutron_source_density_m3_s, 'fusion_total_axial_neutron_rate_s': total_axial_neutron_rate_s, 'fusion_total_fusion_power_density_W_m3': combined.total_fusion_power_density_W_m3, 'fusion_channel_axial_neutron_rate_s': channel_axial_neutron_rate_s, 'fusion_channel_axial_neutron_rate_density_m3_s': channel_axial_neutron_rate_density_m3_s, 'fusion_component_axial_grid_id': 'full_device_population', 'fusion_component_axial_profile_count': len(channel_axial_neutron_rate_s), 'fusion_component_axial_profiles_complete': len(channel_axial_neutron_rate_s) == len(component_tuple), 'fusion_component_axial_neutron_rate_identity_relative_error': component_axial_neutron_rate_identity_error, 'fusion_component_axial_neutron_rate_identity_relative_tolerance': 1e-12, 'fusion_component_axial_neutron_rate_identity_check_passed': bool(component_axial_neutron_rate_identity_error <= 1e-12), 'total_neutron_rate_s': combined.total_neutron_rate_s, 'total_fusion_power_W': combined.total_fusion_power_W, 'total_charged_product_power_W': combined.total_charged_product_power_W}
   
    record_runtime(runtime_profile, "fusion_postprocessing", perf_counter() - fusion_postprocessing_started)
   
    metadata = {**metadata, "runtime_fusion_profile": finalize_runtime_profile(runtime_profile, total_s=perf_counter() - fusion_runtime_started)}
    
    return FusionStageResult(components=component_tuple, combined_profile=combined, populations_by_id=populations, metadata=metadata)
