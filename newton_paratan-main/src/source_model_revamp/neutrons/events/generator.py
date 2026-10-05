"""
Correlated D D and D T neutron event bank construction

Physical component and axial rates are fixed by the fusion calculation
Within each positive rate stratum, sampled reactant pairs are importance weighted by `sigma * g`, evaluated center of mass angles are sampled, and the correlated neutron state is transformed to the lab frame
"""
from __future__ import annotations
from dataclasses import dataclass
import numpy as np
from source_model_revamp.nuclear_data.endf_b_viii1 import load_evaluated_fusion_angular_distribution
from source_model_revamp.neutrons.events.angular_sampling import EvaluatedFusionAngularSampler
from source_model_revamp.neutrons.events.event_kinematics import correlated_neutron_lab_events, directions_from_axis_mu_phi, equivalent_deuteron_lab_energy_eV, projectile_direction_in_cm, reactant_pair_invariant_s_J2, residual_mass_shell_relative_error
from source_model_revamp.neutrons.events.pair_sampling import sample_reactant_pairs
from source_model_revamp.neutrons.events.spatial_sampling import sample_axisymmetric_cell_positions
from source_model_revamp.radial_profiles import KOTELNIKOV_PARABOLIC_FLUX_K1, canonical_radial_profile_model
from source_model_revamp.neutrons.events.types import CorrelatedNeutronEventBank, CorrelatedNeutronEventComponentSpec

_ANGULAR_COVERAGE_EXCLUDED_WEIGHT_TOLERANCE = 1.0e-10

@dataclass(frozen=True)
class _Stratum:
    """One component and axial cell event stratum with fixed physical neutron rate and allocated sample count"""
    component_index: int
    axial_index: int
    physical_rate_s: float
    event_count: int

def _allocate_stratified_event_counts(physical_rates_s: np.ndarray, total_event_count: int) -> np.ndarray:
    """
    Allocate a fixed event count across positive physical rate strata
    
    Each positive stratum receives at least one event and remaining events are apportioned by physical rate using largest remainder allocation
    """
    rates = np.asarray(physical_rates_s, dtype=float)
    count = int(total_event_count)
    if rates.ndim != 1 or np.any(~np.isfinite(rates)) or np.any(rates < 0.0):
        raise ValueError("physical_rates_s must be a finite nonnegative 1D array")
    positive = rates > 0.0
    positive_count = int(np.count_nonzero(positive))
    if positive_count == 0:
        raise ValueError("at least one event stratum must have positive rate")
    if count < positive_count:
        raise ValueError("correlated_event_count must be at least the number of positive component axial strata")
    allocation = np.zeros(rates.size, dtype=np.int64)
    # Preserve every positive physical rate stratum before distributing the remaining samples
    allocation[positive] = 1
    remaining = count - positive_count
    if remaining == 0:
        return allocation
    positive_rates = rates[positive]
    ideal = remaining * positive_rates / float(np.sum(positive_rates))
    base = np.floor(ideal).astype(np.int64)
    allocation[positive] += base
    leftover = remaining - int(np.sum(base))
    if leftover > 0:
        positive_indices = np.flatnonzero(positive)
        remainder = ideal - base
        order = np.argsort(-remainder, kind="stable")
        allocation[positive_indices[order[:leftover]]] += 1
    if int(np.sum(allocation)) != count:
        raise RuntimeError("stratified event allocation does not close")
    
    return allocation

def _strata_from_components(components: tuple[CorrelatedNeutronEventComponentSpec, ...], cell_volumes_m3: np.ndarray, event_count: int) -> tuple[_Stratum, ...]:
    """Build positive event strata from component rate density profiles and axial cell volumes"""
    volumes = np.asarray(cell_volumes_m3, dtype=float)
    flat_rates = []
    flat_indices = []
    for component_index, component in enumerate(components):
        if component.physical_rate_density_m3_s.shape != volumes.shape:
            raise ValueError("component rate profile does not match cell volumes")
        rates = component.physical_rate_density_m3_s * volumes
        for axial_index, rate in enumerate(rates):
            flat_rates.append(float(rate))
            flat_indices.append((component_index, axial_index))
    rate_array = np.asarray(flat_rates, dtype=float)
    allocation = _allocate_stratified_event_counts(rate_array, event_count)

    return tuple(_Stratum(component_index=component_index, axial_index=axial_index, physical_rate_s=float(rate_array[index]), event_count=int(allocation[index])) for index, (component_index, axial_index) in enumerate(flat_indices) if allocation[index] > 0)

def build_correlated_neutron_event_bank(components: tuple[CorrelatedNeutronEventComponentSpec, ...], *, z_edges_m: np.ndarray, cell_volumes_m3: np.ndarray, radial_inner_radius_m_by_z: np.ndarray, radial_outer_radius_m_by_z: np.ndarray, event_count: int, random_seed: int, radial_source_profile_model: str = KOTELNIKOV_PARABOLIC_FLUX_K1, source_connected_radial_rho_limit_edges: np.ndarray | None = None, radial_outer_radius_m_by_edge: np.ndarray | None = None) -> CorrelatedNeutronEventBank:
    """
    Generate a stratified correlated lab neutron event bank
    
    Each stratum carries its fixed physical rate fraction while sampled pair weights determine only the conditional event shape inside that stratum
    Events outside evaluated angular energy coverage retain zero probability weight for audit when their excluded reaction weight is below the configured tolerance
    """
    if not components:
        raise ValueError("at least one neutron producing component is required")
    volumes = np.asarray(cell_volumes_m3, dtype=float)
    if volumes.ndim != 1 or volumes.size == 0:
        raise ValueError("cell_volumes_m3 must be a nonempty 1D array")
    if np.any(~np.isfinite(volumes)) or np.any(volumes <= 0.0):
        raise ValueError("cell_volumes_m3 must be positive and finite")
    for component in components:
        if component.physical_rate_density_m3_s.shape != volumes.shape:
            raise ValueError("all components must use the shared axial grid")
    radial_model = canonical_radial_profile_model(radial_source_profile_model)
    z_edges = np.asarray(z_edges_m, dtype=float)
    radial_outer = np.asarray(radial_outer_radius_m_by_z, dtype=float)
    if radial_outer.shape != volumes.shape or np.any(~np.isfinite(radial_outer)) or np.any(radial_outer <= 0.0):
        raise ValueError("correlated event radial outer centers must match the axial cells and be positive")
    radial_outer_edges = None if radial_outer_radius_m_by_edge is None else np.asarray(radial_outer_radius_m_by_edge, dtype=float)
    if radial_outer_edges is not None and (radial_outer_edges.shape != z_edges.shape or np.any(~np.isfinite(radial_outer_edges)) or np.any(radial_outer_edges <= 0.0)):
        raise ValueError("correlated event radial outer edges must match the axial edges and be positive")
    if source_connected_radial_rho_limit_edges is None:
        radial_rho_edges = np.ones(z_edges.size, dtype=float)
    else:
        radial_rho_edges = np.asarray(source_connected_radial_rho_limit_edges, dtype=float)
    if radial_rho_edges.shape != z_edges.shape:
        raise ValueError("correlated event radial support must match the axial edge grid")
    seed = int(random_seed)
    if seed < 0:
        raise ValueError("random_seed must be nonnegative")
    rng = np.random.default_rng(seed)
    strata = _strata_from_components(components, volumes, int(event_count))
    total_rate = float(sum(stratum.physical_rate_s for stratum in strata))
    if total_rate <= 0.0:
        raise ValueError("correlated event bank has zero physical neutron rate")
    angular_samplers = {key: EvaluatedFusionAngularSampler(load_evaluated_fusion_angular_distribution(key)) for key in {component.reaction.key for component in components}}
    angular_validation = {key: dict(sampler.distribution.runtime_validation) for key, sampler in angular_samplers.items()}
    angular_validation_passed = bool(angular_validation and all(
            validation.get("angular_law_validation_passed") is True
            and validation.get("angular_normalization_measured") is True
            and float(validation["maximum_angular_normalization_error"])
            <= float(validation["angular_normalization_tolerance"])
            and float(validation["minimum_angular_probability_density_per_mu"])
            >= float(validation["minimum_angular_probability_tolerance"])
            for validation in angular_validation.values()
        ))
    if not angular_validation_passed:
        raise RuntimeError("ENDF B VIII.1 angular data validation is unavailable")

    positions_parts = []
    directions_parts = []
    energy_parts = []
    weight_parts = []
    reaction_parts = []
    component_parts = []
    axial_parts = []
    center_energy_parts = []
    incident_energy_parts = []
    mu_parts = []
    velocity_a_parts = []
    velocity_b_parts = []
    stratum_diagnostics = []
    maximum_mass_shell_error = 0.0
    total_excluded_importance = 0.0
    total_importance = 0.0

    for stratum in strata:
        component = components[stratum.component_index]
        sampled = sample_reactant_pairs(component, stratum.axial_index, stratum.event_count, rng)
        invariant = reactant_pair_invariant_s_J2(component.reaction_kinematics, sampled.velocity_a_m_s, sampled.velocity_b_m_s)
        incident_energy = equivalent_deuteron_lab_energy_eV(component.reaction_kinematics, invariant)
        angular_distribution = angular_samplers[component.reaction.key].distribution
        lower_energy, upper_energy = angular_distribution.incident_energy_bounds_eV
        in_coverage = ((incident_energy >= lower_energy) & (incident_energy <= upper_energy))
        importance = sampled.sigma_g_m3_s.copy()
        excluded_importance = float(np.sum(importance[~in_coverage]))
        total_stratum_importance = float(np.sum(importance))
        total_excluded_importance += excluded_importance
        total_importance += total_stratum_importance
        if total_stratum_importance <= 0.0:
            raise RuntimeError(f"event stratum {component.label} z={stratum.axial_index} sampled zero reaction weight")
        excluded_fraction = excluded_importance / total_stratum_importance
        if excluded_fraction > _ANGULAR_COVERAGE_EXCLUDED_WEIGHT_TOLERANCE:
            raise RuntimeError(f"event stratum {component.label} has material reaction weight outside the ENDF angular range")
        importance[~in_coverage] = 0.0
        represented_importance = float(np.sum(importance))
        if represented_importance <= 0.0:
            raise RuntimeError("ENDF angular coverage removed all sampled reaction weight")
        # Physical stratum rate and conditional sampled shape remain separate normalizations
        normalized_local = importance / represented_importance
        normalized_weight = (stratum.physical_rate_s / total_rate) * normalized_local
        mu = np.zeros(stratum.event_count, dtype=float)
        positive_weight = importance > 0.0
        if np.any(positive_weight):
            mu[positive_weight] = angular_samplers[component.reaction.key].sample_mu(incident_energy[positive_weight], rng)
        projectile_axis = projectile_direction_in_cm(component.reaction_kinematics, sampled.velocity_a_m_s, sampled.velocity_b_m_s)
        cm_direction = directions_from_axis_mu_phi(projectile_axis, mu, rng.uniform(0.0, 2.0 * np.pi, stratum.event_count))
        lab_energy, lab_direction = correlated_neutron_lab_events(component.reaction_kinematics, sampled.velocity_a_m_s, sampled.velocity_b_m_s, cm_direction)
        mass_shell_error = residual_mass_shell_relative_error(component.reaction_kinematics, sampled.velocity_a_m_s, sampled.velocity_b_m_s, lab_energy, lab_direction)
        maximum_mass_shell_error = max(maximum_mass_shell_error, float(np.max(mass_shell_error)))
        axial_indices = np.full(stratum.event_count, stratum.axial_index, dtype=np.int64)
        positions = sample_axisymmetric_cell_positions(
            axial_indices,
            z_edges,
            radial_inner_radius_m_by_z,
            radial_outer,
            rng,
            radial_profile_model=radial_model,
            source_connected_radial_rho_limit_edges=radial_rho_edges,
            radial_outer_radius_m_by_edge=radial_outer_edges,
        )
        canonical_rate_density = float(component.physical_rate_density_m3_s[stratum.axial_index])
        raw_relative_error = (sampled.raw_rate_density_estimate_m3_s - canonical_rate_density) / max(canonical_rate_density, np.finfo(float).tiny)
        stratum_diagnostics.append({
                "component_label": component.label,
                "reaction": component.reaction.key,
                "axial_index": stratum.axial_index,
                "event_count": stratum.event_count,
                "physical_rate_s": stratum.physical_rate_s,
                "canonical_rate_density_m3_s": canonical_rate_density,
                "raw_importance_rate_density_estimate_m3_s": (sampled.raw_rate_density_estimate_m3_s),
                "raw_importance_to_canonical_relative_error": raw_relative_error,
                "effective_sample_size": sampled.effective_sample_size,
                "effective_sample_fraction": (sampled.effective_sample_size / stratum.event_count),
                "angular_coverage_excluded_importance_fraction": excluded_fraction,
                "projectile_assignment_model": component.projectile_assignment_model,
                "projectile_exchange_count": sampled.projectile_exchange_count,
                "projectile_exchange_fraction": (sampled.projectile_exchange_count / stratum.event_count),
            })
        positions_parts.append(positions)
        directions_parts.append(lab_direction)
        energy_parts.append(lab_energy)
        weight_parts.append(normalized_weight)
        reaction_parts.append(np.full(stratum.event_count, component.reaction.key, dtype="U8"))
        component_parts.append(np.full(stratum.event_count, component.label, dtype="U96"))
        axial_parts.append(axial_indices)
        center_energy_parts.append(sampled.center_of_mass_energy_J)
        incident_energy_parts.append(incident_energy)
        mu_parts.append(mu)
        velocity_a_parts.append(sampled.velocity_a_m_s)
        velocity_b_parts.append(sampled.velocity_b_m_s)
    normalized_weights = np.concatenate(weight_parts)
    normalized_weights /= float(np.sum(normalized_weights))
    energies = np.concatenate(energy_parts)
    directions = np.concatenate(directions_parts)
    weighted_mean_energy = float(np.sum(normalized_weights * energies))
    centered_energy = energies - weighted_mean_energy
    weighted_direction_mean = np.sum(normalized_weights[:, None] * directions, axis=0)
    energy_direction_covariance = np.sum(normalized_weights[:, None] * centered_energy[:, None] * (directions - weighted_direction_mean[None, :]), axis=0)
    all_positions = np.concatenate(positions_parts)
    all_axial = np.concatenate(axial_parts)
    radial_distance = np.sqrt(all_positions[:, 0] ** 2 + all_positions[:, 1] ** 2)
    cell_width = z_edges[all_axial + 1] - z_edges[all_axial]
    axial_fraction = np.divide(all_positions[:, 2] - z_edges[all_axial], cell_width, out=np.zeros_like(radial_distance), where=cell_width > 0.0)
    sampled_outer = radial_outer[all_axial] if radial_outer_edges is None else radial_outer_edges[all_axial] + axial_fraction * (radial_outer_edges[all_axial + 1] - radial_outer_edges[all_axial])
    sampled_rho = np.divide(radial_distance, sampled_outer, out=np.zeros_like(radial_distance), where=sampled_outer > 0.0)
    sampled_rho_limit = radial_rho_edges[all_axial] + axial_fraction * (radial_rho_edges[all_axial + 1] - radial_rho_edges[all_axial])
    maximum_radial_support_excess = float(np.max(np.maximum(sampled_rho - sampled_rho_limit, 0.0), initial=0.0))
    metadata = {
        "event_model": "dress_style_stratified_importance_sampling",
        "radial_source_profile_model": radial_model,
        "radial_source_profile_role": "prescribed_neutron_birth_position_closure_after_authoritative_one_dimensional_fusion_rate",
        "radial_source_profile_is_radially_self_consistent_FBIS": False,
        "radial_coordinate_definition": "rho_equals_r_over_nominal_a_of_z_and_psi_tilde_equals_rho_squared",
        "radial_source_shape_definition": "G_of_rho_equals_2_times_one_minus_rho_squared",
        "radial_probability_density_definition": "p_of_rho_equals_4_rho_times_one_minus_rho_squared",
        "radial_geometry_conditioned_source_shape_definition": "G_of_rho_divided_by_F_of_rho_max_inside_source_connected_support",
        "radial_profile_cross_section_average_normalized": True,
        "radial_profile_reference": "Kotelnikov_et_al_Journal_of_Plasma_Physics_2025_Eq3_11_k1",
        "radial_profile_reference_role": "published_pressure_profile_used_as_a_prescribed_reduced_neutron_spatial_source_closure_not_a_published_neutron_emissivity_law",
        "radial_profile_unconstrained_fraction_inside_half_radius": 0.4375,
        "radial_profile_unconstrained_mean_rho_squared": 1.0 / 3.0,
        "radial_source_connected_geometry_truncation_active": bool(np.any(radial_rho_edges < 1.0 - 1.0e-12)),
        "radial_source_connected_rho_limit_edges": radial_rho_edges,
        "radial_sampling_uses_flux_tube_radius_interpolated_from_axial_edges": bool(radial_outer_edges is not None),
        "radial_sample_maximum_support_excess": maximum_radial_support_excess,
        "radial_sample_support_check_passed": bool(maximum_radial_support_excess <= 2.0e-12),
        "angular_model": "endf_b_viii1_cm",
        "angular_distribution_role": "normalized_conditional_cm_probability",
        "total_cross_section_model": "bosch_hale_table_iv",
        "reactant_velocity_sampling_model": ("canonical_speed_pitch_cell_centers_with_uniform_absolute_gyroangle"),
        "relative_gyroangle_sampling_model": ("uniform_sampling_of_canonical_midpoint_quadrature_nodes"),
        "event_count": int(normalized_weights.size),
        "positive_weight_event_count": int(np.count_nonzero(normalized_weights > 0.0)),
        "zero_weight_event_count": int(np.count_nonzero(normalized_weights == 0.0)),
        "random_seed": seed,
        "physical_total_rate_s": total_rate,
        "normalized_weight_sum": float(np.sum(normalized_weights)),
        "event_probability_normalization_measured": True,
        "event_probability_normalization_passed": bool(np.isclose(np.sum(normalized_weights), 1.0, rtol=0.0, atol=2.0e-12)),
        "stratum_count": len(strata),
        "stratum_diagnostics": stratum_diagnostics,
        "maximum_residual_mass_shell_relative_error": maximum_mass_shell_error,
        "angular_coverage_excluded_importance_fraction": (total_excluded_importance / max(total_importance, np.finfo(float).tiny)),
        "angular_coverage_excluded_importance_tolerance": (_ANGULAR_COVERAGE_EXCLUDED_WEIGHT_TOLERANCE),
        "angular_coverage_measured": True,
        "angular_coverage_passed": bool(total_excluded_importance / max(total_importance, np.finfo(float).tiny) <= _ANGULAR_COVERAGE_EXCLUDED_WEIGHT_TOLERANCE),
        "angular_law_validation_passed": True,
        "endf_b_viii1_provenance_valid": True,
        "endf_b_viii1_reaction_validation": angular_validation,
        "weighted_mean_energy_J": weighted_mean_energy,
        "weighted_mean_direction": weighted_direction_mean,
        "weighted_energy_direction_covariance_J": energy_direction_covariance,
        "event_component_and_axial_rates_stratified_exactly": True,
        "canonical_component_axial_rate_source": ("bosch_hale_table_iv_deterministic_pair_kernel"),
        "conditional_importance_weight_normalization_role": ("sample_joint_event_shape_within_each_fixed_physical_rate_stratum"),
        "posterior_physical_rate_normalization_applied": False,
        "conditional_event_weights_normalized": True,
        "zero_weight_events_retained_for_audit_only": True,
    }

    return CorrelatedNeutronEventBank(
        positions_m=all_positions,
        directions=directions,
        energies_J=energies,
        normalized_weights=normalized_weights,
        physical_total_rate_s=total_rate,
        reaction_keys=np.concatenate(reaction_parts),
        component_labels=np.concatenate(component_parts),
        axial_cell_indices=all_axial,
        center_of_mass_energy_J=np.concatenate(center_energy_parts),
        equivalent_deuteron_lab_energy_eV=np.concatenate(incident_energy_parts),
        cm_emission_mu=np.concatenate(mu_parts),
        reactant_a_velocity_m_s=np.concatenate(velocity_a_parts),
        reactant_b_velocity_m_s=np.concatenate(velocity_b_parts),
        metadata=metadata,
    )

def event_bank_source_rate_matrix_s(event_bank: CorrelatedNeutronEventBank, energy_edges_J: np.ndarray, axial_cell_count: int, *, allow_out_of_range: bool) -> tuple[np.ndarray, float]:
    """Histogram physical event rate weights by axial cell and neutron energy bin"""
    edges = np.asarray(energy_edges_J, dtype=float)
    if edges.ndim != 1 or edges.size < 2 or np.any(np.diff(edges) <= 0.0):
        raise ValueError("energy_edges_J must be strictly increasing")
    axial_count = int(axial_cell_count)
    if axial_count < 1:
        raise ValueError("axial_cell_count must be positive")
    energy = event_bank.energies_J
    axial = event_bank.axial_cell_indices
    physical_weight = event_bank.physical_rate_weights_s
    inside = (energy >= edges[0]) & (energy < edges[-1])
    inside |= energy == edges[-1]
    out_of_range_rate = float(np.sum(physical_weight[~inside]))
    if (out_of_range_rate > 1.0e-12 * max(event_bank.physical_total_rate_s, np.finfo(float).tiny) and not allow_out_of_range):
        raise ValueError("correlated neutron events lie outside the configured energy range")
    energy_bin = np.searchsorted(edges, energy[inside], side="right") - 1
    energy_bin[energy[inside] == edges[-1]] = edges.size - 2
    matrix = np.zeros((axial_count, edges.size - 1), dtype=float)
    np.add.at(matrix, (axial[inside], energy_bin), physical_weight[inside])

    return matrix, out_of_range_rate