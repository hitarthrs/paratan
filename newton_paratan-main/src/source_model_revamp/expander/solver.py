"""Coupled Eq 63 ion mapping and Eq 70 quasineutral expander solve"""
from __future__ import annotations
from dataclasses import dataclass
import numpy as np
from scipy.optimize import least_squares
from source_model_revamp.constants import ELECTRON_CHARGE_C, KEV_TO_J
from source_model_revamp.expander.mapping import SourceConnectedDensityPlan, map_source_connected_branch, prepare_source_connected_density_plan, source_connected_density_profile
from source_model_revamp.expander.types import ExpanderBranchState, ExpanderPotentialState, ExpanderSystemState
from source_model_revamp.fbis.species import DEUTERON, TRITON, IonSpecies
from source_model_revamp.fbis.velocity_space_grid import PitchGrid, SpeedGrid, gyrotropic_velocity_cell_volumes
from source_model_revamp.fbis.modal.electrostatic.eq70 import _electron_density_fraction_eq70
from source_model_revamp.fusion.populations import FAST_D, FAST_T
from source_model_revamp.geometry.plasma_boundaries import PlasmaFacingSurface, source_connected_radial_support
from source_model_revamp.integration.config.root import SourceModelRunConfig
from source_model_revamp.integration.pipeline_types import GeometryStageResult, KineticStageResult

@dataclass(frozen=True)
class _ExpanderPath:
    """Store geometry from one mirror throat to its selected terminal surface"""
    side: str
    cell_indices: np.ndarray
    z_edges_m: np.ndarray
    segment_lengths_m: np.ndarray
    B_tilde_edges: np.ndarray
    cell_volumes_m3: np.ndarray
    nominal_outer_radius_edges_m: np.ndarray
    radial_rho_limit_edges: np.ndarray
    radial_survival_fraction_edges: np.ndarray
    radial_material_surface_id_by_segment: tuple[str | None, ...]
    radial_profile_model: str
    surface: PlasmaFacingSurface
    terminal_surface_aligned_to_population_grid: bool
    terminal_cell_axial_fraction: float

@dataclass(frozen=True)
class _CenterNodalPathGeometry:
    """Store the synthetic edge and center path used for the Eq 70 solve"""
    z_nodes_m: np.ndarray
    B_tilde_nodes: np.ndarray
    radial_rho_limit_nodes: np.ndarray
    radial_material_surface_id_by_segment: tuple[str | None, ...]
    synthetic_cell_indices: np.ndarray
    parent_cell_volumes_m3: np.ndarray
    segment_lengths_m: np.ndarray

@dataclass(frozen=True)
class _BranchInput:
    """Store directed Eq 63 throat source data for one species and one side"""
    side: str
    population_id: str
    component_id: str
    population_kind: str
    species: IonSpecies
    source_model: str
    particle_rate_s: np.ndarray
    invariant_total_energy_J: np.ndarray
    mu_B0_energy_J: np.ndarray
    target_speed_grid: SpeedGrid
    target_pitch_grid: PitchGrid

def _relative_error(value: float, reference: float) -> float:
    """Return the relative scalar difference using the reference magnitude"""
    return abs(float(value) - float(reference)) / max(abs(float(reference)), np.finfo(float).tiny)

def _electron_temperature_J(kinetic: KineticStageResult) -> float:
    """Return the active electron temperature in J from kinetic metadata"""
    for key in ("operating_point_electron_temperature_J", "electron_temperature_solved_J"):
        value = kinetic.metadata.get(key)
        if value is not None:
            temperature = float(value)
            if np.isfinite(temperature) and temperature > 0.0:
                return temperature
            
    for key in ("electron_temperature_solved_keV", "collision_state_electron_temperature_keV", "electron_temperature_fixed_input_keV", "electron_temperature_keV"):
        value = kinetic.metadata.get(key)
        if value is not None:
            temperature = float(value) * KEV_TO_J
            if np.isfinite(temperature) and temperature > 0.0:
                return temperature
    raise ValueError("electron temperature is unavailable for the expander solve")

def _path_to_selected_terminal(geometry: GeometryStageResult, side: str, radial_profile_model: str) -> _ExpanderPath:
    """Build the outward path from one throat to its selected floating axial terminal"""
    if geometry.full_device_z_edges_m is None or geometry.full_device_cell_volumes_m3 is None or geometry.B_T_function is None or geometry.plasma_boundaries is None:
        raise ValueError("full device geometry direct field and plasma boundaries are required for the expander solve")
    branch = str(side).strip().lower()
    if branch not in {"left", "right"}:
        raise ValueError("expander side must be left or right")
    full_edges = np.asarray(geometry.full_device_z_edges_m, dtype=float)
    full_volumes = np.asarray(geometry.full_device_cell_volumes_m3, dtype=float)
    throat_z = float(geometry.z_edges_m[0] if branch == "left" else geometry.z_edges_m[-1])
    throat_edge = int(np.argmin(np.abs(full_edges - throat_z)))
    tolerance = 256.0 * np.finfo(float).eps * max(float(np.max(np.abs(full_edges))), abs(throat_z), 1.0)
    if abs(float(full_edges[throat_edge]) - throat_z) > tolerance:
        raise ValueError(f"{branch} throat is not an exact full device grid edge")
    surface = geometry.plasma_boundaries.left_selected_loss_surface if branch == "left" else geometry.plasma_boundaries.right_selected_loss_surface
    if not surface.is_terminal or surface.electrical_model != "floating":
        raise ValueError(f"{branch} selected terminal surface must be a floating terminal material interface")
    if surface.geometry != "axial_plane" or surface.z_m is None:
        raise ValueError(f"{branch} selected terminal surface must be an axial plane for the reduced expander path")
    terminal_z = float(surface.z_m)
    terminal_edge = int(np.argmin(np.abs(full_edges - terminal_z)))
    if abs(float(full_edges[terminal_edge]) - terminal_z) > tolerance:
        raise ValueError(f"{branch} selected terminal surface is not an exact full device grid edge")
    if branch == "left":
        if terminal_edge >= throat_edge:
            raise ValueError("left terminal surface must lie outward of the left throat")
        path_edges = full_edges[terminal_edge:throat_edge + 1][::-1].copy()
        path_cells = np.arange(throat_edge - 1, terminal_edge - 1, -1, dtype=int)
    else:
        if terminal_edge <= throat_edge:
            raise ValueError("right terminal surface must lie outward of the right throat")
        path_edges = full_edges[throat_edge:terminal_edge + 1].copy()
        path_cells = np.arange(throat_edge, terminal_edge, dtype=int)
    path_B = np.asarray([float(geometry.B_T_function(float(value))) / float(geometry.B0_T) for value in path_edges], dtype=float)
    # Paraxial flux conservation gives the nominal radius proportional to 1 / sqrt(B_tilde)
    nominal_radius = float(geometry.plasma_radius_m) / np.sqrt(path_B)
    segment_lengths = np.hypot(np.diff(path_edges), np.diff(nominal_radius))
    if np.any(segment_lengths <= 0.0):
        raise ValueError(f"{branch} expander path contains a zero length segment")
    radial_support = source_connected_radial_support(path_z_m=path_edges, nominal_outer_radius_m=nominal_radius, surfaces=geometry.plasma_boundaries.candidate_surfaces, side=branch, terminal_surface_id=surface.surface_id, radial_profile_model=radial_profile_model)
  
    return _ExpanderPath(
        side=branch,
        cell_indices=path_cells,
        z_edges_m=path_edges,
        segment_lengths_m=segment_lengths,
        B_tilde_edges=path_B,
        cell_volumes_m3=full_volumes[path_cells],
        nominal_outer_radius_edges_m=nominal_radius,
        radial_rho_limit_edges=np.asarray(radial_support.rho_limit_edges, dtype=float),
        radial_survival_fraction_edges=np.asarray(radial_support.survival_probability_edges, dtype=float),
        radial_material_surface_id_by_segment=tuple(radial_support.segment_material_surface_ids),
        radial_profile_model=str(radial_support.radial_profile_model),
        surface=surface,
        terminal_surface_aligned_to_population_grid=True,
        terminal_cell_axial_fraction=1.0,
    )

def _truncated_exponential_flux_quantiles(a: np.ndarray, node_count: int) -> np.ndarray:
    """Return midpoint quantiles of p(x) ∝ exp(−a x) on 0 <= x <= 1"""
    values = np.asarray(a, dtype=float)
    count = int(node_count)
    if count < 2:
        raise ValueError("expander fast throat pitch quantile node count must be at least two")
    if np.any(~np.isfinite(values)) or np.any(values < 0.0):
        raise ValueError("Eq 63 truncated exponential parameter must be finite and nonnegative")
    probability = (np.arange(count, dtype=float) + 0.5) / float(count)
    probability = np.broadcast_to(probability[None, :], (values.size, count))
    result = np.empty_like(probability)
    small = values < 1.0e-8
    if np.any(small):
        result[small, :] = probability[small, :]
    regular = ~small
    if np.any(regular):
        normalization = -np.expm1(-values[regular])
        result[regular, :] = -np.log1p(-normalization[:, None] * probability[regular, :]) / values[regular, None]
  
    return np.clip(result, 0.0, 1.0)

def _fast_branch_inputs(config: SourceModelRunConfig, kinetic: KineticStageResult) -> tuple[list[_BranchInput], dict[str, float], dict[str, object]]:
    """Convert Eq 63 one end speed spectra into directed invariant throat samples"""
    system = kinetic.fast_ion_system_state
    if system is None:
        return [], {}, {}
    local_grids = dict(kinetic.local_speed_grid_by_species or {})
    inputs: list[_BranchInput] = []
    prompt_rates: dict[str, float] = {}
    audit: dict[str, object] = {}
    node_count = int(config.kinetic_electrostatic.expander_fast_throat_pitch_quantile_nodes)
    identities = ((DEUTERON, FAST_D), (TRITON, FAST_T))
    for species, population_id in identities:
        state = system.species_states.get(species.species_id)
        if state is None:
            continue
        if state.modal_result is None:
            prompt = float(state.prompt_only_particle_loss_rate_s)
            if prompt > 0.0:
                prompt_rates[f"prompt_{species.species_id}"] = prompt
            continue
        modal = state.modal_result
        prompt = float(modal.metadata.get("modal_prompt_loss_birth_rate_s", 0.0))
        if prompt > 0.0:
            prompt_rates[f"prompt_{species.species_id}"] = prompt
        speed = np.asarray(modal.speed_grid.centers_m_s, dtype=float)
        invariant_energy = 0.5 * float(species.mass_kg) * speed**2
        lambda_boundary = float(modal.basis.lambda_boundary)
        if not np.isfinite(lambda_boundary) or lambda_boundary <= 0.0:
            raise ValueError(f"{species.species_id} magnetic loss boundary is invalid")
        mirror_ratio = 1.0 / lambda_boundary
        target_speed = local_grids.get(species.species_id) or (modal.local_reconstruction.local_speed_grid if modal.local_reconstruction is not None and modal.local_reconstruction.local_speed_grid is not None else modal.speed_grid)
        if modal.pitch_grid is None:
            raise ValueError("active fast ion modal result requires a pitch grid for expander mapping")
        species_audit: dict[str, object] = {}
        for side, key, temperature_key in (
            ("left", "modal_eq63_left_throat_rate_v_lambda_s", "lost_ion_parallel_temperature_left_J"),
            ("right", "modal_eq63_right_throat_rate_v_lambda_s", "lost_ion_parallel_temperature_right_J"),
        ):
            directed = modal.metadata.get(key)
            temperature = modal.metadata.get(temperature_key)
            if directed is None or temperature is None:
                raise ValueError(f"{species.species_id} {side} Eq 63 throat state is unavailable")
            directed_array = np.asarray(directed, dtype=float)
            if directed_array.ndim != 2 or directed_array.shape[0] != speed.size:
                raise ValueError(f"{species.species_id} {side} Eq 63 throat rate has an incompatible shape")
            if np.any(~np.isfinite(directed_array)) or np.any(directed_array < 0.0):
                raise ValueError(f"{species.species_id} {side} Eq 63 throat rate must be finite and nonnegative")
            rate_v = np.sum(directed_array, axis=1)
            T_L = float(temperature)
            if not np.isfinite(T_L) or T_L <= 0.0:
                raise ValueError(f"{species.species_id} {side} Eq 63 parallel temperature must be positive and finite")
            # Eq 63 uses x = 1 − R_M Λ with flux density proportional to exp(−U x / T_L)
            a = invariant_energy / T_L
            x = _truncated_exponential_flux_quantiles(a, node_count)
            Lambda = (1.0 - x) / mirror_ratio
            rate = np.broadcast_to(rate_v[:, None] / float(node_count), x.shape).copy()
            invariant = np.broadcast_to(invariant_energy[:, None], x.shape).copy()
            mu_B0 = invariant * Lambda
            represented_rate_v = np.sum(rate, axis=1)
            target_total = float(np.sum(rate_v))
            spectrum_error = float(np.sum(np.abs(represented_rate_v - rate_v)) / max(target_total, np.finfo(float).tiny))
            if spectrum_error > 4096.0 * np.finfo(float).eps:
                raise ValueError(f"{species.species_id} {side} Eq 63 quantile source failed to preserve the Eq 61 speed spectrum")
            inputs.append(_BranchInput(
                side=side,
                population_id=population_id,
                component_id=f"fast_{species.species_id}",
                population_kind="fast",
                species=species,
                source_model="Egedal_Eq63_Baldwin_TL_exact_truncated_exponential_flux_quantiles",
                particle_rate_s=rate,
                invariant_total_energy_J=invariant,
                mu_B0_energy_J=mu_B0,
                target_speed_grid=target_speed,
                target_pitch_grid=modal.pitch_grid,
            ))
            species_audit[side] = {
                "parallel_temperature_J": T_L,
                "parallel_temperature_model": modal.metadata.get("lost_ion_parallel_temperature_model"),
                "one_end_Eq61_represented_rate_s": target_total,
                "quantile_represented_rate_s": float(np.sum(rate)),
                "speed_spectrum_relative_error": spectrum_error,
                "pitch_quantile_node_count_per_speed_bin": node_count,
                "maximum_empirical_CDF_midpoint_error_bound": 0.5 / float(node_count),
                "source_pitch_variable": "x_equals_one_minus_RM_Lambda",
                "source_flux_PDF": "proportional_to_exp_minus_U_x_over_T_L_on_zero_to_one",
                "ordinary_lambda_cell_midpoints_used": False,
            }
        audit[population_id] = species_audit

    return inputs, prompt_rates, audit


def _initial_potential(path: _ExpanderPath, throat_drop_J: float, wall_drop_J: float) -> np.ndarray:
    """Build the linear throat to wall potential initial guess on path edges"""
    coordinate = np.concatenate(([0.0], np.cumsum(path.segment_lengths_m)))
    fraction = coordinate / max(float(coordinate[-1]), np.finfo(float).tiny)

    return float(throat_drop_J) + fraction * (float(wall_drop_J) - float(throat_drop_J))

def _map_inputs(*, inputs: list[_BranchInput], path: _ExpanderPath, potential_edges_J: np.ndarray, terminal_wall_potential_drop_J: float, geometry: GeometryStageResult, tolerance: float) -> dict[str, ExpanderBranchState]:
    """Map every directed source for one side on a supplied edge potential"""
    mapped: dict[str, ExpanderBranchState] = {}
    full_count = int(np.asarray(geometry.full_device_z_centers_m).size)
    for item in inputs:
        if item.side != path.side:
            continue
        key = f"{item.population_id}:{item.side}"
        mapped[key] = map_source_connected_branch(side=item.side, population_id=item.population_id, species_id=item.species.species_id, source_model=item.source_model, particle_mass_kg=item.species.mass_kg, particle_charge_number=item.species.charge_number, source_particle_rate_s=item.particle_rate_s, invariant_total_energy_J=item.invariant_total_energy_J, mu_B0_energy_J=item.mu_B0_energy_J, full_device_cell_indices_outward=path.cell_indices, z_edges_outward_m=path.z_edges_m, path_segment_lengths_m=path.segment_lengths_m, B_tilde_edges_outward=path.B_tilde_edges, potential_drop_edges_outward_J=potential_edges_J, terminal_potential_drop_J=float(terminal_wall_potential_drop_J), full_device_cell_volumes_m3=path.cell_volumes_m3, full_device_cell_count=full_count, target_speed_grid=item.target_speed_grid, target_pitch_grid=item.target_pitch_grid, terminal_surface_id=path.surface.surface_id, terminal_surface_particle_role=path.surface.particle_role, terminal_surface_aligned_to_population_grid=path.terminal_surface_aligned_to_population_grid, terminal_cell_axial_fraction=path.terminal_cell_axial_fraction, fixed_boundary_nonterminal_fraction_tolerance=tolerance, radial_profile_model=path.radial_profile_model, source_connected_radial_rho_limit_edges=path.radial_rho_limit_edges, radial_material_surface_id_by_segment=path.radial_material_surface_id_by_segment)
    
    return mapped

def _relative_density_residual(ion_density_m3: np.ndarray, electron_density_m3: np.ndarray) -> np.ndarray:
    """Return (n_i − n_e) / max(n_i, n_e) elementwise"""
    ion = np.asarray(ion_density_m3, dtype=float)
    electron = np.asarray(electron_density_m3, dtype=float)
    scale = np.maximum(ion, electron)
  
    return np.divide(ion - electron, scale, out=np.zeros_like(scale), where=scale > 0.0)

def _center_nodal_original_edges(path: _ExpanderPath, center_potential_J: np.ndarray, throat_J: float, wall_J: float) -> np.ndarray:
    """Interpolate cell center potentials to path edges with fixed throat and wall values"""
    centers = np.asarray(center_potential_J, dtype=float)
    if centers.size != path.cell_indices.size:
        raise ValueError("center nodal expander potential must match the path cell count")
    center_z = 0.5 * (np.asarray(path.z_edges_m[:-1], dtype=float) + np.asarray(path.z_edges_m[1:], dtype=float))
    edges = np.empty(centers.size + 1, dtype=float)
    edges[0] = float(throat_J)
    edges[-1] = float(wall_J)
    for index in range(1, centers.size):
        z0 = float(center_z[index - 1])
        z1 = float(center_z[index])
        edge_z = float(path.z_edges_m[index])
        fraction = (edge_z - z0) / (z1 - z0)
        edges[index] = float(centers[index - 1] + fraction * (centers[index] - centers[index - 1]))
 
    return edges

def _prepare_center_nodal_path_geometry(path: _ExpanderPath) -> _CenterNodalPathGeometry:
    """Split each physical expander cell into two synthetic segments around its center"""
    count = path.cell_indices.size
    z_nodes = np.empty(2 * count + 1, dtype=float)
    B_nodes = np.empty_like(z_nodes)
    z_nodes[0::2] = np.asarray(path.z_edges_m, dtype=float)
    z_nodes[1::2] = 0.5 * (np.asarray(path.z_edges_m[:-1], dtype=float) + np.asarray(path.z_edges_m[1:], dtype=float))
    B_nodes[0::2] = np.asarray(path.B_tilde_edges, dtype=float)
    B_nodes[1::2] = 0.5 * (np.asarray(path.B_tilde_edges[:-1], dtype=float) + np.asarray(path.B_tilde_edges[1:], dtype=float))
    rho_nodes = np.empty_like(z_nodes)
    rho_nodes[0::2] = np.asarray(path.radial_rho_limit_edges, dtype=float)
    rho_nodes[1::2] = 0.5 * (np.asarray(path.radial_rho_limit_edges[:-1], dtype=float) + np.asarray(path.radial_rho_limit_edges[1:], dtype=float))
    synthetic_surface_ids: list[str | None] = []
    for surface_id in path.radial_material_surface_id_by_segment:
        synthetic_surface_ids.extend((surface_id, surface_id))
   
    return _CenterNodalPathGeometry(
        z_nodes_m=z_nodes,
        B_tilde_nodes=B_nodes,
        radial_rho_limit_nodes=rho_nodes,
        radial_material_surface_id_by_segment=tuple(synthetic_surface_ids),
        synthetic_cell_indices=np.arange(2 * count, dtype=int),
        parent_cell_volumes_m3=np.repeat(np.asarray(path.cell_volumes_m3, dtype=float), 2),
        segment_lengths_m=np.repeat(0.5 * np.asarray(path.segment_lengths_m, dtype=float), 2),
    )

def _center_nodal_potential_nodes(path: _ExpanderPath, path_geometry: _CenterNodalPathGeometry, center_potential_J: np.ndarray, throat_J: float, wall_J: float) -> tuple[np.ndarray, np.ndarray]:
    """Interleave fixed edge values and independent physical cell center potentials"""
    centers = np.asarray(center_potential_J, dtype=float)
    original_edges = _center_nodal_original_edges(path, centers, throat_J, wall_J)
    potential_nodes = np.empty_like(path_geometry.z_nodes_m)
    potential_nodes[0::2] = original_edges
    potential_nodes[1::2] = centers
  
    return original_edges, potential_nodes

def _center_nodal_augmented_path(path: _ExpanderPath, path_geometry: _CenterNodalPathGeometry, center_potential_J: np.ndarray, throat_J: float, wall_J: float):
    """Return synthetic path arrays used by the center nodal branch mapper"""
    original_edges, potential_nodes = _center_nodal_potential_nodes(path, path_geometry, center_potential_J, throat_J, wall_J)
   
    return original_edges, path_geometry.z_nodes_m, path_geometry.B_tilde_nodes, potential_nodes, path_geometry.synthetic_cell_indices, path_geometry.parent_cell_volumes_m3, path_geometry.segment_lengths_m

def _collapse_center_nodal_branch(branch: ExpanderBranchState, path: _ExpanderPath, original_edges_J: np.ndarray, geometry: GeometryStageResult) -> ExpanderBranchState:
    """Collapse two synthetic segment contributions back onto each physical cell"""
    actual_count = int(np.asarray(geometry.full_device_z_centers_m).size)
    synthetic_indices = np.asarray(branch.full_device_cell_indices, dtype=int)
    expected = 2 * path.cell_indices.size
    if synthetic_indices.size != expected:
        raise ValueError("center nodal synthetic branch does not contain two segments per expander cell")
    synthetic_distribution = np.asarray(branch.full_device_distribution_z_v_pitch, dtype=float)[synthetic_indices]
    synthetic_density = np.asarray(branch.full_device_density_m3, dtype=float)[synthetic_indices]
    distribution = np.zeros((actual_count,) + synthetic_distribution.shape[1:], dtype=float)
    density = np.zeros(actual_count, dtype=float)
    for parent_index, full_index in enumerate(np.asarray(path.cell_indices, dtype=int)):
        distribution[full_index] = synthetic_distribution[2 * parent_index] + synthetic_distribution[2 * parent_index + 1]
        density[full_index] = synthetic_density[2 * parent_index] + synthetic_density[2 * parent_index + 1]
   
    return ExpanderBranchState(
        side=branch.side,
        population_id=branch.population_id,
        species_id=branch.species_id,
        source_model=branch.source_model,
        terminal_surface_id=branch.terminal_surface_id,
        terminal_surface_particle_role=branch.terminal_surface_particle_role,
        terminal_surface_aligned_to_population_grid=branch.terminal_surface_aligned_to_population_grid,
        terminal_cell_axial_fraction=branch.terminal_cell_axial_fraction,
        full_device_cell_indices=np.asarray(path.cell_indices, dtype=int),
        z_edges_outward_m=np.asarray(path.z_edges_m, dtype=float),
        z_centers_outward_m=0.5 * (np.asarray(path.z_edges_m[:-1], dtype=float) + np.asarray(path.z_edges_m[1:], dtype=float)),
        path_segment_lengths_m=np.asarray(path.segment_lengths_m, dtype=float),
        B_tilde_edges_outward=np.asarray(path.B_tilde_edges, dtype=float),
        potential_drop_edges_outward_J=np.asarray(original_edges_J, dtype=float),
        terminal_potential_drop_J=float(branch.terminal_potential_drop_J),
        trajectory_status_source=np.asarray(branch.trajectory_status_source),
        local_well_topology_source=np.asarray(branch.local_well_topology_source, dtype=bool),
        full_device_distribution_z_v_pitch=distribution,
        full_device_density_m3=density,
        mapped_cell_volumes_m3=np.asarray(path.cell_volumes_m3, dtype=float),
        velocity_cell_measure_m3_s3=np.asarray(branch.velocity_cell_measure_m3_s3, dtype=float),
        throat_particle_rate_s=float(branch.throat_particle_rate_s),
        terminal_particle_rate_s=float(branch.terminal_particle_rate_s),
        reflected_particle_rate_s=float(branch.reflected_particle_rate_s),
        local_well_topology_incident_rate_s=float(branch.local_well_topology_incident_rate_s),
        unclassified_particle_rate_s=float(branch.unclassified_particle_rate_s),
        throat_midplane_kinetic_power_W=float(branch.throat_midplane_kinetic_power_W),
        terminal_midplane_kinetic_power_W=float(branch.terminal_midplane_kinetic_power_W),
        reflected_midplane_kinetic_power_W=float(branch.reflected_midplane_kinetic_power_W),
        unclassified_midplane_kinetic_power_W=float(branch.unclassified_midplane_kinetic_power_W),
        terminal_wall_power_W=float(branch.terminal_wall_power_W),
        terminal_electrostatic_power_W=float(branch.terminal_electrostatic_power_W),
        represented_inventory_particles=float(branch.represented_inventory_particles),
        mapped_inventory_particles=float(branch.mapped_inventory_particles),
        speed_domain_overflow_inventory_particles=float(branch.speed_domain_overflow_inventory_particles),
        pitch_domain_overflow_inventory_particles=float(branch.pitch_domain_overflow_inventory_particles),
        particle_balance_relative_error=float(branch.particle_balance_relative_error),
        energy_partition_relative_error=float(branch.energy_partition_relative_error),
        inventory_mapping_relative_error=float(branch.inventory_mapping_relative_error),
        finite_transit_inventory_assigned=bool(branch.finite_transit_inventory_assigned),
        fixed_boundary_applicable=bool(branch.fixed_boundary_applicable),
        failure_reason=branch.failure_reason,
        radial_profile_model=branch.radial_profile_model,
        source_connected_radial_rho_limit_edges=np.asarray(path.radial_rho_limit_edges, dtype=float),
        source_connected_radial_survival_fraction_edges=np.asarray(path.radial_survival_fraction_edges, dtype=float),
        material_particle_rate_s_by_surface=dict(branch.material_particle_rate_s_by_surface),
        material_midplane_kinetic_power_W_by_surface=dict(branch.material_midplane_kinetic_power_W_by_surface),
        material_wall_power_W_by_surface=dict(branch.material_wall_power_W_by_surface),
    )

def _map_inputs_center_nodal(*, inputs: list[_BranchInput], path: _ExpanderPath, path_geometry: _CenterNodalPathGeometry, center_potential_J: np.ndarray, throat_drop_J: float, wall_drop_J: float, geometry: GeometryStageResult, tolerance: float) -> tuple[dict[str, ExpanderBranchState], np.ndarray]:
    """Map directed sources on the center nodal potential and restore physical cells"""
    original_edges, z_nodes, B_nodes, potential_nodes, synthetic_indices, parent_volumes, segment_lengths = _center_nodal_augmented_path(path, path_geometry, center_potential_J, throat_drop_J, wall_drop_J)
    mapped: dict[str, ExpanderBranchState] = {}
    for item in inputs:
        if item.side != path.side:
            continue
        synthetic = map_source_connected_branch(side=item.side, population_id=item.population_id, species_id=item.species.species_id, source_model=item.source_model, particle_mass_kg=item.species.mass_kg, particle_charge_number=item.species.charge_number, source_particle_rate_s=item.particle_rate_s, invariant_total_energy_J=item.invariant_total_energy_J, mu_B0_energy_J=item.mu_B0_energy_J, full_device_cell_indices_outward=synthetic_indices, z_edges_outward_m=z_nodes, path_segment_lengths_m=segment_lengths, B_tilde_edges_outward=B_nodes, potential_drop_edges_outward_J=potential_nodes, terminal_potential_drop_J=float(wall_drop_J), full_device_cell_volumes_m3=parent_volumes, full_device_cell_count=int(synthetic_indices.size), target_speed_grid=item.target_speed_grid, target_pitch_grid=item.target_pitch_grid, terminal_surface_id=path.surface.surface_id, terminal_surface_particle_role=path.surface.particle_role, terminal_surface_aligned_to_population_grid=path.terminal_surface_aligned_to_population_grid, terminal_cell_axial_fraction=path.terminal_cell_axial_fraction, fixed_boundary_nonterminal_fraction_tolerance=tolerance, radial_profile_model=path.radial_profile_model, source_connected_radial_rho_limit_edges=path_geometry.radial_rho_limit_nodes, radial_material_surface_id_by_segment=path_geometry.radial_material_surface_id_by_segment)
        mapped[f"{item.population_id}:{item.side}"] = _collapse_center_nodal_branch(synthetic, path, original_edges, geometry)
   
    return mapped, original_edges

def _center_nodal_density_state(*, config: SourceModelRunConfig, geometry: GeometryStageResult, kinetic: KineticStageResult, path: _ExpanderPath, path_geometry: _CenterNodalPathGeometry, inputs: list[_BranchInput], density_plans: dict[str, SourceConnectedDensityPlan], center_potential_J: np.ndarray, throat_drop_J: float, wall_drop_J: float, build_mapped_state: bool) -> tuple[dict[str, ExpanderBranchState], np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Evaluate ion and Eq 70 electron densities for one center nodal potential trial"""
    original_edges, potential_nodes = _center_nodal_potential_nodes(path, path_geometry, center_potential_J, throat_drop_J, wall_drop_J)
    ion = np.zeros(path.cell_indices.size, dtype=float)
    mapped: dict[str, ExpanderBranchState] = {}
    if build_mapped_state:
        mapped, original_edges = _map_inputs_center_nodal(inputs=inputs, path=path, path_geometry=path_geometry, center_potential_J=center_potential_J, throat_drop_J=throat_drop_J, wall_drop_J=wall_drop_J, geometry=geometry, tolerance=float(config.kinetic_electrostatic.expander_fast_nonterminal_fraction_tolerance))
        charge_by_population = {item.population_id: float(item.species.charge_number) for item in inputs}
        for branch in mapped.values():
            ion += np.asarray(branch.full_device_density_m3, dtype=float)[path.cell_indices] * charge_by_population[branch.population_id]
    else:
        # Density only evaluation avoids constructing full branch states during least squares
        for item in inputs:
            if item.side != path.side:
                continue
            key = f"{item.population_id}:{item.side}"
            synthetic_density = source_connected_density_profile(density_plans[key], potential_nodes, float(wall_drop_J))
            ion += float(item.species.charge_number) * np.sum(synthetic_density.reshape(path.cell_indices.size, 2), axis=1)
    system = kinetic.fast_ion_system_state
    if system is None or system.shared_electrostatic_profile is None:
        raise ValueError("center nodal expander Eq 70 solve requires the shared electrostatic profile")
    n0 = float(system.shared_electrostatic_profile.electron_parent_maxwellian_n0_m3)
    temperature = _electron_temperature_J(kinetic)
    electron = n0 * np.asarray([_electron_density_fraction_eq70(float(value), float(wall_drop_J), temperature) for value in np.asarray(center_potential_J, dtype=float)], dtype=float)
    residual = _relative_density_residual(ion, electron)
   
    return mapped, ion, electron, residual, original_edges

def _side_eq70_center_nodal_solve(*, config: SourceModelRunConfig, geometry: GeometryStageResult, kinetic: KineticStageResult, path: _ExpanderPath, inputs: list[_BranchInput], throat_drop_J: float, wall_drop_J: float) -> tuple[ExpanderPotentialState, dict[str, ExpanderBranchState]]:
    """Solve one side potential from cell center Eq 70 quasineutrality residuals"""
    controls = config.kinetic_electrostatic
    temperature = _electron_temperature_J(kinetic)
    initial_edges = _initial_potential(path, throat_drop_J, wall_drop_J)
    initial_centers = 0.5 * (initial_edges[:-1] + initial_edges[1:])
    scale_J = max(float(temperature), np.finfo(float).tiny)
    upper = float(wall_drop_J) / scale_J if wall_drop_J > 0.0 else 0.0
    evaluation_count = 0
    path_geometry = _prepare_center_nodal_path_geometry(path)
    density_plans: dict[str, SourceConnectedDensityPlan] = {}
    for item in inputs:
        if item.side != path.side:
            continue
        key = f"{item.population_id}:{item.side}"
        density_plans[key] = prepare_source_connected_density_plan(side=item.side, particle_mass_kg=item.species.mass_kg, particle_charge_number=item.species.charge_number, source_particle_rate_s=item.particle_rate_s, invariant_total_energy_J=item.invariant_total_energy_J, mu_B0_energy_J=item.mu_B0_energy_J, path_segment_lengths_m=path_geometry.segment_lengths_m, B_tilde_edges_outward=path_geometry.B_tilde_nodes, cell_volumes_m3=path_geometry.parent_cell_volumes_m3, target_speed_grid=item.target_speed_grid, target_pitch_grid=item.target_pitch_grid, radial_profile_model=path.radial_profile_model, source_connected_radial_rho_limit_edges=path_geometry.radial_rho_limit_nodes)

    def evaluate_residual(values: np.ndarray) -> np.ndarray:
        """Evaluate the dimensionless Eq 70 residual for one least squares trial"""
        nonlocal evaluation_count
        centers = np.asarray(values, dtype=float) * scale_J
        state = _center_nodal_density_state(config=config, geometry=geometry, kinetic=kinetic, path=path, path_geometry=path_geometry, inputs=inputs, density_plans=density_plans, center_potential_J=centers, throat_drop_J=throat_drop_J, wall_drop_J=wall_drop_J, build_mapped_state=False)
        evaluation_count += 1
        return np.asarray(state[3], dtype=float)

    if initial_centers.size and upper > 0.0:
        x0 = np.clip(initial_centers / scale_J, 0.0, upper)
        solver_tolerance = 1.0e-12
        result = least_squares(evaluate_residual, x0, bounds=(np.zeros_like(x0), np.full_like(x0, upper)), xtol=solver_tolerance, ftol=solver_tolerance, gtol=solver_tolerance, max_nfev=max(int(controls.expander_potential_iterations), 1))
        centers = np.asarray(result.x, dtype=float) * scale_J
        solver_success = bool(result.success)
        solver_message = str(result.message)
    else:
        centers = np.zeros_like(initial_centers) if initial_centers.size else initial_centers
        solver_success = True
        solver_message = "zero_wall_barrier" if initial_centers.size else "no_expander_cells"
    state = _center_nodal_density_state(config=config, geometry=geometry, kinetic=kinetic, path=path, path_geometry=path_geometry, inputs=inputs, density_plans=density_plans, center_potential_J=centers, throat_drop_J=throat_drop_J, wall_drop_J=wall_drop_J, build_mapped_state=True)
    mapped, ion, electron, residual, original_edges = state
    error = float(np.max(np.abs(residual))) if residual.size else 0.0
    tolerance = float(controls.expander_potential_relative_tolerance)
    converged = bool(error <= tolerance)
    failure_reason = None if converged else ("expander_center_nodal_Eq70_not_converged" if solver_success else "expander_center_nodal_Eq70_solver_failed")
    state = ExpanderPotentialState(side=path.side, full_device_cell_indices=path.cell_indices, z_edges_outward_m=path.z_edges_m, potential_drop_edges_outward_J=original_edges, terminal_wall_potential_drop_J=float(wall_drop_J), wall_sheath_potential_jump_J=float(wall_drop_J) - float(original_edges[-1]), ion_charge_density_m3=ion, electron_density_m3=electron, quasineutrality_relative_error=error, iterations=max(int(evaluation_count), 1), converged=converged, failure_reason=failure_reason, potential_solver_failure_reason=None if converged or solver_success else solver_message, potential_drop_centers_outward_J=np.asarray(centers, dtype=float))
    
    return state, mapped

def _full_device_radial_support_edges(geometry: GeometryStageResult, paths: dict[str, _ExpanderPath]) -> tuple[np.ndarray, np.ndarray]:
    """Assemble radial support and survival probability on the full device edge grid"""
    if geometry.full_device_z_edges_m is None:
        raise ValueError("full device radial support requires the full device edge grid")
    full_edges = np.asarray(geometry.full_device_z_edges_m, dtype=float)
    rho = np.zeros(full_edges.size, dtype=float)
    throat_left = float(geometry.z_edges_m[0])
    throat_right = float(geometry.z_edges_m[-1])
    central = (full_edges >= throat_left) & (full_edges <= throat_right)
    rho[central] = 1.0
    scale = max(float(np.max(np.abs(full_edges))), 1.0)
    tolerance = 256.0 * np.finfo(float).eps * scale
    for path in paths.values():
        for z_value, rho_value in zip(path.z_edges_m, path.radial_rho_limit_edges, strict=True):
            index = int(np.argmin(np.abs(full_edges - float(z_value))))
            if abs(float(full_edges[index]) - float(z_value)) > tolerance:
                raise ValueError("expander radial support edge is not aligned to the full device grid")
            rho[index] = max(rho[index], float(rho_value))
    survival = np.zeros_like(rho)
    for path in paths.values():
        model = path.radial_profile_model
        break
    else:
        return rho, survival
    from source_model_revamp.radial_profiles import radial_probability_cdf
    survival = np.asarray(radial_probability_cdf(rho, model), dtype=float)
  
    return rho, survival

def _sum_surface_maps(branches: tuple[ExpanderBranchState, ...], attribute: str) -> dict[str, float]:
    """Sum one material surface diagnostic map across expander branches"""
    totals: dict[str, float] = {}
    for branch in branches:
        values = getattr(branch, attribute)
        for surface_id, value in values.items():
            totals[str(surface_id)] = totals.get(str(surface_id), 0.0) + float(value)
  
    return totals

def solve_expander_system(config: SourceModelRunConfig, geometry: GeometryStageResult, kinetic: KineticStageResult) -> ExpanderSystemState:
    """Solve both expanders from Eq 63 ions, the Eq 68 wall barrier, and Eq 70 density closure"""
    system = kinetic.fast_ion_system_state
    if system is None or system.shared_current_balance is None or system.shared_electrostatic_profile is None:
        raise ValueError("Pass 13 expander solve requires the shared Pass 11 and Eq 70 states")
    if geometry.full_device_z_centers_m is None:
        raise ValueError("Pass 13 expander solve requires the full device population grid")
    fast_inputs, prompt_rates, fast_source_audit = _fast_branch_inputs(config, kinetic)
    inputs = fast_inputs
    if not inputs:
        raise ValueError("Pass 13 expander solve has no directed throat populations")
    # Eq 68 supplies the shared wall barrier before the side Eq 70 potential solves
    wall_drop = float(system.shared_current_balance.electron_balance.barrier_energy_J)
    profile = system.shared_electrostatic_profile
    paths = {side: _path_to_selected_terminal(geometry, side, config.kinetic_electrostatic.expander_radial_profile_model) for side in ("left", "right")}
    potential_states: dict[str, ExpanderPotentialState] = {}
    branches: dict[str, ExpanderBranchState] = {}
    for side, throat_drop in (("left", float(profile.throat_potential_left_energy_J)), ("right", float(profile.throat_potential_right_energy_J))):
        potential_state, side_branches = _side_eq70_center_nodal_solve(config=config, geometry=geometry, kinetic=kinetic, path=paths[side], inputs=inputs, throat_drop_J=throat_drop, wall_drop_J=wall_drop)
        potential_states[side] = potential_state
        branches.update(side_branches)
    full_count = int(np.asarray(geometry.full_device_z_centers_m).size)
    radial_rho_full_edges, radial_survival_full_edges = _full_device_radial_support_edges(geometry, paths)
    branch_tuple = tuple(branches.values())
    material_rate_by_surface = _sum_surface_maps(branch_tuple, "material_particle_rate_s_by_surface")
    material_midplane_power_by_surface = _sum_surface_maps(branch_tuple, "material_midplane_kinetic_power_W_by_surface")
    material_wall_power_by_surface = _sum_surface_maps(branch_tuple, "material_wall_power_W_by_surface")
    material_rate_by_population_surface = {population_id: _sum_surface_maps(tuple(branch for branch in branch_tuple if branch.population_id == population_id), "material_particle_rate_s_by_surface") for population_id in sorted({branch.population_id for branch in branch_tuple})}
    distribution_by_population: dict[str, np.ndarray] = {}
    for item in inputs:
        shape = (full_count, item.target_speed_grid.centers_m_s.size, item.target_pitch_grid.centers.size)
        distribution_by_population.setdefault(item.population_id, np.zeros(shape, dtype=float))
    for branch in branches.values():
        distribution_by_population[branch.population_id] += np.asarray(branch.full_device_distribution_z_v_pitch, dtype=float)
    full_device_volumes = np.asarray(geometry.full_device_cell_volumes_m3, dtype=float)
    grid_by_population: dict[str, tuple[SpeedGrid, PitchGrid]] = {}
    for item in inputs:
        grid = grid_by_population.get(item.population_id)
        if grid is None:
            grid_by_population[item.population_id] = (item.target_speed_grid, item.target_pitch_grid)
        elif grid[0] is not item.target_speed_grid or grid[1] is not item.target_pitch_grid:
            if not np.array_equal(grid[0].faces_m_s, item.target_speed_grid.faces_m_s) or not np.array_equal(grid[1].faces, item.target_pitch_grid.faces):
                raise ValueError(f"{item.population_id} expander branches use inconsistent velocity grids")
    population_density_error: dict[str, float] = {}
    population_inventory_error: dict[str, float] = {}
    population_inventory_from_distribution: dict[str, float] = {}
    population_inventory_from_branches: dict[str, float] = {}
    for population_id, distribution in distribution_by_population.items():
        speed_grid, pitch_grid = grid_by_population[population_id]
        velocity_measure = gyrotropic_velocity_cell_volumes(speed_grid, pitch_grid)
        density_from_distribution = np.sum(distribution * velocity_measure[None, :, :], axis=(1, 2))
        density_from_branches = sum((np.asarray(branch.full_device_density_m3, dtype=float) for branch in branches.values() if branch.population_id == population_id), start=np.zeros(full_count, dtype=float))
        density_scale = max(float(np.max(density_from_distribution)) if density_from_distribution.size else 0.0, float(np.max(density_from_branches)) if density_from_branches.size else 0.0, np.finfo(float).tiny)
        population_density_error[population_id] = float(np.max(np.abs(density_from_distribution - density_from_branches)) / density_scale)
        distribution_inventory = float(np.sum(density_from_distribution * full_device_volumes))
        branch_inventory = float(sum(branch.mapped_inventory_particles for branch in branches.values() if branch.population_id == population_id))
        population_inventory_from_distribution[population_id] = distribution_inventory
        population_inventory_from_branches[population_id] = branch_inventory
        population_inventory_error[population_id] = _relative_error(distribution_inventory, branch_inventory)
    terminal_rates: dict[str, float] = {}
    terminal_midplane_power: dict[str, float] = {}
    terminal_wall_power: dict[str, float] = {}
    component_by_population = {item.population_id: item.component_id for item in inputs}
    charge_by_population = {item.population_id: float(item.species.charge_number) for item in inputs}
    for branch in branches.values():
        component = component_by_population[branch.population_id]
        terminal_rates[component] = terminal_rates.get(component, 0.0) + float(branch.terminal_particle_rate_s)
        terminal_midplane_power[component] = terminal_midplane_power.get(component, 0.0) + float(branch.terminal_midplane_kinetic_power_W)
        terminal_wall_power[component] = terminal_wall_power.get(component, 0.0) + float(branch.terminal_wall_power_W)
    for component in ('fast_deuterium', 'fast_tritium'):
        terminal_rates.setdefault(component, 0.0)
        terminal_midplane_power.setdefault(component, 0.0)
        terminal_wall_power.setdefault(component, 0.0)
    total_throat_mapped = float(sum(branch.throat_particle_rate_s for branch in branches.values()))
    total_terminal = float(sum(branch.terminal_particle_rate_s for branch in branches.values()))
    total_reflected = float(sum(branch.reflected_particle_rate_s for branch in branches.values()))
    total_unclassified = float(sum(branch.unclassified_particle_rate_s for branch in branches.values()))
    prompt_unresolved = float(sum(prompt_rates.values()))
    fast_branches = tuple(branch for branch in branches.values() if branch.population_id in {FAST_D, FAST_T})
    fast_throat = float(sum(branch.throat_particle_rate_s for branch in fast_branches))
    fast_nonterminal = float(sum(branch.reflected_particle_rate_s + branch.unclassified_particle_rate_s for branch in fast_branches))
    prompt_unresolved_fraction = prompt_unresolved / max(fast_throat + prompt_unresolved, np.finfo(float).tiny)
    prompt_population_mapped = prompt_unresolved <= np.finfo(float).tiny
    fast_fraction = fast_nonterminal / max(fast_throat, np.finfo(float).tiny)
    overflow_fraction = max((branch.speed_domain_overflow_inventory_particles + branch.pitch_domain_overflow_inventory_particles) / max(branch.represented_inventory_particles, np.finfo(float).tiny) for branch in branches.values())
    terminal_current = ELECTRON_CHARGE_C * sum(float(charge_by_population[population_id]) * float(terminal_rates[component_by_population[population_id]]) for population_id in sorted(component_by_population))
    terminal_electrostatic_power = float(sum(branch.terminal_electrostatic_power_W for branch in branches.values()))
    pass11_target = float(system.shared_current_balance.ion_current_loss.total_ion_current_A)
    current_difference = _relative_error(terminal_current, pass11_target)
    population_mapping_conservation = bool(all(value <= 1.0e-12 for value in population_density_error.values()) and all(value <= 1.0e-12 for value in population_inventory_error.values()))
    branch_conservation = bool(all(branch.particle_balance_relative_error <= 1.0e-10 and branch.energy_partition_relative_error <= 1.0e-10 and branch.inventory_mapping_relative_error <= 1.0e-8 for branch in branches.values()) and population_mapping_conservation)
    fast_tolerance = float(config.kinetic_electrostatic.expander_fast_nonterminal_fraction_tolerance)
    unclosed_tolerance = float(config.kinetic_electrostatic.expander_unclosed_population_fraction_tolerance)
    fixed_fast = bool(fast_fraction <= fast_tolerance and all(branch.fixed_boundary_applicable for branch in fast_branches))
    velocity_mapping_complete = overflow_fraction <= unclosed_tolerance
    potential_converged = bool(all(state.converged for state in potential_states.values()))
    terminal_surface_grid_aligned = bool(all(path.terminal_surface_aligned_to_population_grid for path in paths.values()))
    terminal_equivalent = bool(current_difference <= float(config.kinetic_electrostatic.total_current_balance_relative_tolerance) and prompt_population_mapped)
    population_classification_passed = bool(branch_conservation and fixed_fast and velocity_mapping_complete and terminal_surface_grid_aligned and prompt_population_mapped)
    qualified = bool(population_classification_passed and potential_converged)
    failure_reasons: list[str] = []
    if not branch_conservation:
        failure_reasons.append("expander_particle_energy_or_inventory_identity_failed")
    if not fixed_fast:
        failure_reasons.append("fast_ion_fixed_boundary_inapplicable_due_to_material_expander_return")
    if not velocity_mapping_complete:
        failure_reasons.append("expander_velocity_domain_overflow")
    if not potential_converged:
        failure_reasons.append("expander_quasineutral_potential_not_converged")
    if not prompt_population_mapped:
        failure_reasons.append("prompt_ion_expander_distribution_unavailable")
    if not terminal_surface_grid_aligned:
        failure_reasons.append("expander_terminal_surface_not_aligned_to_population_grid")
    metadata: dict[str, object] = {'expander_model': 'Egedal_Eq63_source_connected_ions_with_Eq68_current_and_Eq70_density', 'expander_fast_throat_source_model': 'Egedal_Eq63_Baldwin_TL_exact_truncated_exponential_flux_quantiles', 'expander_fast_throat_pitch_quantile_nodes': int(config.kinetic_electrostatic.expander_fast_throat_pitch_quantile_nodes), 'expander_fast_throat_source_audit_by_population': fast_source_audit, 'expander_fast_throat_source_uses_ordinary_lambda_midpoints': False, 'expander_potential_model': 'direct_coupled_Egedal_Eq70_center_nodal_quasineutrality', 'expander_potential_spatial_discretization': 'independent_cell_center_nodes_piecewise_linear_between_fixed_throat_and_wall', 'expander_initial_potential_model': 'linear_throat_to_shared_wall_barrier_initial_guess_only', 'expander_geometry_model': 'paraxial_flux_tube_with_radial_material_accessibility_and_arc_length', 'expander_radial_profile_model': str(config.kinetic_electrostatic.expander_radial_profile_model), 'expander_radial_profile_role': 'prescribed_normalized_flux_profile_for_Eq63_exhaust_material_routing_not_radially_self_consistent_FBIS', 'expander_radial_coordinate_definition': 'rho_equals_r_over_nominal_a_of_z_and_psi_tilde_equals_rho_squared', 'expander_radial_source_shape_definition': 'G_of_rho_equals_2_times_one_minus_rho_squared', 'expander_radial_profile_reference': 'Kotelnikov_et_al_Journal_of_Plasma_Physics_2025_Eq3_11_k1', 'expander_radial_profile_reference_role': 'published_pressure_profile_used_as_a_prescribed_reduced_Eq63_exhaust_weighting_closure_not_a_radially_self_consistent_FBiS_solution', 'expander_source_connected_radial_rho_limit_edges_by_side': {side: np.asarray(path.radial_rho_limit_edges, dtype=float) for side, path in paths.items()}, 'expander_source_connected_radial_survival_fraction_edges_by_side': {side: np.asarray(path.radial_survival_fraction_edges, dtype=float) for side, path in paths.items()}, 'expander_source_connected_radial_rho_limit_edges_full_device': radial_rho_full_edges, 'expander_source_connected_radial_survival_fraction_edges_full_device': radial_survival_full_edges, 'expander_radial_material_surface_id_by_segment_by_side': {side: tuple(path.radial_material_surface_id_by_segment) for side, path in paths.items()}, 'expander_material_particle_rate_s_by_surface': material_rate_by_surface, 'expander_material_midplane_kinetic_power_W_by_surface': material_midplane_power_by_surface, 'expander_material_wall_power_W_by_surface': material_wall_power_by_surface, 'expander_material_particle_rate_s_by_population_surface': material_rate_by_population_surface, 'expander_selected_final_terminal_cap_by_side': {side: path.surface.surface_id for side, path in paths.items()}, 'expander_outer_flux_tube_first_material_hit_is_whole_population_terminal': False, 'expander_source_disconnected_population_closure_available': False, 'expander_source_disconnected_population_closure_role': 'ion_reflection_local_well_and_unclassified_population_only', 'expander_fast_ion_source_connectivity_model': 'causal_first_turning_point', 'expander_source_disconnected_local_well_density_included': False, 'expander_source_disconnected_population_fraction_tolerance': unclosed_tolerance, 'expander_Eq68_current_closure_active': True, 'expander_Eq70_density_closure_active': True, 'expander_center_nodal_potential_active': True, 'expander_center_nodal_quasineutrality_error_by_side': {side: float(state.quasineutrality_relative_error) for side, state in potential_states.items()}, 'expander_center_nodal_potential_J_by_side': {side: np.asarray(state.potential_drop_centers_outward_J, dtype=float) for side, state in potential_states.items()}, 'expander_potential_converged': potential_converged, 'expander_branch_conservation_passed': branch_conservation, 'expander_population_distribution_density_identity_relative_error_by_population': population_density_error, 'expander_population_distribution_inventory_particles_by_population': population_inventory_from_distribution, 'expander_population_branch_inventory_particles_by_population': population_inventory_from_branches, 'expander_population_distribution_inventory_identity_relative_error_by_population': population_inventory_error, 'expander_population_distribution_conservation_passed': population_mapping_conservation, 'expander_population_classification_passed': population_classification_passed, 'expander_source_connected_population_complete': population_classification_passed, 'expander_fixed_fast_boundary_applicable': fixed_fast, 'expander_fast_nonterminal_fraction_tolerance': fast_tolerance, 'expander_velocity_domain_overflow_inventory_fraction_max': overflow_fraction, 'expander_velocity_mapping_complete': velocity_mapping_complete, 'expander_terminal_current_equivalent_to_Pass11_target': terminal_equivalent, 'expander_terminal_to_Pass11_current_comparison_role': 'diagnostic_initial_throat_target_only', 'expander_terminal_current_A': terminal_current, 'expander_Pass11_target_current_A': pass11_target, 'expander_terminal_to_Pass11_current_relative_difference': current_difference, 'expander_electron_source_model': 'Egedal_Eq70_truncated_Maxwellian_density_with_shared_Eq68_wall_barrier', 'expander_electron_terminal_occupancy_model': 'Eq68_sets_terminal_current_and_wall_barrier_Eq70_sets_local_density', 'expander_potential_failure_reason_by_side': {side: state.potential_solver_failure_reason for side, state in potential_states.items()}, 'expander_total_throat_particle_rate_s': total_throat_mapped + prompt_unresolved, 'expander_total_terminal_particle_rate_s': total_terminal, 'expander_total_reflected_particle_rate_s': total_reflected, 'expander_total_unclassified_particle_rate_s': total_unclassified, 'expander_prompt_unresolved_particle_rate_s': prompt_unresolved, 'expander_prompt_unresolved_fraction': prompt_unresolved_fraction, 'expander_prompt_population_mapping_available': prompt_population_mapped, 'expander_prompt_population_mapping_qualified': prompt_population_mapped, 'expander_prompt_population_mapping_failure_reason': None if prompt_population_mapped else 'prompt_loss_is_total_only_without_side_resolved_throat_distribution', 'expander_fast_nonterminal_fraction': fast_fraction, 'expander_terminal_electrostatic_power_W': terminal_electrostatic_power, 'expander_terminal_surface_by_side': {side: paths[side].surface.surface_id for side in paths}, 'expander_terminal_surface_semantics': 'selected_final_axial_cap_for_remaining_source_connected_population_material_loss_also_occurs_on_intermediate_radial_surfaces', 'expander_terminal_surface_role_by_side': {side: paths[side].surface.particle_role for side in paths}, 'expander_terminal_surface_aligned_to_population_grid': terminal_surface_grid_aligned, 'expander_terminal_surface_aligned_to_population_grid_by_side': {side: bool(paths[side].terminal_surface_aligned_to_population_grid) for side in paths}, 'expander_terminal_cell_axial_fraction_by_side': {side: float(paths[side].terminal_cell_axial_fraction) for side in paths}, 'expander_terminal_surface_alignment_required_for_population_export': True, 'expander_terminal_particle_rate_s_by_component': terminal_rates, 'expander_terminal_midplane_kinetic_power_W_by_component': terminal_midplane_power, 'expander_terminal_wall_power_W_by_component': terminal_wall_power, 'expander_failure_reasons': tuple(failure_reasons), 'expander_qualified': qualified}

    return ExpanderSystemState(branches_by_population_side=branches, potential_by_side=potential_states, full_device_distribution_by_population=distribution_by_population, terminal_particle_rate_s_by_component=terminal_rates, terminal_midplane_kinetic_power_W_by_component=terminal_midplane_power, terminal_wall_power_W_by_component=terminal_wall_power, total_throat_particle_rate_s=total_throat_mapped + prompt_unresolved, total_terminal_particle_rate_s=total_terminal, total_reflected_particle_rate_s=total_reflected, total_unclassified_particle_rate_s=total_unclassified, fast_nonterminal_fraction=fast_fraction, prompt_unresolved_particle_rate_s=prompt_unresolved, terminal_current_A=terminal_current, terminal_electrostatic_power_W=terminal_electrostatic_power, fixed_fast_boundary_applicable=fixed_fast, terminal_current_equivalent_to_pass11_target=terminal_equivalent, potential_converged=potential_converged, status='qualified' if qualified else 'diagnostic_unqualified', metadata=metadata)