"""Source connected guiding center mapping from one mirror throat through an expander"""
from __future__ import annotations
from collections.abc import Sequence
from dataclasses import dataclass
import numpy as np
from source_model_revamp.expander.types import ExpanderBranchState
from source_model_revamp.fbis.velocity_space_grid import PitchGrid, SpeedGrid, gyrotropic_velocity_cell_volumes
from source_model_revamp.radial_profiles import KOTELNIKOV_PARABOLIC_FLUX_K1, canonical_radial_profile_model, radial_probability_cdf, radial_survival_average

_TERMINAL_ROLES = {"absorbing", "collector", "end_ring"}

@dataclass(frozen=True)
class SourceConnectedDensityPlan:
    """Cache invariant source and grid data for repeated Eq 70 density evaluations"""
    side: str
    particle_mass_kg: float
    particle_charge_number: float
    flat_rate_s: np.ndarray
    flat_invariant_total_energy_J: np.ndarray
    flat_mu_B0_energy_J: np.ndarray
    path_segment_lengths_m: np.ndarray
    B_tilde_edges_outward: np.ndarray
    cell_volumes_m3: np.ndarray
    target_speed_faces_m_s: np.ndarray
    target_pitch_faces: np.ndarray
    velocity_cell_measure_m3_s3: np.ndarray
    radial_profile_model: str
    source_connected_radial_rho_limit_edges: np.ndarray
    source_connected_radial_survival_fraction_edges: np.ndarray

def _validated_radial_support(edge_count: int, radial_profile_model: str, rho_limit_edges: np.ndarray | None) -> tuple[str, np.ndarray, np.ndarray]:
    """Validate outward radial support and return the associated survival probability"""
    model = canonical_radial_profile_model(radial_profile_model)
    if rho_limit_edges is None:
        rho = np.ones(int(edge_count), dtype=float)
    else:
        rho = np.asarray(rho_limit_edges, dtype=float)
    if rho.shape != (int(edge_count),):
        raise ValueError("source connected radial rho limits must match the path edges")
    if np.any(~np.isfinite(rho)) or np.any(rho < 0.0) or np.any(rho > 1.0):
        raise ValueError("source connected radial rho limits must lie inside zero to one")
    if np.any(np.diff(rho) > 256.0 * np.finfo(float).eps):
        raise ValueError("source connected radial rho limit must not increase outward")
    if rho.size and rho[0] < 1.0 - 1.0e-10:
        raise ValueError("source connected radial support must begin with the full throat profile")
    survival = np.asarray(radial_probability_cdf(rho, model), dtype=float)
   
    return model, rho, survival

def prepare_source_connected_density_plan(*, side: str, particle_mass_kg: float, particle_charge_number: float, source_particle_rate_s: np.ndarray, invariant_total_energy_J: np.ndarray, mu_B0_energy_J: np.ndarray, path_segment_lengths_m: np.ndarray, B_tilde_edges_outward: np.ndarray, cell_volumes_m3: np.ndarray, target_speed_grid: SpeedGrid, target_pitch_grid: PitchGrid, radial_profile_model: str = KOTELNIKOV_PARABOLIC_FLUX_K1, source_connected_radial_rho_limit_edges: np.ndarray | None = None) -> SourceConnectedDensityPlan:
    """Flatten source invariants and cache grid measures for repeated density evaluations"""
    branch = str(side).strip().lower()
    if branch not in {"left", "right"}:
        raise ValueError("side must be left or right")
    mass = float(particle_mass_kg)
    charge = float(particle_charge_number)
    if not np.isfinite(mass) or mass <= 0.0 or not np.isfinite(charge) or charge <= 0.0:
        raise ValueError("particle mass and charge must be positive and finite")
    rate = np.asarray(source_particle_rate_s, dtype=float)
    invariant = np.asarray(invariant_total_energy_J, dtype=float)
    mu_B0 = np.asarray(mu_B0_energy_J, dtype=float)
    if rate.shape != invariant.shape or rate.shape != mu_B0.shape or rate.size == 0:
        raise ValueError("source rate and invariant arrays must have matching nonempty shapes")
    if np.any(~np.isfinite(rate)) or np.any(rate < 0.0) or np.any(~np.isfinite(invariant)) or np.any(~np.isfinite(mu_B0)) or np.any(mu_B0 < 0.0):
        raise ValueError("source arrays must be finite with nonnegative rate and mu B0")
    segment_lengths = np.asarray(path_segment_lengths_m, dtype=float)
    field_edges = np.asarray(B_tilde_edges_outward, dtype=float)
    volumes = np.asarray(cell_volumes_m3, dtype=float)
    if segment_lengths.ndim != 1 or field_edges.shape != (segment_lengths.size + 1,) or volumes.shape != segment_lengths.shape:
        raise ValueError("density plan path arrays have inconsistent shapes")
    if np.any(~np.isfinite(segment_lengths)) or np.any(segment_lengths <= 0.0) or np.any(~np.isfinite(field_edges)) or np.any(field_edges <= 0.0) or np.any(~np.isfinite(volumes)) or np.any(volumes <= 0.0):
        raise ValueError("density plan path arrays must be finite and positive")
    model, rho, survival = _validated_radial_support(field_edges.size, radial_profile_model, source_connected_radial_rho_limit_edges)
    
    return SourceConnectedDensityPlan(
        side=branch,
        particle_mass_kg=mass,
        particle_charge_number=charge,
        flat_rate_s=rate.reshape(-1),
        flat_invariant_total_energy_J=invariant.reshape(-1),
        flat_mu_B0_energy_J=mu_B0.reshape(-1),
        path_segment_lengths_m=segment_lengths,
        B_tilde_edges_outward=field_edges,
        cell_volumes_m3=volumes,
        target_speed_faces_m_s=np.asarray(target_speed_grid.faces_m_s, dtype=float),
        target_pitch_faces=np.asarray(target_pitch_grid.faces, dtype=float),
        velocity_cell_measure_m3_s3=gyrotropic_velocity_cell_volumes(target_speed_grid, target_pitch_grid),
        radial_profile_model=model,
        source_connected_radial_rho_limit_edges=rho,
        source_connected_radial_survival_fraction_edges=survival,
    )

def _kinematic_state(*, rate: np.ndarray, invariant: np.ndarray, mu_B0: np.ndarray, charge: float, potential_edges: np.ndarray, field_edges: np.ndarray, terminal_potential: float, rho_limit_edges: np.ndarray, radial_profile_model: str) -> dict[str, np.ndarray | float]:
    """Classify throat states using K_parallel = U + Z ΔΦ − (μ B0) B_tilde"""
    parallel_edges = invariant[:, None] + charge * potential_edges[None, :] - mu_B0[:, None] * field_edges[None, :]
    terminal_parallel = invariant + charge * terminal_potential - mu_B0 * float(field_edges[-1])
    classification_parallel = np.concatenate((parallel_edges, terminal_parallel[:, None]), axis=1)
    scale = max(float(np.max(np.abs(classification_parallel))) if classification_parallel.size else 0.0, float(np.max(np.abs(invariant))) if invariant.size else 0.0, np.finfo(float).tiny)
    tolerance = 512.0 * np.finfo(float).eps * scale
    parallel_edges = np.where(np.abs(parallel_edges) <= tolerance, 0.0, parallel_edges)
    classification_parallel = np.where(np.abs(classification_parallel) <= tolerance, 0.0, classification_parallel)
    active = rate > 0.0
    source_connected = active & np.isfinite(parallel_edges[:, 0]) & (parallel_edges[:, 0] >= -tolerance)
    # Negative parallel energy blocks causal travel outward from the throat
    blocked = classification_parallel < -tolerance
    edge_accessible_from_throat = np.logical_and.accumulate(parallel_edges >= -tolerance, axis=1)
    terminal_kinematic = source_connected & ~np.any(blocked, axis=1)
    reflected_kinematic = source_connected & ~terminal_kinematic
    unclassified = active & ~source_connected
    first_block = np.argmax(blocked, axis=1)
    has_block = np.any(blocked, axis=1)
    coordinate = np.arange(classification_parallel.shape[1])[None, :]
    # Positive parallel energy after the first block marks local well topology
    positive_after_block = np.any((coordinate > first_block[:, None]) & (classification_parallel >= -tolerance), axis=1)
    local_well_topology = reflected_kinematic & has_block & positive_after_block
    final_survival = float(radial_probability_cdf(float(rho_limit_edges[-1]), radial_profile_model))
    turn_rho = np.full(rate.shape, float(rho_limit_edges[-1]), dtype=float)
    bad_edges = parallel_edges < -tolerance
    reflected_with_bad_edge = reflected_kinematic & np.any(bad_edges, axis=1)
    for row in np.flatnonzero(reflected_with_bad_edge):
        bad_index = int(np.flatnonzero(bad_edges[row])[0])
        if bad_index <= 0:
            turn_rho[row] = float(rho_limit_edges[0])
            continue
        E0 = float(parallel_edges[row, bad_index - 1])
        E1 = float(parallel_edges[row, bad_index])
        denominator = E0 - E1
        fraction = 0.0 if denominator <= 0.0 else float(np.clip(E0 / denominator, 0.0, 1.0))
        turn_rho[row] = float(rho_limit_edges[bad_index - 1] + fraction * (rho_limit_edges[bad_index] - rho_limit_edges[bad_index - 1]))
    turn_survival = np.asarray(radial_probability_cdf(turn_rho, radial_profile_model), dtype=float)
    turn_survival[terminal_kinematic] = final_survival
    turn_survival[unclassified] = 0.0

    return {
        "parallel_edges": parallel_edges,
        "tolerance": float(tolerance),
        "source_connected": source_connected,
        "edge_accessible_from_throat": edge_accessible_from_throat,
        "terminal_kinematic": terminal_kinematic,
        "reflected_kinematic": reflected_kinematic,
        "unclassified": unclassified,
        "local_well_topology": local_well_topology,
        "turn_survival": turn_survival,
    }

def _segment_visit_fraction(E0: np.ndarray, E1: np.ndarray, starts_accessible: np.ndarray, tolerance: float) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return the fraction of a path segment reached before K_parallel crosses zero"""
    full_segment = starts_accessible & (E1 >= -tolerance)
    crossing_segment = starts_accessible & (E0 > tolerance) & (E1 < -tolerance)
    segment_fraction = np.zeros(E0.shape, dtype=float)
    segment_fraction[full_segment] = 1.0
    segment_fraction[crossing_segment] = np.clip(E0[crossing_segment] / (E0[crossing_segment] - E1[crossing_segment]), 0.0, 1.0)
   
    return full_segment, crossing_segment, segment_fraction

def source_connected_density_profile(plan: SourceConnectedDensityPlan, potential_drop_edges_outward_J: np.ndarray, terminal_potential_drop_J: float) -> np.ndarray:
    """Map radial weighted finite transit inventory to axial density for one trial potential"""
    potential_edges = np.asarray(potential_drop_edges_outward_J, dtype=float)
    if potential_edges.shape != plan.B_tilde_edges_outward.shape or np.any(~np.isfinite(potential_edges)):
        raise ValueError("potential edges must match the density plan")
    terminal_potential = float(terminal_potential_drop_J)
    if not np.isfinite(terminal_potential):
        raise ValueError("terminal potential must be finite")
    rate = plan.flat_rate_s
    invariant = plan.flat_invariant_total_energy_J
    mu_B0 = plan.flat_mu_B0_energy_J
    field_edges = plan.B_tilde_edges_outward
    state = _kinematic_state(rate=rate, invariant=invariant, mu_B0=mu_B0, charge=plan.particle_charge_number, potential_edges=potential_edges, field_edges=field_edges, terminal_potential=terminal_potential, rho_limit_edges=plan.source_connected_radial_rho_limit_edges, radial_profile_model=plan.radial_profile_model)
    parallel_edges = np.asarray(state["parallel_edges"], dtype=float)
    source_connected = np.asarray(state["source_connected"], dtype=bool)
    edge_accessible_from_throat = np.asarray(state["edge_accessible_from_throat"], dtype=bool)
    reflected = np.asarray(state["reflected_kinematic"], dtype=bool)
    turn_survival = np.asarray(state["turn_survival"], dtype=float)
    tolerance = float(state["tolerance"])
    density = np.zeros(plan.path_segment_lengths_m.size, dtype=float)
    outgoing_sign = -1.0 if plan.side == "left" else 1.0
    sqrt_2m = np.sqrt(2.0 * plan.particle_mass_kg)
    speed_faces = plan.target_speed_faces_m_s
    pitch_faces = plan.target_pitch_faces
    speed_count = speed_faces.size - 1
    pitch_count = pitch_faces.size - 1
    speed_tolerance = 64.0 * np.finfo(float).eps * max(float(speed_faces[-1]), 1.0)
    pitch_tolerance = 64.0 * np.finfo(float).eps
    rho_edges = plan.source_connected_radial_rho_limit_edges
    for local_cell in range(plan.path_segment_lengths_m.size):
        E0 = parallel_edges[:, local_cell]
        E1 = parallel_edges[:, local_cell + 1]
        starts_accessible = source_connected & edge_accessible_from_throat[:, local_cell]
        full_segment, _, segment_fraction = _segment_visit_fraction(E0, E1, starts_accessible, tolerance)
        visited = segment_fraction > 0.0
        if not np.any(visited):
            continue
        end_energy = np.where(full_segment, np.maximum(E1, 0.0), 0.0)
        start_energy = np.maximum(E0, 0.0)
        segment_length = plan.path_segment_lengths_m[local_cell] * segment_fraction
        # Linear K_parallel along a segment gives this finite transit time integral
        denominator = np.sqrt(start_energy) + np.sqrt(end_energy)
        transit_time = np.divide(sqrt_2m * segment_length, denominator, out=np.zeros_like(rate), where=visited & (denominator > 0.0))
        # Map the visited inventory to the target velocity cell at the segment midpoint
        midpoint_fraction = 0.5 * segment_fraction
        midpoint_field = field_edges[local_cell] + midpoint_fraction * (field_edges[local_cell + 1] - field_edges[local_cell])
        midpoint_potential = potential_edges[local_cell] + midpoint_fraction * (potential_edges[local_cell + 1] - potential_edges[local_cell])
        midpoint_parallel_energy = np.maximum(invariant + plan.particle_charge_number * midpoint_potential - mu_B0 * midpoint_field, 0.0)
        midpoint_perpendicular_energy = np.maximum(mu_B0 * midpoint_field, 0.0)
        speed = np.sqrt(2.0 * (midpoint_parallel_energy + midpoint_perpendicular_energy) / plan.particle_mass_kg)
        parallel_speed = np.sqrt(2.0 * midpoint_parallel_energy / plan.particle_mass_kg)
        pitch_magnitude = np.divide(parallel_speed, speed, out=np.zeros_like(speed), where=speed > 0.0)
        outward_pitch = outgoing_sign * pitch_magnitude
        speed_index = np.searchsorted(speed_faces, speed, side="right") - 1
        speed_index[np.isclose(speed, speed_faces[-1], rtol=0.0, atol=speed_tolerance)] = speed_count - 1
        outward_pitch_index = np.searchsorted(pitch_faces, outward_pitch, side="right") - 1
        outward_pitch_index[np.isclose(outward_pitch, pitch_faces[-1], rtol=0.0, atol=pitch_tolerance)] = pitch_count - 1
        speed_valid = (speed_index >= 0) & (speed_index < speed_count)
        outward_pitch_valid = (outward_pitch_index >= 0) & (outward_pitch_index < pitch_count)
        survival_average = radial_survival_average(float(rho_edges[local_cell]), float(rho_edges[local_cell + 1]), segment_fraction, plan.radial_profile_model)
        outward_inventory = rate * transit_time * survival_average
        cell_distribution = np.zeros((speed_count, pitch_count), dtype=float)
        outward_active = visited & speed_valid & outward_pitch_valid
        if np.any(outward_active):
            np.add.at(cell_distribution, (speed_index[outward_active], outward_pitch_index[outward_active]), outward_inventory[outward_active] / (float(plan.cell_volumes_m3[local_cell]) * plan.velocity_cell_measure_m3_s3[speed_index[outward_active], outward_pitch_index[outward_active]]),)
        returning_pitch = -outward_pitch
        returning_pitch_index = np.searchsorted(pitch_faces, returning_pitch, side="right") - 1
        returning_pitch_index[np.isclose(returning_pitch, pitch_faces[-1], rtol=0.0, atol=pitch_tolerance)] = pitch_count - 1
        returning_pitch_valid = (returning_pitch_index >= 0) & (returning_pitch_index < pitch_count)
        returning_inventory = rate * transit_time * turn_survival
        returning_active = visited & reflected & speed_valid & returning_pitch_valid
        if np.any(returning_active):
            np.add.at(cell_distribution, (speed_index[returning_active], returning_pitch_index[returning_active]), returning_inventory[returning_active] / (float(plan.cell_volumes_m3[local_cell]) * plan.velocity_cell_measure_m3_s3[speed_index[returning_active], returning_pitch_index[returning_active]]),)
        density[local_cell] = float(np.sum(cell_distribution * plan.velocity_cell_measure_m3_s3))
  
    return density

def _relative_error(value: float, reference: float) -> float:
    """Return the relative scalar difference using the reference magnitude"""
    return abs(float(value) - float(reference)) / max(abs(float(reference)), np.finfo(float).tiny)

def map_source_connected_branch(*, side: str, population_id: str, species_id: str, source_model: str, particle_mass_kg: float, particle_charge_number: float, source_particle_rate_s: np.ndarray, invariant_total_energy_J: np.ndarray, mu_B0_energy_J: np.ndarray, full_device_cell_indices_outward: np.ndarray, z_edges_outward_m: np.ndarray, path_segment_lengths_m: np.ndarray | None, B_tilde_edges_outward: np.ndarray, potential_drop_edges_outward_J: np.ndarray, terminal_potential_drop_J: float | None = None, full_device_cell_volumes_m3: np.ndarray, full_device_cell_count: int, target_speed_grid: SpeedGrid, target_pitch_grid: PitchGrid, terminal_surface_id: str, terminal_surface_particle_role: str, terminal_surface_aligned_to_population_grid: bool, terminal_cell_axial_fraction: float, fixed_boundary_nonterminal_fraction_tolerance: float, radial_profile_model: str = KOTELNIKOV_PARABOLIC_FLUX_K1, source_connected_radial_rho_limit_edges: np.ndarray | None = None, radial_material_surface_id_by_segment: Sequence[str | None] | None = None) -> ExpanderBranchState:
    """Map one directed throat flux with conserved U, μ B0, and radial material accessibility"""
    branch = str(side).strip().lower()
    if branch not in {"left", "right"}:
        raise ValueError("side must be left or right")
    role = str(terminal_surface_particle_role).strip().lower()
    if role not in _TERMINAL_ROLES:
        raise ValueError("terminal surface role must be absorbing, collector, or end_ring")
    mass = float(particle_mass_kg)
    charge = float(particle_charge_number)
    if not np.isfinite(mass) or mass <= 0.0:
        raise ValueError("particle_mass_kg must be positive and finite")
    if not np.isfinite(charge) or charge <= 0.0:
        raise ValueError("particle_charge_number must be positive and finite")
    rate = np.asarray(source_particle_rate_s, dtype=float)
    invariant = np.asarray(invariant_total_energy_J, dtype=float)
    mu_B0 = np.asarray(mu_B0_energy_J, dtype=float)
    if rate.shape != invariant.shape or rate.shape != mu_B0.shape or rate.size == 0:
        raise ValueError("source rate invariant energy and mu B0 arrays must have matching nonempty shapes")
    if np.any(~np.isfinite(rate)) or np.any(rate < 0.0) or np.any(~np.isfinite(invariant)) or np.any(~np.isfinite(mu_B0)) or np.any(mu_B0 < 0.0):
        raise ValueError("source arrays must be finite with nonnegative rate and mu B0")
    indices = np.asarray(full_device_cell_indices_outward, dtype=int)
    edges = np.asarray(z_edges_outward_m, dtype=float)
    field_edges = np.asarray(B_tilde_edges_outward, dtype=float)
    potential_edges = np.asarray(potential_drop_edges_outward_J, dtype=float)
    terminal_potential = float(potential_edges[-1]) if terminal_potential_drop_J is None else float(terminal_potential_drop_J)
    terminal_fraction = float(terminal_cell_axial_fraction)
    volumes = np.asarray(full_device_cell_volumes_m3, dtype=float)
    if indices.ndim != 1 or indices.size == 0 or edges.size != indices.size + 1 or field_edges.shape != edges.shape or potential_edges.shape != edges.shape or volumes.shape != indices.shape:
        raise ValueError("expander path arrays have inconsistent shapes")
    if np.any(indices < 0) or np.any(indices >= int(full_device_cell_count)) or np.any(~np.isfinite(edges)) or np.any(~np.isfinite(field_edges)) or np.any(field_edges <= 0.0) or np.any(~np.isfinite(potential_edges)) or not np.isfinite(terminal_potential) or np.any(~np.isfinite(volumes)) or np.any(volumes <= 0.0):
        raise ValueError("expander path arrays and terminal potential contain invalid values")
    if not np.isfinite(terminal_fraction) or terminal_fraction <= 0.0 or terminal_fraction > 1.0:
        raise ValueError("terminal_cell_axial_fraction must be finite in the interval zero to one")
    axial_width = np.abs(np.diff(edges))
    segment_lengths = axial_width if path_segment_lengths_m is None else np.asarray(path_segment_lengths_m, dtype=float)
    if segment_lengths.shape != indices.shape or np.any(~np.isfinite(segment_lengths)) or np.any(segment_lengths <= 0.0):
        raise ValueError("expander path segment lengths must be finite and positive")
    if np.any(axial_width <= 0.0):
        raise ValueError("expander path edges must be strictly ordered outward")
    model, rho_edges, survival_edges = _validated_radial_support(edges.size, radial_profile_model, source_connected_radial_rho_limit_edges)
    if radial_material_surface_id_by_segment is None:
        material_surface_ids = tuple(None for _ in indices)
    else:
        material_surface_ids = tuple(None if value is None else str(value) for value in radial_material_surface_id_by_segment)
    if len(material_surface_ids) != indices.size:
        raise ValueError("radial material surface IDs must match the expander path segments")
    for index, surface_id in enumerate(material_surface_ids):
        radial_loss = float(survival_edges[index] - survival_edges[index + 1])
        if radial_loss > 256.0 * np.finfo(float).eps and surface_id is None:
            raise ValueError("radial material loss must be assigned to a material surface")

    source_shape = rate.shape
    flat_rate = rate.reshape(-1)
    flat_invariant = invariant.reshape(-1)
    flat_mu_B0 = mu_B0.reshape(-1)
    state = _kinematic_state(rate=flat_rate, invariant=flat_invariant, mu_B0=flat_mu_B0, charge=charge, potential_edges=potential_edges, field_edges=field_edges, terminal_potential=terminal_potential, rho_limit_edges=rho_edges, radial_profile_model=model)
    parallel_edges = np.asarray(state["parallel_edges"], dtype=float)
    source_connected = np.asarray(state["source_connected"], dtype=bool)
    edge_accessible_from_throat = np.asarray(state["edge_accessible_from_throat"], dtype=bool)
    terminal_kinematic = np.asarray(state["terminal_kinematic"], dtype=bool)
    reflected_kinematic = np.asarray(state["reflected_kinematic"], dtype=bool)
    unclassified = np.asarray(state["unclassified"], dtype=bool)
    local_well_topology = np.asarray(state["local_well_topology"], dtype=bool)
    turn_survival = np.asarray(state["turn_survival"], dtype=float)
    tolerance = float(state["tolerance"])
    status = np.full(flat_rate.shape, "inactive", dtype="U16")
    status[terminal_kinematic] = "terminal"
    status[reflected_kinematic] = "reflected"
    status[unclassified] = "unclassified"
    material_rate_by_surface: dict[str, float] = {}
    material_midplane_power_by_surface: dict[str, float] = {}
    material_wall_power_by_surface: dict[str, float] = {}

    def add_material_loss(surface_id: str, fraction: np.ndarray) -> None:
        """Accumulate rate and power deposited on one material surface"""
        active_fraction = np.asarray(fraction, dtype=float)
        if np.any(active_fraction < -256.0 * np.finfo(float).eps):
            raise ValueError("material loss fraction must be nonnegative")
        active_fraction = np.clip(active_fraction, 0.0, 1.0)
        local_rate = float(np.sum(flat_rate * active_fraction))
        local_midplane = float(np.sum(flat_rate * active_fraction * flat_invariant))
        wall_energy = flat_invariant + charge * terminal_potential
        local_wall = float(np.sum(flat_rate * active_fraction * wall_energy))
        material_rate_by_surface[surface_id] = material_rate_by_surface.get(surface_id, 0.0) + local_rate
        material_midplane_power_by_surface[surface_id] = material_midplane_power_by_surface.get(surface_id, 0.0) + local_midplane
        material_wall_power_by_surface[surface_id] = material_wall_power_by_surface.get(surface_id, 0.0) + local_wall

    full_shape = (int(full_device_cell_count), target_speed_grid.centers_m_s.size, target_pitch_grid.centers.size)
    distribution = np.zeros(full_shape, dtype=float)
    represented_inventory = 0.0
    mapped_inventory = 0.0
    speed_overflow_inventory = 0.0
    pitch_overflow_inventory = 0.0
    velocity_measure = gyrotropic_velocity_cell_volumes(target_speed_grid, target_pitch_grid)
    outgoing_sign = -1.0 if branch == "left" else 1.0
    sqrt_2m = np.sqrt(2.0 * mass)
    speed_tolerance = 64.0 * np.finfo(float).eps * max(float(target_speed_grid.faces_m_s[-1]), 1.0)
    pitch_tolerance = 64.0 * np.finfo(float).eps

    for local_cell, global_cell in enumerate(indices):
        E0 = parallel_edges[:, local_cell]
        E1 = parallel_edges[:, local_cell + 1]
        starts_accessible = source_connected & edge_accessible_from_throat[:, local_cell]
        full_segment, _, segment_fraction = _segment_visit_fraction(E0, E1, starts_accessible, tolerance)
        visited = segment_fraction > 0.0
        if not np.any(visited):
            continue
        rho_start = float(rho_edges[local_cell])
        rho_end = float(rho_edges[local_cell + 1])
        rho_stop = rho_start + segment_fraction * (rho_end - rho_start)
        survival_start = float(radial_probability_cdf(rho_start, model))
        survival_stop = np.asarray(radial_probability_cdf(rho_stop, model), dtype=float)
        radial_loss_fraction = np.where(visited, np.maximum(survival_start - survival_stop, 0.0), 0.0)
        if np.any(radial_loss_fraction > 256.0 * np.finfo(float).eps):
            surface_id = material_surface_ids[local_cell]
            if surface_id is None:
                raise ValueError("visited radial material loss does not have a material surface")
            add_material_loss(surface_id, radial_loss_fraction)

        end_energy = np.where(full_segment, np.maximum(E1, 0.0), 0.0)
        start_energy = np.maximum(E0, 0.0)
        segment_length = segment_lengths[local_cell] * segment_fraction
        # Linear K_parallel along a segment gives this finite transit time integral
        denominator = np.sqrt(start_energy) + np.sqrt(end_energy)
        transit_time = np.divide(sqrt_2m * segment_length, denominator, out=np.zeros_like(flat_rate), where=visited & (denominator > 0.0))
        # Radial survival weights the source connected inventory that remains inside material limits
        survival_average = radial_survival_average(rho_start, rho_end, segment_fraction, model)
        outward_inventory = flat_rate * transit_time * survival_average
        returning_inventory = flat_rate * transit_time * turn_survival
        outward_inventory[~visited] = 0.0
        returning_inventory[~(visited & reflected_kinematic)] = 0.0
        total_inventory = outward_inventory + returning_inventory
        represented_inventory += float(np.sum(total_inventory))

        midpoint_fraction = 0.5 * segment_fraction
        midpoint_field = field_edges[local_cell] + midpoint_fraction * (field_edges[local_cell + 1] - field_edges[local_cell])
        midpoint_potential = potential_edges[local_cell] + midpoint_fraction * (potential_edges[local_cell + 1] - potential_edges[local_cell])
        midpoint_parallel_energy = np.maximum(flat_invariant + charge * midpoint_potential - flat_mu_B0 * midpoint_field, 0.0)
        midpoint_perpendicular_energy = np.maximum(flat_mu_B0 * midpoint_field, 0.0)
        speed = np.sqrt(2.0 * (midpoint_parallel_energy + midpoint_perpendicular_energy) / mass)
        parallel_speed = np.sqrt(2.0 * midpoint_parallel_energy / mass)
        pitch_magnitude = np.divide(parallel_speed, speed, out=np.zeros_like(speed), where=speed > 0.0)
        outward_pitch = outgoing_sign * pitch_magnitude
        speed_index = np.searchsorted(target_speed_grid.faces_m_s, speed, side="right") - 1
        speed_index[np.isclose(speed, target_speed_grid.faces_m_s[-1], rtol=0.0, atol=speed_tolerance)] = target_speed_grid.centers_m_s.size - 1
        outward_pitch_index = np.searchsorted(target_pitch_grid.faces, outward_pitch, side="right") - 1
        outward_pitch_index[np.isclose(outward_pitch, target_pitch_grid.faces[-1], rtol=0.0, atol=pitch_tolerance)] = target_pitch_grid.centers.size - 1
        speed_valid = (speed_index >= 0) & (speed_index < target_speed_grid.centers_m_s.size)
        outward_pitch_valid = (outward_pitch_index >= 0) & (outward_pitch_index < target_pitch_grid.centers.size)
        outward_active = visited & speed_valid & outward_pitch_valid & (outward_inventory > 0.0)
        if np.any(outward_active):
            np.add.at(distribution[global_cell], (speed_index[outward_active], outward_pitch_index[outward_active]), outward_inventory[outward_active] / (float(volumes[local_cell]) * velocity_measure[speed_index[outward_active], outward_pitch_index[outward_active]]),)
            mapped_inventory += float(np.sum(outward_inventory[outward_active]))
        speed_overflow_inventory += float(np.sum(outward_inventory[visited & ~speed_valid]))
        pitch_overflow_inventory += float(np.sum(outward_inventory[visited & speed_valid & ~outward_pitch_valid]))

        returning = visited & reflected_kinematic & (returning_inventory > 0.0)
        if np.any(returning):
            returning_pitch = -outward_pitch
            returning_pitch_index = np.searchsorted(target_pitch_grid.faces, returning_pitch, side="right") - 1
            returning_pitch_index[np.isclose(returning_pitch, target_pitch_grid.faces[-1], rtol=0.0, atol=pitch_tolerance)] = target_pitch_grid.centers.size - 1
            returning_pitch_valid = (returning_pitch_index >= 0) & (returning_pitch_index < target_pitch_grid.centers.size)
            returning_active = returning & speed_valid & returning_pitch_valid
            if np.any(returning_active):
                np.add.at(distribution[global_cell], (speed_index[returning_active], returning_pitch_index[returning_active]), returning_inventory[returning_active] / (float(volumes[local_cell]) * velocity_measure[speed_index[returning_active], returning_pitch_index[returning_active]]),)
                mapped_inventory += float(np.sum(returning_inventory[returning_active]))
            speed_overflow_inventory += float(np.sum(returning_inventory[returning & ~speed_valid]))
            pitch_overflow_inventory += float(np.sum(returning_inventory[returning & speed_valid & ~returning_pitch_valid]))

    final_survival = float(survival_edges[-1])
    terminal_cap_fraction = np.where(terminal_kinematic, final_survival, 0.0)
    if np.any(terminal_cap_fraction > 0.0):
        add_material_loss(str(terminal_surface_id), terminal_cap_fraction)

    terminal_rate = float(sum(material_rate_by_surface.values()))
    terminal_midplane_power = float(sum(material_midplane_power_by_surface.values()))
    terminal_wall_power = float(sum(material_wall_power_by_surface.values()))
    reflected_rate = float(np.sum(flat_rate * turn_survival * reflected_kinematic))
    unclassified_rate = float(np.sum(flat_rate[unclassified]))
    local_well_rate = float(np.sum(flat_rate * turn_survival * local_well_topology))
    throat_rate = float(np.sum(flat_rate))
    throat_power = float(np.sum(flat_rate * flat_invariant))
    reflected_midplane_power = float(np.sum(flat_rate * flat_invariant * turn_survival * reflected_kinematic))
    unclassified_midplane_power = float(np.sum(flat_rate[unclassified] * flat_invariant[unclassified]))
    terminal_electrostatic_power = terminal_wall_power - terminal_midplane_power

    density = np.sum(distribution * velocity_measure[None, :, :], axis=(1, 2))
    integrated_inventory = float(np.sum(density[indices] * volumes))
    partition_rate = terminal_rate + reflected_rate + unclassified_rate
    partition_power = terminal_midplane_power + reflected_midplane_power + unclassified_midplane_power
    nonterminal_fraction = (reflected_rate + unclassified_rate) / max(throat_rate, np.finfo(float).tiny)
    overflow_inventory = speed_overflow_inventory + pitch_overflow_inventory
    boundary_failure_reasons: list[str] = []
    mapping_failure_reasons: list[str] = []
    if nonterminal_fraction > float(fixed_boundary_nonterminal_fraction_tolerance):
        boundary_failure_reasons.append("material_nonterminal_fraction")
    if overflow_inventory > 1.0e-10 * max(represented_inventory, np.finfo(float).tiny):
        mapping_failure_reasons.append("velocity_domain_overflow")
    if not bool(terminal_surface_aligned_to_population_grid):
        mapping_failure_reasons.append("terminal_surface_not_aligned_to_population_grid")
    if local_well_rate > float(fixed_boundary_nonterminal_fraction_tolerance) * max(throat_rate, np.finfo(float).tiny):
        boundary_failure_reasons.append("source_disconnected_local_well_topology")
    failure_reasons = boundary_failure_reasons + mapping_failure_reasons

    return ExpanderBranchState(
        side=branch,
        population_id=str(population_id),
        species_id=str(species_id),
        source_model=str(source_model),
        terminal_surface_id=str(terminal_surface_id),
        terminal_surface_particle_role=role,
        terminal_surface_aligned_to_population_grid=bool(terminal_surface_aligned_to_population_grid),
        terminal_cell_axial_fraction=terminal_fraction,
        full_device_cell_indices=indices,
        z_edges_outward_m=edges,
        z_centers_outward_m=0.5 * (edges[:-1] + edges[1:]),
        path_segment_lengths_m=segment_lengths,
        B_tilde_edges_outward=field_edges,
        potential_drop_edges_outward_J=potential_edges,
        terminal_potential_drop_J=terminal_potential,
        trajectory_status_source=status.reshape(source_shape),
        local_well_topology_source=local_well_topology.reshape(source_shape),
        full_device_distribution_z_v_pitch=distribution,
        full_device_density_m3=density,
        mapped_cell_volumes_m3=volumes,
        velocity_cell_measure_m3_s3=velocity_measure,
        throat_particle_rate_s=throat_rate,
        terminal_particle_rate_s=terminal_rate,
        reflected_particle_rate_s=reflected_rate,
        local_well_topology_incident_rate_s=local_well_rate,
        unclassified_particle_rate_s=unclassified_rate,
        throat_midplane_kinetic_power_W=throat_power,
        terminal_midplane_kinetic_power_W=terminal_midplane_power,
        reflected_midplane_kinetic_power_W=reflected_midplane_power,
        unclassified_midplane_kinetic_power_W=unclassified_midplane_power,
        terminal_wall_power_W=terminal_wall_power,
        terminal_electrostatic_power_W=terminal_electrostatic_power,
        represented_inventory_particles=represented_inventory,
        mapped_inventory_particles=integrated_inventory,
        speed_domain_overflow_inventory_particles=speed_overflow_inventory,
        pitch_domain_overflow_inventory_particles=pitch_overflow_inventory,
        particle_balance_relative_error=_relative_error(partition_rate, throat_rate),
        energy_partition_relative_error=_relative_error(partition_power, throat_power),
        inventory_mapping_relative_error=_relative_error(integrated_inventory + overflow_inventory, represented_inventory),
        finite_transit_inventory_assigned=True,
        fixed_boundary_applicable=not boundary_failure_reasons,
        failure_reason=None if not failure_reasons else ";".join(failure_reasons),
        radial_profile_model=model,
        source_connected_radial_rho_limit_edges=rho_edges,
        source_connected_radial_survival_fraction_edges=survival_edges,
        material_particle_rate_s_by_surface=material_rate_by_surface,
        material_midplane_kinetic_power_W_by_surface=material_midplane_power_by_surface,
        material_wall_power_W_by_surface=material_wall_power_by_surface,
    )