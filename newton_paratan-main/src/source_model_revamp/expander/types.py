"""Validated state containers for source connected expander populations"""
from __future__ import annotations
from collections.abc import Mapping
from dataclasses import dataclass, field
import numpy as np
from source_model_revamp.constants import ELECTRON_CHARGE_C
from source_model_revamp.radial_profiles import KOTELNIKOV_PARABOLIC_FLUX_K1, canonical_radial_profile_model, radial_probability_cdf

def _relative_error(value: float, reference: float) -> float:
    """Return the relative scalar difference using the reference magnitude"""
    return abs(float(value) - float(reference)) / max(abs(float(reference)), np.finfo(float).tiny)

def _require_close(name: str, value: float, reference: float, tolerance: float) -> None:
    """Require two scalar diagnostics to agree within a relative tolerance"""
    if _relative_error(value, reference) > float(tolerance):
        raise ValueError(f"{name} is inconsistent")

@dataclass(frozen=True)
class ExpanderBranchState:
    """Store one directed throat population mapped through one expander

    Path arrays run outward from the throat to the terminal surface
    The full device distribution has shape (n_z, n_speed, n_pitch)
    """
    side: str
    population_id: str
    species_id: str
    source_model: str
    terminal_surface_id: str
    terminal_surface_particle_role: str
    terminal_surface_aligned_to_population_grid: bool
    terminal_cell_axial_fraction: float
    full_device_cell_indices: np.ndarray
    z_edges_outward_m: np.ndarray
    z_centers_outward_m: np.ndarray
    path_segment_lengths_m: np.ndarray
    B_tilde_edges_outward: np.ndarray
    potential_drop_edges_outward_J: np.ndarray
    terminal_potential_drop_J: float
    trajectory_status_source: np.ndarray
    local_well_topology_source: np.ndarray
    full_device_distribution_z_v_pitch: np.ndarray
    full_device_density_m3: np.ndarray
    mapped_cell_volumes_m3: np.ndarray
    velocity_cell_measure_m3_s3: np.ndarray
    throat_particle_rate_s: float
    terminal_particle_rate_s: float
    reflected_particle_rate_s: float
    local_well_topology_incident_rate_s: float
    unclassified_particle_rate_s: float
    throat_midplane_kinetic_power_W: float
    terminal_midplane_kinetic_power_W: float
    reflected_midplane_kinetic_power_W: float
    unclassified_midplane_kinetic_power_W: float
    terminal_wall_power_W: float
    terminal_electrostatic_power_W: float
    represented_inventory_particles: float
    mapped_inventory_particles: float
    speed_domain_overflow_inventory_particles: float
    pitch_domain_overflow_inventory_particles: float
    particle_balance_relative_error: float
    energy_partition_relative_error: float
    inventory_mapping_relative_error: float
    finite_transit_inventory_assigned: bool
    fixed_boundary_applicable: bool
    failure_reason: str | None
    radial_profile_model: str = KOTELNIKOV_PARABOLIC_FLUX_K1
    source_connected_radial_rho_limit_edges: np.ndarray = field(default_factory=lambda: np.asarray([], dtype=float))
    source_connected_radial_survival_fraction_edges: np.ndarray = field(default_factory=lambda: np.asarray([], dtype=float))
    material_particle_rate_s_by_surface: Mapping[str, float] = field(default_factory=dict)
    material_midplane_kinetic_power_W_by_surface: Mapping[str, float] = field(default_factory=dict)
    material_wall_power_W_by_surface: Mapping[str, float] = field(default_factory=dict)

    def __post_init__(self) -> None:
        """Validate branch geometry, distributions, rate partitions, power partitions, and inventories"""
        side = str(self.side).strip().lower()
        if side not in {"left", "right"}:
            raise ValueError("expander branch side must be left or right")
        indices = np.asarray(self.full_device_cell_indices, dtype=int)
        edges = np.asarray(self.z_edges_outward_m, dtype=float)
        centers = np.asarray(self.z_centers_outward_m, dtype=float)
        lengths = np.asarray(self.path_segment_lengths_m, dtype=float)
        field = np.asarray(self.B_tilde_edges_outward, dtype=float)
        potential = np.asarray(self.potential_drop_edges_outward_J, dtype=float)
        terminal_potential = float(self.terminal_potential_drop_J)
        terminal_fraction = float(self.terminal_cell_axial_fraction)
        status = np.asarray(self.trajectory_status_source)
        well = np.asarray(self.local_well_topology_source, dtype=bool)
        distribution = np.asarray(self.full_device_distribution_z_v_pitch, dtype=float)
        density = np.asarray(self.full_device_density_m3, dtype=float)
        volumes = np.asarray(self.mapped_cell_volumes_m3, dtype=float)
        velocity_measure = np.asarray(self.velocity_cell_measure_m3_s3, dtype=float)
        radial_model = canonical_radial_profile_model(self.radial_profile_model)
        rho = np.asarray(self.source_connected_radial_rho_limit_edges, dtype=float)
        survival = np.asarray(self.source_connected_radial_survival_fraction_edges, dtype=float)
        if rho.size == 0:
            rho = np.ones(edges.shape, dtype=float)
        if survival.size == 0:
            survival = np.asarray(radial_probability_cdf(rho, radial_model), dtype=float)
        material_rate = {str(key): float(value) for key, value in self.material_particle_rate_s_by_surface.items()}
        material_midplane = {str(key): float(value) for key, value in self.material_midplane_kinetic_power_W_by_surface.items()}
        material_wall = {str(key): float(value) for key, value in self.material_wall_power_W_by_surface.items()}
        if indices.ndim != 1 or edges.size != indices.size + 1 or centers.shape != indices.shape or lengths.shape != indices.shape:
            raise ValueError("expander branch path arrays have inconsistent shapes")
        if field.shape != edges.shape or potential.shape != edges.shape:
            raise ValueError("expander branch field and potential must match the path edges")
        if rho.shape != edges.shape or survival.shape != edges.shape:
            raise ValueError("expander branch radial support must match the path edges")
        if np.any(~np.isfinite(rho)) or np.any(rho < 0.0) or np.any(rho > 1.0) or np.any(np.diff(rho) > 256.0 * np.finfo(float).eps):
            raise ValueError("expander branch radial rho support must be finite bounded and nonincreasing")
        expected_survival = np.asarray(radial_probability_cdf(rho, radial_model), dtype=float)
        if not np.allclose(survival, expected_survival, rtol=1.0e-13, atol=1.0e-15):
            raise ValueError("expander branch radial survival does not match the radial profile")
        if set(material_rate) != set(material_midplane) or set(material_rate) != set(material_wall):
            raise ValueError("expander branch material surface maps must have identical keys")
        if any(not np.isfinite(value) or value < 0.0 for mapping in (material_rate, material_midplane, material_wall) for value in mapping.values()):
            raise ValueError("expander branch material surface values must be finite and nonnegative")
        if np.any(~np.isfinite(edges)) or np.any(~np.isfinite(centers)) or np.any(~np.isfinite(lengths)) or np.any(lengths <= 0.0) or np.any(~np.isfinite(field)) or np.any(field <= 0.0) or np.any(~np.isfinite(potential)) or not np.isfinite(terminal_potential):
            raise ValueError("expander branch path arrays and terminal potential must be finite with positive lengths and field")
        if not np.isfinite(terminal_fraction) or terminal_fraction <= 0.0 or terminal_fraction > 1.0:
            raise ValueError("expander terminal cell axial fraction must be finite in the interval zero to one")
        if status.shape != well.shape:
            raise ValueError("expander branch source classification arrays must match")
        if distribution.ndim != 3 or density.ndim != 1 or distribution.shape[0] != density.size:
            raise ValueError("expander branch full device distribution and density have inconsistent shapes")
        if volumes.shape != indices.shape or velocity_measure.shape != distribution.shape[1:]:
            raise ValueError("expander branch cell volumes and velocity measure have inconsistent shapes")
        if np.any(~np.isfinite(distribution)) or np.any(distribution < 0.0) or np.any(~np.isfinite(density)) or np.any(density < 0.0):
            raise ValueError("expander branch distributions must be finite and nonnegative")
        if np.any(~np.isfinite(volumes)) or np.any(volumes <= 0.0) or np.any(~np.isfinite(velocity_measure)) or np.any(velocity_measure <= 0.0):
            raise ValueError("expander branch cell volumes and velocity measure must be finite and positive")
        if np.any(indices < 0) or np.any(indices >= density.size) or np.unique(indices).size != indices.size:
            raise ValueError("expander branch cell indices must be unique and inside the full device distribution")
        nonnegative = (
            self.throat_particle_rate_s,
            self.terminal_particle_rate_s,
            self.reflected_particle_rate_s,
            self.local_well_topology_incident_rate_s,
            self.unclassified_particle_rate_s,
            self.throat_midplane_kinetic_power_W,
            self.terminal_midplane_kinetic_power_W,
            self.reflected_midplane_kinetic_power_W,
            self.unclassified_midplane_kinetic_power_W,
            self.terminal_wall_power_W,
            self.terminal_electrostatic_power_W,
            self.represented_inventory_particles,
            self.mapped_inventory_particles,
            self.speed_domain_overflow_inventory_particles,
            self.pitch_domain_overflow_inventory_particles,
            self.particle_balance_relative_error,
            self.energy_partition_relative_error,
            self.inventory_mapping_relative_error,
        )
        if any(not np.isfinite(float(value)) or float(value) < 0.0 for value in nonnegative):
            raise ValueError("expander branch rates powers inventories and errors must be finite and nonnegative")
        particle_partition = float(self.terminal_particle_rate_s) + float(self.reflected_particle_rate_s) + float(self.unclassified_particle_rate_s)
        energy_partition = float(self.terminal_midplane_kinetic_power_W) + float(self.reflected_midplane_kinetic_power_W) + float(self.unclassified_midplane_kinetic_power_W)
        inventory_partition = float(self.mapped_inventory_particles) + float(self.speed_domain_overflow_inventory_particles) + float(self.pitch_domain_overflow_inventory_particles)
        density_from_distribution = np.sum(distribution * velocity_measure[None, :, :], axis=(1, 2))
        density_scale = max(float(np.max(density)) if density.size else 0.0, float(np.max(density_from_distribution)) if density_from_distribution.size else 0.0, 1.0)
        if not np.allclose(density, density_from_distribution, rtol=1.0e-12, atol=2048.0 * np.finfo(float).eps * density_scale):
            raise ValueError("expander branch density does not equal the integrated distribution")
        mapped_inventory_from_density = float(np.sum(density[indices] * volumes))
        if material_rate:
            _require_close("expander branch material particle rate", sum(material_rate.values()), self.terminal_particle_rate_s, 1.0e-10)
            _require_close("expander branch material midplane power", sum(material_midplane.values()), self.terminal_midplane_kinetic_power_W, 1.0e-10)
            _require_close("expander branch material wall power", sum(material_wall.values()), self.terminal_wall_power_W, 1.0e-10)
        _require_close("expander branch particle partition", particle_partition, self.throat_particle_rate_s, 1.0e-10)
        _require_close("expander branch energy partition", energy_partition, self.throat_midplane_kinetic_power_W, 1.0e-10)
        _require_close("expander branch inventory partition", inventory_partition, self.represented_inventory_particles, 1.0e-8)
        _require_close("expander branch mapped inventory", mapped_inventory_from_density, self.mapped_inventory_particles, 1.0e-12)
        _require_close("expander branch terminal wall power", self.terminal_midplane_kinetic_power_W + self.terminal_electrostatic_power_W, self.terminal_wall_power_W, 1.0e-10)
        _require_close("expander branch particle error", self.particle_balance_relative_error, _relative_error(particle_partition, self.throat_particle_rate_s), 1.0e-10)
        _require_close("expander branch energy error", self.energy_partition_relative_error, _relative_error(energy_partition, self.throat_midplane_kinetic_power_W), 1.0e-10)
        _require_close("expander branch inventory error", self.inventory_mapping_relative_error, _relative_error(inventory_partition, self.represented_inventory_particles), 1.0e-8)
        if float(self.local_well_topology_incident_rate_s) > float(self.reflected_particle_rate_s) + 1.0e-10 * max(float(self.throat_particle_rate_s), 1.0):
            raise ValueError("expander branch local well rate cannot exceed the reflected rate")
        if not self.finite_transit_inventory_assigned and float(self.represented_inventory_particles) > 0.0:
            raise ValueError("expander branch represented inventory requires finite transit assignment")
        object.__setattr__(self, "side", side)
        object.__setattr__(self, "full_device_cell_indices", indices)
        object.__setattr__(self, "z_edges_outward_m", edges)
        object.__setattr__(self, "z_centers_outward_m", centers)
        object.__setattr__(self, "path_segment_lengths_m", lengths)
        object.__setattr__(self, "B_tilde_edges_outward", field)
        object.__setattr__(self, "potential_drop_edges_outward_J", potential)
        object.__setattr__(self, "terminal_cell_axial_fraction", terminal_fraction)
        object.__setattr__(self, "trajectory_status_source", status)
        object.__setattr__(self, "local_well_topology_source", well)
        object.__setattr__(self, "full_device_distribution_z_v_pitch", distribution)
        object.__setattr__(self, "full_device_density_m3", density)
        object.__setattr__(self, "mapped_cell_volumes_m3", volumes)
        object.__setattr__(self, "velocity_cell_measure_m3_s3", velocity_measure)
        object.__setattr__(self, "radial_profile_model", radial_model)
        object.__setattr__(self, "source_connected_radial_rho_limit_edges", rho)
        object.__setattr__(self, "source_connected_radial_survival_fraction_edges", survival)
        object.__setattr__(self, "material_particle_rate_s_by_surface", material_rate)
        object.__setattr__(self, "material_midplane_kinetic_power_W_by_surface", material_midplane)
        object.__setattr__(self, "material_wall_power_W_by_surface", material_wall)

@dataclass(frozen=True)
class ExpanderPotentialState:
    """Store one side Eq 70 expander potential and quasineutral density state"""
    side: str
    full_device_cell_indices: np.ndarray
    z_edges_outward_m: np.ndarray
    potential_drop_edges_outward_J: np.ndarray
    terminal_wall_potential_drop_J: float
    wall_sheath_potential_jump_J: float
    ion_charge_density_m3: np.ndarray
    electron_density_m3: np.ndarray
    quasineutrality_relative_error: float
    iterations: int
    converged: bool
    failure_reason: str | None
    potential_solver_failure_reason: str | None = None
    potential_drop_centers_outward_J: np.ndarray | None = None

    def __post_init__(self) -> None:
        """Validate the side potential grid, density closure, sheath jump, and convergence state"""
        side = str(self.side).strip().lower()
        if side not in {"left", "right"}:
            raise ValueError("expander potential side must be left or right")
        indices = np.asarray(self.full_device_cell_indices, dtype=int)
        edges = np.asarray(self.z_edges_outward_m, dtype=float)
        potential = np.asarray(self.potential_drop_edges_outward_J, dtype=float)
        terminal_wall_potential = float(self.terminal_wall_potential_drop_J)
        sheath_jump = float(self.wall_sheath_potential_jump_J)
        ion = np.asarray(self.ion_charge_density_m3, dtype=float)
        electron = np.asarray(self.electron_density_m3, dtype=float)
        potential_centers = 0.5 * (potential[:-1] + potential[1:]) if self.potential_drop_centers_outward_J is None else np.asarray(self.potential_drop_centers_outward_J, dtype=float)
        if indices.ndim != 1 or edges.size != indices.size + 1 or potential.shape != edges.shape or potential_centers.shape != indices.shape or ion.shape != indices.shape or electron.shape != indices.shape:
            raise ValueError("expander potential arrays have inconsistent shapes")
        if np.any(~np.isfinite(edges)) or np.any(~np.isfinite(potential)) or np.any(~np.isfinite(potential_centers)) or not np.isfinite(terminal_wall_potential) or not np.isfinite(sheath_jump) or np.any(~np.isfinite(ion)) or np.any(ion < 0.0) or np.any(~np.isfinite(electron)) or np.any(electron < 0.0):
            raise ValueError("expander potential arrays and wall sheath values must be finite with nonnegative densities")
        if not np.isfinite(float(self.quasineutrality_relative_error)) or float(self.quasineutrality_relative_error) < 0.0:
            raise ValueError("expander quasineutrality error must be finite and nonnegative")
        if int(self.iterations) < 1:
            raise ValueError("expander potential iterations must be positive")
        if bool(self.converged) and self.failure_reason is not None:
            raise ValueError("converged expander potential cannot retain a failure reason")
        sheath_scale = max(abs(terminal_wall_potential), abs(float(potential[-1])), 1.0)
        if abs(sheath_jump - (terminal_wall_potential - float(potential[-1]))) > 256.0 * np.finfo(float).eps * sheath_scale:
            raise ValueError("expander wall sheath jump must equal the terminal wall potential minus the plasma side potential")
        object.__setattr__(self, "side", side)
        object.__setattr__(self, "full_device_cell_indices", indices)
        object.__setattr__(self, "z_edges_outward_m", edges)
        object.__setattr__(self, "potential_drop_edges_outward_J", potential)
        object.__setattr__(self, "potential_drop_centers_outward_J", potential_centers)
        object.__setattr__(self, "ion_charge_density_m3", ion)
        object.__setattr__(self, "electron_density_m3", electron)

@dataclass(frozen=True)
class ExpanderSystemState:
    """Store both expander sides, terminal totals, mapped populations, and closure status"""
    branches_by_population_side: Mapping[str, ExpanderBranchState]
    potential_by_side: Mapping[str, ExpanderPotentialState]
    full_device_distribution_by_population: Mapping[str, np.ndarray]
    terminal_particle_rate_s_by_component: Mapping[str, float]
    terminal_midplane_kinetic_power_W_by_component: Mapping[str, float]
    terminal_wall_power_W_by_component: Mapping[str, float]
    total_throat_particle_rate_s: float
    total_terminal_particle_rate_s: float
    total_reflected_particle_rate_s: float
    total_unclassified_particle_rate_s: float
    fast_nonterminal_fraction: float
    prompt_unresolved_particle_rate_s: float
    terminal_current_A: float
    terminal_electrostatic_power_W: float
    fixed_fast_boundary_applicable: bool
    terminal_current_equivalent_to_pass11_target: bool
    potential_converged: bool
    status: str
    metadata: Mapping[str, object]

    def __post_init__(self) -> None:
        """Validate branch sums, terminal component identities, current, and system status"""
        branches = dict(self.branches_by_population_side)
        potential = dict(self.potential_by_side)
        distributions = {key: np.asarray(value, dtype=float) for key, value in self.full_device_distribution_by_population.items()}
        if set(potential) != {"left", "right"}:
            raise ValueError("expander system requires left and right potential states")
        if not branches:
            raise ValueError("expander system requires at least one mapped branch")
        population_identity = {'fast_D': ('deuterium', 'fast_deuterium'), 'fast_T': ('tritium', 'fast_tritium')}
        for key, branch in branches.items():
            if key != f"{branch.population_id}:{branch.side}":
                raise ValueError("expander branch keys must match population and side")
            if branch.population_id not in population_identity or branch.species_id != population_identity[branch.population_id][0]:
                raise ValueError("expander branch population and species identity is inconsistent")
            side_potential = potential[branch.side]
            if not np.array_equal(branch.full_device_cell_indices, side_potential.full_device_cell_indices) or not np.allclose(branch.z_edges_outward_m, side_potential.z_edges_outward_m, rtol=0.0, atol=0.0) or not np.allclose(branch.potential_drop_edges_outward_J, side_potential.potential_drop_edges_outward_J, rtol=1.0e-13, atol=0.0):
                raise ValueError("expander branch path does not match the side potential state")
        for key, value in distributions.items():
            if value.ndim != 3 or np.any(~np.isfinite(value)) or np.any(value < 0.0):
                raise ValueError(f"expander distribution {key!r} must be finite and nonnegative")
        population_ids = {branch.population_id for branch in branches.values()}
        if set(distributions) != population_ids:
            raise ValueError("expander system distributions must match the mapped populations")
        for population_id in population_ids:
            branch_sum = sum((np.asarray(branch.full_device_distribution_z_v_pitch, dtype=float) for branch in branches.values() if branch.population_id == population_id), start=np.zeros_like(distributions[population_id]))
            if not np.allclose(distributions[population_id], branch_sum, rtol=1.0e-12, atol=0.0):
                raise ValueError("expander system population distribution does not equal the branch sum")
        rate_map = {str(key): float(value) for key, value in self.terminal_particle_rate_s_by_component.items()}
        midplane_map = {str(key): float(value) for key, value in self.terminal_midplane_kinetic_power_W_by_component.items()}
        wall_map = {str(key): float(value) for key, value in self.terminal_wall_power_W_by_component.items()}
        if set(rate_map) != set(midplane_map) or set(rate_map) != set(wall_map):
            raise ValueError("expander terminal component maps must have identical keys")
        if any(not np.isfinite(value) or value < 0.0 for mapping in (rate_map, midplane_map, wall_map) for value in mapping.values()):
            raise ValueError("expander terminal component values must be finite and nonnegative")
        expected_components = {value[1] for value in population_identity.values()}
        if set(rate_map) != expected_components:
            raise ValueError("expander terminal component maps must contain the D and T fast components")
        nonnegative = (self.total_throat_particle_rate_s, self.total_terminal_particle_rate_s, self.total_reflected_particle_rate_s, self.total_unclassified_particle_rate_s, self.fast_nonterminal_fraction, self.prompt_unresolved_particle_rate_s, self.terminal_current_A, self.terminal_electrostatic_power_W)
        if any(not np.isfinite(float(value)) or float(value) < 0.0 for value in nonnegative):
            raise ValueError("expander system totals must be finite and nonnegative")
        if float(self.fast_nonterminal_fraction) > 1.0:
            raise ValueError("expander nonterminal fractions must not exceed one")
        branch_throat = float(sum(branch.throat_particle_rate_s for branch in branches.values()))
        branch_terminal = float(sum(branch.terminal_particle_rate_s for branch in branches.values()))
        branch_reflected = float(sum(branch.reflected_particle_rate_s for branch in branches.values()))
        branch_unclassified = float(sum(branch.unclassified_particle_rate_s for branch in branches.values()))
        branch_midplane = float(sum(branch.terminal_midplane_kinetic_power_W for branch in branches.values()))
        branch_wall = float(sum(branch.terminal_wall_power_W for branch in branches.values()))
        branch_electrostatic = float(sum(branch.terminal_electrostatic_power_W for branch in branches.values()))
        expected_rate_map = {component: 0.0 for component in expected_components}
        expected_midplane_map = {component: 0.0 for component in expected_components}
        expected_wall_map = {component: 0.0 for component in expected_components}
        for branch in branches.values():
            component = population_identity[branch.population_id][1]
            expected_rate_map[component] += float(branch.terminal_particle_rate_s)
            expected_midplane_map[component] += float(branch.terminal_midplane_kinetic_power_W)
            expected_wall_map[component] += float(branch.terminal_wall_power_W)
        for component in expected_components:
            _require_close(f"expander system {component} terminal rate", rate_map[component], expected_rate_map[component], 1.0e-10)
            _require_close(f"expander system {component} terminal midplane power", midplane_map[component], expected_midplane_map[component], 1.0e-10)
            _require_close(f"expander system {component} terminal wall power", wall_map[component], expected_wall_map[component], 1.0e-10)
        _require_close("expander system throat total", branch_throat + float(self.prompt_unresolved_particle_rate_s), self.total_throat_particle_rate_s, 1.0e-10)
        _require_close("expander system terminal total", branch_terminal, self.total_terminal_particle_rate_s, 1.0e-10)
        _require_close("expander system reflected total", branch_reflected, self.total_reflected_particle_rate_s, 1.0e-10)
        _require_close("expander system unclassified total", branch_unclassified, self.total_unclassified_particle_rate_s, 1.0e-10)
        _require_close("expander system terminal rate map", sum(rate_map.values()), self.total_terminal_particle_rate_s, 1.0e-10)
        _require_close("expander system terminal midplane power map", sum(midplane_map.values()), branch_midplane, 1.0e-10)
        _require_close("expander system terminal wall power map", sum(wall_map.values()), branch_wall, 1.0e-10)
        _require_close("expander system terminal electrostatic power", branch_electrostatic, self.terminal_electrostatic_power_W, 1.0e-10)
        _require_close("expander system terminal current", ELECTRON_CHARGE_C * sum(rate_map.values()), self.terminal_current_A, 1.0e-10)
        _require_close("expander system particle partition", self.total_terminal_particle_rate_s + self.total_reflected_particle_rate_s + self.total_unclassified_particle_rate_s + self.prompt_unresolved_particle_rate_s, self.total_throat_particle_rate_s, 1.0e-10)
        fast_branches = tuple(branch for branch in branches.values() if branch.population_id in {"fast_D", "fast_T"})
        fast_throat = float(sum(branch.throat_particle_rate_s for branch in fast_branches))
        fast_nonterminal = float(sum(branch.reflected_particle_rate_s + branch.unclassified_particle_rate_s for branch in fast_branches))
        _require_close("expander system fast nonterminal fraction", fast_nonterminal / max(fast_throat, np.finfo(float).tiny), self.fast_nonterminal_fraction, 1.0e-10)
        if self.fixed_fast_boundary_applicable and any(not branch.fixed_boundary_applicable for branch in fast_branches):
            raise ValueError("expander system cannot claim a fixed fast boundary when a fast branch fails")
        if bool(self.potential_converged) != all(state.converged for state in potential.values()):
            raise ValueError("expander system potential convergence does not match the side states")
        if self.status not in {"qualified", "diagnostic_unqualified"}:
            raise ValueError("expander system status is invalid")
        if self.status == "qualified" and (not self.fixed_fast_boundary_applicable or not self.potential_converged):
            raise ValueError("qualified expander system requires the fast boundary and potential checks")
        qualified_metadata = self.metadata.get("expander_qualified")
        if qualified_metadata is not None and bool(qualified_metadata) != (self.status == "qualified"):
            raise ValueError("expander system status does not match its qualification metadata")
        object.__setattr__(self, "branches_by_population_side", branches)
        object.__setattr__(self, "potential_by_side", potential)
        object.__setattr__(self, "full_device_distribution_by_population", distributions)
        object.__setattr__(self, "terminal_particle_rate_s_by_component", rate_map)
        object.__setattr__(self, "terminal_midplane_kinetic_power_W_by_component", midplane_map)
        object.__setattr__(self, "terminal_wall_power_W_by_component", wall_map)
        object.__setattr__(self, "metadata", dict(self.metadata))