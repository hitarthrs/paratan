"""Direct consistency validators for already calculated source model states"""
from __future__ import annotations
from collections.abc import Mapping
from typing import Any
import numpy as np
from source_model_revamp.fbis.velocity_space_grid import _distribution_nonnegativity_metrics

NEUTRON_RAW_TO_FUSION_RELATIVE_TOLERANCE = 1.0e-3
NEUTRON_OUT_OF_RANGE_FRACTION_TOLERANCE = 1.0e-14

class SourceConstructionError(ValueError):
    """Structured error raised when the correlated source lacks required finite physical state"""
    def __init__(self, code: str, reason: str, *, quantity: str, value: object = None) -> None:
        """Store a stable error code, reason, quantity, and optional offending value"""
        detail = f"{code}: {quantity}: {reason}"
        if value is not None:
            detail += f"; value={value}"
        super().__init__(detail)
        self.code = str(code)
        self.reason = str(reason)
        self.quantity = str(quantity)
        self.value = value

def _array(name: str, value: object, *, ndim: int | None = None, nonnegative: bool = False, positive: bool = False) -> np.ndarray:
    """Convert one value to a finite numeric array and apply optional dimensionality and sign requirements"""
    result = np.asarray(value, dtype=float)
    if ndim is not None and result.ndim != ndim:
        raise ValueError(f"{name} must have {ndim} dimensions")
    if result.size == 0:
        raise ValueError(f"{name} must not be empty")
    if np.any(~np.isfinite(result)):
        raise ValueError(f"{name} contains a nonfinite value")
    if nonnegative and np.any(result < 0.0):
        raise ValueError(f"{name} contains a negative value")
    if positive and np.any(result <= 0.0):
        raise ValueError(f"{name} contains a nonpositive value")
  
    return result

def validate_geometry_state(geometry: object) -> None:
    """Validate confined geometry topology, positive field and volume quantities, and required plasma boundaries"""
    edges = _array("geometry.z_edges_m", getattr(geometry, "z_edges_m"), ndim=1)
    centers = _array("geometry.z_centers_m", getattr(geometry, "z_centers_m"), ndim=1)
    field = _array("geometry.B_tilde_centers", getattr(geometry, "B_tilde_centers"), ndim=1, positive=True)
    volumes = _array("geometry.cell_volumes_m3", getattr(geometry, "cell_volumes_m3"), ndim=1, positive=True)
    if edges.size != centers.size + 1 or centers.shape != field.shape or centers.shape != volumes.shape:
        raise ValueError("geometry grid arrays have incompatible shapes")
    if np.any(np.diff(edges) <= 0.0):
        raise ValueError("geometry.z_edges_m is not strictly increasing")
    for name in ("B0_T", "mirror_ratio", "half_length_m", "plasma_radius_m", "midplane_area_m2", "volume_m3"):
        value = float(getattr(geometry, name))
        if not np.isfinite(value) or value <= 0.0:
            raise ValueError(f"geometry.{name} must be positive and finite")
    if getattr(geometry, "plasma_boundaries", None) is None:
        raise ValueError("geometry.plasma_boundaries is missing")

def validate_density_state(geometry: object, kinetic: object, beam: object, *, identity_tolerance: float) -> None:
    """
    Validate startup or reference profiles, solved density arrays, beam target density support, and exact magnetic midplane quasineutrality
    
    All density profiles are checked on the confined geometry grid
    """
    profiles = getattr(geometry, "background_profiles", None)
    if profiles is None:
        raise ValueError("geometry.background_profiles is missing")
    density_arrays = [_array(f"background_profiles.{name}", getattr(profiles, name), ndim=1, nonnegative=True) for name in ("shape", "deuterium_density_m3", "tritium_density_m3")]
    temperature = _array("background_profiles.ion_temperature_keV", profiles.ion_temperature_keV, ndim=1, nonnegative=True)
    if len({array.shape for array in (*density_arrays, temperature)}) != 1:
        raise ValueError("seed or reference profile arrays have incompatible shapes")
    maintained = density_arrays[1] + density_arrays[2] > 0.0
    if np.any(maintained & (temperature <= 0.0)):
        raise ValueError("seed or reference temperature must be positive where density is present")
    midplane_shape = float(profiles.normalized_shape_at(profiles.magnetic_midplane_z_m))
    if not np.isfinite(midplane_shape) or midplane_shape != 1.0:
        raise ValueError("seed or reference profile is not normalized at the magnetic midplane")
    symmetry_error = float(profiles.symmetry_relative_error)
    symmetry_tolerance = float(profiles.symmetry_tolerance)
    if not np.isfinite(symmetry_error) or not np.isfinite(symmetry_tolerance) or symmetry_error < 0.0 or symmetry_error > symmetry_tolerance:
        raise ValueError("seed or reference profile symmetry is invalid")
    density_state = getattr(kinetic, "operating_point_density_state", None)
    if density_state is None:
        raise ValueError("kinetic operating point density state is missing")
    for name in ('fast_deuterium_cell_density_m3', 'fast_tritium_cell_density_m3', 'total_positive_charge_cell_density_m3', 'electron_cell_density_m3'):
        values = _array(f"operating_point_density_state.{name}", getattr(density_state, name), ndim=1, nonnegative=True)
        if values.shape != np.asarray(geometry.z_centers_m).shape:
            raise ValueError(f"operating_point_density_state.{name} has an incompatible shape")
    for name in ("electron_midplane_density_m3", "electron_collision_density_m3", "electron_parent_maxwellian_n0_m3"):
        value = float(getattr(density_state, name))
        if not np.isfinite(value) or value <= 0.0:
            raise ValueError(f"operating_point_density_state.{name} must be positive and finite")
    beam_metadata = getattr(beam, "metadata")
    target = _array("beam_target_density_profile_m3", beam_metadata.get("beam_target_density_profile_m3"), ndim=1, nonnegative=True)
    target_D = _array("beam_thermal_deuterium_target_density_profile_m3", beam_metadata.get("beam_thermal_deuterium_target_density_profile_m3"), ndim=1, nonnegative=True)
    target_T = _array("beam_thermal_tritium_target_density_profile_m3", beam_metadata.get("beam_thermal_tritium_target_density_profile_m3"), ndim=1, nonnegative=True)
    defined = np.asarray(beam_metadata.get("beam_target_density_profile_defined_mask"), dtype=bool)
    support = np.asarray(beam_metadata.get("beam_target_density_profile_support_mask"), dtype=bool)
    if any(values.shape != target.shape for values in (target_D, target_T, defined, support)) or target.shape != np.asarray(geometry.z_centers_m).shape:
        raise ValueError("beam target density profile arrays have incompatible shapes")
    if np.any(support & ~defined):
        raise ValueError("beam target density support contains undefined values")
    residual = float(density_state.exact_midplane_quasineutrality_relative_error)
    tolerance = float(identity_tolerance)
    if not np.isfinite(residual) or not np.isfinite(tolerance) or tolerance < 0.0 or abs(residual) > tolerance:
        raise ValueError("exact midplane quasineutrality residual is invalid: "f"residual={residual:.16e}, tolerance={tolerance:.16e}")

def validate_kinetic_state(kinetic: object) -> None:
    """Reject materially negative invariant, local, density, or full device kinetic arrays using the shared distribution roundoff criterion"""
    arrays: dict[str, object] = {"final_distribution_v_lambda": getattr(kinetic, "final_distribution_v_lambda")}
    for name in ("local_distribution_z_v_lambda", "local_distribution_z_v_pitch", "local_density_m3", "full_device_local_distribution_z_v_pitch"):
        value = getattr(kinetic, name, None)
        if value is not None:
            arrays[name] = value
    for mapping_name in ("final_distribution_v_lambda_by_species", "local_distribution_z_v_lambda_by_species", "local_distribution_z_v_pitch_by_species", "local_density_m3_by_species", "full_device_local_distribution_z_v_pitch_by_species"):
        mapping = getattr(kinetic, mapping_name, None)
        if mapping is not None:
            arrays.update({f"{mapping_name}.{key}": value for key, value in mapping.items()})
    for name, value in arrays.items():
        _, _, _, _, significant = _distribution_nonnegativity_metrics(value, name=name)
        if significant:
            raise ValueError(f"{name} contains a materially negative distribution or moment")

def validate_electrostatic_state(kinetic: object) -> None:
    """
    Validate the shared Eq 70 electrostatic profile and exact left throat, magnetic midplane, and right throat nodes
    
    Floating reference solutions must also provide a converged scalar closure and a valid exact midplane root
    """
    system = getattr(kinetic, "fast_ion_system_state", None)
    profile = None if system is None else system.shared_electrostatic_profile
    if profile is None:
        raise ValueError("shared electrostatic profile is missing")
    potential = _array("electrostatic_profile.potential_energy_J", profile.potential_energy_J, ndim=1, nonnegative=True)
    node_zeta = _array("electrostatic_profile.electrostatic_node_zeta", profile.electrostatic_node_zeta, ndim=1)
    node_potential = _array("electrostatic_profile.electrostatic_node_potential_energy_J", profile.electrostatic_node_potential_energy_J, ndim=1)
    node_electron = _array("electrostatic_profile.electrostatic_node_electron_density_m3", profile.electrostatic_node_electron_density_m3, ndim=1, nonnegative=True)
    node_ion = _array("electrostatic_profile.electrostatic_node_ion_density_m3", profile.electrostatic_node_ion_density_m3, ndim=1, nonnegative=True)
    if potential.shape != np.asarray(getattr(kinetic, "local_density_m3")).shape:
        raise ValueError("electrostatic potential and kinetic density grids are incompatible")
    if len({values.shape for values in (node_zeta, node_potential, node_electron, node_ion)}) != 1:
        raise ValueError("electrostatic node arrays have incompatible shapes")
    if node_zeta.size < 3:
        raise ValueError("electrostatic profile does not contain the required Eq 70 nodes")
    indices = tuple(int(getattr(profile, name)) for name in ("electrostatic_left_throat_index", "electrostatic_midplane_index", "electrostatic_right_throat_index"))
    if not (0 <= indices[0] < indices[1] < indices[2] < node_zeta.size):
        raise ValueError("electrostatic profile is missing ordered throat and midplane nodes")
    if not all(getattr(profile, name) is True for name in ("electrostatic_left_throat_node_exact", "electrostatic_midplane_node_exact", "electrostatic_right_throat_node_exact")):
        raise ValueError("electrostatic profile required exact throat or magnetic midplane nodes are invalid")
    floating = bool(getattr(profile, "floating_nonnegative_reference_active", False))
    if floating:
        if getattr(profile, "floating_scalar_closure_converged", False) is not True or getattr(profile, "direct_exact_midplane_root_valid", False) is not True:
            raise ValueError("floating electrostatic reference or exact magnetic midplane root is invalid")
        exact_midplane = float(getattr(profile, "exact_midplane_potential_energy_J", np.nan))
        if not np.isfinite(exact_midplane) or exact_midplane < 0.0:
            raise ValueError("floating electrostatic exact magnetic midplane drop is invalid")
        if float(node_potential[indices[1]]) != exact_midplane:
            raise ValueError("electrostatic exact magnetic midplane node is inconsistent with the solved drop")
    elif getattr(profile, "electrostatic_midplane_gauge_exact", False) is not True:
        raise ValueError("legacy electrostatic profile requires the exact magnetic midplane gauge")

def validate_current_identities(metadata: Mapping[str, Any]) -> None:
    """Validate current component sums and applicable left plus right side current identities"""
    components = metadata.get("total_current_balance_component_currents_A")
    if isinstance(components, Mapping) and components:
        values = _array("total_current_balance_component_currents_A", tuple(components.values()), ndim=1)
        reported = float(metadata.get("total_current_balance_component_current_sum_A", np.nan))
        calculated = float(sum(float(value) for value in values))
        if not np.isfinite(reported) or reported != calculated:
            raise ValueError("total current component sum is internally inconsistent")
    for prefix in ("total_current_balance", "terminal_current_balance"):
        if metadata.get(f"{prefix}_side_sum_identity_applicable") is not True:
            continue
        left = float(metadata.get(f"{prefix}_side_resolved_left_ion_current_A", metadata.get("terminal_current_left_ion_current_A", np.nan)))
        right = float(metadata.get(f"{prefix}_side_resolved_right_ion_current_A", metadata.get("terminal_current_right_ion_current_A", np.nan)))
        reported = metadata.get(f"{prefix}_side_resolved_current_sum_A", metadata.get("terminal_current_side_sum_A"))
        if reported is None:
            continue
        reported = float(reported)
        if not all(np.isfinite(value) for value in (left, right, reported)) or left + right != reported:
            raise ValueError(f"{prefix} left plus right current sum is internally inconsistent")

def validate_expander_state(expander: object) -> None:
    """Validate expander branch distribution shapes and exact particle and kinetic power classification sums"""
    system = getattr(expander, "system_state", expander)
    if system is None:
        raise ValueError("expander system state is missing")
    branches = getattr(system, "branches_by_population_side", None)
    if not isinstance(branches, Mapping) or not branches:
        raise ValueError("expander branch state is missing")
    for key, branch in branches.items():
        distribution = _array(f"expander.{key}.full_device_distribution_z_v_pitch", branch.full_device_distribution_z_v_pitch, ndim=3, nonnegative=True)
        density = _array(f"expander.{key}.full_device_density_m3", branch.full_device_density_m3, ndim=1, nonnegative=True)
        if distribution.shape[0] != density.size:
            raise ValueError(f"expander.{key} branch arrays have incompatible shapes")
        for name in ("throat_particle_rate_s", "terminal_particle_rate_s", "reflected_particle_rate_s", "unclassified_particle_rate_s", "throat_midplane_kinetic_power_W", "terminal_midplane_kinetic_power_W", "reflected_midplane_kinetic_power_W", "unclassified_midplane_kinetic_power_W"):
            value = float(getattr(branch, name))
            if not np.isfinite(value) or value < 0.0:
                raise ValueError(f"expander.{key}.{name} must be finite and nonnegative")
        particle_sum = float(branch.terminal_particle_rate_s) + float(branch.reflected_particle_rate_s) + float(branch.unclassified_particle_rate_s)
        if not np.isclose(particle_sum, float(branch.throat_particle_rate_s), rtol=1.0e-10, atol=0.0):
            raise ValueError(f"expander.{key} particle classification is inconsistent")
        energy_sum = float(branch.terminal_midplane_kinetic_power_W) + float(branch.reflected_midplane_kinetic_power_W) + float(branch.unclassified_midplane_kinetic_power_W)
        if not np.isclose(energy_sum, float(branch.throat_midplane_kinetic_power_W), rtol=1.0e-10, atol=0.0):
            raise ValueError(f"expander.{key} energy classification is inconsistent")

def validate_fusion_neutron_state(fusion: object, neutrons: object, *, relative_tolerance: float) -> None:
    """Validate finite nonnegative fusion and neutron rates and require total and raw quadrature rate conservation within tolerance"""
    profile = getattr(fusion, "combined_profile")
    fusion_rates = _array("fusion.total_neutron_source_density_m3_s", profile.total_neutron_source_density_m3_s, ndim=1, nonnegative=True)
    neutron_rates = _array("neutrons.source_rate_matrix_s", getattr(neutrons, "source_rate_matrix_s"), ndim=2, nonnegative=True)
    fusion_total = float(profile.total_neutron_rate_s)
    neutron_total = float(getattr(neutrons, "total_neutron_rate_s"))
    if not np.isfinite(fusion_total) or fusion_total < 0.0 or not np.isfinite(neutron_total) or neutron_total < 0.0:
        raise ValueError("fusion or neutron total rate is invalid")
    residual = abs(neutron_total - fusion_total) / max(abs(fusion_total), np.finfo(float).tiny)
    if residual > float(relative_tolerance):
        raise ValueError("fusion and neutron total rates are inconsistent")
    raw_residual = getattr(neutrons, "metadata", {}).get("neutron_raw_to_fusion_conservative_max_relative_error")
    if raw_residual is not None and (not np.isfinite(float(raw_residual)) or abs(float(raw_residual)) > float(relative_tolerance)):
        raise ValueError("raw neutron quadrature and fusion rates are inconsistent")
    if fusion_rates.size == 0 or neutron_rates.size == 0:
        raise ValueError("fusion or neutron rate support is empty")

def validate_neutron_domain(metadata: Mapping[str, Any], *, fraction_tolerance: float) -> None:
    """Require strict neutron energy domain handling and bound the measured out of range raw rate fraction"""
    if metadata.get("neutron_allow_out_of_range") is not False:
        raise ValueError("neutron energy domain is not configured strictly")
    fraction = float(metadata.get("neutron_out_of_range_fraction_of_raw", np.nan))
    fractions = [fraction]
    components = metadata.get("neutron_component_raw_diagnostics")
    if isinstance(components, (tuple, list)):
        for component in components:
            if not isinstance(component, Mapping):
                raise ValueError("neutron out of range component evidence is malformed")
            raw = float(component.get("raw_total_rate_s", np.nan))
            outside = float(component.get("out_of_range_rate_s", np.nan))
            if not np.isfinite(raw) or not np.isfinite(outside) or raw < 0.0 or outside < 0.0:
                raise ValueError("neutron out of range component evidence is invalid")
            fractions.append(outside / max(raw, np.finfo(float).tiny))
    if any(not np.isfinite(value) or value < 0.0 for value in fractions) or max(fractions) > float(fraction_tolerance):
        raise ValueError("neutron out of range quadrature fraction is invalid")

def validate_correlated_event_state(event_bank: object, *, expected_physical_rate_s: float | None = None, relative_tolerance: float | None = None, mass_shell_tolerance: float | None = None) -> None:
    """Validate complete correlated neutron event arrays, unit directions, normalized probabilities, physical rate agreement, and optional relativistic mass shell evidence"""
    positions = _array("event_bank.positions_m", getattr(event_bank, "positions_m"), ndim=2)
    directions = _array("event_bank.directions", getattr(event_bank, "directions"), ndim=2)
    energies = _array("event_bank.energies_J", getattr(event_bank, "energies_J"), ndim=1, nonnegative=True)
    weights = _array("event_bank.normalized_weights", getattr(event_bank, "normalized_weights"), ndim=1, nonnegative=True)
    count = positions.shape[0]
    if positions.shape != (count, 3) or directions.shape != (count, 3) or energies.shape != (count,) or weights.shape != (count,):
        raise ValueError("correlated event arrays have incompatible shapes")
    for name in ("reaction_keys", "component_labels"):
        labels = np.asarray(getattr(event_bank, name))
        if labels.shape != (count,) or any(not str(value) for value in labels):
            raise ValueError(f"correlated event {name} is incomplete")
    axial = np.asarray(getattr(event_bank, "axial_cell_indices"), dtype=int)
    if axial.shape != (count,) or np.any(axial < 0):
        raise ValueError("correlated event axial cell indices are invalid")
    for name in ("center_of_mass_energy_J", "equivalent_deuteron_lab_energy_eV"):
        audit_energy = _array(f"event_bank.{name}", getattr(event_bank, name), ndim=1, nonnegative=True)
        if audit_energy.shape != (count,):
            raise ValueError(f"correlated event {name} has an incompatible shape")
    mu = _array("event_bank.cm_emission_mu", getattr(event_bank, "cm_emission_mu"), ndim=1)
    if mu.shape != (count,) or np.any(mu < -1.0) or np.any(mu > 1.0):
        raise ValueError("correlated event CM emission cosine is invalid")
    for name in ("reactant_a_velocity_m_s", "reactant_b_velocity_m_s"):
        velocity = _array(f"event_bank.{name}", getattr(event_bank, name), ndim=2)
        if velocity.shape != (count, 3):
            raise ValueError(f"correlated event {name} has an incompatible shape")
    if not np.allclose(np.linalg.norm(directions, axis=1), 1.0, rtol=0.0, atol=2.0e-12):
        raise ValueError("correlated event directions are not unit vectors")
    # Event probabilities are normalized separately from the retained physical neutron rate
    if not np.isclose(np.sum(weights), 1.0, rtol=0.0, atol=2.0e-12):
        raise ValueError("correlated event probabilities do not sum to one")
    physical_rate = float(getattr(event_bank, "physical_total_rate_s"))
    if not np.isfinite(physical_rate) or physical_rate <= 0.0:
        raise ValueError("correlated event physical rate is not positive and finite")
    if expected_physical_rate_s is not None:
        if relative_tolerance is None:
            raise ValueError("correlated event rate comparison requires the existing tolerance")
        expected = float(expected_physical_rate_s)
        residual = abs(physical_rate - expected) / max(abs(expected), np.finfo(float).tiny)
        if not np.isfinite(expected) or expected <= 0.0 or residual > float(relative_tolerance):
            raise ValueError("correlated event physical rate is inconsistent")
    event_metadata = getattr(event_bank, "metadata", None)
    if isinstance(event_metadata, Mapping) and "maximum_residual_mass_shell_relative_error" in event_metadata:
        if mass_shell_tolerance is None:
            raise ValueError("correlated event kinematics comparison requires the existing tolerance")
        mass_shell_error = float(event_metadata["maximum_residual_mass_shell_relative_error"])
        if not np.isfinite(mass_shell_error) or mass_shell_error < 0.0 or mass_shell_error > float(mass_shell_tolerance):
            raise ValueError("correlated event relativistic two body kinematics are invalid")

def validate_openmc_file_source(file_metadata: object, *, probability_tolerance: float) -> None:
    """Validate measured OpenMC HDF5 structure, roundtrip state, probability normalization, particle identity, correlation preservation, and physical rate"""
    for name in ("hdf5_structure_valid", "roundtrip_valid", "roundtrip_position_valid", "roundtrip_direction_valid", "roundtrip_energy_valid", "roundtrip_probability_weights_valid", "roundtrip_particle_identity_valid", "joint_lab_energy_direction_correlation_preserved"):
        if getattr(file_metadata, name, None) is not True:
            raise ValueError(f"OpenMC file source {name} is invalid")
    probability_sum = float(getattr(file_metadata, "event_probability_weight_sum", np.nan))
    if not np.isfinite(probability_sum) or not np.isclose(probability_sum, 1.0, rtol=0.0, atol=float(probability_tolerance)):
        raise ValueError("OpenMC file source event probabilities are invalid")
    physical_rate = float(getattr(file_metadata, "physical_total_rate_s", np.nan))
    if not np.isfinite(physical_rate) or physical_rate <= 0.0:
        raise ValueError("OpenMC file source physical rate is invalid")

def validate_source_constructibility(event_bank: object | None = None, *, physical_total_rate_s: object = None, weights: object = None, positions_m: object = None, directions: object = None, energies_J: object = None) -> None:
    """
    Validate the minimum correlated source state before export
    
    A positive physical rate, nonempty compatible event arrays, finite nonnegative energies and probabilities, and at least one positive event weight are required
    """
    if event_bank is not None:
        physical_total_rate_s = getattr(event_bank, "physical_total_rate_s", None)
        weights = getattr(event_bank, "normalized_weights", None)
        positions_m = getattr(event_bank, "positions_m", None)
        directions = getattr(event_bank, "directions", None)
        energies_J = getattr(event_bank, "energies_J", None)
    elif all(value is None for value in (physical_total_rate_s, weights, positions_m, directions, energies_J)):
        raise SourceConstructionError("missing_event_bank", "correlated event bank is required", quantity="correlated_event_bank")
    missing = tuple(name for name, value in (("physical_total_rate_s", physical_total_rate_s), ("normalized_weights", weights), ("positions_m", positions_m), ("directions", directions), ("energies_J", energies_J)) if value is None)
    if missing:
        raise SourceConstructionError("missing_source_quantity", "required correlated source data is missing", quantity=missing[0], value=missing)
    try:
        rate = float(physical_total_rate_s)
        position = np.asarray(positions_m, dtype=float)
        direction = np.asarray(directions, dtype=float)
        energy = np.asarray(energies_J, dtype=float)
        probability = np.asarray(weights, dtype=float)
    except (TypeError, ValueError) as exc:
        raise SourceConstructionError("nonnumeric_source_data", "source arrays are not numerical", quantity="correlated_event_bank") from exc
    if not np.isfinite(rate) or rate <= 0.0:
        raise SourceConstructionError("invalid_total_rate", "must be positive and finite", quantity="physical_total_rate_s", value=rate)
    if probability.ndim != 1 or probability.size == 0:
        raise SourceConstructionError("empty_event_support", "event support must be a nonempty vector", quantity="normalized_weights", value=probability.shape)
    count = probability.size
    if position.shape != (count, 3) or direction.shape != (count, 3) or energy.shape != (count,):
        raise SourceConstructionError("incompatible_event_shapes", "source event arrays have incompatible shapes", quantity="correlated_event_bank", value={"positions_m": position.shape, "directions": direction.shape, "energies_J": energy.shape, "normalized_weights": probability.shape})
    for name, value in (("positions_m", position), ("directions", direction), ("energies_J", energy), ("weights", probability)):
        if np.any(~np.isfinite(value)):
            raise SourceConstructionError("nonfinite_event_data", "contains a nonfinite value", quantity=name)
    if np.any(energy < 0.0) or np.any(probability < 0.0):
        raise SourceConstructionError("negative_event_data", "source energies and weights must be nonnegative", quantity="energies_J_or_normalized_weights")
    if not np.any(probability > 0.0):
        raise SourceConstructionError("no_positive_event_weight", "at least one positive event weight is required", quantity="normalized_weights")
    total = float(np.sum(probability))
    if not np.isfinite(total) or total <= 0.0:
        raise SourceConstructionError("nonnormalizable_event_weights", "event weights are not normalizable", quantity="normalized_weights", value=total)

__all__ = [
    "SourceConstructionError",
    "validate_correlated_event_state",
    "validate_current_identities",
    "validate_density_state",
    "validate_electrostatic_state",
    "validate_expander_state",
    "validate_fusion_neutron_state",
    "validate_geometry_state",
    "validate_kinetic_state",
    "validate_neutron_domain",
    "validate_openmc_file_source",
    "validate_source_constructibility",
]