"""
Deterministic neutron spectrum formation from local gyrotropic reactant distributions

Reactant cells are weighted by `f d³v`, total reaction weighting uses `σ(E_cm) |v_a − v_b|`, and neutron lab energies are obtained from relativistic two body kinematics
"""
from __future__ import annotations
from dataclasses import dataclass
from collections.abc import Callable
import numpy as np
from numpy.typing import ArrayLike
from source_model_revamp.fbis.velocity_space_grid import SpeedGrid, PitchGrid, gyrotropic_velocity_cell_volumes, require_gyrotropic_distribution_shape
from source_model_revamp.fusion.cross_sections import reduced_mass_kg
from source_model_revamp.neutrons.kinematics import TwoBodyReactionKinematics, center_of_mass_kinetic_energy_J, isotropic_emission_directions, neutron_lab_events_relativistic_batch, neutron_lab_kinetic_energy_endpoints_relativistic_batch, velocity_vector_from_speed_pitch_gyro
from source_model_revamp.pipeline_errors import NeutronEnergyDomainError

ANGULAR_MODEL_ISOTROPIC_CM = "isotropic_cm"
_MAX_PAIR_STATES_PER_BATCH = 131072

@dataclass(frozen=True)
class NeutronEnergySpectrum:
    """
    Energy binned neutron rate density spectrum
    
    `bin_rate_density_m3_s` stores the rate density integrated over each energy bin
    `raw_total_rate_density_m3_s` retains the rate before any optional target normalization
    """
    energy_edges_J: np.ndarray
    bin_rate_density_m3_s: np.ndarray
    total_rate_density_m3_s: float
    raw_total_rate_density_m3_s: float
    out_of_range_rate_density_m3_s: float
    angular_model: str = ANGULAR_MODEL_ISOTROPIC_CM

    @property
    def energy_centers_J(self) -> np.ndarray:
        """Return neutron energy bin centers in J"""
        return 0.5 * (self.energy_edges_J[:-1] + self.energy_edges_J[1:])
    @property
    def spectral_density_m3_s_J(self) -> np.ndarray:
        """Return `dR / dE` in m⁻³ s⁻¹ J⁻¹ from the energy integrated bin rates"""
        widths = np.diff(self.energy_edges_J)
        return np.divide(self.bin_rate_density_m3_s, widths, out=np.zeros_like(self.bin_rate_density_m3_s), where=widths > 0.0)

@dataclass(frozen=True)
class NeutronAngularEnergySpectrum:
    """
    Joint neutron energy and lab pitch cosine rate density spectrum
    
    `energy_mu_rate_density_m3_s` has shape `(n_energy, n_mu)` and stores the rate density integrated over each two dimensional bin
    """
    energy_edges_J: np.ndarray
    mu_edges: np.ndarray
    energy_mu_rate_density_m3_s: np.ndarray
    total_rate_density_m3_s: float
    raw_total_rate_density_m3_s: float
    out_of_range_rate_density_m3_s: float
    angular_model: str = ANGULAR_MODEL_ISOTROPIC_CM

    @property
    def energy_centers_J(self) -> np.ndarray:
        """Return neutron energy bin centers in J"""
        return 0.5 * (self.energy_edges_J[:-1] + self.energy_edges_J[1:])
    @property
    def mu_centers(self) -> np.ndarray:
        """Return lab pitch cosine bin centers"""
        return 0.5 * (self.mu_edges[:-1] + self.mu_edges[1:])
    @property
    def energy_bin_rate_density_m3_s(self) -> np.ndarray:
        """Return the neutron rate density in each energy bin after summing lab pitch cosine"""
        return np.sum(self.energy_mu_rate_density_m3_s, axis=1)
    @property
    def mu_bin_rate_density_m3_s(self) -> np.ndarray:
        """Return the neutron rate density in each lab pitch cosine bin after summing energy"""
        return np.sum(self.energy_mu_rate_density_m3_s, axis=0)

def validate_bin_edges(edges: ArrayLike, name: str) -> np.ndarray:
    """Validate finite strictly increasing histogram edges"""
    arr = np.asarray(edges, dtype=float)
    if arr.ndim != 1:
        raise ValueError(f"{name} must be 1D")
    if arr.size < 2:
        raise ValueError(f"{name} must have at least two entries")
    if not np.all(np.isfinite(arr)):
        raise ValueError(f"{name} must contain only finite values")
    if np.any(np.diff(arr) <= 0.0):
        raise ValueError(f"{name} must be strictly increasing")
    
    return arr

def _normalize_binned_rates(binned: np.ndarray, raw_total: float, normalize_to_rate_density_m3_s: float | None) -> tuple[np.ndarray, float]:
    """Optionally rescale represented bins to a requested nonnegative total rate density"""
    if normalize_to_rate_density_m3_s is None:
        return binned, raw_total
    target = float(normalize_to_rate_density_m3_s)
    if not np.isfinite(target) or target < 0.0:
        raise ValueError("normalize_to_rate_density_m3_s must be finite and nonnegative")
    if target == 0.0:
        return np.zeros_like(binned), 0.0
    represented_total = float(np.sum(binned))
    if represented_total <= 0.0:
        raise ValueError("cannot normalize an empty represented spectrum to a positive target rate")
    
    return binned * (target / represented_total), target

def _histogram_energy_mu_events(energies_J: np.ndarray, mu_values: np.ndarray, weights_m3_s: np.ndarray, energy_edges_J: np.ndarray, mu_edges: np.ndarray, *, allow_out_of_range: bool, normalize_to_rate_density_m3_s: float | None) -> NeutronAngularEnergySpectrum:
    """Histogram weighted neutron events in energy and lab pitch cosine while retaining raw and out of range rates"""
    raw_total = float(np.sum(weights_m3_s))
    inside_e = (energies_J >= energy_edges_J[0]) & (energies_J < energy_edges_J[-1])
    inside_e |= energies_J == energy_edges_J[-1]
    inside_mu = (mu_values >= mu_edges[0]) & (mu_values < mu_edges[-1])
    inside_mu |= mu_values == mu_edges[-1]
    inside = inside_e & inside_mu
    out_of_range = float(np.sum(weights_m3_s[~inside]))
    if out_of_range > max(1.0e-30, 1.0e-12 * max(raw_total, 1.0)) and not allow_out_of_range:
        raise NeutronEnergyDomainError("neutron energy/mu bins do not cover all weighted events, " f"out of range rate density is {out_of_range:.6e} m^-3 s^-1")
    hist, _, _ = np.histogram2d(energies_J[inside], mu_values[inside], bins=(energy_edges_J, mu_edges), weights=weights_m3_s[inside])
    normalized, total = _normalize_binned_rates(hist.astype(float), raw_total, normalize_to_rate_density_m3_s)

    return NeutronAngularEnergySpectrum(energy_edges_J=energy_edges_J, mu_edges=mu_edges, energy_mu_rate_density_m3_s=normalized, total_rate_density_m3_s=float(total), raw_total_rate_density_m3_s=raw_total, out_of_range_rate_density_m3_s=out_of_range)

def _validate_reactants(speed_grid_a: SpeedGrid, pitch_grid_a: PitchGrid, distribution_a_v_xi: ArrayLike, speed_grid_b: SpeedGrid, pitch_grid_b: PitchGrid, distribution_b_v_xi: ArrayLike, num_gyroangle_points: int, num_emission_directions: int) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, int, np.ndarray]:
    """
    Validate two local gyrotropic reactant distributions and build `f d³v` cell weights
    
    The returned emission directions are deterministic isotropic center of momentum directions
    """
    f_a = require_gyrotropic_distribution_shape(distribution_a_v_xi, speed_grid_a, pitch_grid_a, name="distribution_a_v_xi", allow_negative=False)
    f_b = require_gyrotropic_distribution_shape(distribution_b_v_xi, speed_grid_b, pitch_grid_b, name="distribution_b_v_xi", allow_negative=False)
    nphi = int(num_gyroangle_points)
    if nphi < 1:
        raise ValueError("num_gyroangle_points must be at least one")
    emission_dirs = isotropic_emission_directions(num_emission_directions)
    weights_a = f_a * gyrotropic_velocity_cell_volumes(speed_grid_a, pitch_grid_a)
    weights_b = f_b * gyrotropic_velocity_cell_volumes(speed_grid_b, pitch_grid_b)

    return f_a, f_b, weights_a, weights_b, nphi, emission_dirs

def _require_gyrotropic_distribution_stack(values: ArrayLike, speed_grid: SpeedGrid, pitch_grid: PitchGrid, name: str) -> np.ndarray:
    """Validate one or more local gyrotropic distributions with shape `(n_z, n_speed, n_pitch)`"""
    array = np.asarray(values, dtype=float)
    expected = (speed_grid.centers_m_s.size, pitch_grid.centers.size)
    if array.ndim == 2:
        array = array[None, :, :]
    if array.ndim != 3 or array.shape[1:] != expected:
        raise ValueError(f"{name} must have shape (n_z, {expected[0]}, {expected[1]})")
    if not np.all(np.isfinite(array)) or np.any(array < 0.0):
        raise ValueError(f"{name} must be finite and nonnegative")

    return array

def _velocity_vectors_for_grid(speed_grid: SpeedGrid, pitch_grid: PitchGrid, gyro_angle_rad: float = 0.0) -> np.ndarray:
    """Return local velocity vectors in speed major and pitch minor order"""
    speeds = np.repeat(np.asarray(speed_grid.centers_m_s, dtype=float), pitch_grid.centers.size)
    pitches = np.tile(np.asarray(pitch_grid.centers, dtype=float), speed_grid.centers_m_s.size)
    perpendicular = speeds * np.sqrt(np.maximum(1.0 - pitches**2, 0.0))
    angle = float(gyro_angle_rad)

    return np.column_stack((perpendicular * np.cos(angle), perpendicular * np.sin(angle), speeds * pitches))

def _velocity_vectors_for_grid_and_gyroangles(speed_grid: SpeedGrid, pitch_grid: PitchGrid, gyro_angles_rad: np.ndarray) -> np.ndarray:
    """Return local velocity vectors with reactant cell major and relative gyro angle minor ordering"""
    base = _velocity_vectors_for_grid(speed_grid, pitch_grid, 0.0)
    speeds_perpendicular = np.linalg.norm(base[:, :2], axis=1)
    parallel = base[:, 2]
    angles = np.asarray(gyro_angles_rad, dtype=float)
    perpendicular = np.repeat(speeds_perpendicular, angles.size)
    repeated_angles = np.tile(angles, base.shape[0])

    return np.column_stack((perpendicular * np.cos(repeated_angles), perpendicular * np.sin(repeated_angles), np.repeat(parallel, angles.size))).reshape(base.shape[0], angles.size, 3)

def _cross_section_array(cross_section_function_m2: Callable[[ArrayLike], ArrayLike], energy_J: np.ndarray) -> np.ndarray:
    """Evaluate a total cross section callable on an energy array with scalar fallback"""
    try:
        sigma = np.asarray(cross_section_function_m2(energy_J), dtype=float)
        if sigma.ndim == 0:
            sigma = np.full(energy_J.shape, float(sigma), dtype=float)
        else:
            sigma = np.broadcast_to(sigma, energy_J.shape).astype(float, copy=False)
    except (TypeError, ValueError):
        sigma = np.asarray([cross_section_function_m2(float(value)) for value in energy_J], dtype=float)
    if sigma.shape != energy_J.shape or np.any(~np.isfinite(sigma)) or np.any(sigma < 0.0):
        raise ValueError("cross_section_function_m2 must return finite nonnegative values")

    return sigma

def neutron_energy_spectrum_matrix_from_gyrotropic_distributions(speed_grid_a: SpeedGrid, pitch_grid_a: PitchGrid, distributions_a_z_v_xi: ArrayLike, speed_grid_b: SpeedGrid, pitch_grid_b: PitchGrid, distributions_b_z_v_xi: ArrayLike, reaction_kinematics: TwoBodyReactionKinematics, cross_section_function_m2: Callable[[ArrayLike], ArrayLike], energy_edges_J: ArrayLike, *, identical_population: bool = False, num_gyroangle_points: int = 8, num_emission_directions: int = 24, normalize_to_rate_density_m3_s: ArrayLike | None = None, allow_out_of_range: bool = False, max_pair_states_per_batch: int = _MAX_PAIR_STATES_PER_BATCH, exact_isotropic_energy_marginal: bool = False) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Build neutron energy spectra for multiple axial gyrotropic reactant pairs
    
    The raw pair kernel is `f_a d³v_a f_b d³v_b σ(E_cm) |v_a − v_b|` with a `1 / 2` factor for identical populations
    When `exact_isotropic_energy_marginal` is true, isotropic center of momentum emission is integrated exactly as a uniform lab energy interval for each reactant pair
    Otherwise deterministic center of momentum directions approximate the energy marginal
    """
    stack_a = _require_gyrotropic_distribution_stack(distributions_a_z_v_xi, speed_grid_a, pitch_grid_a, "distributions_a_z_v_xi")
    stack_b = _require_gyrotropic_distribution_stack(distributions_b_z_v_xi, speed_grid_b, pitch_grid_b, "distributions_b_z_v_xi")
    if stack_a.shape[0] != stack_b.shape[0]:
        raise ValueError("reactant distribution stacks must have the same axial size")
    n_z = stack_a.shape[0]
    energy_edges = validate_bin_edges(energy_edges_J, "energy_edges_J")
    n_energy = energy_edges.size - 1
    nphi = int(num_gyroangle_points)
    if nphi < 1:
        raise ValueError("num_gyroangle_points must be at least one")
    nemission = int(num_emission_directions)
    if nemission < 1:
        raise ValueError("num_emission_directions must be at least one")
    use_exact_energy_marginal = bool(exact_isotropic_energy_marginal)
    emission_dirs = None
    if not use_exact_energy_marginal:
        emission_dirs = isotropic_emission_directions(nemission)
        nemission = emission_dirs.shape[0]
    max_states = int(max_pair_states_per_batch)
    if max_states < 1:
        raise ValueError("max_pair_states_per_batch must be positive")
    if max_states < nphi:
        raise ValueError("max_pair_states_per_batch must be at least num_gyroangle_points")

    velocity_weights_a = gyrotropic_velocity_cell_volumes(speed_grid_a, pitch_grid_a).reshape(-1)
    velocity_weights_b = gyrotropic_velocity_cell_volumes(speed_grid_b, pitch_grid_b).reshape(-1)
    weighted_a = stack_a.reshape(n_z, -1) * velocity_weights_a[None, :]
    weighted_b = stack_b.reshape(n_z, -1) * velocity_weights_b[None, :]
    active_a = np.flatnonzero(np.any(weighted_a != 0.0, axis=0))
    active_b = np.flatnonzero(np.any(weighted_b != 0.0, axis=0))
    histograms = np.zeros((n_z, n_energy), dtype=float)
    raw_totals = np.zeros(n_z, dtype=float)
    out_of_range_totals = np.zeros(n_z, dtype=float)
    if active_a.size == 0 or active_b.size == 0:
        if normalize_to_rate_density_m3_s is not None:
            targets = np.asarray(normalize_to_rate_density_m3_s, dtype=float)
            if targets.ndim == 0:
                targets = np.full(n_z, float(targets), dtype=float)
            if targets.shape != (n_z,) or np.any(~np.isfinite(targets)) or np.any(targets < 0.0):
                raise ValueError("normalize_to_rate_density_m3_s must contain one finite nonnegative value per axial plane")
            if np.any(targets > 0.0):
                raise ValueError("cannot normalize a zero raw spectrum to a positive target rate")
        return histograms, raw_totals, out_of_range_totals

    vectors_a = _velocity_vectors_for_grid(speed_grid_a, pitch_grid_a)
    gyro_angles = (np.arange(nphi, dtype=float) + 0.5) * (2.0 * np.pi / nphi)
    vectors_b_by_gyro = _velocity_vectors_for_grid_and_gyroangles(speed_grid_b, pitch_grid_b, gyro_angles)
    # Identical reactant populations require the 1 / 2 pair counting factor
    pair_factor = 0.5 if bool(identical_population) else 1.0
    reduced_mass = reduced_mass_kg(reaction_kinematics.reactant_a_mass_kg, reaction_kinematics.reactant_b_mass_kg)
    b_chunk_size = min(active_b.size, max(1, int(np.sqrt(max_states / nphi))))
    a_chunk_size = min(active_a.size, max(1, max_states // (b_chunk_size * nphi)))

    for a_start in range(0, active_a.size, a_chunk_size):
        a_indices = active_a[a_start:a_start + a_chunk_size]
        for b_start in range(0, active_b.size, b_chunk_size):
            b_indices = active_b[b_start:b_start + b_chunk_size]
            b_vectors = vectors_b_by_gyro[b_indices].reshape(-1, 3)
            states_per_a = b_vectors.shape[0]
            v_a = np.repeat(vectors_a[a_indices], states_per_a, axis=0)
            v_b = np.tile(b_vectors, (a_indices.size, 1))
            a_local_indices = np.repeat(np.arange(a_indices.size), states_per_a)
            b_local_indices = np.tile(np.repeat(np.arange(b_indices.size), nphi), a_indices.size)
            relative_velocity = v_a - v_b
            relative_speed2 = np.einsum("ij,ij->i", relative_velocity, relative_velocity)
            relative_speed = np.sqrt(np.maximum(relative_speed2, 0.0))
            center_of_mass_energy = 0.5 * reduced_mass * relative_speed2
            sigma = _cross_section_array(cross_section_function_m2, center_of_mass_energy)
            sigma_g = sigma * relative_speed
            contributing = sigma_g > 0.0
            if not np.any(contributing):
                continue
            v_a = v_a[contributing]
            v_b = v_b[contributing]
            a_local_indices = a_local_indices[contributing]
            b_local_indices = b_local_indices[contributing]
            sigma_g = sigma_g[contributing]
            pair_prefactor = pair_factor * sigma_g / nphi

            def contract_pair_kernel(kernel_a_b: np.ndarray) -> np.ndarray:
                """Contract one pair kernel with the active reactant cell weights at every axial plane"""
                contracted_b = kernel_a_b @ weighted_b[:, b_indices].T
                return np.sum(weighted_a[:, a_indices].T * contracted_b, axis=0)

            pair_kernel = np.zeros((a_indices.size, b_indices.size), dtype=float)
            np.add.at(pair_kernel, (a_local_indices, b_local_indices), pair_prefactor)
            raw_totals += contract_pair_kernel(pair_kernel)
            if use_exact_energy_marginal:
                # Isotropic CM emission maps each reactant pair to a uniform lab energy interval
                energy_min, energy_max = (neutron_lab_kinetic_energy_endpoints_relativistic_batch(reaction_kinematics, v_a, v_b))
                interval_width = energy_max - energy_min
                continuous = interval_width > 0.0
                delta = ~continuous
                has_continuous = bool(np.any(continuous))
                delta_inside = delta & (energy_min >= energy_edges[0]) & (energy_min <= energy_edges[-1])
                has_delta_inside = bool(np.any(delta_inside))
                delta_bins = np.full(energy_min.shape, -1, dtype=int)
                if has_delta_inside:
                    delta_bins[delta_inside] = (np.searchsorted(energy_edges, energy_min[delta_inside], side="right",) - 1)
                    delta_bins[delta_inside & (energy_min == energy_edges[-1])] = n_energy - 1
                represented_fraction = np.zeros(energy_min.shape, dtype=float)
                pair_kernel_entries = a_indices.size * b_indices.size
                bins_per_block = max(1, min(n_energy, 4_000_000 // max(pair_kernel_entries, 1)))
                minimum_pair_energy = float(np.min(energy_min))
                first_candidate_bin = int(
                    np.clip(np.searchsorted(energy_edges, minimum_pair_energy, side="right")- 1, 0, n_energy,))
                if minimum_pair_energy == energy_edges[-1]:
                    first_candidate_bin = n_energy - 1
                candidate_stop = int(np.clip(np.searchsorted(energy_edges, float(np.max(energy_max)), side="right", ), 0, n_energy,))
                for energy_start in range(first_candidate_bin, candidate_stop, bins_per_block):
                    energy_stop = min(energy_start + bins_per_block, candidate_stop)
                    energy_pair_kernels = np.zeros((energy_stop - energy_start, a_indices.size, b_indices.size,), dtype=float)
                    for energy_bin in range(energy_start, energy_stop):
                        fraction = np.zeros(energy_min.shape, dtype=float)
                        if has_continuous:
                            overlap = np.minimum(energy_max[continuous], energy_edges[energy_bin + 1]) - np.maximum(energy_min[continuous], energy_edges[energy_bin])
                            fraction[continuous] = np.maximum(overlap, 0.0) / ( interval_width[continuous])
                        if has_delta_inside:
                            fraction[delta_inside & (delta_bins == energy_bin)] = 1.0
                        np.clip(fraction, 0.0, 1.0, out=fraction)
                        represented_fraction += fraction
                        contributing_to_bin = fraction > 0.0
                        if np.any(contributing_to_bin):
                            np.add.at(energy_pair_kernels[energy_bin - energy_start], (a_local_indices[contributing_to_bin], b_local_indices[contributing_to_bin]),pair_prefactor[contributing_to_bin] * fraction[contributing_to_bin])
                    contracted_b = (energy_pair_kernels.reshape((energy_stop - energy_start) * a_indices.size, b_indices.size) @ weighted_b[:, b_indices].T).reshape( energy_stop - energy_start, a_indices.size, n_z)
                    histograms[:, energy_start:energy_stop] += np.sum(contracted_b * weighted_a[:, a_indices].T[None, :, :], axis=1).T
                out_of_range_fraction = 1.0 - np.clip(represented_fraction, 0.0, 1.0)
                out_of_range_kernel = np.zeros_like(pair_kernel)
                np.add.at(out_of_range_kernel, (a_local_indices, b_local_indices), pair_prefactor * out_of_range_fraction)
                out_of_range_totals += contract_pair_kernel(out_of_range_kernel)
            else:
                # Deterministic CM directions approximate the energy marginal when exact interval integration is disabled
                energies, _ = neutron_lab_events_relativistic_batch(reaction_kinematics, v_a, v_b, emission_dirs, include_lab_directions=False)
                inside = (energies >= energy_edges[0]) & (energies < energy_edges[-1])
                inside |= energies == energy_edges[-1]
                out_of_range_fraction = 1.0 - np.mean(inside, axis=1)
                out_of_range_kernel = np.zeros_like(pair_kernel)
                np.add.at(out_of_range_kernel, (a_local_indices, b_local_indices), pair_prefactor * out_of_range_fraction)
                out_of_range_totals += contract_pair_kernel(out_of_range_kernel)
                bin_indices = (np.searchsorted(energy_edges, energies, side="right") - 1)
                bin_indices[energies == energy_edges[-1]] = n_energy - 1
                for emission_index in range(nemission):
                    valid = inside[:, emission_index]
                    if not np.any(valid):
                        continue
                    bins = bin_indices[:, emission_index]
                    for energy_bin in np.unique(bins[valid]):
                        selected = valid & (bins == energy_bin)
                        energy_kernel = np.zeros_like(pair_kernel)
                        np.add.at(energy_kernel, (a_local_indices[selected], b_local_indices[selected],), pair_prefactor[selected] / nemission)
                        histograms[:, int(energy_bin)] += contract_pair_kernel(energy_kernel)

    for z_index in range(n_z):
        threshold = max(1.0e-30, 1.0e-12 * max(raw_totals[z_index], 1.0))
        if out_of_range_totals[z_index] > threshold and not allow_out_of_range:
            raise NeutronEnergyDomainError("neutron energy_edges_J do not cover all weighted events" f"out of range rate density is {out_of_range_totals[z_index]:.6e} m^-3 s^-1")

    if normalize_to_rate_density_m3_s is not None:
        targets = np.asarray(normalize_to_rate_density_m3_s, dtype=float)
        if targets.ndim == 0:
            targets = np.full(n_z, float(targets), dtype=float)
        if targets.shape != (n_z,) or np.any(~np.isfinite(targets)) or np.any(targets < 0.0):
            raise ValueError("normalize_to_rate_density_m3_s must contain one finite nonnegative value per axial plane")
        for z_index in range(n_z):
            histograms[z_index], _ = _normalize_binned_rates(histograms[z_index], raw_totals[z_index], float(targets[z_index]))

    return histograms, raw_totals, out_of_range_totals

def _event_arrays_from_gyrotropic_distributions(speed_grid_a: SpeedGrid, pitch_grid_a: PitchGrid, distribution_a_v_xi: ArrayLike, speed_grid_b: SpeedGrid, pitch_grid_b: PitchGrid, distribution_b_v_xi: ArrayLike, reaction_kinematics: TwoBodyReactionKinematics, cross_section_function_m2: Callable[[ArrayLike], ArrayLike], *, identical_population: bool, num_gyroangle_points: int, num_emission_directions: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Expand one pair of local gyrotropic distributions into weighted deterministic neutron events
    
    The event path samples relative gyro angle and isotropic center of momentum direction and returns lab neutron energy, lab `μ = cos(θ)`, and rate density weight
    """
    _, _, weights_a, weights_b, nphi, emission_dirs = _validate_reactants(speed_grid_a, pitch_grid_a, distribution_a_v_xi, speed_grid_b, pitch_grid_b, distribution_b_v_xi, num_gyroangle_points, num_emission_directions)
    pair_factor = 0.5 if bool(identical_population) else 1.0
    mu_reduced = reduced_mass_kg(reaction_kinematics.reactant_a_mass_kg, reaction_kinematics.reactant_b_mass_kg)
    gyro_angles = (np.arange(nphi, dtype=float) + 0.5) * (2.0 * np.pi / nphi)
    energies: list[float] = []
    mu_values: list[float] = []
    weights: list[float] = []
    nemission = emission_dirs.shape[0]
    for i_v, v_a in enumerate(speed_grid_a.centers_m_s):
        for i_xi, xi_a in enumerate(pitch_grid_a.centers):
            weight_a = float(weights_a[i_v, i_xi])
            if weight_a == 0.0:
                continue
            vvec_a = velocity_vector_from_speed_pitch_gyro(float(v_a), float(xi_a), 0.0)
            for j_v, v_b in enumerate(speed_grid_b.centers_m_s):
                for j_xi, xi_b in enumerate(pitch_grid_b.centers):
                    base_pair_weight = weight_a * float(weights_b[j_v, j_xi])
                    if base_pair_weight == 0.0:
                        continue
                    for gyro_angle in gyro_angles:
                        vvec_b = velocity_vector_from_speed_pitch_gyro(float(v_b), float(xi_b), float(gyro_angle))
                        g_energy = center_of_mass_kinetic_energy_J(reaction_kinematics.reactant_a_mass_kg, vvec_a, reaction_kinematics.reactant_b_mass_kg, vvec_b)
                        sigma = float(np.asarray(cross_section_function_m2(g_energy), dtype=float))
                        if sigma < 0.0 or not np.isfinite(sigma):
                            raise ValueError("cross_section_function_m2 must return finite nonnegative values")
                        g = float(np.sqrt(max(2.0 * g_energy / mu_reduced, 0.0)))
                        rate_weight_for_pair = pair_factor * base_pair_weight * sigma * g / nphi
                        if rate_weight_for_pair == 0.0:
                            continue
                        event_weight = rate_weight_for_pair / nemission
                        event_energies, lab_directions = neutron_lab_events_relativistic_batch(reaction_kinematics, vvec_a, vvec_b, emission_dirs, include_lab_directions=True)
                        energies.extend(event_energies[0].tolist())
                        mu_values.extend(lab_directions[0, :, 2].tolist())
                        weights.extend(np.full(nemission, event_weight, dtype=float).tolist())
    if not weights:
        return np.zeros(0), np.zeros(0), np.zeros(0)
    
    return np.asarray(energies, dtype=float), np.asarray(mu_values, dtype=float), np.asarray(weights, dtype=float)

def neutron_energy_spectrum_from_gyrotropic_distributions(speed_grid_a: SpeedGrid, pitch_grid_a: PitchGrid, distribution_a_v_xi: ArrayLike, speed_grid_b: SpeedGrid, pitch_grid_b: PitchGrid, distribution_b_v_xi: ArrayLike, reaction_kinematics: TwoBodyReactionKinematics, cross_section_function_m2: Callable[[ArrayLike], ArrayLike], energy_edges_J: ArrayLike, *, identical_population: bool = False, num_gyroangle_points: int = 8, num_emission_directions: int = 24, normalize_to_rate_density_m3_s: float | None = None, allow_out_of_range: bool = False, exact_isotropic_energy_marginal: bool = False) -> NeutronEnergySpectrum:
    """Build one energy binned neutron spectrum from two local gyrotropic reactant distributions"""
    energy_edges = validate_bin_edges(energy_edges_J, "energy_edges_J")
    target = None if normalize_to_rate_density_m3_s is None else np.asarray([normalize_to_rate_density_m3_s], dtype=float)
    matrix, raw_totals, out_of_range_totals = neutron_energy_spectrum_matrix_from_gyrotropic_distributions(speed_grid_a, pitch_grid_a, distribution_a_v_xi, speed_grid_b, pitch_grid_b, distribution_b_v_xi, reaction_kinematics, cross_section_function_m2, energy_edges, identical_population=identical_population, num_gyroangle_points=num_gyroangle_points, num_emission_directions=num_emission_directions, normalize_to_rate_density_m3_s=target, allow_out_of_range=allow_out_of_range, exact_isotropic_energy_marginal=exact_isotropic_energy_marginal)
    total = float(raw_totals[0] if normalize_to_rate_density_m3_s is None else normalize_to_rate_density_m3_s)

    return NeutronEnergySpectrum(energy_edges_J=energy_edges, bin_rate_density_m3_s=matrix[0], total_rate_density_m3_s=total, raw_total_rate_density_m3_s=float(raw_totals[0]), out_of_range_rate_density_m3_s=float(out_of_range_totals[0]))

def neutron_angular_energy_spectrum_from_gyrotropic_distributions(speed_grid_a: SpeedGrid, pitch_grid_a: PitchGrid, distribution_a_v_xi: ArrayLike, speed_grid_b: SpeedGrid, pitch_grid_b: PitchGrid, distribution_b_v_xi: ArrayLike, reaction_kinematics: TwoBodyReactionKinematics, cross_section_function_m2: Callable[[ArrayLike], ArrayLike], energy_edges_J: ArrayLike, mu_edges: ArrayLike, *, identical_population: bool = False, num_gyroangle_points: int = 8, num_emission_directions: int = 24, normalize_to_rate_density_m3_s: float | None = None, allow_out_of_range: bool = False) -> NeutronAngularEnergySpectrum:
    """Build one joint energy and lab pitch cosine neutron spectrum from two local gyrotropic reactant distributions"""
    energy_edges = validate_bin_edges(energy_edges_J, "energy_edges_J")
    mu = validate_bin_edges(mu_edges, "mu_edges")
    if mu[0] < -1.0 or mu[-1] > 1.0:
        raise ValueError("mu_edges must lie inside [-1, 1]")
    energies, mu_values, weights = _event_arrays_from_gyrotropic_distributions(speed_grid_a, pitch_grid_a, distribution_a_v_xi, speed_grid_b, pitch_grid_b, distribution_b_v_xi, reaction_kinematics, cross_section_function_m2, identical_population=identical_population, num_gyroangle_points=num_gyroangle_points, num_emission_directions=num_emission_directions)

    return _histogram_energy_mu_events(energies, mu_values, weights, energy_edges, mu, allow_out_of_range=allow_out_of_range, normalize_to_rate_density_m3_s=normalize_to_rate_density_m3_s)
