"""
Axial wrappers for deterministic neutron spectrum construction

The helpers accept invariant `F(v, Λ)` or local gyrotropic `f(z, v, ξ)` reactant populations and return energy binned neutron rate density by axial cell
"""
from __future__ import annotations
from dataclasses import dataclass
from collections.abc import Callable
import numpy as np
from numpy.typing import ArrayLike
from source_model_revamp.fbis.mirror_invariant_grid import LambdaGrid, require_mirror_distribution_shape
from source_model_revamp.fbis.mirror_pitch_remap import conservative_local_pitch_distribution_from_lambda_distribution
from source_model_revamp.fbis.velocity_space_grid import PitchGrid, SpeedGrid
from source_model_revamp.fusion.cross_sections import bosch_hale_cross_section_m2_from_J
from source_model_revamp.fusion.reactions import FusionReaction
from source_model_revamp.neutrons.kinematics import TwoBodyReactionKinematics, neutron_kinematics_for_reaction
from source_model_revamp.neutrons.spectrum import NeutronEnergySpectrum, neutron_angular_energy_spectrum_from_gyrotropic_distributions, neutron_energy_spectrum_matrix_from_gyrotropic_distributions, validate_bin_edges

@dataclass(frozen=True)
class AxialNeutronSpectrumMatrix:
    """
    Axial and energy neutron source matrix
    
    `energy_bin_rate_density_m3_s` has shape `(n_z, n_energy)` and stores the rate density integrated over each energy bin
    Optional cell volumes convert rate density into physical rate in s⁻¹
    """
    coordinate: np.ndarray
    energy_edges_J: np.ndarray
    energy_bin_rate_density_m3_s: np.ndarray
    cell_volumes_m3: np.ndarray | None = None
    angular_model: str = "isotropic_cm"
    raw_total_rate_density_m3_s: np.ndarray | None = None
    out_of_range_rate_density_m3_s: np.ndarray | None = None
    posterior_rate_normalization_applied: bool = False

    @property
    def num_axial_bins(self) -> int:
        """Return the number of axial cells"""
        return int(self.energy_bin_rate_density_m3_s.shape[0])

    @property
    def num_energy_bins(self) -> int:
        """Return the number of neutron energy bins"""
        return int(self.energy_bin_rate_density_m3_s.shape[1])

    @property
    def energy_centers_J(self) -> np.ndarray:
        """Return neutron energy bin centers in J"""
        return 0.5 * (self.energy_edges_J[:-1] + self.energy_edges_J[1:])

    @property
    def axial_rate_density_m3_s(self) -> np.ndarray:
        """Return the energy integrated neutron rate density in each axial cell"""
        return np.sum(self.energy_bin_rate_density_m3_s, axis=1)

    @property
    def source_rate_matrix_s(self) -> np.ndarray | None:
        """Return the physical neutron rate in each axial and energy bin when cell volumes are available"""
        if self.cell_volumes_m3 is None:
            return None
        return self.energy_bin_rate_density_m3_s * self.cell_volumes_m3[:, None]

    @property
    def axial_bin_rates_s(self) -> np.ndarray | None:
        """Return the neutron rate in each axial cell when cell volumes are available"""
        matrix = self.source_rate_matrix_s
        if matrix is None:
            return None
        return np.sum(matrix, axis=1)

    @property
    def energy_bin_rates_s(self) -> np.ndarray | None:
        """Return the volume integrated neutron rate in each energy bin when cell volumes are available"""
        matrix = self.source_rate_matrix_s
        if matrix is None:
            return None
        return np.sum(matrix, axis=0)

    @property
    def total_source_rate_s(self) -> float | None:
        """Return the volume and energy integrated neutron rate when cell volumes are available"""
        matrix = self.source_rate_matrix_s
        if matrix is None:
            return None
        return float(np.sum(matrix))

    @property
    def raw_axial_bin_rates_s(self) -> np.ndarray | None:
        """Return raw quadrature rates integrated over each axial cell before optional target normalization"""
        if self.cell_volumes_m3 is None or self.raw_total_rate_density_m3_s is None:
            return None
        return (np.asarray(self.raw_total_rate_density_m3_s, dtype=float) * self.cell_volumes_m3)

    @property
    def raw_total_source_rate_s(self) -> float | None:
        """Return the volume integrated raw quadrature rate before optional target normalization"""
        rates = self.raw_axial_bin_rates_s
        if rates is None:
            return None
        return float(np.sum(rates))

    @property
    def out_of_range_axial_bin_rates_s(self) -> np.ndarray | None:
        """Return raw quadrature rate outside the configured energy range in each axial cell"""
        if self.cell_volumes_m3 is None or self.out_of_range_rate_density_m3_s is None:
            return None
        return (np.asarray(self.out_of_range_rate_density_m3_s, dtype=float) * self.cell_volumes_m3)

    @property
    def out_of_range_total_source_rate_s(self) -> float | None:
        """Return the volume integrated raw rate outside the configured energy range"""
        rates = self.out_of_range_axial_bin_rates_s
        if rates is None:
            return None
        return float(np.sum(rates))

    def global_energy_spectrum(self) -> NeutronEnergySpectrum:
        """Return the globally summed or volume integrated NeutronEnergySpectrum with raw and out of range diagnostics"""
        if (self.raw_total_rate_density_m3_s is None or self.out_of_range_rate_density_m3_s is None):
            raise ValueError("raw and out of range quadrature diagnostics are required for a global neutron spectrum")
        if self.source_rate_matrix_s is None:
            bin_rates = np.sum(self.energy_bin_rate_density_m3_s, axis=0)
            raw_total = float(np.sum(self.raw_total_rate_density_m3_s))
            out_of_range = float(np.sum(self.out_of_range_rate_density_m3_s))
        else:
            bin_rates = np.sum(self.source_rate_matrix_s, axis=0)
            raw_total = float(self.raw_total_source_rate_s)
            out_of_range = float(self.out_of_range_total_source_rate_s)
        total = float(np.sum(bin_rates))

        return NeutronEnergySpectrum(energy_edges_J=self.energy_edges_J, bin_rate_density_m3_s=bin_rates, total_rate_density_m3_s=total, raw_total_rate_density_m3_s=raw_total, out_of_range_rate_density_m3_s=out_of_range, angular_model=self.angular_model)

def _as_1d(values: ArrayLike, name: str) -> np.ndarray:
    """Validate one finite nonempty one dimensional array"""
    arr = np.asarray(values, dtype=float)
    if arr.ndim != 1:
        raise ValueError(f"{name} must be 1D")
    if arr.size == 0:
        raise ValueError(f"{name} must not be empty")
    if not np.all(np.isfinite(arr)):
        raise ValueError(f"{name} must contain only finite values")
    
    return arr

def _broadcast_profile(values: ArrayLike, size: int, name: str) -> np.ndarray:
    """Broadcast a scalar or validate a one dimensional profile with the requested size"""
    arr = np.asarray(values, dtype=float)
    if arr.ndim == 0:
        out = np.full(size, float(arr), dtype=float)
    else:
        out = _as_1d(arr, name)
        if out.size != size:
            raise ValueError(f"{name} must have length {size}, got {out.size}")
    if not np.all(np.isfinite(out)):
        raise ValueError(f"{name} must contain only finite values")
    
    return out

def _volumes_or_none(cell_volumes_m3: ArrayLike | None, size: int) -> np.ndarray | None:
    """Validate one positive cell volume per coordinate when volumes are supplied"""
    if cell_volumes_m3 is None:
        return None
    volumes = _as_1d(cell_volumes_m3, "cell_volumes_m3")
    if volumes.size != size:
        raise ValueError("cell_volumes_m3 must have one value per coordinate")
    if np.any(volumes <= 0.0):
        raise ValueError("cell_volumes_m3 must be positive")
    
    return volumes

def _target_profile_or_none(target_rate_density_m3_s: ArrayLike | None, size: int) -> np.ndarray | None:
    """Validate an optional nonnegative target rate density profile"""
    if target_rate_density_m3_s is None:
        return None
    target = _broadcast_profile(target_rate_density_m3_s, size, "target_rate_density_m3_s")
    if np.any(target < 0.0):
        raise ValueError("target_rate_density_m3_s must be nonnegative")
    
    return target

def _reaction_kinematics_and_cross_section(reaction: str | FusionReaction | TwoBodyReactionKinematics, cross_section_function_m2: Callable[[ArrayLike], ArrayLike] | None) -> tuple[TwoBodyReactionKinematics, Callable[[ArrayLike], ArrayLike]]:
    """Resolve neutron branch kinematics and use the supplied cross section or the Bosch and Hale total cross section"""
    if isinstance(reaction, TwoBodyReactionKinematics):
        kin = reaction
    else:
        key = reaction.key if isinstance(reaction, FusionReaction) else str(reaction)
        kin = neutron_kinematics_for_reaction(key)
    if cross_section_function_m2 is None:
        cross_section_function_m2 = lambda E_J: bosch_hale_cross_section_m2_from_J(kin.reaction.key, E_J)

    return kin, cross_section_function_m2

def local_lambda_distribution_neutron_spectrum_matrix(coordinate: ArrayLike, speed_grid_a: SpeedGrid, lambda_grid_a: LambdaGrid, distribution_a_v_lambda: ArrayLike, speed_grid_b: SpeedGrid, lambda_grid_b: LambdaGrid, distribution_b_v_lambda: ArrayLike, pitch_grid: PitchGrid, B_tilde_profile: ArrayLike, reaction: str | FusionReaction | TwoBodyReactionKinematics, energy_edges_J: ArrayLike, *, cross_section_function_m2: Callable[[ArrayLike], ArrayLike] | None = None, cell_volumes_m3: ArrayLike | None = None, target_rate_density_m3_s: ArrayLike | None = None, identical_population: bool = False, num_gyroangle_points: int = 8, num_emission_directions: int = 24, allow_out_of_range: bool = False, max_pair_states_per_batch: int = 131072, exact_isotropic_energy_marginal: bool = True) -> AxialNeutronSpectrumMatrix:
    """
    Build axial neutron spectra from two invariant `F(v, Λ)` distributions
    
    Each invariant distribution is conservatively mapped to local pitch `ξ` at the supplied `B_tilde_profile` before the shared gyrotropic pair spectrum kernel is evaluated
    """
    z = _as_1d(coordinate, "coordinate")
    B_profile = _broadcast_profile(B_tilde_profile, z.size, "B_tilde_profile")
    if np.any(B_profile <= 0.0):
        raise ValueError("B_tilde_profile must be positive")
    dist_a = require_mirror_distribution_shape(distribution_a_v_lambda, speed_grid_a, lambda_grid_a, name="distribution_a_v_lambda")
    dist_b = require_mirror_distribution_shape(distribution_b_v_lambda, speed_grid_b, lambda_grid_b, name="distribution_b_v_lambda")
    if np.any(dist_a < 0.0) or np.any(dist_b < 0.0):
        raise ValueError("Lambda space distributions must be nonnegative")
    energy_edges = validate_bin_edges(energy_edges_J, "energy_edges_J")
    volumes = _volumes_or_none(cell_volumes_m3, z.size)
    target = _target_profile_or_none(target_rate_density_m3_s, z.size)
    kin, sigma = _reaction_kinematics_and_cross_section(reaction, cross_section_function_m2)
    # Map the invariant Λ distributions to the local pitch grid at each B_tilde value
    local_a = np.stack([conservative_local_pitch_distribution_from_lambda_distribution(lambda_grid_a, dist_a, pitch_grid, float(B_tilde)) for B_tilde in B_profile])
    local_b = np.stack([conservative_local_pitch_distribution_from_lambda_distribution(lambda_grid_b, dist_b, pitch_grid, float(B_tilde)) for B_tilde in B_profile])
    matrix, raw_totals, out_of_range_totals = neutron_energy_spectrum_matrix_from_gyrotropic_distributions(speed_grid_a, pitch_grid, local_a, speed_grid_b, pitch_grid, local_b, kin, sigma, energy_edges, identical_population=identical_population, num_gyroangle_points=num_gyroangle_points, num_emission_directions=num_emission_directions, normalize_to_rate_density_m3_s=target, allow_out_of_range=allow_out_of_range, max_pair_states_per_batch=max_pair_states_per_batch, exact_isotropic_energy_marginal=exact_isotropic_energy_marginal)

    return AxialNeutronSpectrumMatrix(
        coordinate=z,
        energy_edges_J=energy_edges,
        energy_bin_rate_density_m3_s=matrix,
        cell_volumes_m3=volumes,
        raw_total_rate_density_m3_s=raw_totals,
        out_of_range_rate_density_m3_s=out_of_range_totals,
        posterior_rate_normalization_applied=target is not None,
    )

def _require_local_pitch_distribution(local_distribution_z_v_pitch: ArrayLike, speed_grid: SpeedGrid, pitch_grid: PitchGrid, name: str = "local_distribution_z_v_pitch") -> np.ndarray:
    """Validate a local gyrotropic distribution `f(z, v, ξ)` with shape `(n_z, n_speed, n_pitch)`"""
    arr = np.asarray(local_distribution_z_v_pitch, dtype=float)
    expected_tail = (speed_grid.centers_m_s.size, pitch_grid.centers.size)
    if arr.ndim != 3 or arr.shape[1:] != expected_tail:
        raise ValueError(f"{name} must have shape (n_z, n_speed, n_pitch); got {arr.shape}, expected (*, {expected_tail[0]}, {expected_tail[1]})")
    if not np.all(np.isfinite(arr)):
        raise ValueError(f"{name} must contain only finite values")
    if np.any(arr < 0.0):
        raise ValueError(f"{name} must be nonnegative")
    
    return arr

def local_gyrotropic_pair_neutron_spectrum_matrix(coordinate: ArrayLike, speed_grid_a: SpeedGrid, pitch_grid_a: PitchGrid, local_distribution_a_z_v_pitch: ArrayLike, speed_grid_b: SpeedGrid, pitch_grid_b: PitchGrid, local_distribution_b_z_v_pitch: ArrayLike, reaction: str | FusionReaction | TwoBodyReactionKinematics, energy_edges_J: ArrayLike, *, cross_section_function_m2: Callable[[ArrayLike], ArrayLike] | None = None, cell_volumes_m3: ArrayLike | None = None, target_rate_density_m3_s: ArrayLike | None = None, identical_population: bool = False, num_gyroangle_points: int = 8, num_emission_directions: int = 24, allow_out_of_range: bool = False, max_pair_states_per_batch: int = 131072, exact_isotropic_energy_marginal: bool = True) -> AxialNeutronSpectrumMatrix:
    """
    Build axial neutron spectra from two represented local gyrotropic populations
    
    Both local distributions must provide one `(v, ξ)` plane per axial coordinate
    """
    z = _as_1d(coordinate, "coordinate")
    local_a = _require_local_pitch_distribution(local_distribution_a_z_v_pitch, speed_grid_a, pitch_grid_a, "local_distribution_a_z_v_pitch")
    local_b = _require_local_pitch_distribution(local_distribution_b_z_v_pitch, speed_grid_b, pitch_grid_b, "local_distribution_b_z_v_pitch")
    if local_a.shape[0] != z.size or local_b.shape[0] != z.size:
        raise ValueError("local reactant distributions must have one z plane per coordinate")
    energy_edges = validate_bin_edges(energy_edges_J, "energy_edges_J")
    volumes = _volumes_or_none(cell_volumes_m3, z.size)
    target = _target_profile_or_none(target_rate_density_m3_s, z.size)
    kin, sigma = _reaction_kinematics_and_cross_section(reaction, cross_section_function_m2)
    matrix, raw_totals, out_of_range_totals = neutron_energy_spectrum_matrix_from_gyrotropic_distributions(speed_grid_a, pitch_grid_a, local_a, speed_grid_b, pitch_grid_b, local_b, kin, sigma, energy_edges, identical_population=identical_population, num_gyroangle_points=num_gyroangle_points, num_emission_directions=num_emission_directions, normalize_to_rate_density_m3_s=target, allow_out_of_range=allow_out_of_range, max_pair_states_per_batch=max_pair_states_per_batch, exact_isotropic_energy_marginal=exact_isotropic_energy_marginal)

    return AxialNeutronSpectrumMatrix(coordinate=z, energy_edges_J=energy_edges, energy_bin_rate_density_m3_s=matrix, cell_volumes_m3=volumes, raw_total_rate_density_m3_s=raw_totals, out_of_range_rate_density_m3_s=out_of_range_totals, posterior_rate_normalization_applied=target is not None)

def local_fast_pitch_self_neutron_spectrum_matrix(coordinate: ArrayLike, speed_grid: SpeedGrid, pitch_grid: PitchGrid, local_distribution_z_v_pitch: ArrayLike, reaction: str | FusionReaction | TwoBodyReactionKinematics, energy_edges_J: ArrayLike, *, cross_section_function_m2: Callable[[ArrayLike], ArrayLike] | None = None, cell_volumes_m3: ArrayLike | None = None, target_rate_density_m3_s: ArrayLike | None = None, identical_population: bool = True, num_gyroangle_points: int = 8, num_emission_directions: int = 24, allow_out_of_range: bool = False, max_pair_states_per_batch: int = 131072, exact_isotropic_energy_marginal: bool = True) -> AxialNeutronSpectrumMatrix:
    """
    Build axial neutron spectra for one local fast population reacting with itself
    
    The identical population factor is enabled by default so unordered reactant pairs receive the required `1 / 2` counting factor
    """
    z = _as_1d(coordinate, "coordinate")
    local = _require_local_pitch_distribution(local_distribution_z_v_pitch, speed_grid, pitch_grid)
    if local.shape[0] != z.size:
        raise ValueError("local_distribution_z_v_pitch must have one z plane per coordinate")
    energy_edges = validate_bin_edges(energy_edges_J, "energy_edges_J")
    volumes = _volumes_or_none(cell_volumes_m3, z.size)
    target = _target_profile_or_none(target_rate_density_m3_s, z.size)
    kin, sigma = _reaction_kinematics_and_cross_section(reaction, cross_section_function_m2)
    matrix, raw_totals, out_of_range_totals = neutron_energy_spectrum_matrix_from_gyrotropic_distributions(speed_grid, pitch_grid, local, speed_grid, pitch_grid, local, kin, sigma, energy_edges, identical_population=identical_population, num_gyroangle_points=num_gyroangle_points, num_emission_directions=num_emission_directions, normalize_to_rate_density_m3_s=target, allow_out_of_range=allow_out_of_range, max_pair_states_per_batch=max_pair_states_per_batch, exact_isotropic_energy_marginal=exact_isotropic_energy_marginal)
   
    return AxialNeutronSpectrumMatrix(
        coordinate=z,
        energy_edges_J=energy_edges,
        energy_bin_rate_density_m3_s=matrix,
        cell_volumes_m3=volumes,
        raw_total_rate_density_m3_s=raw_totals,
        out_of_range_rate_density_m3_s=out_of_range_totals,
        posterior_rate_normalization_applied=target is not None,
    )

