"""Construct Egedal Eq 61 through Eq 64 ion loss states

This module evaluates the production fixed magnetic boundary closure, calculates the throat parallel temperature state, matches each directed Eq 63 branch to its one end Eq 61 spectrum, and retains reference loss reconstructions for diagnostics
"""
from __future__ import annotations
from collections.abc import Mapping
from dataclasses import dataclass
import numpy as np
from source_model_revamp.fbis.collision_parameters import FBISCollisionParameterState
from source_model_revamp.fbis.mirror_invariant_grid import LambdaGrid
from source_model_revamp.fbis.modal.eq59.collision_operator_state import Eq59CollisionOperatorState
from source_model_revamp.fbis.modal.models import COLD_ION_EQ14, EGEDAL_2022_PUBLISHED_ELECTROSTATIC_FBIS_APPROXIMATION, HOT_ION_ROSENBLUTH_EQ59, MAGNETIC_ONLY_FBIS, ZERO_POTENTIAL_MAGNETIC_REFERENCE, electrostatic_feedback_model as canonical_electrostatic_feedback_model, velocity_solution_model as canonical_velocity_solution_model
from source_model_revamp.fbis.modal.utils import _EPS
from source_model_revamp.fbis.velocity_space_grid import PitchGrid, SpeedGrid

EGEDAL_HOT_FIXED_MAGNETIC_BOUNDARY = "egedal_hot_fixed_magnetic_boundary"
CALCULATED_THROAT_PARALLEL_TEMPERATURE = "calculated_throat_parallel_temperature"
BALDWIN_1972_THROAT_DENSITY_CLOSURE = "baldwin_1972_throat_density_closure"
BALDWIN_1972_A1 = 6.5
BALDWIN_1972_M_MINUS_HALF_FROM_A1 = BALDWIN_1972_A1 / (np.sqrt(2.0) * np.pi)

def _eq63_pitch_scattering_state(*, speed_grid: SpeedGrid, collision_state: FBISCollisionParameterState, collision_operator_state: Eq59CollisionOperatorState | None) -> tuple[np.ndarray, float]:
    """Return the speed dependent Eq 63 pitch scattering coefficient and slowing time
    
    The cold reference uses the reduced scalar coefficient while the Eq 59 path uses the solved ion pitch scattering array
    """
    v = np.asarray(speed_grid.centers_m_s, dtype=float)
    if collision_operator_state is None:
        return np.full_like(v, float(collision_state.beta_m) * float(collision_state.critical_velocity_m_s) ** 3), float(collision_state.spitzer_slowing_down_time_s)
    state_speed = np.asarray(collision_operator_state.speed_grid.centers_m_s, dtype=float)
    pitch = np.asarray(collision_operator_state.ion_pitch_scattering_velocity_cubed_m3_s3, dtype=float)
    tau_s = float(collision_operator_state.spitzer_slowing_down_time_s)
    if state_speed.shape != v.shape or not np.array_equal(state_speed, v) or pitch.shape != v.shape or np.any(~np.isfinite(pitch)) or np.any(pitch < 0.0):
        raise ValueError("Eq 63 collision operator pitch array must match the speed grid")
    if not np.isclose(tau_s, float(collision_state.spitzer_slowing_down_time_s), rtol=4.0e-15, atol=0.0):
        raise ValueError("Eq 63 collision operator slowing time must match the collision state")
  
    return pitch, tau_s

class ThroatParallelTemperatureUnavailableError(ValueError):
    """The represented throat state cannot define the selected T_L closure"""

@dataclass(frozen=True)
class ThroatParallelTemperatureResult:
    """Parallel temperature state for one directed throat branch
    
    Temperatures use J and velocity moments use SI units
    """
    side: str
    eq63_parallel_temperature_J: float
    pressure_parallel_temperature_J: float | None
    outgoing_number_measure: float
    outgoing_mean_parallel_velocity_m_s: float
    outgoing_mean_parallel_velocity_squared_m2_s2: float
    source_axial_cell_index: int
    calculation_model: str
    diagnostics: Mapping[str, object] | None = None

@dataclass(frozen=True)
class FixedBoundaryEq63Result:
    """Directed Eq 63 branch matched independently to one end of the Eq 61 speed spectrum
    
    Throat arrays have shape (n_speed, n_lambda)
    """
    side: str
    parallel_temperature_J: float
    H_U: np.ndarray
    throat_distribution_v_lambda: np.ndarray
    throat_rate_v_lambda_s: np.ndarray
    target_rate_v_s: np.ndarray
    reconstructed_rate_v_s: np.ndarray
    target_total_rate_s: float
    reconstructed_total_rate_s: float
    relative_rate_error: float
    forced_posterior_normalization_applied: bool
    matching_model: str

def ion_loss_model(value: str | None) -> str:
    """Return the canonical production ion loss closure identifier"""
    model = EGEDAL_HOT_FIXED_MAGNETIC_BOUNDARY if value is None else str(value).strip().lower()
    if model != EGEDAL_HOT_FIXED_MAGNETIC_BOUNDARY:
        raise ValueError(f"Unsupported ion loss closure model {value!r}")
    
    return model

def _weighted_quantile(values: np.ndarray, weights: np.ndarray, quantile: float) -> float | None:
    """Return a weighted quantile of finite samples with positive weights"""
    x = np.asarray(values, dtype=float)
    w = np.asarray(weights, dtype=float)
    mask = np.isfinite(x) & np.isfinite(w) & (w > 0.0)
    if not np.any(mask):
        return None
    x = x[mask]
    w = w[mask]
    order = np.argsort(x)
    x = x[order]
    w = w[order]
    cumulative = np.cumsum(w)
    target = float(quantile) * float(cumulative[-1])
    index = min(int(np.searchsorted(cumulative, target, side="left")), x.size - 1)
  
    return float(x[index])

def calculate_baldwin_throat_parallel_temperature( *, side: str, speed_grid: SpeedGrid, one_end_eq61_rate_v_s: np.ndarray, dfdlambda_boundary_v: np.ndarray, mirror_ratio: float, midplane_area_m2: float, half_length_m: float, geometry_factor_G: float, particle_mass_kg: float, collision_state: FBISCollisionParameterState, collision_operator_state: Eq59CollisionOperatorState | None, volume_geometry_integral: float | None = None, volume_m3: float | None = None) -> ThroatParallelTemperatureResult:
    """Calculate the production Eq 63 parallel temperature from the one end Eq 61 rate and encoded throat density closure
    
    The result is inferred independently for the requested left or right branch and retains identity diagnostics against the Eq 61 differential rate
    """
    branch = str(side).strip().lower()
    if branch not in {"left", "right"}:
        raise ValueError("side must be left or right")
    v = np.asarray(speed_grid.centers_m_s, dtype=float)
    shell = np.asarray(speed_grid.shell_volumes_m3_s3, dtype=float)
    rate = np.asarray(one_end_eq61_rate_v_s, dtype=float)
    slope = np.asarray(dfdlambda_boundary_v, dtype=float)
    if rate.shape != v.shape or slope.shape != v.shape:
        raise ValueError("Baldwin throat closure arrays must match the invariant speed grid")
    if np.any(~np.isfinite(rate)) or np.any(rate < 0.0):
        raise ValueError("Eq 61 one end rate must be finite and nonnegative")
    if np.any(~np.isfinite(slope)) or np.any(slope < 0.0):
        raise ValueError("Eq 61 boundary slope must be finite and nonnegative")
    R = float(mirror_ratio)
    area0 = float(midplane_area_m2)
    length = float(half_length_m)
    G = float(geometry_factor_G)
    mass = float(particle_mass_kg)
    if not np.isfinite(R) or R <= 1.0:
        raise ValueError("mirror_ratio must be finite and greater than one")
    if not np.isfinite(area0) or area0 <= 0.0:
        raise ValueError("midplane_area_m2 must be positive and finite")
    if not np.isfinite(length) or length <= 0.0:
        raise ValueError("half_length_m must be positive and finite")
    if not np.isfinite(G) or G <= 0.0:
        raise ValueError("geometry_factor_G must be positive and finite")
    if not np.isfinite(mass) or mass <= 0.0:
        raise ValueError("particle_mass_kg must be positive and finite")
    total_rate = float(np.sum(rate))
    if not np.isfinite(total_rate) or total_rate <= 0.0:
        raise ThroatParallelTemperatureUnavailableError(f"{branch} Eq 61 one end rate is not positive")
    mu_c = float(np.sqrt(1.0 - 1.0 / R))
    pitch, tau_s = _eq63_pitch_scattering_state(speed_grid=speed_grid, collision_state=collision_state, collision_operator_state=collision_operator_state)
    gamma_tilde = np.divide(pitch, tau_s * v**2, out=np.zeros_like(v), where=v > 0.0)
    lambda_tilde = np.divide(gamma_tilde * length * G, v**2 * mu_c**2 * R, out=np.zeros_like(v), where=v > 0.0)
    if np.any(~np.isfinite(lambda_tilde)) or np.any(lambda_tilde < 0.0):
        raise ThroatParallelTemperatureUnavailableError(f"{branch} Baldwin boundary layer width is invalid")
    active = rate > 0.0
    if np.any(active & (lambda_tilde <= 0.0)):
        raise ThroatParallelTemperatureUnavailableError(f"{branch} positive Eq 61 rate has zero Baldwin boundary layer width")
    mean_v_parallel = np.zeros_like(v)
    mean_v_parallel[active] = ( np.sqrt(2.0) * v[active] * np.sqrt(R * mu_c) * np.power(lambda_tilde[active], 0.25) / float(BALDWIN_1972_M_MINUS_HALF_FROM_A1))
    throat_area = area0 / R
    density_by_speed = np.divide(rate, throat_area * mean_v_parallel, out=np.zeros_like(rate), where=mean_v_parallel > 0.0)
    total_density = float(np.sum(density_by_speed))
    if not np.isfinite(total_density) or total_density <= 0.0:
        raise ThroatParallelTemperatureUnavailableError(f"{branch} Baldwin throat density is not positive")
    mean_parallel_total = total_rate / (throat_area * total_density)
    T_L = 0.5 * np.pi * mass * mean_parallel_total**2
    if not np.isfinite(T_L) or T_L <= 0.0:
        raise ThroatParallelTemperatureUnavailableError(f"{branch} Baldwin throat density closure did not produce a positive T_L")
    gaussian_second_moment = T_L / mass
    G0 = 2.0 * mu_c * slope
    egedal_differential = 4.0 * np.pi * area0 * gamma_tilde * length * (1.0 / R) * G * v * slope
    baldwin_differential = 2.0 * np.pi * area0 * lambda_tilde * G0 * mu_c * v**3
    differential_scale = np.maximum(np.abs(egedal_differential), np.finfo(float).tiny)
    differential_identity = float(np.max(np.abs(baldwin_differential - egedal_differential) / differential_scale))
    one_end_rate_from_mapping = None
    one_end_rate_mapping_error = None
    if volume_geometry_integral is not None and volume_m3 is not None:
        volume_geometry = float(volume_geometry_integral)
        volume = float(volume_m3)
        if np.isfinite(volume_geometry) and volume_geometry > 0.0 and np.isfinite(volume) and volume > 0.0:
            gamma_from_lambda = lambda_tilde * v**2 * mu_c**2 * R / (length * G)
            flux_density = np.divide(gamma_from_lambda * slope * (1.0 / R) * G, volume_geometry * v, out=np.zeros_like(v), where=v > 0.0)
            mapped_rate = flux_density * shell * volume
            one_end_rate_from_mapping = float(np.sum(mapped_rate))
            one_end_rate_mapping_error = abs(one_end_rate_from_mapping - total_rate) / max(total_rate, np.finfo(float).tiny)
    sqrt_lambda = np.sqrt(lambda_tilde)
    diagnostics: dict[str, object] = {
        "reference": "Baldwin_et_al_Nuclear_Fusion_12_307_1972",
        "closure_equations": "Baldwin_Eq66_Eq68_Eq69_mapped_to_Egedal_Eq61_Eq63",
        "A1_published": float(BALDWIN_1972_A1),
        "mirror_loss_pitch_cosine_mu_c": mu_c,
        "Baldwin_total_throat_density_m3": total_density,
        "Baldwin_total_mean_parallel_velocity_m_s": mean_parallel_total,
        "Baldwin_rate_weighted_sqrt_lambda_q10": _weighted_quantile(sqrt_lambda, rate, 0.10),
        "Baldwin_rate_weighted_sqrt_lambda_q50": _weighted_quantile(sqrt_lambda, rate, 0.50),
        "Baldwin_rate_weighted_sqrt_lambda_q90": _weighted_quantile(sqrt_lambda, rate, 0.90),
        "Baldwin_maximum_sqrt_lambda": float(np.max(sqrt_lambda[active])) if np.any(active) else 0.0,
        "Baldwin_Eq68_vs_Egedal_Eq61_maximum_differential_relative_error": differential_identity,
        "Baldwin_Eq61_one_end_rate_from_lambda_mapping_s": one_end_rate_from_mapping,
        "Baldwin_Eq61_one_end_rate_mapping_relative_error": one_end_rate_mapping_error,
        "Baldwin_full_anisotropic_collision_operator_claimed": False,
        "Baldwin_fixed_magnetic_loss_boundary_used": True,
        "Baldwin_electrostatic_moving_loss_boundary_used": False,
        "Baldwin_low_energy_asymptotic_limitation_retained": True,
    }

    return ThroatParallelTemperatureResult(
        side=branch,
        eq63_parallel_temperature_J=float(T_L),
        pressure_parallel_temperature_J=None,
        outgoing_number_measure=total_density,
        outgoing_mean_parallel_velocity_m_s=float(mean_parallel_total),
        outgoing_mean_parallel_velocity_squared_m2_s2=float(gaussian_second_moment),
        source_axial_cell_index=-1,
        calculation_model=BALDWIN_1972_THROAT_DENSITY_CLOSURE,
        diagnostics=diagnostics,
    )

def _integrated_eq63_shape_over_x(*, total_energy_J: np.ndarray, x_lower: np.ndarray, x_upper: np.ndarray, parallel_temperature_J: float) -> np.ndarray:
    """Integrate the Eq 63 exponential exactly over x equals one minus mirror ratio times Λ cell intervals"""
    U = np.asarray(total_energy_J, dtype=float)[:, None]
    lo = np.asarray(x_lower, dtype=float)[None, :]
    hi = np.asarray(x_upper, dtype=float)[None, :]
    width = np.maximum(hi - lo, 0.0)
    out = np.zeros((U.shape[0], lo.shape[1]), dtype=float)
    active = width > 0.0
    if not np.any(active):
        return out
    T = float(parallel_temperature_J)
    a = U / T
    small = np.abs(a) < 1.0e-10
    exponential_difference = np.exp(np.clip(-a * lo, -700.0, 0.0)) - np.exp(np.clip(-a * hi, -700.0, 0.0))
    regular_value = np.divide(exponential_difference, a, out=np.zeros_like(exponential_difference), where=~small)
    out = np.where(small, width, regular_value)
    out[:, ~active.ravel()] = 0.0
   
    return np.maximum(out, 0.0)

def build_fixed_boundary_eq63_branch(*, side: str, speed_grid: SpeedGrid, lambda_grid: LambdaGrid, one_end_eq61_rate_v_s: np.ndarray, mirror_ratio: float, midplane_area_m2: float, geometry_factor_G: float, parallel_temperature_J: float, particle_mass_kg: float) -> FixedBoundaryEq63Result:
    """Build one directed Eq 63 throat branch by exact Eq 64 flux matching
    
    The analytic loss cone cell integral determines H(U) independently for each invariant speed cell so the reconstructed one end speed spectrum matches Eq 61 without posterior normalization
    """
    branch = str(side).strip().lower()
    if branch not in {"left", "right"}:
        raise ValueError("side must be left or right")
    rate_v = np.asarray(one_end_eq61_rate_v_s, dtype=float)
    speed_faces = np.asarray(speed_grid.faces_m_s, dtype=float)
    speed_centers = np.asarray(speed_grid.centers_m_s, dtype=float)
    lambda_faces = np.asarray(lambda_grid.faces, dtype=float)
    if rate_v.shape != speed_centers.shape:
        raise ValueError("one_end_eq61_rate_v_s must match the invariant speed grid")
    if np.any(~np.isfinite(rate_v)) or np.any(rate_v < 0.0):
        raise ValueError("one_end_eq61_rate_v_s must be finite and nonnegative")
    R = float(mirror_ratio)
    area = float(midplane_area_m2) / R
    G = float(geometry_factor_G)
    T_L = float(parallel_temperature_J)
    mass = float(particle_mass_kg)
    if not np.isfinite(R) or R <= 1.0:
        raise ValueError("mirror_ratio must be finite and greater than one")
    if not np.isfinite(area) or area <= 0.0:
        raise ValueError("midplane_area_m2 must define a positive throat area")
    if not np.isfinite(G) or G <= 0.0:
        raise ValueError("geometry_factor_G must be positive and finite")
    if not np.isfinite(T_L) or T_L <= 0.0:
        raise ValueError("parallel_temperature_J must be positive and finite")
    if not np.isfinite(mass) or mass <= 0.0:
        raise ValueError("particle_mass_kg must be positive and finite")
    open_lower = np.maximum(lambda_faces[:-1], 0.0)
    open_upper = np.minimum(lambda_faces[1:], 1.0 / R)
    open_width = np.maximum(open_upper - open_lower, 0.0)
    x_lower = np.maximum(1.0 - R * open_upper, 0.0)
    x_upper = np.maximum(1.0 - R * open_lower, 0.0)
    delta_x = np.maximum(x_upper - x_lower, 0.0)
    kinetic = 0.5 * mass * speed_centers**2
    integrated_shape = _integrated_eq63_shape_over_x(total_energy_J=kinetic, x_lower=x_lower, x_upper=x_upper, parallel_temperature_J=T_L)
    speed_flux_measure = 0.25 * np.pi * (speed_faces[1:] ** 4 - speed_faces[:-1] ** 4)
    # Fix H times G analytically from the one end Eq 61 rate in each invariant speed cell
    denominator = speed_flux_measure * area * np.sum(integrated_shape, axis=1)
    H_times_G = np.divide(rate_v, denominator, out=np.zeros_like(rate_v), where=denominator > 0.0)
    rate_matrix = H_times_G[:, None] * speed_flux_measure[:, None] * area * integrated_shape
    cell_average_shape = np.divide(integrated_shape, delta_x[None, :], out=np.zeros_like(integrated_shape), where=delta_x[None, :] > 0.0)
    throat_distribution = H_times_G[:, None] * cell_average_shape
    throat_distribution[:, open_width <= 0.0] = 0.0
    reconstructed_rate_v = np.sum(rate_matrix, axis=1)
    target_total = float(np.sum(rate_v))
    reconstructed_total = float(np.sum(reconstructed_rate_v))
    relative_error = abs(reconstructed_total - target_total) / max(abs(target_total), np.finfo(float).tiny)
    speed_error = float(np.sum(np.abs(reconstructed_rate_v - rate_v)) / max(target_total, np.finfo(float).tiny))
    if speed_error > 4096.0 * np.finfo(float).eps:
        raise ValueError("continuous Eq 64 flux matching failed to preserve the Eq 61 spectrum")
    H_U = H_times_G / G
  
    return FixedBoundaryEq63Result(
        side=branch,
        parallel_temperature_J=T_L,
        H_U=H_U,
        throat_distribution_v_lambda=np.maximum(throat_distribution, 0.0),
        throat_rate_v_lambda_s=np.maximum(rate_matrix, 0.0),
        target_rate_v_s=rate_v,
        reconstructed_rate_v_s=reconstructed_rate_v,
        target_total_rate_s=target_total,
        reconstructed_total_rate_s=reconstructed_total,
        relative_rate_error=float(relative_error),
        forced_posterior_normalization_applied=False,
        matching_model="analytic_continuous_Eq64_loss_cone_flux_integral_per_invariant_speed_cell",
    )

def _cumulative_G_profile(zeta_faces: np.ndarray, B_tilde_midpoints: np.ndarray, mirror_ratio: float) -> np.ndarray:
    """Return midpoint cumulative magnetic Eq 62 geometry values from the left boundary"""
    widths = np.diff(np.asarray(zeta_faces, dtype=float))
    B = np.asarray(B_tilde_midpoints, dtype=float)
    integrand = (1.0 / B) * np.sqrt(np.maximum(1.0 - B / float(mirror_ratio), 0.0))
    cumulative_faces = np.concatenate([[0.0], np.cumsum(widths * integrand)])
   
    return 0.5 * (cumulative_faces[:-1] + cumulative_faces[1:])

def _lost_ion_distribution_eq63(*, speed_grid: SpeedGrid, lambda_grid: LambdaGrid, dfdlambda_boundary_v: np.ndarray, zeta_faces: np.ndarray, B_tilde_midpoints: np.ndarray, mirror_ratio: float, collision_state: FBISCollisionParameterState, half_length_m: float, lost_ion_parallel_temperature_J: float, particle_mass_kg: float, collision_operator_state: Eq59CollisionOperatorState | None = None) -> tuple[np.ndarray, np.ndarray]:
    """Evaluate the magnetic reference Eq 63 distribution on the axial, speed, and Λ grid"""
    Gz = _cumulative_G_profile(zeta_faces, B_tilde_midpoints, mirror_ratio)
    lost = _lost_ion_distribution_eq63_for_G(speed_grid=speed_grid, lambda_grid=lambda_grid, dfdlambda_boundary_v=dfdlambda_boundary_v, mirror_ratio=mirror_ratio, collision_state=collision_state, half_length_m=half_length_m, lost_ion_parallel_temperature_J=lost_ion_parallel_temperature_J, geometry_G_values=Gz, particle_mass_kg=particle_mass_kg, collision_operator_state=collision_operator_state)
  
    return lost, Gz

def _lost_ion_distribution_eq63_for_G(*, speed_grid: SpeedGrid, lambda_grid: LambdaGrid, dfdlambda_boundary_v: np.ndarray, mirror_ratio: float, collision_state: FBISCollisionParameterState, half_length_m: float, lost_ion_parallel_temperature_J: float, geometry_G_values: np.ndarray, particle_mass_kg: float, collision_operator_state: Eq59CollisionOperatorState | None = None) -> np.ndarray:
    """Evaluate the reference Eq 63 distribution at explicit geometry factor values"""
    T_L = max(float(lost_ion_parallel_temperature_J), _EPS)
    v = np.asarray(speed_grid.centers_m_s, dtype=float)
    Lambda = np.asarray(lambda_grid.centers, dtype=float)
    faces = np.asarray(lambda_grid.faces, dtype=float)
    R = float(mirror_ratio)
    vLt2 = 2.0 * T_L / float(particle_mass_kg)
    pitch, tau_s = _eq63_pitch_scattering_state(speed_grid=speed_grid, collision_state=collision_state, collision_operator_state=collision_operator_state)
    H_U = 4.0 * pitch
    H_U /= tau_s * np.maximum(v**2, _EPS) * vLt2
    H_U *= float(half_length_m) * np.maximum(dfdlambda_boundary_v, 0.0)
    Gz = np.asarray(geometry_G_values, dtype=float)
    if Gz.ndim != 1 or Gz.size < 1 or np.any(~np.isfinite(Gz)) or np.any(Gz < 0.0):
        raise ValueError("geometry_G_values must be a nonnegative finite vector")
    lost = np.zeros((Gz.size, v.size, Lambda.size), dtype=float)
    lambda_loss_edge = 1.0 / R
    kinetic = 0.5 * float(particle_mass_kg) * v**2
    for iz, G in enumerate(Gz):
        for il in range(Lambda.size):
            lo = max(float(faces[il]), 0.0)
            hi = min(float(faces[il + 1]), lambda_loss_edge)
            if hi <= lo:
                continue
            lam_eff = 0.5 * (lo + hi)
            exponent = -kinetic * max(1.0 - lam_eff * R, 0.0) / T_L
            lost[iz, :, il] = H_U * np.exp(np.clip(exponent, -700.0, 0.0)) * float(G)
  
    return np.maximum(lost, 0.0)

def _lost_ion_distribution_eq63_exact_right_throat(*, speed_grid: SpeedGrid, lambda_grid: LambdaGrid, dfdlambda_boundary_v: np.ndarray, zeta_faces: np.ndarray, B_tilde_midpoints: np.ndarray, mirror_ratio: float, collision_state: FBISCollisionParameterState, half_length_m: float, lost_ion_parallel_temperature_J: float, particle_mass_kg: float, collision_operator_state: Eq59CollisionOperatorState | None = None) -> np.ndarray:
    """Return the reference Eq 63 distribution evaluated at the exact right throat geometry factor"""
    faces = np.asarray(zeta_faces, dtype=float)
    field = np.asarray(B_tilde_midpoints, dtype=float)
    if faces.ndim != 1 or field.ndim != 1 or faces.size != field.size + 1:
        raise ValueError("zeta_faces and B_tilde_midpoints have inconsistent shapes")
    widths = np.diff(faces)
    G_right_face = float(np.sum(widths * (1.0 / field) * np.sqrt(np.maximum(1.0 - field / float(mirror_ratio), 0.0))))
  
    return _lost_ion_distribution_eq63_for_G(
        speed_grid=speed_grid,
        lambda_grid=lambda_grid,
        dfdlambda_boundary_v=dfdlambda_boundary_v,
        mirror_ratio=mirror_ratio,
        collision_state=collision_state,
        half_length_m=half_length_m,
        lost_ion_parallel_temperature_J=lost_ion_parallel_temperature_J,
        geometry_G_values=np.asarray([G_right_face], dtype=float),
        particle_mass_kg=particle_mass_kg,
        collision_operator_state=collision_operator_state,
    )[0]
