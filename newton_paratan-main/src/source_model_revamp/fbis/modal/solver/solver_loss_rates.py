"""Eq 61 boundary loss reconstruction and modal sink diagnostics"""
from __future__ import annotations
import json
import numpy as np
from source_model_revamp.fbis.collision_parameters import FBISCollisionParameterState
from source_model_revamp.fbis.mirror_invariant_grid import LambdaGrid
from source_model_revamp.fbis.velocity_space_grid import SpeedGrid
from source_model_revamp.fbis.modal.lost_ion_distribution import EGEDAL_HOT_FIXED_MAGNETIC_BOUNDARY
from source_model_revamp.fbis.modal.types import ModalPhysicalDistributionError
from source_model_revamp.fbis.modal.eq59.collision_operator_state import Eq59CollisionOperatorState
from source_model_revamp.fbis.modal.utils import _EPS
from source_model_revamp.fbis.modal.types import ModalFBISBasis

def _raw_directed_eq63_throat_rate(*, speed_grid: SpeedGrid, lambda_grid: LambdaGrid, lost_ion_distribution_z_v_lambda: np.ndarray, mirror_ratio: float, midplane_area_m2: float) -> np.ndarray:
    """
    Integrate the outbound F_L+ cell flux at one magnetic throat

    Eq 63 is the directed + branch, treating the stored cell value as a piecewise constant phase space density, its exact one sign cell flux is
        (π / 4) * (v_hi^4 - v_lo^4) * R_M * Δ_Λ_open

    This follows from integrating v_parallel d^3v with ξ = sqrt(1 - R_M * Λ)
    The square root singularity in the Λ Jacobian cancels the parallel speed factor
    """
    lost = np.asarray(lost_ion_distribution_z_v_lambda, dtype=float)
    expected = (np.asarray(speed_grid.centers_m_s, dtype=float).size, np.asarray(lambda_grid.centers, dtype=float).size,)
    if lost.ndim == 2 and lost.shape == expected:
        throat_distribution = lost
    elif lost.ndim == 3 and lost.shape[1:] == expected and lost.shape[0] >= 1:
        throat_distribution = lost[-1]
    else:
        raise ValueError("Eq. 63 lost distribution has inconsistent shape")
    R = float(mirror_ratio)
    area = float(midplane_area_m2) / R
    if not np.isfinite(R) or R <= 1.0 or not np.isfinite(area) or area <= 0.0:
        raise ValueError("Eq. 63 throat integration requires physical geometry")
    if np.any(~np.isfinite(throat_distribution)):
        raise ValueError("Eq. 63 throat distribution must be finite")
    distribution_scale = max(float(np.max(np.abs(throat_distribution))) if throat_distribution.size else 0.0, np.finfo(float).tiny)
    negative_tolerance = 128.0 * np.finfo(float).eps * distribution_scale
    if float(np.min(throat_distribution)) < -negative_tolerance:
        raise ValueError("Eq 63 throat distribution contains material negative values")
    speed_faces = np.asarray(speed_grid.faces_m_s, dtype=float)
    lambda_faces = np.asarray(lambda_grid.faces, dtype=float)
    speed_flux_measure = (0.25 * np.pi * (speed_faces[1:] ** 4 - speed_faces[:-1] ** 4))
    open_lambda_lower = np.maximum(lambda_faces[:-1], 0.0)
    open_lambda_upper = np.minimum(lambda_faces[1:], 1.0 / R)
    open_lambda_width = np.maximum(open_lambda_upper - open_lambda_lower, 0.0)
    one_sign_flux_measure = (speed_flux_measure[:, None] * (R * open_lambda_width)[None, :])

    return (np.maximum(throat_distribution, 0.0) * one_sign_flux_measure * area)

def _modal_mode_eta_integrals(basis: ModalFBISBasis) -> np.ndarray:
    """Integrate each physical I_j(η) over the confined η interval used by the modal density measure"""
    physical = basis.physical_basis

    return np.trapezoid(np.asarray(physical.eigenfunctions, dtype=float), np.asarray(physical.eta_grid, dtype=float), axis=1)

def _conservative_eq61_mode_slopes_dI_dlambda(basis: ModalFBISBasis) -> np.ndarray:
    """Return weak form Eq 61 dI_j/dΛ slopes consistent with the resolved modal eigenvalue sink"""
    geometry_factor = float(basis.geometry_factor_G)
    volume_geometry = float(basis.volume_geometry_integral)
    mirror_ratio = 1.0 / float(basis.lambda_boundary)
    if not np.isfinite(geometry_factor) or geometry_factor <= 0.0:
        raise ValueError("Eq 61 geometry factor must be positive and finite")
    if not np.isfinite(volume_geometry) or volume_geometry <= 0.0:
        raise ValueError("Eq 61 volume geometry integral must be positive and finite")
    eigenvalues = np.asarray(basis.physical_basis.eigenvalues, dtype=float)
    eta_integrals = _modal_mode_eta_integrals(basis)
    slopes = (eigenvalues * eta_integrals * mirror_ratio * volume_geometry / (2.0 * geometry_factor))
    if np.any(~np.isfinite(slopes)):
        raise ValueError("weak form Eq 61 boundary slopes must be finite")

    return slopes



def _eq61_boundary_reconstruction_diagnostics(*, speed_grid: SpeedGrid, basis: ModalFBISBasis, modal_f_j_v: np.ndarray, collision_state: FBISCollisionParameterState, fast_ion_mass_kg: float, collision_operator_state: Eq59CollisionOperatorState | None = None, species_id: str | None = None) -> dict[str, object]:
    """Compare weak form and cell average Eq 61 boundary reconstructions including sign and integrated flux metrics"""
    modal = np.asarray(modal_f_j_v, dtype=float)
    weak_slopes = _conservative_eq61_mode_slopes_dI_dlambda(basis)
    direct_slopes = np.asarray(basis.mode_slopes_dI_dlambda_at_boundary, dtype=float)
    if modal.ndim != 2 or modal.shape[0] != weak_slopes.size or direct_slopes.shape != weak_slopes.shape:
        raise ValueError("Eq 61 diagnostics require modal coefficients matching both boundary slope representations")
    # Weak form slopes preserve the resolved modal sink while cell average slopes diagnose boundary resolution
    weak_raw = weak_slopes @ modal
    direct_raw = direct_slopes @ modal
    speed = np.asarray(speed_grid.centers_m_s, dtype=float)
    if weak_raw.shape != speed.shape or direct_raw.shape != speed.shape:
        raise ValueError("Eq 61 diagnostics require one reconstructed boundary slope per speed cell")

    def sign_metrics(values: np.ndarray) -> dict[str, object]:
        """Summarize negative reconstructed slopes relative to floating point roundoff"""
        scale = max(float(np.max(np.abs(values))) if values.size else 0.0, np.finfo(float).tiny)
        tolerance = 128.0 * np.finfo(float).eps * scale
        material_negative = values < -tolerance
        negative = values < 0.0
        minimum_index = int(np.argmin(values)) if values.size else None
        return {
            "minimum": float(np.min(values)) if values.size else 0.0,
            "maximum_absolute": scale,
            "minimum_over_maximum_absolute": float(np.min(values) / scale) if values.size else 0.0,
            "roundoff_negative_tolerance": float(tolerance),
            "negative_cell_count": int(np.count_nonzero(negative)),
            "material_negative_cell_count": int(np.count_nonzero(material_negative)),
            "material_negative_cell_fraction": float(np.count_nonzero(material_negative) / values.size) if values.size else 0.0,
            "most_negative_speed_index": minimum_index,
            "most_negative_speed_m_s": None if minimum_index is None else float(speed[minimum_index]),
        }

    weak_metrics = sign_metrics(weak_raw)
    direct_metrics = sign_metrics(direct_raw)
    slope_scale = max(float(np.max(np.abs(direct_slopes))) if direct_slopes.size else 0.0, np.finfo(float).tiny)
    reconstruction_scale = max(float(np.linalg.norm(direct_raw)), np.finfo(float).tiny)
    diagnostics: dict[str, object] = {
        "schema": "eq61_boundary_reconstruction_diagnostic_v1",
        "species_id": None if species_id is None else str(species_id),
        "n_modes": int(modal.shape[0]),
        "n_speed_cells": int(modal.shape[1]),
        "n_eta_grid": int(np.asarray(basis.physical_basis.eta_grid, dtype=float).size),
        "weak_form": weak_metrics,
        "direct_cell_average": direct_metrics,
        "mode_slope_maximum_relative_difference": float(np.max(np.abs(weak_slopes - direct_slopes)) / slope_scale) if weak_slopes.size else 0.0,
        "boundary_reconstruction_l2_relative_difference": float(np.linalg.norm(weak_raw - direct_raw) / reconstruction_scale),
        "boundary_reconstruction_sign_disagreement_cell_count": int(np.count_nonzero(np.signbit(weak_raw) != np.signbit(direct_raw))),
        "integrated_flux_metrics_available": False,
        "integrated_flux_metrics_failure_reason": None,
    }

    try:
        if np.any(~np.isfinite(speed)) or np.any(speed <= 0.0):
            raise ValueError("speed grid must be positive and finite")
        if collision_operator_state is None:
            pitch = np.full_like(speed, float(collision_state.beta_m) * float(collision_state.critical_velocity_m_s) ** 3)
            tau_s = float(collision_state.spitzer_slowing_down_time_s)
        else:
            state_speed = np.asarray(collision_operator_state.speed_grid.centers_m_s, dtype=float)
            pitch = np.asarray(collision_operator_state.ion_pitch_scattering_velocity_cubed_m3_s3, dtype=float)
            tau_s = float(collision_operator_state.spitzer_slowing_down_time_s)
            if state_speed.shape != speed.shape or not np.array_equal(state_speed, speed) or pitch.shape != speed.shape or np.any(~np.isfinite(pitch)) or np.any(pitch < 0.0):
                raise ValueError("collision operator pitch array does not match the speed grid")
        if not np.isfinite(tau_s) or tau_s <= 0.0:
            raise ValueError("slowing time must be positive and finite")
        mirror_ratio = 1.0 / float(basis.lambda_boundary)
        factor = pitch / (tau_s * speed**2) * (1.0 / mirror_ratio) * float(basis.geometry_factor_G) / (float(basis.volume_geometry_integral) * speed)
        if np.any(~np.isfinite(factor)) or np.any(factor < 0.0):
            raise ValueError("Eq 61 transport factor must be finite and nonnegative")
        shell = np.asarray(speed_grid.shell_volumes_m3_s3, dtype=float)
        kinetic = 0.5 * float(fast_ion_mass_kg) * speed**2
        if shell.shape != speed.shape or np.any(~np.isfinite(shell)) or np.any(shell <= 0.0):
            raise ValueError("speed shell measure must be positive and finite")

        def integrated(values: np.ndarray, metrics: dict[str, object]) -> dict[str, object]:
            """Integrate positive and material negative Eq 61 flux contributions"""
            raw_flux = factor * values
            material_negative = values < -float(metrics["roundoff_negative_tolerance"])
            positive_flux = np.where(raw_flux > 0.0, raw_flux, 0.0)
            negative_flux = np.where(material_negative, -raw_flux, 0.0)
            positive_particle = float(np.sum(positive_flux * shell))
            negative_particle = float(np.sum(negative_flux * shell))
            positive_power = float(np.sum(positive_flux * kinetic * shell))
            negative_power = float(np.sum(negative_flux * kinetic * shell))
            return {
                "positive_particle_rate_density_m3_s": positive_particle,
                "material_negative_particle_rate_density_m3_s_magnitude": negative_particle,
                "material_negative_particle_fraction_of_absolute": float(negative_particle / max(positive_particle + negative_particle, np.finfo(float).tiny)),
                "material_negative_to_positive_particle_rate_ratio": float(negative_particle / max(positive_particle, np.finfo(float).tiny)),
                "positive_kinetic_power_density_W_m3": positive_power,
                "material_negative_kinetic_power_density_W_m3_magnitude": negative_power,
                "material_negative_kinetic_power_fraction_of_absolute": float(negative_power / max(positive_power + negative_power, np.finfo(float).tiny)),
                "material_negative_to_positive_kinetic_power_ratio": float(negative_power / max(positive_power, np.finfo(float).tiny)),
            }

        diagnostics["weak_form"] = {**weak_metrics, **integrated(weak_raw, weak_metrics)}
        diagnostics["direct_cell_average"] = {**direct_metrics, **integrated(direct_raw, direct_metrics)}
        diagnostics["integrated_flux_metrics_available"] = True
    except Exception as exc:
        diagnostics["integrated_flux_metrics_failure_reason"] = f"{type(exc).__name__}: {exc}"

    return diagnostics

def _ion_boundary_flux_from_modal(*, speed_grid: SpeedGrid, basis: ModalFBISBasis, modal_f_j_v: np.ndarray, collision_state: FBISCollisionParameterState, fast_ion_mass_kg: float, collision_operator_state: Eq59CollisionOperatorState | None = None, species_id: str | None = None) -> tuple[np.ndarray, np.ndarray, float, float, float]:
    """Reconstruct one end Eq 61 dF/dΛ and its particle and kinetic power loss densities"""
    slopes = _conservative_eq61_mode_slopes_dI_dlambda(basis)
    dfdlambda_raw = slopes @ np.asarray(modal_f_j_v, dtype=float)
    slope_scale = max(float(np.max(np.abs(dfdlambda_raw))) if dfdlambda_raw.size else 0.0, np.finfo(float).tiny)
    negative_tolerance = 128.0 * np.finfo(float).eps * slope_scale
    if np.any(dfdlambda_raw < -negative_tolerance):
        diagnostics = _eq61_boundary_reconstruction_diagnostics(
            speed_grid=speed_grid,
            basis=basis,
            modal_f_j_v=modal_f_j_v,
            collision_state=collision_state,
            fast_ion_mass_kg=fast_ion_mass_kg,
            collision_operator_state=collision_operator_state,
            species_id=species_id,
        )
        raise ModalPhysicalDistributionError(
            "weak form Eq 61 reconstruction produced a negative physical boundary flux; "
            "eq61_diagnostic_json=" + json.dumps(diagnostics, sort_keys=True, separators=(",", ":"), allow_nan=False)
        )
    dfdlambda = np.where(dfdlambda_raw < 0.0, 0.0, dfdlambda_raw)
    v = np.asarray(speed_grid.centers_m_s, dtype=float)
    if collision_operator_state is None:
        pitch = np.full_like(v, float(collision_state.beta_m) * float(collision_state.critical_velocity_m_s) ** 3)
        tau_s = float(collision_state.spitzer_slowing_down_time_s)
    else:
        state_speed = np.asarray(collision_operator_state.speed_grid.centers_m_s, dtype=float)
        pitch = np.asarray(collision_operator_state.ion_pitch_scattering_velocity_cubed_m3_s3, dtype=float)
        tau_s = float(collision_operator_state.spitzer_slowing_down_time_s)
        if state_speed.shape != v.shape or not np.array_equal(state_speed, v) or pitch.shape != v.shape or np.any(~np.isfinite(pitch)) or np.any(pitch < 0.0):
            raise ValueError("Eq 61 collision operator pitch array must match the speed grid")
        if not np.isclose(tau_s, float(collision_state.spitzer_slowing_down_time_s), rtol=4.0e-15, atol=0.0):
            raise ValueError("Eq 61 collision operator slowing time must match the collision state")
    gamma = pitch / (tau_s * v**2)
    mirror_ratio = 1.0 / basis.lambda_boundary
    flux = gamma * dfdlambda * (1.0 / mirror_ratio) * basis.geometry_factor_G / (basis.volume_geometry_integral * v)
    flux = np.maximum(np.where(np.isfinite(flux), flux, 0.0), 0.0)
    shell = np.asarray(speed_grid.shell_volumes_m3_s3, dtype=float)
    kinetic = 0.5 * float(fast_ion_mass_kg) * v**2
    particle_rate_density = float(np.sum(flux * shell))
    kinetic_power_density = float(np.sum(flux * kinetic * shell))
    mean_energy = kinetic_power_density / particle_rate_density if particle_rate_density > 0.0 else 0.0

    return dfdlambda, flux, particle_rate_density, kinetic_power_density, mean_energy

def _modal_eigenvalue_sink_audit(*, speed_grid: SpeedGrid, basis: ModalFBISBasis, modal_f_j_v: np.ndarray, collision_state: FBISCollisionParameterState, fast_ion_mass_kg: float, collision_operator_state: Eq59CollisionOperatorState | None = None, eigenvalues_by_speed: np.ndarray | None = None) -> dict[str, object]:
    """Integrate the modal eigenvalue sink by mode for comparison with the Eq 61 boundary loss"""
    f = np.asarray(modal_f_j_v, dtype=float)
    v = np.asarray(speed_grid.centers_m_s, dtype=float)
    if f.ndim != 2 or f.shape[1] != v.size:
        raise ValueError("modal_f_j_v must have shape (n_modes, n_speed)")
    magnetic_lambdas = np.asarray(basis.physical_basis.eigenvalues, dtype=float)
    if f.shape[0] != magnetic_lambdas.size:
        raise ValueError("modal_f_j_v mode count must match the physical eigenbasis")
    if eigenvalues_by_speed is None:
        lambdas = np.broadcast_to(magnetic_lambdas[:, None], f.shape)
    else:
        lambdas = np.asarray(eigenvalues_by_speed, dtype=float)
        if lambdas.shape != f.shape or np.any(~np.isfinite(lambdas)) or np.any(lambdas <= 0.0):
            raise ValueError("eigenvalues_by_speed must be positive, finite, and match modal_f_j_v")
    mode_eta_integrals = _modal_mode_eta_integrals(basis)
    if collision_operator_state is None:
        pitch = np.full_like(v, float(collision_state.beta_m) * float(collision_state.critical_velocity_m_s) ** 3)
        tau_s = float(collision_state.spitzer_slowing_down_time_s)
        pitch_model = "legacy_scalar_beta_m_vc3_cold_ion_reference"
    else:
        state_speed = np.asarray(collision_operator_state.speed_grid.centers_m_s, dtype=float)
        pitch = np.asarray(collision_operator_state.ion_pitch_scattering_velocity_cubed_m3_s3, dtype=float)
        tau_s = float(collision_operator_state.spitzer_slowing_down_time_s)
        if state_speed.shape != v.shape or not np.array_equal(state_speed, v) or pitch.shape != v.shape or np.any(~np.isfinite(pitch)) or np.any(pitch < 0.0):
            raise ValueError("Eq 59 modal sink pitch array must match the speed grid")
        if not np.isclose(tau_s, float(collision_state.spitzer_slowing_down_time_s), rtol=4.0e-15, atol=0.0):
            raise ValueError("Eq 59 modal sink slowing time must match the collision state")
        pitch_model = collision_operator_state.operator_model
    shell = np.asarray(speed_grid.shell_volumes_m3_s3, dtype=float)
    kinetic = 0.5 * float(fast_ion_mass_kg) * v**2
    loss_frequency = pitch / (tau_s * np.maximum(v, _EPS) ** 3)
    mode_particle_rate_density = np.zeros(f.shape[0], dtype=float)
    mode_midplane_power_density = np.zeros(f.shape[0], dtype=float)
    for j, lam in enumerate(lambdas):
        mode_weight = float(mode_eta_integrals[j])
        rate_integrand = shell * loss_frequency * lam * f[j, :] * mode_weight
        rate_integrand = np.where(np.isfinite(rate_integrand), rate_integrand, 0.0)
        mode_particle_rate_density[j] = float(np.sum(rate_integrand))
        mode_midplane_power_density[j] = float(np.sum(rate_integrand * kinetic))
    total_particle_rate_density = float(np.sum(mode_particle_rate_density))
    total_midplane_power_density = float(np.sum(mode_midplane_power_density))
    first_mode_particle_rate_density = float(mode_particle_rate_density[0]) if mode_particle_rate_density.size else 0.0
    first_mode_midplane_power_density = float(mode_midplane_power_density[0]) if mode_midplane_power_density.size else 0.0
    mean_energy_J = total_midplane_power_density / total_particle_rate_density if total_particle_rate_density > 0.0 else 0.0
    first_mode_mean_energy_J = first_mode_midplane_power_density / first_mode_particle_rate_density if first_mode_particle_rate_density > 0.0 else 0.0

    return {
        "pitch_scattering_model": pitch_model,
        "mode_eta_integrals": mode_eta_integrals,
        "mode_particle_rate_density_m3_s": mode_particle_rate_density,
        "mode_midplane_power_density_W_m3": mode_midplane_power_density,
        "total_particle_rate_density_m3_s": total_particle_rate_density,
        "total_midplane_power_density_W_m3": total_midplane_power_density,
        "mean_energy_J": float(mean_energy_J),
        "first_mode_particle_rate_density_m3_s": first_mode_particle_rate_density,
        "first_mode_midplane_power_density_W_m3": first_mode_midplane_power_density,
        "first_mode_mean_energy_J": float(first_mode_mean_energy_J),
    }

def _select_global_loss_rate(*, model: str, source_particle_rate_s: float, eigenvalue_sink_particle_rate_s: float, eigenvalue_sink_midplane_power_W: float, boundary_single_particle_rate_s: float, boundary_single_midplane_power_W: float, boundary_mean_loss_energy_J: float, volume_m3: float, boundary_flux_density_v: np.ndarray) -> dict[str, object]:
    """Select the configured symmetric two ended Eq 61 device loss representation and scale the boundary spectrum consistently"""
    normalized_model = str(model).strip().lower()
    allowed = {"egedal_eq61_two_end_total_device", EGEDAL_HOT_FIXED_MAGNETIC_BOUNDARY, "egedal_2022_published_electrostatic_modal_sink"}
    if normalized_model not in allowed:
        raise ValueError("the modal backend has one production ion loss convention: egedal_eq61_two_end_total_device")
    published_feedback = normalized_model == "egedal_2022_published_electrostatic_modal_sink"
    fixed_boundary = normalized_model == EGEDAL_HOT_FIXED_MAGNETIC_BOUNDARY
    requested_rate = (float(eigenvalue_sink_particle_rate_s) if published_feedback else 2.0 * float(boundary_single_particle_rate_s))
    active_model = ("egedal_2022_published_electrostatic_modal_sink" if published_feedback else (EGEDAL_HOT_FIXED_MAGNETIC_BOUNDARY if fixed_boundary else "egedal_eq61_two_end_total_device"))
    end_count = 2.0
    end_count_model = ("folded_eq59_sink_counts_both_parallel_branches" if published_feedback else "two_symmetric_loss_cone_surfaces_Egedal_eq61_and_f_L_plus_are_one_end")
    if not np.isfinite(requested_rate) or requested_rate < 0.0:
        raise ValueError("the two-ended Eq. 61 particle-loss rate must be finite and nonnegative")
    if boundary_single_particle_rate_s > 0.0:
        scale_from_single = requested_rate / float(boundary_single_particle_rate_s)
    else:
        scale_from_single = 0.0
    active_flux = np.asarray(boundary_flux_density_v, dtype=float) * scale_from_single
    active_midplane_power_W = (float(eigenvalue_sink_midplane_power_W) if published_feedback else 2.0 * float(boundary_single_midplane_power_W))
    active_rate_density_m3_s = requested_rate / max(float(volume_m3), _EPS)

    return {
        "model": active_model,
        "end_count": end_count,
        "end_count_model": end_count_model,
        "particle_loss_rate_s": float(requested_rate),
        "midplane_power_loss_W": float(active_midplane_power_W),
        "loss_rate_density_m3_s": float(active_rate_density_m3_s),
        "boundary_flux_density_v": active_flux,
        "scale_from_single_boundary_flux": float(scale_from_single),
    }
