"""Project attenuated beam births onto the retained physical modal basis

Physical deposited births are partitioned into prompt magnetic losses and confined births before projection
The signed modal projection is audited without forced rate normalization
"""
from __future__ import annotations
from dataclasses import dataclass
from typing import Any
import numpy as np
from scipy.integrate import trapezoid
from scipy.interpolate import PchipInterpolator
from source_model_revamp.fbis.modal.types import ModalFBISBasis
from source_model_revamp.scattering.physical_eigenbasis import PhysicalEigenbasis

NO_CONFINED_FBIS_SOURCE = "no_confined_fbis_source"
PARTIALLY_CONFINED_BEAM_SOURCE = "partially_confined_beam_source"
CONFINED_BEAM_SOURCE = "confined_beam_source"

def eigenfunction_gram_matrix(physical_basis) -> np.ndarray:
    """Return the retained physical eigenfunction Gram matrix over η"""
    eta = np.asarray(physical_basis.eta_grid, dtype=float)
    eigenfunctions = np.asarray(physical_basis.eigenfunctions, dtype=float)

    return np.asarray([[float(np.trapezoid(left * right, eta)) for right in eigenfunctions] for left in eigenfunctions ],dtype=float,)

@dataclass(frozen=True)
class ModalSourceComponentAudit:
    """Physical birth partition and signed modal projection audit for one beam energy component
    
    Rates use particles s^-1 and modal source arrays follow the retained mode axis
    """
    total_deposited_birth_rate_s: float
    quadrature_total_birth_rate_s: float
    prompt_magnetic_loss_birth_rate_s: float
    confined_physical_birth_rate_s: float
    modal_projected_confined_rate_s: float
    unrepresented_confined_truncation_residual_s: float
    total_partition_residual_s: float
    source_state: str
    pitch_total_birth_rates_s: np.ndarray
    pitch_prompt_birth_rates_s: np.ndarray
    pitch_confined_birth_rates_s: np.ndarray
    modal_rhs: np.ndarray
    modal_source_coefficients: np.ndarray
    retained_mode_counts: np.ndarray
    retained_mode_coefficients: tuple[np.ndarray, ...]
    retained_mode_projected_rates_s: np.ndarray
    retained_mode_relative_discrepancies: np.ndarray
    retained_mode_successive_relative_changes: np.ndarray
    retained_mode_negative_reconstruction_fractions: np.ndarray

@dataclass(frozen=True)
class ModalSourceProjectionAudit:
    """Combined physical birth partition and modal truncation audit across beam components"""
    component_audits: tuple[ModalSourceComponentAudit, ...]
    total_deposited_birth_rate_s: float
    prompt_magnetic_loss_birth_rate_s: float
    confined_physical_birth_rate_s: float
    modal_projected_confined_rate_s: float
    unrepresented_confined_truncation_residual_s: float
    total_partition_residual_s: float
    source_state: str
    retained_mode_counts: np.ndarray
    retained_mode_coefficients: tuple[np.ndarray, ...]
    retained_mode_projected_rates_s: np.ndarray
    retained_mode_relative_discrepancies: np.ndarray
    retained_mode_successive_relative_changes: np.ndarray
    retained_mode_negative_reconstruction_fractions: np.ndarray

def _normalized_pitch_weights(component: Any, sample_count: int) -> np.ndarray:
    """Return nonnegative pitch sample weights normalized to unit sum"""
    weights = np.asarray(getattr(component.spatial_source, "pitch_angle_weights", np.ones(1, dtype=float)), dtype=float)
    if weights.shape != (sample_count,):
        raise ValueError("component pitch weights must match sampled Lambda profiles")
    if np.any(~np.isfinite(weights)) or np.any(weights < 0.0):
        raise ValueError("component pitch weights must be finite and nonnegative")
    weight_sum = float(np.sum(weights))
    if weight_sum <= 0.0:
        raise ValueError("component pitch weights must have positive sum")
    
    return weights / weight_sum

def _source_state(total_rate_s: float, confined_rate_s: float, prompt_rate_s: float) -> str:
    """Classify a deposited source as unconfined, fully confined, or partially confined"""
    scale = max(abs(float(total_rate_s)), 1.0)
    tolerance = 64.0 * np.finfo(float).eps * scale
    if abs(float(confined_rate_s)) <= tolerance:
        return NO_CONFINED_FBIS_SOURCE
    if abs(float(prompt_rate_s)) <= tolerance:
        return CONFINED_BEAM_SOURCE
    
    return PARTIALLY_CONFINED_BEAM_SOURCE

def _component_sample_accounting(*, component: Any, basis: ModalFBISBasis) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Return sampled Λ, mapped η, physical birth rates, confinement mask, and normalized pitch weights
    
    Sample rate arrays have shape (n_pitch, n_z) and are formed from birth rate density in m^-3 s^-1 times physical cell volume
    """
    q = np.asarray(component.spatial_source.axial_birth_rate_density_m3_s, dtype=float)
    volumes = np.asarray(component.attenuation.cell_volumes_m3, dtype=float)
    Lambda_birth = np.asarray(component.spatial_source.Lambda_birth_profile, dtype=float)
    Lambda_samples = np.asarray(getattr(component.spatial_source, "Lambda_birth_sample_profiles", Lambda_birth[None, :]), dtype=float)
    if q.shape != volumes.shape or q.shape != Lambda_birth.shape:
        raise ValueError("component source profile arrays have inconsistent shapes")
    if q.ndim != 1 or np.any(~np.isfinite(q)) or np.any(q < 0.0):
        raise ValueError("component birth rate density must be finite and nonnegative")
    if np.any(~np.isfinite(volumes)) or np.any(volumes <= 0.0):
        raise ValueError("component attenuation volumes must be positive and finite")
    if Lambda_samples.ndim != 2 or Lambda_samples.shape[1] != q.size:
        raise ValueError("component sampled Lambda profiles must have shape (n_pitch, n_z)")
    if np.any(~np.isfinite(Lambda_samples)) or np.any((Lambda_samples < 0.0) | (Lambda_samples > 1.0)):
        raise ValueError("component sampled Lambda profiles must be finite and inside [0, 1]")
    pitch_weights = _normalized_pitch_weights(component, Lambda_samples.shape[0])
    eta_birth_raw = np.asarray(basis.eta_lambda_map.eta_of_lambda(Lambda_samples), dtype=float)
    eta_tolerance = 1.0e-10
    sample_cell_rates_s = pitch_weights[:, None] * q[None, :] * volumes[None, :]
    # Only physical births mapped inside the retained η interval contribute to the confined modal source
    confined_mask = ((sample_cell_rates_s > 0.0) & np.isfinite(eta_birth_raw) & (eta_birth_raw >= -eta_tolerance) & (eta_birth_raw <= basis.eta_boundary + eta_tolerance))
    eta_birth = np.clip(eta_birth_raw, 0.0, basis.eta_boundary)

    return Lambda_samples, eta_birth, sample_cell_rates_s, confined_mask, pitch_weights

def _modal_projected_rate_s(*, coefficients: np.ndarray, physical_basis: PhysicalEigenbasis, volume_m3: float) -> float:
    """Integrate the signed modal source with 4pi * v^2 dv deta"""
    coeff = np.asarray(coefficients, dtype=float)
    if coeff.shape != np.asarray(physical_basis.eigenvalues).shape:
        raise ValueError("modal source coefficients must match the retained physical modes")
    mode_eta_integrals = np.asarray(trapezoid(physical_basis.eigenfunctions, physical_basis.eta_grid, axis=1), dtype=float)

    return float(4.0 * np.pi * float(volume_m3) * np.sum(coeff * mode_eta_integrals))

def _successive_relative_changes(values: np.ndarray) -> np.ndarray:
    """Return signed successive changes scaled by adjacent magnitudes"""
    sequence = np.asarray(values, dtype=float)
    changes = np.zeros_like(sequence)
    if sequence.size > 1:
        scale = np.maximum.reduce((np.abs(sequence[1:]), np.abs(sequence[:-1]), np.full(sequence.size - 1, np.finfo(float).tiny)))
        changes[1:] = (sequence[1:] - sequence[:-1]) / scale

    return changes

def _negative_reconstruction_l1_fraction(*, coefficients: np.ndarray, physical_basis: PhysicalEigenbasis) -> float:
    """Return integral(max(-S, 0)) / integral(abs(S)) on the eta grid"""
    coeff = np.asarray(coefficients, dtype=float)
    eigenfunctions = np.asarray(physical_basis.eigenfunctions, dtype=float)
    if coeff.ndim != 1 or coeff.size > eigenfunctions.shape[0]:
        raise ValueError("coefficient prefix must fit the retained physical basis")
    reconstructed = np.asarray(coeff @ eigenfunctions[: coeff.size], dtype=float)
    absolute_l1 = float(trapezoid(np.abs(reconstructed), physical_basis.eta_grid))
    if absolute_l1 <= np.finfo(float).tiny:
        return 0.0
    negative_l1 = float(trapezoid(np.maximum(-reconstructed, 0.0), physical_basis.eta_grid))

    return float(np.clip(negative_l1 / absolute_l1, 0.0, 1.0))

def _retained_mode_projection_sequence(*, rhs: np.ndarray, physical_basis: PhysicalEigenbasis, volume_m3: float, confined_rate_s: float) -> tuple[np.ndarray, tuple[np.ndarray, ...], np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Resolve every nested retained mode subspace without rescaling the physical source"""
    full_rhs = np.asarray(rhs, dtype=float)
    gram = eigenfunction_gram_matrix(physical_basis)
    mode_integrals = np.asarray(trapezoid(physical_basis.eigenfunctions, physical_basis.eta_grid, axis=1), dtype=float)
    counts = np.arange(1, full_rhs.size + 1, dtype=int)
    rates = np.empty(counts.size, dtype=float)
    coefficient_prefixes: list[np.ndarray] = []
    negative_fractions = np.empty(counts.size, dtype=float)
    for index, count in enumerate(counts):
        # Solve every retained prefix independently because the numerical physical modes are not assumed orthonormal
        coefficients = np.linalg.solve(gram[:count, :count], full_rhs[:count])
        coefficient_prefixes.append(np.asarray(coefficients, dtype=float))
        rates[index] = float(4.0 * np.pi * float(volume_m3) * np.sum(coefficients * mode_integrals[:count]))
        negative_fractions[index] = _negative_reconstruction_l1_fraction(coefficients=coefficients, physical_basis=physical_basis)
    if confined_rate_s > 0.0:
        discrepancies = (rates - float(confined_rate_s)) / float(confined_rate_s)
    else:
        discrepancies = np.zeros_like(rates)
    changes = _successive_relative_changes(rates)

    return (counts, tuple(coefficient_prefixes), rates, discrepancies, changes, negative_fractions)

def _component_modal_rhs_from_spatial_birth(*, component, basis: ModalFBISBasis, volume_m3: float) -> tuple[np.ndarray, int, float, float, float, float]:
    """Project an attenuated spatial birth profile into modal RHS b_j

    The returned right hand side has one value per retained mode and units particles m^-3 s^-1
    """
    physical = basis.physical_basis
    Lambda_samples, eta_birth, spatial_sample_weights, active, _ = (_component_sample_accounting(component=component, basis=basis))
    rhs = np.zeros(physical.eigenvalues.size, dtype=float)
    lambda_mean = float("nan")
    eta_mean = float("nan")
    first_mode_mean = float("nan")
    confined_birth_rate_s = 0.0
    if np.any(active):
        weights = spatial_sample_weights[active]
        weight_sum = float(np.sum(weights))
        confined_birth_rate_s = weight_sum
        if weight_sum > 0.0:
            lambda_mean = float(np.sum(weights * Lambda_samples[active]) / weight_sum)
            eta_mean = float(np.sum(weights * eta_birth[active]) / weight_sum)
        interpolators = [PchipInterpolator(physical.eta_grid, physical.eigenfunctions[j], extrapolate=False) for j in range(physical.eigenvalues.size)]
        for j, interp in enumerate(interpolators):
            I_values = np.asarray(interp(eta_birth[active]), dtype=float)
            rhs[j] = float(np.sum(weights * I_values) / (volume_m3 * 4.0 * np.pi))
            if j == 0 and weight_sum > 0.0:
                first_mode_mean = float(np.sum(weights * I_values) / weight_sum)

    return (rhs, int(np.count_nonzero(np.any(active, axis=0))), lambda_mean, eta_mean, first_mode_mean, confined_birth_rate_s)

def _solve_modal_source_coefficients(rhs: np.ndarray, physical_basis: PhysicalEigenbasis) -> np.ndarray:
    """Solve the physical eigenfunction Gram system for signed modal source coefficients"""
    G = eigenfunction_gram_matrix(physical_basis)
  
    return np.linalg.solve(G, rhs)

def audit_component_modal_source_projection(*, component: Any, basis: ModalFBISBasis, volume_m3: float) -> ModalSourceComponentAudit:
    """Audit one component physical partition and unnormalized signed modal projection"""
    volume = float(volume_m3)
    if not np.isfinite(volume) or volume <= 0.0:
        raise ValueError("volume_m3 must be positive and finite")
    _, _, sample_cell_rates_s, confined_mask, _ = _component_sample_accounting(component=component, basis=basis)
    pitch_total = np.sum(sample_cell_rates_s, axis=1)
    pitch_confined = np.sum(np.where(confined_mask, sample_cell_rates_s, 0.0), axis=1)
    pitch_prompt = pitch_total - pitch_confined
    if np.any(pitch_prompt < -64.0 * np.finfo(float).eps * np.maximum(pitch_total, 1.0)):
        raise ValueError("confined pitch sample birth rate exceeds its deposited birth rate")
    pitch_prompt = np.maximum(pitch_prompt, 0.0)
    quadrature_total = float(np.sum(pitch_total))
    deposited_total = float(component.birth_rate_s)
    if not np.isfinite(deposited_total) or deposited_total < 0.0:
        raise ValueError("component deposited birth rate must be finite and nonnegative")
    quadrature_error = quadrature_total - deposited_total
    quadrature_tolerance = 1.0e-12 * max(abs(deposited_total), abs(quadrature_total), 1.0)
    if abs(quadrature_error) > quadrature_tolerance:
        raise ValueError("pitch quadrature total does not reproduce the deposited component birth rate")
    confined = float(np.sum(pitch_confined))
    prompt = deposited_total - confined
    if prompt < -quadrature_tolerance:
        raise ValueError("confined component birth exceeds the deposited component birth rate")
    if abs(prompt) <= quadrature_tolerance:
        prompt = 0.0
        confined = deposited_total
    elif abs(confined) <= quadrature_tolerance:
        confined = 0.0
        prompt = deposited_total
    rhs, _, _, _, _, rhs_confined = _component_modal_rhs_from_spatial_birth(component=component, basis=basis, volume_m3=volume)
    if not np.isclose(rhs_confined, confined, rtol=1.0e-12, atol=quadrature_tolerance):
        raise RuntimeError("modal RHS and physical source partition use inconsistent confined rates")
    coefficients = _solve_modal_source_coefficients(rhs, basis.physical_basis)
    projected = _modal_projected_rate_s(coefficients=coefficients, physical_basis=basis.physical_basis, volume_m3=volume)
    (mode_counts, mode_coefficients, mode_rates, mode_discrepancies, mode_changes, mode_negative_fractions) = (_retained_mode_projection_sequence(rhs=rhs, physical_basis=basis.physical_basis, volume_m3=volume, confined_rate_s=confined))
    partition_residual = deposited_total - prompt - confined

    return ModalSourceComponentAudit(
        total_deposited_birth_rate_s=deposited_total,
        quadrature_total_birth_rate_s=quadrature_total,
        prompt_magnetic_loss_birth_rate_s=prompt,
        confined_physical_birth_rate_s=confined,
        modal_projected_confined_rate_s=projected,
        unrepresented_confined_truncation_residual_s=confined - projected,
        total_partition_residual_s=partition_residual,
        source_state=_source_state(deposited_total, confined, prompt),
        pitch_total_birth_rates_s=np.asarray(pitch_total, dtype=float),
        pitch_prompt_birth_rates_s=np.asarray(pitch_prompt, dtype=float),
        pitch_confined_birth_rates_s=np.asarray(pitch_confined, dtype=float),
        modal_rhs=np.asarray(rhs, dtype=float),
        modal_source_coefficients=np.asarray(coefficients, dtype=float),
        retained_mode_counts=mode_counts,
        retained_mode_coefficients=mode_coefficients,
        retained_mode_projected_rates_s=mode_rates,
        retained_mode_relative_discrepancies=mode_discrepancies,
        retained_mode_successive_relative_changes=mode_changes,
        retained_mode_negative_reconstruction_fractions=mode_negative_fractions,
    )

def audit_attenuated_modal_source_projection(*, attenuated_source: Any, basis: ModalFBISBasis, volume_m3: float) -> ModalSourceProjectionAudit:
    """Combine component audits without altering any physical source rate or modal coefficient"""
    components = tuple(audit_component_modal_source_projection(component=component, basis=basis, volume_m3=volume_m3) for component in attenuated_source.component_sources)
    total = float(sum(item.total_deposited_birth_rate_s for item in components))
    prompt = float(sum(item.prompt_magnetic_loss_birth_rate_s for item in components))
    confined = float(sum(item.confined_physical_birth_rate_s for item in components))
    projected = float(sum(item.modal_projected_confined_rate_s for item in components))
    if components:
        mode_counts = components[0].retained_mode_counts.copy()
        if any(not np.array_equal(item.retained_mode_counts, mode_counts) for item in components[1:]):
            raise RuntimeError("source components use inconsistent retained-mode sequences")
        mode_coefficients = tuple(np.sum([item.retained_mode_coefficients[index] for item in components], axis=0) for index in range(mode_counts.size))
        mode_rates = np.sum([item.retained_mode_projected_rates_s for item in components], axis=0)
        mode_negative_fractions = np.asarray([_negative_reconstruction_l1_fraction(coefficients=coefficients, physical_basis=basis.physical_basis) for coefficients in mode_coefficients], dtype=float,)
    else:
        mode_counts = np.empty(0, dtype=int)
        mode_coefficients = ()
        mode_rates = np.empty(0, dtype=float)
        mode_negative_fractions = np.empty(0, dtype=float)
    if confined > 0.0:
        mode_discrepancies = (mode_rates - confined) / confined
    else:
        mode_discrepancies = np.zeros_like(mode_rates)
    mode_changes = _successive_relative_changes(mode_rates)

    return ModalSourceProjectionAudit(
        component_audits=components,
        total_deposited_birth_rate_s=total,
        prompt_magnetic_loss_birth_rate_s=prompt,
        confined_physical_birth_rate_s=confined,
        modal_projected_confined_rate_s=projected,
        unrepresented_confined_truncation_residual_s=confined - projected,
        total_partition_residual_s=total - prompt - confined,
        source_state=_source_state(total, confined, prompt),
        retained_mode_counts=mode_counts,
        retained_mode_coefficients=mode_coefficients,
        retained_mode_projected_rates_s=np.asarray(mode_rates, dtype=float),
        retained_mode_relative_discrepancies=np.asarray(mode_discrepancies, dtype=float),
        retained_mode_successive_relative_changes=np.asarray(mode_changes, dtype=float),
        retained_mode_negative_reconstruction_fractions=mode_negative_fractions,
    )

def modal_source_projection_audit_metadata(audit: ModalSourceProjectionAudit) -> dict[str, Any]:
    """Return stable metadata for source partition, folded phase space measure, and retained mode truncation"""
    confined = float(audit.confined_physical_birth_rate_s)
    projected = float(audit.modal_projected_confined_rate_s)
    relative_error = ((projected - confined) / confined if confined > 0.0 else 0.0)
    return {
        "kinetic_source_state": audit.source_state,
        "modal_source_projection_reference_population": "confined_physical_beam_births_only",
        "modal_source_phase_space_measure": "4_pi_v_squared_dv_deta_folded_even_parallel_signs",
        "modal_source_projection_rhs_units": "particles_per_cubic_meter_per_second",
        "modal_source_projection_coefficient_units": "particles_per_cubic_meter_per_second",
        "modal_source_rate_units": "particles_per_second",
        "modal_source_volume_weighting": "axial_birth_rate_density_times_attenuation_cell_volume_divided_by_total_flux_tube_volume",
        "modal_source_component_weighting": "independent_physical_energy_component_birth_rates",
        "modal_source_pitch_quadrature_weighting": "normalized_pitch_sample_weights_shared_with_the_physical_beam_source",
        "modal_source_end_counting": "not_applicable_birth_source_is_not_an_end_loss_rate",
        "modal_source_folded_parallel_sign_count": 2,
        "modal_source_extra_parallel_sign_factor_applied": False,
        "modal_source_projection_forced_normalization_applied": False,
        "modal_total_deposited_birth_rate_s": audit.total_deposited_birth_rate_s,
        "modal_prompt_magnetic_loss_birth_rate_s": audit.prompt_magnetic_loss_birth_rate_s,
        "modal_confined_physical_birth_rate_s": audit.confined_physical_birth_rate_s,
        "modal_projected_confined_source_rate_s": audit.modal_projected_confined_rate_s,
        "modal_unrepresented_confined_truncation_residual_s": audit.unrepresented_confined_truncation_residual_s,
        "modal_source_total_partition_residual_s": audit.total_partition_residual_s,
        "modal_source_projection_relative_error_against_confined_birth": relative_error,
        "modal_component_source_states": [item.source_state for item in audit.component_audits],
        "modal_component_total_deposited_birth_rates_s": [item.total_deposited_birth_rate_s for item in audit.component_audits],
        "modal_component_quadrature_total_birth_rates_s": [item.quadrature_total_birth_rate_s for item in audit.component_audits],
        "modal_component_prompt_magnetic_loss_birth_rates_s": [item.prompt_magnetic_loss_birth_rate_s for item in audit.component_audits],
        "modal_component_confined_physical_birth_rates_s": [item.confined_physical_birth_rate_s for item in audit.component_audits],
        "modal_component_projected_confined_source_rates_s": [item.modal_projected_confined_rate_s for item in audit.component_audits],
        "modal_component_unrepresented_confined_truncation_residuals_s": [item.unrepresented_confined_truncation_residual_s for item in audit.component_audits],
        "modal_component_total_partition_residuals_s": [item.total_partition_residual_s for item in audit.component_audits],
        "modal_component_pitch_total_birth_rates_s": [item.pitch_total_birth_rates_s for item in audit.component_audits],
        "modal_component_pitch_prompt_birth_rates_s": [item.pitch_prompt_birth_rates_s for item in audit.component_audits],
        "modal_component_pitch_confined_birth_rates_s": [item.pitch_confined_birth_rates_s for item in audit.component_audits],
        "modal_component_signed_source_coefficients": [item.modal_source_coefficients for item in audit.component_audits],
        "modal_source_projection_retained_mode_counts": audit.retained_mode_counts,
        "modal_source_projection_retained_mode_signed_coefficients": audit.retained_mode_coefficients,
        "modal_source_projection_retained_mode_component_signed_coefficients": [item.retained_mode_coefficients for item in audit.component_audits],
        "modal_source_projection_retained_mode_rates_s": audit.retained_mode_projected_rates_s,
        "modal_source_projection_retained_mode_component_rates_s": [item.retained_mode_projected_rates_s for item in audit.component_audits],
        "modal_source_projection_retained_mode_relative_discrepancies": audit.retained_mode_relative_discrepancies,
        "modal_source_projection_retained_mode_component_relative_discrepancies": [item.retained_mode_relative_discrepancies for item in audit.component_audits],
        "modal_source_projection_retained_mode_successive_relative_changes": audit.retained_mode_successive_relative_changes,
        "modal_source_projection_retained_mode_component_successive_relative_changes": [item.retained_mode_successive_relative_changes for item in audit.component_audits],
        "modal_source_projection_negative_reconstruction_fraction_definition": "integral_maximum_of_negative_reconstruction_and_zero_deta_divided_by_integral_absolute_reconstruction_deta",
        "modal_source_projection_retained_mode_negative_reconstruction_fractions": audit.retained_mode_negative_reconstruction_fractions,
        "modal_source_projection_retained_mode_component_negative_reconstruction_fractions": [item.retained_mode_negative_reconstruction_fractions for item in audit.component_audits],
        "modal_source_projection_max_pointwise_error_status": "not_applicable_without_a_resolved_analytic_source_profile",
    }
