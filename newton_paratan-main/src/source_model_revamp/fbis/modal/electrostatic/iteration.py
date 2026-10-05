"""Single species Eq 70 to Eq 72 quasineutral potential iteration"""
from __future__ import annotations
import numpy as np
from scipy.optimize import brentq
from source_model_revamp.constants import ELECTRON_CHARGE_C
from source_model_revamp.fbis.mirror_invariant_grid import LambdaGrid
from source_model_revamp.fbis.velocity_space_grid import PitchGrid, SpeedGrid
from source_model_revamp.fbis.modal.local.mapping import _base_distribution_interpolator
from source_model_revamp.fbis.modal.local.pitch import _density_from_local_pitch, _phi_corrected_local_pitch_distribution
from source_model_revamp.fbis.modal.local.speed import derive_local_physical_speed_grid
from source_model_revamp.fbis.modal.types import ModalElectrostaticProfile
from source_model_revamp.fbis.modal.utils import _EPS
from source_model_revamp.fbis.modal.electrostatic.eq70 import _clip_potential_roundoff_only, _electron_density_fraction_eq70, _energy_roundoff_tolerance, _quasineutrality_residual_metrics
from source_model_revamp.fbis.modal.electrostatic.nodes import _exact_electrostatic_node_profiles, _symmetric_average, _symmetric_input, _midplane_density_value, _enforce_midplane_reference
from source_model_revamp.fbis.modal.electrostatic.accessibility import _closed_eq71_modal_inventory_fraction, _effective_potential_throat_boundary_diagnostics

def _classify_iteration_history(history: list[float], tolerance: float) -> tuple[bool, bool, bool, str]:
    """Classify a residual history as converged, stagnant, oscillatory, or divergent"""
    if not history or not np.all(np.isfinite(history)):
        return False, False, True, "nonfinite_quasineutrality_residual"
    if history[-1] <= tolerance:
        return False, False, False, "converged"
    window_size = min(6, len(history))
    window = np.asarray(history[-window_size:], dtype=float)
    stagnated = False
    oscillatory = False
    diverged = False
    if window.size >= 5:
        spread = float(np.max(window) - np.min(window))
        stagnated = spread <= 0.01 * max(float(np.max(window)), tolerance, _EPS)
        changes = np.diff(window)
        nonzero = changes[np.abs(changes) > 1.0e-14 * max(float(np.max(window)), 1.0)]
        if nonzero.size >= 4:
            oscillatory = bool(np.count_nonzero(nonzero[:-1] * nonzero[1:] < 0.0) >= nonzero.size - 2)
    minimum = max(float(np.min(history)), tolerance, _EPS)
    diverged = bool(history[-1] > 2.0 * minimum and history[-1] > 1.2 * history[0])
    if diverged:
        reason = "quasineutrality_iteration_diverged"
    elif oscillatory:
        reason = "quasineutrality_iteration_oscillatory"
    elif stagnated:
        reason = "quasineutrality_iteration_stagnated"
    else:
        reason = "quasineutrality_tolerance_not_met"

    return stagnated, oscillatory, diverged, reason

def _solve_phi_profile_quasineutrality(*, speed_grid: SpeedGrid, lambda_grid: LambdaGrid, pitch_grid: PitchGrid, base_distribution_v_lambda: np.ndarray, zeta: np.ndarray, B_tilde: np.ndarray, cell_volumes_m3: np.ndarray, mirror_ratio: float, wall_barrier_energy_J: float, electron_temperature_J: float, particle_mass_kg: float, iterations: int, relaxation: float, relative_tolerance: float, electron_midplane_density_m3: float | None = None, electron_collision_density_m3: float | None = None, background_positive_charge_density_m3: np.ndarray | None = None, background_midplane_positive_charge_density_m3: float | None = None, background_left_throat_positive_charge_density_m3: float = 0.0, background_right_throat_positive_charge_density_m3: float = 0.0, target_volume_averaged_ion_density_m3: float | None = None, eta_to_local_phase_space_normalization: float = 1.0, local_velocity_quadrature_order: int = 3, low_energy_weight_fraction_tolerance: float = 1.0e-3, invariant_energy_cell_population_weights: np.ndarray | None = None, initial_potential_energy_J: np.ndarray | None = None) -> tuple[ModalElectrostaticProfile, np.ndarray, np.ndarray, np.ndarray]:
    """Iterate Egedal Eqs 70 through 72 for one fast ion species

    Ion reconstruction preserves total energy U and magnetic moment through the local mapping
    The exact magnetic midplane is the zero potential reference in this single species path
    """
    z = np.asarray(zeta, dtype=float)
    B = np.asarray(B_tilde, dtype=float)
    volumes = np.asarray(cell_volumes_m3, dtype=float)
    n_z = z.size
    if B.shape != (n_z,) or volumes.shape != (n_z,):
        raise ValueError("zeta, B_tilde, and cell_volumes must have matching lengths")
    if n_z == 0 or not np.all(np.isfinite(z)):
        raise ValueError("zeta must contain at least one finite axial cell")
    if np.any(~np.isfinite(B)) or np.any(B <= 0.0):
        raise ValueError("B_tilde must contain positive finite values")
    if np.any(~np.isfinite(volumes)) or np.any(volumes <= 0.0):
        raise ValueError("cell_volumes_m3 must contain positive finite values")
    R = float(mirror_ratio)
    if not np.isfinite(R) or R <= 1.0:
        raise ValueError("mirror_ratio must be finite and greater than one")
    magnetic_tolerance = 128.0 * np.finfo(float).eps * max(R, float(np.max(B)))
    if np.any(B > R + magnetic_tolerance):
        raise ValueError("B_tilde exceeds mirror_ratio inside the fitted throat domain")
    wall = float(wall_barrier_energy_J)
    if not np.isfinite(wall) or wall < 0.0:
        raise ValueError("wall_barrier_energy_J must be finite and nonnegative")
    Te = float(electron_temperature_J)
    if not np.isfinite(Te) or Te <= 0.0:
        raise ValueError("electron_temperature_J must be positive and finite")
    target_volume_average = (None if target_volume_averaged_ion_density_m3 is None else float(target_volume_averaged_ion_density_m3))
    if target_volume_average is not None:
        if not np.isfinite(target_volume_average) or target_volume_average < 0.0:
            raise ValueError("target_volume_averaged_ion_density_m3 must be finite and nonnegative")
    background = (np.zeros(n_z, dtype=float) if background_positive_charge_density_m3 is None else np.asarray(background_positive_charge_density_m3, dtype=float))
    if background.shape != (n_z,) or np.any(~np.isfinite(background)) or np.any(background < 0.0):
        raise ValueError("background_positive_charge_density_m3 must be finite and nonnegative with one value per cell")
    background_left = float(background_left_throat_positive_charge_density_m3)
    background_right = float(background_right_throat_positive_charge_density_m3)
    if not np.isfinite(background_left) or background_left < 0.0:
        raise ValueError("background_left_throat_positive_charge_density_m3 must be finite and nonnegative")
    if not np.isfinite(background_right) or background_right < 0.0:
        raise ValueError("background_right_throat_positive_charge_density_m3 must be finite and nonnegative")
    background_midplane_input = None if background_midplane_positive_charge_density_m3 is None else float(background_midplane_positive_charge_density_m3)
    if background_midplane_input is not None and (not np.isfinite(background_midplane_input) or background_midplane_input < 0.0):
        raise ValueError("background_midplane_positive_charge_density_m3 must be finite and nonnegative")
    electron_midplane = (None if electron_midplane_density_m3 is None else float(electron_midplane_density_m3))
    if electron_midplane is not None and (not np.isfinite(electron_midplane) or electron_midplane < 0.0):
        raise ValueError("electron_midplane_density_m3 must be finite and nonnegative")
    electron_collision = (electron_midplane if electron_collision_density_m3 is None else float(electron_collision_density_m3))
    if electron_collision is not None and (not np.isfinite(electron_collision) or electron_collision < 0.0):
        raise ValueError("electron_collision_density_m3 must be finite and nonnegative")
    normalization = float(eta_to_local_phase_space_normalization)
    if not np.isfinite(normalization) or normalization <= 0.0:
        raise ValueError("eta_to_local_phase_space_normalization must be positive and finite")
    local_speed_grid, _ = derive_local_physical_speed_grid(invariant_speed_grid=speed_grid, maximum_potential_drop_magnitude_J=wall, particle_mass_kg=particle_mass_kg)
    base_interpolator = _base_distribution_interpolator(speed_grid, lambda_grid, base_distribution_v_lambda,)
    symmetric = _symmetric_input(z, B, volumes)
    initial_roundoff_clips = 0
    continuation_accepted = False
    if initial_potential_energy_J is not None:
        initial_phi = np.asarray(initial_potential_energy_J, dtype=float)
        if initial_phi.shape != (n_z,):
            raise ValueError("initial_potential_energy_J must have one value per axial cell")
        if np.any(~np.isfinite(initial_phi)):
            raise ValueError("initial_potential_energy_J must contain only finite values")
        continuation_tolerance = _energy_roundoff_tolerance(wall, float(np.max(np.abs(initial_phi))) if initial_phi.size else 0.0)
        continuation_accepted = bool(np.all(initial_phi >= -continuation_tolerance) and np.all(initial_phi <= wall + continuation_tolerance))
        if continuation_accepted:
            phi, initial_roundoff_clips = _clip_potential_roundoff_only(initial_phi, wall, context="initial potential continuation profile")
    if not continuation_accepted:
        if wall <= 0.0:
            phi = np.zeros(n_z, dtype=float)
        else:
            B_midplane = float(np.min(B))
            denominator = R - B_midplane
            if denominator <= magnetic_tolerance:
                raise ValueError("mirror_ratio must exceed the midplane normalized magnetic field")
            phi_raw = wall * (B - B_midplane) / denominator
            phi, initial_roundoff_clips = _clip_potential_roundoff_only(phi_raw, wall, context="initial potential profile")
    if symmetric:
        phi = _symmetric_average(phi)
    phi = _enforce_midplane_reference(z, phi)
    n_iter = max(int(iterations), 1)
    relax = float(relaxation)
    if not np.isfinite(relax) or relax <= 0.0 or relax > 1.0:
        raise ValueError("relaxation must lie in the interval (0, 1]")
    tolerance = float(relative_tolerance)
    if not np.isfinite(tolerance) or tolerance < 0.0:
        raise ValueError("relative_tolerance must be finite and nonnegative")
    low_energy_tolerance = float(low_energy_weight_fraction_tolerance)
    if not np.isfinite(low_energy_tolerance) or low_energy_tolerance < 0.0:
        raise ValueError("low_energy_weight_fraction_tolerance must be finite and nonnegative")
    residual_history: list[float] = []
    n0_history: list[float] = []
    roundoff_clip_count = int(initial_roundoff_clips)
    throat_extrapolation_clip_count = 0
    low_energy_sample_count = 0
    nonpositive_energy_sample_count = 0
    closed_boundary_sample_count = 0
    throat_drop = wall
    throat_drop_left = wall
    throat_drop_right = wall
    direct_throat_roots_valid = True
    background_order = np.argsort(z)
    background_midplane = (float(np.interp(0.0, z[background_order], background[background_order])) if background_midplane_input is None else background_midplane_input)

    def fast_density_at_local_state(*, field_value: float, potential_drop_J: float, throat_drop_J: float) -> float:
        """Evaluate confined fast ion density at one physical B_tilde and potential state

        The common local speed grid includes the maximum kinetic energy gain allowed by the wall barrier
        """
        _, local_pitch_at_node, _ = _phi_corrected_local_pitch_distribution(
            speed_grid=speed_grid,
            local_speed_grid=local_speed_grid,
            lambda_grid=lambda_grid,
            pitch_grid=pitch_grid,
            base_distribution_v_lambda=base_distribution_v_lambda,
            mirror_ratio=R,
            B_tilde=field_value,
            local_potential_drop_magnitude_J=potential_drop_J,
            throat_potential_drop_magnitude_J=throat_drop_J,
            eta_to_local_phase_space_normalization=normalization,
            quadrature_order=int(local_velocity_quadrature_order),
            base_distribution_interpolator=base_interpolator,
            particle_mass_kg=particle_mass_kg,
        )

        return _density_from_local_pitch(local_speed_grid, pitch_grid, local_pitch_at_node)

    def ion_density_for_phi(phi_profile: np.ndarray, throat_potential_override_J: float | None = None) -> tuple[np.ndarray, np.ndarray, np.ndarray, float, dict[str, float | int]]:
        """Reconstruct local f(v, Λ), f(v, ξ), and fast ion density for one axial potential profile"""
        nonlocal throat_drop, throat_drop_left, throat_drop_right
        extrapolation_clips = 0
        if throat_potential_override_J is not None:
            override = float(throat_potential_override_J)
            if not np.isfinite(override) or override < 0.0 or override > wall:
                raise ValueError("throat_potential_override_J must lie inside the physical wall barrier")
            throat_drop = override
            throat_drop_left = override
            throat_drop_right = override
            extrapolation_clips = 0
        local_lambda = np.zeros((n_z, local_speed_grid.centers_m_s.size, lambda_grid.centers.size), dtype=float)
        local_pitch = np.zeros((n_z, local_speed_grid.centers_m_s.size, pitch_grid.centers.size), dtype=float)
        density = np.zeros(n_z, dtype=float)
        diagnostics: dict[str, float | int] = {"roundoff_clip_count": 0, "throat_extrapolation_clip_count": extrapolation_clips, "low_energy_approximation_sample_count": 0, "nonpositive_total_energy_sample_count": 0, "closed_eq71_boundary_sample_count": 0}
        for index in range(n_z):
            local_lambda[index], local_pitch[index], mapping = _phi_corrected_local_pitch_distribution(
                speed_grid=speed_grid,
                local_speed_grid=local_speed_grid,
                lambda_grid=lambda_grid,
                pitch_grid=pitch_grid,
                base_distribution_v_lambda=base_distribution_v_lambda,
                mirror_ratio=mirror_ratio,
                B_tilde=float(B[index]),
                local_potential_drop_magnitude_J=float(phi_profile[index]),
                throat_potential_drop_magnitude_J=float(throat_drop_left if z[index] < 0.0 else throat_drop_right),
                eta_to_local_phase_space_normalization=normalization,
                quadrature_order=int(local_velocity_quadrature_order),
                base_distribution_interpolator=base_interpolator,
                particle_mass_kg=particle_mass_kg,
            )
            density[index] = _density_from_local_pitch(local_speed_grid, pitch_grid, local_pitch[index])
            diagnostics["roundoff_clip_count"] = int(diagnostics["roundoff_clip_count"]) + int(mapping["roundoff_clip_count"])
            diagnostics["low_energy_approximation_sample_count"] = int(diagnostics["low_energy_approximation_sample_count"]) + int(mapping["low_energy_approximation_sample_count"])
            diagnostics["nonpositive_total_energy_sample_count"] = int(diagnostics["nonpositive_total_energy_sample_count"]) + int(mapping["nonpositive_total_energy_sample_count"])
            diagnostics["closed_eq71_boundary_sample_count"] = int(diagnostics["closed_eq71_boundary_sample_count"]) + int(mapping["closed_eq71_boundary_sample_count"])
        if symmetric:
            density = _symmetric_average(density)
            local_lambda = 0.5 * (local_lambda + local_lambda[::-1])
            local_pitch = 0.5 * (local_pitch + local_pitch[::-1])
        scale = normalization

        return density, local_lambda, local_pitch, float(scale), diagnostics

    electron_fraction_midplane = max(_electron_density_fraction_eq70(0.0, wall, Te), 1.0e-300)

    def electron_density_for_state(phi_profile: np.ndarray, ion_density_profile: np.ndarray) -> tuple[float, np.ndarray]:
        """Set the Eq 70 parent Maxwellian amplitude and evaluate n_e(z) in m⁻³"""
        actual_midplane = (max(fast_density_at_local_state(field_value=1.0, potential_drop_J=0.0, throat_drop_J=throat_drop,) + background_midplane, 0.0,) if electron_midplane is None else electron_midplane)
        n0_local = actual_midplane / electron_fraction_midplane
        electron_density_profile = np.array([n0_local * _electron_density_fraction_eq70(value, wall, Te) for value in phi_profile], dtype=float)
        if symmetric:
            electron_density_profile = _symmetric_average(electron_density_profile)

        return float(n0_local), electron_density_profile

    def midplane_residual_reference_density(electron_density_profile: np.ndarray, ion_density_profile: np.ndarray) -> float:
        """Tie relative residual support to the actual magnetic midplane density"""
        if electron_midplane is not None:
            return float(electron_midplane)
        
        return max(float(_midplane_density_value(z, electron_density_profile)), float(_midplane_density_value(z, ion_density_profile)), np.finfo(float).tiny)

    def evaluate_state(phi_profile: np.ndarray, throat_potential_override_J: float | None = None):
        """Evaluate ion mapping, Eq 70 electrons, and quasineutrality residuals for one profile"""
        nonlocal roundoff_clip_count, throat_extrapolation_clip_count
        nonlocal low_energy_sample_count, nonpositive_energy_sample_count, closed_boundary_sample_count
        fast_ion_density, local_lambda, local_pitch, inventory_scale, mapping = ion_density_for_phi(phi_profile, throat_potential_override_J=throat_potential_override_J)
        ion_density = fast_ion_density + background
        n0_local, electron_density = electron_density_for_state(phi_profile, ion_density)
        residual_metrics = _quasineutrality_residual_metrics(electron_density_m3=electron_density, ion_density_m3=ion_density, zeta=z, cell_volumes_m3=volumes, reference_density_m3=midplane_residual_reference_density(electron_density, ion_density))
        relative = np.asarray(residual_metrics["supported_relative_residual"], dtype=float)
        maximum_supported = float(residual_metrics["maximum_supported_relative_residual"])
        maximum_absolute_normalized = float(np.max(np.abs(np.asarray(residual_metrics["absolute_residual_normalized_to_reference"], dtype=float,))))
        total_volume = float(np.sum(volumes))
        integrated_absolute_fraction = float(residual_metrics["volume_integrated_absolute_particle_mismatch"]) / max(float(residual_metrics["reference_density_m3"]) * total_volume, np.finfo(float).tiny)
        maximum = max(maximum_supported, maximum_absolute_normalized, integrated_absolute_fraction)
        roundoff_clip_count += int(mapping["roundoff_clip_count"])
        throat_extrapolation_clip_count += int(mapping["throat_extrapolation_clip_count"])
        low_energy_sample_count += int(mapping["low_energy_approximation_sample_count"])
        nonpositive_energy_sample_count += int(mapping["nonpositive_total_energy_sample_count"])
        closed_boundary_sample_count += int(mapping["closed_eq71_boundary_sample_count"])

        return fast_ion_density, ion_density, local_lambda, local_pitch, inventory_scale, n0_local, electron_density, relative, maximum

    fast_ion_density, ion_density, local_lambda, local_pitch, inventory_scale, n0, electron_density, relative_residual, max_relative = evaluate_state(phi, throat_potential_override_J=wall)
    residual_history.append(max_relative)
    n0_history.append(n0)
    completed_iterations = 0

    def direct_throat_potential(background_throat_m3: float) -> tuple[float, bool]:
        """Solve one throat Eq 70 density root on the bounded interval from zero to wall"""
        if wall <= 0.0:
            return 0.0, True

        def residual(candidate: float) -> float:
            """Return n_e − n_i at one candidate throat potential energy"""
            fast_density = fast_density_at_local_state(field_value=R, potential_drop_J=candidate, throat_drop_J=candidate)
            if electron_midplane is None:
                candidate_midplane_fast_density = fast_density_at_local_state(field_value=1.0, potential_drop_J=0.0, throat_drop_J=candidate)
                candidate_n0 = (candidate_midplane_fast_density + background_midplane) / electron_fraction_midplane
            else:
                candidate_n0 = float(electron_midplane) / electron_fraction_midplane
            electron = candidate_n0 * _electron_density_fraction_eq70(candidate, wall, Te)

            return electron - fast_density - background_throat_m3

        low = residual(0.0)
        high = residual(wall)
        residual_scale = max(abs(low), abs(high), electron_midplane or 0.0, 1.0)
        residual_tolerance = tolerance * residual_scale
        if abs(low) <= residual_tolerance:
            return 0.0, True
        if abs(high) <= residual_tolerance:
            return wall, True
        if low * high < 0.0:
            return (float(brentq(residual, 0.0, wall, xtol=_energy_roundoff_tolerance(wall, Te), rtol=8.0 * np.finfo(float).eps)), True,)
        return (0.0 if abs(low) <= abs(high) else wall), False

    # Alternate direct throat roots with relaxed axial Eq 70 updates
    for iteration in range(n_iter):
        if max_relative <= tolerance:
            break
        completed_iterations = iteration + 1
        throat_drop_left, left_throat_valid = direct_throat_potential(background_left)
        if symmetric and np.isclose(background_left, background_right):
            throat_drop_right = throat_drop_left
            right_throat_valid = left_throat_valid
        else:
            throat_drop_right, right_throat_valid = direct_throat_potential(background_right)
        throat_drop = 0.5 * (throat_drop_left + throat_drop_right)
        direct_throat_roots_valid = bool(left_throat_valid and right_throat_valid)
        fast_ion_density, ion_density, local_lambda, local_pitch, inventory_scale, n0, electron_density, relative_residual, max_relative = evaluate_state(phi)
        new_phi = np.zeros_like(phi)
        for index in range(n_z):
            if wall <= 0.0 or np.isclose(z[index], 0.0, rtol=0.0, atol=1.0e-14):
                continue
            target = float(max(ion_density[index], 0.0))

            def residual(candidate: float) -> float:
                """Return the cell Eq 70 electron density minus the current ion target"""
                return n0 * _electron_density_fraction_eq70(candidate, wall, Te) - target

            residual_low = residual(0.0)
            residual_high = residual(wall)
            if residual_low == 0.0:
                root = 0.0
            elif residual_high == 0.0:
                root = wall
            elif residual_low * residual_high < 0.0:
                root = float(brentq(residual, 0.0, wall, xtol=_energy_roundoff_tolerance(wall, Te), rtol=8.0 * np.finfo(float).eps))
            else:
                root = 0.0 if abs(residual_low) <= abs(residual_high) else wall
            new_phi[index] = root

        if symmetric:
            new_phi = _symmetric_average(new_phi)
        phi = (1.0 - relax) * phi + relax * new_phi
        phi, update_roundoff_clips = _clip_potential_roundoff_only(phi, wall, context="updated potential profile")
        roundoff_clip_count += int(update_roundoff_clips)
        phi = _enforce_midplane_reference(z, phi)
        fast_ion_density, ion_density, local_lambda, local_pitch, inventory_scale, n0, electron_density, relative_residual, max_relative = evaluate_state(phi)
        residual_history.append(max_relative)
        n0_history.append(n0)
        stagnated, oscillatory, diverged, _ = _classify_iteration_history(residual_history, tolerance)
        if diverged or (stagnated and len(residual_history) >= 8):
            break
        if oscillatory and len(residual_history) >= 12 and max_relative > tolerance:
            break

    throat_drop_left, left_throat_valid = direct_throat_potential(background_left)
    if symmetric and np.isclose(background_left, background_right):
        throat_drop_right = throat_drop_left
        right_throat_valid = left_throat_valid
    else:
        throat_drop_right, right_throat_valid = direct_throat_potential(background_right)
    throat_drop = 0.5 * (throat_drop_left + throat_drop_right)
    direct_throat_roots_valid = bool(left_throat_valid and right_throat_valid)
    fast_ion_density, ion_density, local_lambda, local_pitch, inventory_scale, n0, electron_density, relative_residual, max_relative = evaluate_state(phi)
    residual_metrics = _quasineutrality_residual_metrics(electron_density_m3=electron_density, ion_density_m3=ion_density, zeta=z, cell_volumes_m3=volumes, reference_density_m3=midplane_residual_reference_density(electron_density, ion_density),)
    absolute_residual = np.asarray(residual_metrics["signed_absolute_density_residual_m3"], dtype=float)
    relative_residual = np.asarray(residual_metrics["supported_relative_residual"], dtype=float)
    max_relative = float(residual_metrics["maximum_supported_relative_residual"])
    maximum_absolute_normalized = float(residual_metrics["maximum_absolute_residual_normalized_to_reference"])
    volume_integrated_signed_fraction = float(residual_metrics["volume_integrated_signed_particle_mismatch_fraction"])
    volume_integrated_absolute_fraction = float(residual_metrics["volume_integrated_absolute_particle_mismatch_fraction"])
    final_iteration_residual = max(max_relative, maximum_absolute_normalized, volume_integrated_absolute_fraction)
    if not residual_history or not np.isclose(final_iteration_residual, residual_history[-1], rtol=0.0, atol=0.0):
        residual_history.append(final_iteration_residual)
        n0_history.append(n0)
    total_volume_m3 = float(np.sum(volumes))
    stagnated, oscillatory, diverged, iteration_failure_reason = (_classify_iteration_history(residual_history, tolerance))
    effective_potential = _effective_potential_throat_boundary_diagnostics(speed_grid=speed_grid, lambda_grid=lambda_grid, base_distribution_v_lambda=base_distribution_v_lambda, zeta=z, B_tilde=B, potential_drop_magnitude_J=phi, mirror_ratio=R, throat_potential_left_energy_J=throat_drop_left, throat_potential_right_energy_J=throat_drop_right, particle_mass_kg=particle_mass_kg)
    closed_inventory = _closed_eq71_modal_inventory_fraction(speed_grid=speed_grid, mirror_ratio=R, throat_potential_drop_magnitude_J=throat_drop, invariant_energy_cell_population_weights=invariant_energy_cell_population_weights, particle_mass_kg=particle_mass_kg)
    low_energy_valid = bool(closed_inventory["assessed"] and float(closed_inventory["fraction"]) <= low_energy_tolerance)
    # Reinsert exact physical anchors so cell centered interpolation cannot define the qualification nodes
    node_profiles = _exact_electrostatic_node_profiles(
        cell_zeta=z,
        cell_B_tilde=B,
        cell_potential_energy_J=phi,
        cell_electron_density_m3=electron_density,
        cell_ion_density_m3=ion_density,
        cell_fast_ion_density_m3=fast_ion_density,
        cell_background_density_m3=background,
        mirror_ratio=R,
        left_throat_potential_energy_J=throat_drop_left,
        right_throat_potential_energy_J=throat_drop_right,
        electron_parent_maxwellian_n0_m3=n0,
        wall_barrier_energy_J=wall,
        electron_temperature_J=Te,
        exact_midplane_fast_density_m3=fast_density_at_local_state(field_value=1.0, potential_drop_J=0.0, throat_drop_J=throat_drop),
        exact_midplane_potential_energy_J=0.0,
        exact_left_throat_fast_density_m3=fast_density_at_local_state(field_value=R, potential_drop_J=throat_drop_left, throat_drop_J=throat_drop_left),
        exact_right_throat_fast_density_m3=fast_density_at_local_state(field_value=R, potential_drop_J=throat_drop_right, throat_drop_J=throat_drop_right),
        exact_midplane_background_density_m3=background_midplane,
        exact_left_throat_background_density_m3=background_left,
        exact_right_throat_background_density_m3=background_right,
    )
    node_relative_residual = np.asarray(node_profiles["relative_residual_profile"], dtype=float)
    node_absolute_residual = np.asarray(node_profiles["absolute_residual_profile_m3"], dtype=float)
    node_electron_density = np.asarray(node_profiles["electron_density_m3"], dtype=float)
    node_ion_density = np.asarray(node_profiles["ion_density_m3"], dtype=float)
    anchor_indices = np.array((int(node_profiles["left_throat_index"]), int(node_profiles["midplane_index"]), int(node_profiles["right_throat_index"])), dtype=int)
    anchor_density_scale = np.maximum(node_electron_density[anchor_indices], node_ion_density[anchor_indices])
    anchor_supported = anchor_density_scale >= float(residual_metrics["density_support_floor_m3"])
    anchor_normalized_absolute = np.abs(node_absolute_residual[anchor_indices]) / float(residual_metrics["reference_density_m3"])
    anchor_convergence_residual = np.where(anchor_supported, np.abs(node_relative_residual[anchor_indices]), anchor_normalized_absolute)
    exact_left_throat_residual = float(anchor_convergence_residual[0])
    exact_midplane_residual = float(anchor_convergence_residual[1])
    exact_right_throat_residual = float(anchor_convergence_residual[2])
    exact_throat_maximum_residual = max(exact_left_throat_residual, exact_right_throat_residual)
    exact_midplane_valid = bool(exact_midplane_residual <= tolerance)
    exact_throat_valid = bool(exact_throat_maximum_residual <= tolerance)
    exact_anchor_maximum_residual = max(exact_midplane_residual, exact_throat_maximum_residual)
    exact_anchor_residual_valid = bool(exact_midplane_valid and exact_throat_valid)
    exact_midplane_electron_density = float(node_electron_density[anchor_indices[1]])
    convergence_failures: list[str] = []
    if max_relative > tolerance:
        convergence_failures.append("maximum_supported_relative_residual_above_tolerance")
    if maximum_absolute_normalized > tolerance:
        convergence_failures.append("maximum_absolute_reference_normalized_residual_above_tolerance")
    if volume_integrated_absolute_fraction > tolerance:
        convergence_failures.append("volume_integrated_absolute_residual_above_tolerance")
    if not exact_midplane_valid:
        convergence_failures.append("exact_midplane_quasineutrality_residual_above_tolerance")
    if diverged:
        convergence_failures.append("quasineutrality_iteration_diverged")
    elif oscillatory and max_relative > tolerance:
        convergence_failures.append("quasineutrality_iteration_oscillatory")
    elif stagnated and max_relative > tolerance:
        convergence_failures.append("quasineutrality_iteration_stagnated")
    elif max_relative > tolerance and iteration_failure_reason not in {"converged", "quasineutrality_tolerance_not_met"}:
        convergence_failures.append(iteration_failure_reason)
    converged = not convergence_failures
    failure_reason = None if converged else ";".join(convergence_failures)
    central_supported_residual = max(float(residual_metrics["central_region_maximum_supported_relative_residual"]), abs(float(node_relative_residual[anchor_indices[1]])) if anchor_supported[1] else 0.0)
    throat_supported_residual = max(float(residual_metrics["throat_region_maximum_supported_relative_residual"]), abs(float(node_relative_residual[anchor_indices[0]])) if anchor_supported[0] else 0.0, abs(float(node_relative_residual[anchor_indices[2]])) if anchor_supported[2] else 0.0)
    fast_ion_inventory_particles = float(np.sum(fast_ion_density * volumes))
    fast_ion_volume_average_density_m3 = fast_ion_inventory_particles / total_volume_m3
    target_absolute_error_m3 = (None if target_volume_average is None else fast_ion_volume_average_density_m3 - target_volume_average)
    target_relative_error = (None if target_volume_average is None or target_volume_average == 0.0 else (fast_ion_volume_average_density_m3 - target_volume_average) / target_volume_average)
    profile = ModalElectrostaticProfile(
        zeta=z,
        B_tilde=B,
        potential_energy_J=phi,
        potential_relative_to_midplane_V=-phi / ELECTRON_CHARGE_C,
        electron_density_m3=electron_density,
        ion_density_m3=ion_density,
        iterations=completed_iterations,
        max_relative_quasineutrality_error=max_relative,
        converged=converged,
        electron_density_scale_n0_m3=float(n0),
        electron_midplane_density_m3=float(electron_midplane if electron_midplane is not None else exact_midplane_electron_density),
        electron_parent_maxwellian_n0_m3=float(n0),
        electron_collision_density_m3=float(electron_collision if electron_collision is not None else _midplane_density_value(z, electron_density)),
        electron_volume_average_density_m3=float(np.sum(electron_density * volumes) / np.sum(volumes)),
        fast_ion_density_m3=fast_ion_density,
        background_positive_charge_density_m3=background,
        electron_density_input_semantics=("derived_from_total_positive_charge_at_exact_magnetic_midplane" if electron_midplane is None else "benchmark_prescribed_actual_total_electron_density_at_magnetic_midplane"),
        electron_density_closure_model=("eq70_parent_amplitude_from_derived_exact_midplane_total_positive_charge" if electron_midplane is None else "eq70_parent_amplitude_from_benchmark_prescribed_midplane_density"),
        relative_residual_profile=relative_residual,
        absolute_residual_profile_m3=absolute_residual,
        residual_history=tuple(float(value) for value in residual_history),
        electron_density_scale_history_m3=tuple(float(value) for value in n0_history),
        stagnated=stagnated,
        oscillatory=oscillatory,
        diverged=diverged,
        failure_reason=failure_reason,
        energy_mapping_model="conserved_total_ion_energy_U_equals_K_plus_ePhi",
        invariant_mapping_model="egedal_eq71_global_throat_boundary_and_eq72_compression",
        density_measure="local_2pi_v_squared_dv_dxi_with_global_bounce_weighted_eta_inventory",
        roundoff_clip_count=int(roundoff_clip_count),
        throat_extrapolation_clip_count=int(throat_extrapolation_clip_count),
        throat_potential_energy_J=float(throat_drop),
        throat_potential_left_energy_J=float(throat_drop_left),
        throat_potential_right_energy_J=float(throat_drop_right),
        low_energy_approximation_sample_count=int(low_energy_sample_count),
        nonpositive_total_energy_sample_count=int(nonpositive_energy_sample_count),
        closed_eq71_boundary_sample_count=int(closed_boundary_sample_count),
        low_energy_approximation_max_weight_fraction=float(closed_inventory["fraction"]),
        low_energy_approximation_valid=low_energy_valid,
        low_energy_approximation_weight_fraction_assessed=bool(closed_inventory["assessed"]),
        low_energy_approximation_weight_model=str(closed_inventory["model"]),
        eq71_closed_interval_intersected_population_cell_count=int(closed_inventory["intersected_population_cell_count"]),
        eq71_closed_interval_intersects_distribution_support=bool(closed_inventory["intersects_distribution_support"]),
        eta_to_local_phase_space_normalization=normalization,
        ion_inventory_normalization_history=(),
        phase_space_measure_conversion_factor=normalization,
        phase_space_measure_conversion_history=(normalization,),
        phase_space_measure_conversion_model=("eq45_to_eq48_eta_measure_to_local_2pi_v_squared_dv_dxi_no_posterior_normalization"),
        fast_ion_inventory_particles=fast_ion_inventory_particles,
        fast_ion_volume_average_density_m3=fast_ion_volume_average_density_m3,
        target_volume_averaged_fast_ion_density_m3=target_volume_average,
        fast_ion_volume_average_density_absolute_error_m3=target_absolute_error_m3,
        fast_ion_volume_average_density_relative_error=target_relative_error,
        symmetric_input=symmetric,
        midplane_reference_V=0.0,
        throat_extrapolation_valid=direct_throat_roots_valid,
        throat_potential_definition="direct_quasineutral_throat_node_solve_no_extrapolation",
        potential_clipping_policy="roundoff_only_otherwise_raise",
        eq71_closed_interval_threshold_J=(float(throat_drop / (R - 1.0)) if throat_drop > 0.0 else 0.0),
        effective_potential_check_assessed=bool(effective_potential["assessed"]),
        effective_potential_throat_boundary_valid=bool(effective_potential["valid"]),
        effective_potential_failure_reason=(None if effective_potential["failure_reason"] is None else str(effective_potential["failure_reason"])),
        effective_potential_active_energy_sample_count=int(effective_potential["active_energy_sample_count"]),
        effective_potential_open_energy_sample_count=int(effective_potential["open_confined_energy_sample_count"]),
        effective_potential_interior_minimum_energy_sample_count=int(effective_potential["interior_minimum_energy_sample_count"]),
        effective_potential_interior_minimum_weight_fraction=float(effective_potential["interior_minimum_weight_fraction"]),
        effective_potential_max_relative_boundary_shortfall=float(effective_potential["max_relative_boundary_shortfall"]),
        effective_potential_max_relative_throat_asymmetry=float(effective_potential["max_relative_throat_asymmetry"]),
        effective_potential_worst_total_energy_J=float(effective_potential["worst_total_energy_J"]),
        effective_potential_worst_interior_location_zeta=(None if effective_potential["worst_interior_location_zeta"] is None else float(effective_potential["worst_interior_location_zeta"])),
        effective_potential_relative_tolerance=float(effective_potential["relative_tolerance"]),
        absolute_residual_normalized_to_reference=np.asarray(residual_metrics["absolute_residual_normalized_to_reference"], dtype=float),
        density_support_mask=np.asarray(residual_metrics["density_support_mask"], dtype=bool),
        density_support_reference_m3=float(residual_metrics["reference_density_m3"]),
        density_support_relative_floor=float(residual_metrics["density_support_relative_floor"]),
        density_support_floor_m3=float(residual_metrics["density_support_floor_m3"]),
        density_support_floor_derivation=str(residual_metrics["density_support_floor_derivation"]),
        density_support_masked_node_count=int(residual_metrics["masked_node_count"]),
        density_support_masked_volume_m3=float(residual_metrics["masked_volume_m3"]),
        density_support_masked_volume_fraction=float(residual_metrics["masked_volume_fraction"]),
        volume_integrated_electron_minus_ion_particles=float(residual_metrics["volume_integrated_electron_minus_ion_particles"]),
        volume_integrated_absolute_particle_mismatch=float(residual_metrics["volume_integrated_absolute_particle_mismatch"]),
        volume_integrated_signed_charge_mismatch_C=float(residual_metrics["volume_integrated_signed_charge_mismatch_C"]),
        volume_integrated_absolute_charge_mismatch_C=float(residual_metrics["volume_integrated_absolute_charge_mismatch_C"]),
        central_region_maximum_supported_relative_residual=central_supported_residual,
        throat_region_maximum_supported_relative_residual=throat_supported_residual,
        maximum_absolute_density_residual_m3=float(residual_metrics["maximum_absolute_density_residual_m3"]),
        maximum_absolute_density_residual_normalized_to_reference=maximum_absolute_normalized,
        volume_integrated_signed_particle_mismatch_fraction=(volume_integrated_signed_fraction),
        volume_integrated_absolute_particle_mismatch_fraction=volume_integrated_absolute_fraction,
        electrostatic_node_zeta=np.asarray(node_profiles["zeta"], dtype=float),
        electrostatic_node_B_tilde=np.asarray(node_profiles["B_tilde"], dtype=float),
        electrostatic_node_potential_energy_J=np.asarray(node_profiles["potential_energy_J"], dtype=float),
        electrostatic_node_potential_relative_to_midplane_V=np.asarray(node_profiles["potential_relative_to_midplane_V"], dtype=float),
        electrostatic_node_electron_density_m3=np.asarray(node_profiles["electron_density_m3"], dtype=float),
        electrostatic_node_ion_density_m3=np.asarray(node_profiles["ion_density_m3"], dtype=float),
        electrostatic_node_fast_ion_density_m3=np.asarray(node_profiles["fast_ion_density_m3"], dtype=float),
        electrostatic_node_background_positive_charge_density_m3=np.asarray(node_profiles["background_positive_charge_density_m3"], dtype=float),
        electrostatic_node_absolute_residual_profile_m3=np.asarray(node_profiles["absolute_residual_profile_m3"], dtype=float),
        electrostatic_node_relative_residual_profile=node_relative_residual,
        electrostatic_left_throat_index=int(node_profiles["left_throat_index"]),
        electrostatic_midplane_index=int(node_profiles["midplane_index"]),
        electrostatic_right_throat_index=int(node_profiles["right_throat_index"]),
        electrostatic_left_throat_node_exact=bool(node_profiles["left_throat_node_exact"]),
        electrostatic_midplane_node_exact=bool(node_profiles["midplane_node_exact"]),
        electrostatic_right_throat_node_exact=bool(node_profiles["right_throat_node_exact"]),
        electrostatic_midplane_gauge_exact=bool(node_profiles["midplane_gauge_exact"]),
        electrostatic_node_model=str(node_profiles["node_model"]),
        exact_anchor_maximum_quasineutrality_residual=exact_anchor_maximum_residual,
        exact_anchor_quasineutrality_valid=exact_anchor_residual_valid,
        exact_anchor_residual_support_model=("supported_relative_else_absolute_normalized_to_actual_midplane_density"),
        local_speed_grid=local_speed_grid,
        exact_midplane_quasineutrality_residual=exact_midplane_residual,
        exact_left_throat_quasineutrality_residual=exact_left_throat_residual,
        exact_right_throat_quasineutrality_residual=exact_right_throat_residual,
        exact_throat_maximum_quasineutrality_residual=exact_throat_maximum_residual,
        exact_midplane_quasineutrality_valid=exact_midplane_valid,
        exact_throat_quasineutrality_valid=exact_throat_valid,
        direct_left_throat_root_valid=bool(left_throat_valid),
        direct_right_throat_root_valid=bool(right_throat_valid),
        direct_throat_roots_valid=direct_throat_roots_valid,
    )

    return profile, local_lambda, local_pitch, fast_ion_density
