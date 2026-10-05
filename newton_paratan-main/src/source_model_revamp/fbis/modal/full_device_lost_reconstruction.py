"""
Reconstruct fixed boundary Eq 63 lost ion branches in local velocity space

The reconstruction maps invariant total energy into local kinetic energy after the electrostatic potential is known and integrates the directed loss cone pitch intervals analytically
"""
from __future__ import annotations
from dataclasses import dataclass
import numpy as np
from scipy.special import erfcx
from source_model_revamp.fbis.modal.local.mapping import _cell_quadrature
from source_model_revamp.fbis.velocity_space_grid import PitchGrid, SpeedGrid

def _G_profiles(*, z_edges_m: np.ndarray, B_tilde_centers: np.ndarray, left_throat_z_m: float, right_throat_z_m: float, half_length_m: float, mirror_ratio: float, target_G_at_throat: float) -> tuple[np.ndarray, np.ndarray]:
    """Build directed cumulative Eq 62 geometry profiles across the central mirror cell
    
    The left and right profiles are scaled so their throat value matches the supplied geometry factor
    Cells outside the central throat interval receive zero contribution
    """
    centers = 0.5 * (z_edges_m[:-1] + z_edges_m[1:])
    widths = np.diff(z_edges_m)
    if not np.isfinite(half_length_m) or half_length_m <= 0.0:
        raise ValueError("half_length_m must be positive and finite")
    if not np.isfinite(mirror_ratio) or mirror_ratio <= 1.0:
        raise ValueError("mirror_ratio must be finite and greater than one")
    if not np.isfinite(target_G_at_throat) or target_G_at_throat <= 0.0:
        raise ValueError("target_G_at_throat must be positive and finite")
    integrand = widths / float(half_length_m) / B_tilde_centers * np.sqrt(np.maximum(1.0 - B_tilde_centers / float(mirror_ratio), 0.0))
    confined = (centers >= left_throat_z_m) & (centers <= right_throat_z_m)
    integrand[~confined] = 0.0
    central_integral = float(np.sum(integrand[confined]))
    if not np.isfinite(central_integral) or central_integral <= 0.0:
        raise ValueError("confined field profile gives a nonpositive Eq 63 geometry integral")
    scale = float(target_G_at_throat) / central_integral
    integrand *= scale
    right = np.zeros_like(centers)
    active_right = centers >= left_throat_z_m
    right[active_right] = np.cumsum(integrand[active_right]) - 0.5 * integrand[active_right]
    left = np.zeros_like(centers)
    active_left = centers <= right_throat_z_m
    reverse = integrand[active_left][::-1]
    left[active_left] = (np.cumsum(reverse) - 0.5 * reverse)[::-1]
  
    return np.maximum(left, 0.0), np.maximum(right, 0.0)

def _loss_cone_lower_abs_pitch(total_energy_J: np.ndarray, local_kinetic_J: np.ndarray, B_tilde: float, mirror_ratio: float) -> np.ndarray:
    """Return the local lower absolute pitch boundary for fixed magnetic Eq 63 losses"""
    U = np.asarray(total_energy_J, dtype=float)
    K = np.asarray(local_kinetic_J, dtype=float)
    ratio = np.divide(float(B_tilde) * U, float(mirror_ratio) * K, out=np.zeros_like(U), where=(U > 0.0) & (K > 0.0))
    result = np.sqrt(np.maximum(1.0 - ratio, 0.0))
    result[U <= 0.0] = 1.0
  
    return result

def _exact_eq63_pitch_integral(*, total_energy_J: np.ndarray, local_kinetic_J: np.ndarray, B_tilde: float, mirror_ratio: float, parallel_temperature_J: float, lower_abs_pitch: np.ndarray, upper_abs_pitch: float) -> np.ndarray:
    """Integrate the Eq 63 pitch exponential analytically over one absolute pitch interval
    
    The returned array follows the broadcast total energy shape and remains stable near zero exponential argument
    """
    U = np.asarray(total_energy_J, dtype=float)
    K = np.asarray(local_kinetic_J, dtype=float)
    lo = np.asarray(lower_abs_pitch, dtype=float)
    hi = float(upper_abs_pitch)
    out = np.zeros_like(U)
    active = (U > 0.0) & (K > 0.0) & (lo < hi)

    if not np.any(active):
        return out
    
    T = float(parallel_temperature_J)
    a = float(mirror_ratio) * K[active] / (float(B_tilde) * T)
    lo_a = lo[active]
    width = hi - lo_a
    C = -U[active] / T + a
    value = np.zeros_like(a)
    small_a = a < 1.0e-10
    if np.any(small_a):
        value[small_a] = np.exp(np.clip(C[small_a], -700.0, 0.0)) * width[small_a]
    regular = ~small_a
    if np.any(regular):
        ar = a[regular]
        lor = lo_a[regular]
        xlo = np.sqrt(ar) * lor
        xhi = np.sqrt(ar) * hi
        lower_exponent = C[regular] - ar * lor**2
        exponent_delta = xlo**2 - xhi**2
        raw = erfcx(xlo) - np.exp(np.clip(exponent_delta, -700.0, 0.0)) * erfcx(xhi)
        scale = np.maximum(np.abs(erfcx(xlo)), np.finfo(float).tiny)
        if np.any(raw < -256.0 * np.finfo(float).eps * scale):
            raise RuntimeError("exact Eq 63 pitch integral produced a material negative value")
        raw = np.maximum(raw, 0.0)
        value[regular] = np.exp(np.clip(lower_exponent, -700.0, 0.0)) * np.sqrt(np.pi) / (2.0 * np.sqrt(ar)) * raw
    out[active] = value
  
    return out

def _eq63_H_index_mapping(*, total_energy_J: np.ndarray, particle_mass_kg: float, invariant_speed_faces_m_s: np.ndarray, H_size: int) -> tuple[np.ndarray, np.ndarray]:
    """Map local total energies to piecewise constant invariant speed cells for H(U) sampling"""
    U = np.asarray(total_energy_J, dtype=float)
    positive = U > 0.0
    indices = np.full(U.shape, -1, dtype=np.intp)
    represented = np.zeros(U.shape, dtype=bool)

    if not np.any(positive):
        return indices, represented
    
    speed = np.zeros_like(U)
    speed[positive] = np.sqrt(2.0 * U[positive] / float(particle_mass_kg))
    indices = np.searchsorted(invariant_speed_faces_m_s, speed, side="right") - 1
    indices[speed == invariant_speed_faces_m_s[-1]] = int(H_size) - 1
    represented = positive & (indices >= 0) & (indices < int(H_size))
  
    return indices, represented

def _sample_H_from_index_mapping(*, indices: np.ndarray, represented: np.ndarray, H_U: np.ndarray) -> np.ndarray:
    """Sample H(U) from precomputed invariant speed cell indices"""
    result = np.zeros(indices.shape, dtype=float)
    result[represented] = H_U[indices[represented]]
   
    return result

def _sample_H_on_total_energy(*, total_energy_J: np.ndarray, particle_mass_kg: float, invariant_speed_faces_m_s: np.ndarray, H_U: np.ndarray) -> np.ndarray:
    """Sample piecewise constant H(U) values at local total energies"""
    U = np.asarray(total_energy_J, dtype=float)
    indices, represented = _eq63_H_index_mapping(total_energy_J=U, particle_mass_kg=particle_mass_kg, invariant_speed_faces_m_s=invariant_speed_faces_m_s, H_size=np.asarray(H_U).size)
   
    return _sample_H_from_index_mapping(indices=indices, represented=represented, H_U=np.asarray(H_U, dtype=float))

def _pitch_direction_domains(pitch_faces: np.ndarray) -> tuple[tuple[float, float], tuple[float, float]]:
    """Return the negative and positive pitch domains represented by the pitch grid"""
    faces = np.asarray(pitch_faces, dtype=float)
    lower = float(faces[0])
    upper = float(faces[-1])
    positive = (max(0.0, lower), max(0.0, upper))
    negative = (max(0.0, -upper), max(0.0, -lower))
 
    return negative, positive

@dataclass(frozen=True)
class _FixedBoundaryEq63ChargeDensityState:
    """Potential independent inputs for exact fixed boundary Eq 63 charge density evaluation"""
    invariant_speed_faces_m_s: np.ndarray
    left_H_U: np.ndarray
    right_H_U: np.ndarray
    local_kinetic_J: np.ndarray
    speed_squared_dv_weight: np.ndarray
    negative_pitch_domain: tuple[float, float]
    positive_pitch_domain: tuple[float, float]
    left_parallel_temperature_J: float
    right_parallel_temperature_J: float
    mirror_ratio: float
    particle_mass_kg: float
    charge_number: float

def _prepare_fixed_boundary_eq63_charge_density_exact_pitch(*, invariant_speed_grid: SpeedGrid, local_speed_grid: SpeedGrid, pitch_grid: PitchGrid, left_H_U: np.ndarray, right_H_U: np.ndarray, left_parallel_temperature_J: float, right_parallel_temperature_J: float, mirror_ratio: float, particle_mass_kg: float, charge_number: float = 1.0, speed_quadrature_order: int = 3) -> _FixedBoundaryEq63ChargeDensityState:
    """Prepare local speed quadrature, H(U) mappings, and directed pitch domains for repeated density evaluation"""
    R = float(mirror_ratio)
    mass = float(particle_mass_kg)
    charge = float(charge_number)
    T_left = float(left_parallel_temperature_J)
    T_right = float(right_parallel_temperature_J)
    if not np.isfinite(R) or R <= 1.0:
        raise ValueError("mirror_ratio must be finite and greater than one")
    if not np.isfinite(mass) or mass <= 0.0:
        raise ValueError("particle_mass_kg must be positive and finite")
    if not np.isfinite(charge) or charge < 0.0:
        raise ValueError("charge_number must be finite and nonnegative")
    if not np.isfinite(T_left) or T_left <= 0.0 or not np.isfinite(T_right) or T_right <= 0.0:
        raise ValueError("Eq 63 parallel temperatures must be positive and finite")
    invariant_faces = np.asarray(invariant_speed_grid.faces_m_s, dtype=float)
    H_left = np.asarray(left_H_U, dtype=float)
    H_right = np.asarray(right_H_U, dtype=float)
    if H_left.shape != np.asarray(invariant_speed_grid.centers_m_s).shape or H_right.shape != H_left.shape:
        raise ValueError("Eq 63 H_U arrays must match the invariant speed grid")
    if np.any(~np.isfinite(H_left)) or np.any(H_left < 0.0) or np.any(~np.isfinite(H_right)) or np.any(H_right < 0.0):
        raise ValueError("Eq 63 H_U arrays must be finite and nonnegative")
    speed_nodes, speed_weights = _cell_quadrature(np.asarray(local_speed_grid.faces_m_s, dtype=float), int(speed_quadrature_order))
    speed = np.asarray(speed_nodes, dtype=float).reshape(-1)
    dv_weight = np.asarray(speed_weights, dtype=float).reshape(-1)
    local_kinetic = 0.5 * mass * speed**2
    negative_domain, positive_domain = _pitch_direction_domains(np.asarray(pitch_grid.faces, dtype=float))
  
    return _FixedBoundaryEq63ChargeDensityState(
        invariant_speed_faces_m_s=invariant_faces,
        left_H_U=H_left,
        right_H_U=H_right,
        local_kinetic_J=local_kinetic,
        speed_squared_dv_weight=speed**2 * dv_weight,
        negative_pitch_domain=negative_domain,
        positive_pitch_domain=positive_domain,
        left_parallel_temperature_J=T_left,
        right_parallel_temperature_J=T_right,
        mirror_ratio=R,
        particle_mass_kg=mass,
        charge_number=charge,
    )

def _fixed_boundary_eq63_charge_density_exact_pitch_from_state(*, state: _FixedBoundaryEq63ChargeDensityState, B_tilde: float, potential_drop_magnitude_J: float, left_geometry_G: float, right_geometry_G: float) -> float:
    """Evaluate charge number weighted Eq 63 number density from prepared inputs
    
    The returned value has units m^-3 and is weighted by the supplied ion charge number rather than by Coulombs
    """
    B = float(B_tilde)
    P = float(potential_drop_magnitude_J)
    G_left = float(left_geometry_G)
    G_right = float(right_geometry_G)
    if not np.isfinite(B) or B <= 0.0:
        raise ValueError("B_tilde must be positive and finite")
    if not np.isfinite(P) or P < 0.0:
        raise ValueError("potential_drop_magnitude_J must be finite and nonnegative")
    if not np.isfinite(G_left) or G_left < 0.0 or not np.isfinite(G_right) or G_right < 0.0:
        raise ValueError("Eq 63 geometry factors must be finite and nonnegative")
    left_active = G_left != 0.0
    right_active = G_right != 0.0
    if not left_active and not right_active:
        return 0.0
    U = state.local_kinetic_J - P
    xi_min = _loss_cone_lower_abs_pitch(U, state.local_kinetic_J, B, state.mirror_ratio)
    indices, represented = _eq63_H_index_mapping(total_energy_J=U, particle_mass_kg=state.particle_mass_kg, invariant_speed_faces_m_s=state.invariant_speed_faces_m_s, H_size=state.left_H_U.size)
    left_term = np.zeros_like(U)
    right_term = np.zeros_like(U)
    shared_pitch_integral = bool(
        left_active
        and right_active
        and state.left_parallel_temperature_J == state.right_parallel_temperature_J
        and state.negative_pitch_domain == state.positive_pitch_domain
    )
    if shared_pitch_integral:
        shared_lower = np.maximum(xi_min, float(state.negative_pitch_domain[0]))
        integral = _exact_eq63_pitch_integral(total_energy_J=U, local_kinetic_J=state.local_kinetic_J, B_tilde=B, mirror_ratio=state.mirror_ratio, parallel_temperature_J=state.left_parallel_temperature_J, lower_abs_pitch=shared_lower, upper_abs_pitch=float(state.negative_pitch_domain[1]))
        sampled_left = _sample_H_from_index_mapping(indices=indices, represented=represented, H_U=state.left_H_U)
        sampled_right = _sample_H_from_index_mapping(indices=indices, represented=represented, H_U=state.right_H_U)
        left_term = sampled_left * G_left * integral
        right_term = sampled_right * G_right * integral
    else:
        if left_active:
            left_lower = np.maximum(xi_min, float(state.negative_pitch_domain[0]))
            sampled_left = _sample_H_from_index_mapping(indices=indices, represented=represented, H_U=state.left_H_U)
            left_integral = _exact_eq63_pitch_integral(total_energy_J=U, local_kinetic_J=state.local_kinetic_J, B_tilde=B, mirror_ratio=state.mirror_ratio, parallel_temperature_J=state.left_parallel_temperature_J, lower_abs_pitch=left_lower, upper_abs_pitch=float(state.negative_pitch_domain[1]))
            left_term = sampled_left * G_left * left_integral
        if right_active:
            right_lower = np.maximum(xi_min, float(state.positive_pitch_domain[0]))
            sampled_right = _sample_H_from_index_mapping(indices=indices, represented=represented, H_U=state.right_H_U)
            right_integral = _exact_eq63_pitch_integral(total_energy_J=U, local_kinetic_J=state.local_kinetic_J, B_tilde=B, mirror_ratio=state.mirror_ratio, parallel_temperature_J=state.right_parallel_temperature_J, lower_abs_pitch=right_lower, upper_abs_pitch=float(state.positive_pitch_domain[1]))
            right_term = sampled_right * G_right * right_integral
    number_density = 2.0 * np.pi * float(np.sum(state.speed_squared_dv_weight * (left_term + right_term)))
  
    return state.charge_number * number_density

def fixed_boundary_eq63_charge_density_exact_pitch(*, invariant_speed_grid: SpeedGrid, local_speed_grid: SpeedGrid, pitch_grid: PitchGrid, left_H_U: np.ndarray, right_H_U: np.ndarray, left_parallel_temperature_J: float, right_parallel_temperature_J: float, B_tilde: float, potential_drop_magnitude_J: float, left_geometry_G: float, right_geometry_G: float, mirror_ratio: float, particle_mass_kg: float, charge_number: float = 1.0, speed_quadrature_order: int = 3) -> float:
    """Evaluate local charge number weighted Eq 63 density with analytic pitch integration"""
    B = float(B_tilde)
    P = float(potential_drop_magnitude_J)
    G_left = float(left_geometry_G)
    G_right = float(right_geometry_G)
    if not np.isfinite(B) or B <= 0.0:
        raise ValueError("B_tilde must be positive and finite")
    if not np.isfinite(P) or P < 0.0:
        raise ValueError("potential_drop_magnitude_J must be finite and nonnegative")
    if not np.isfinite(G_left) or G_left < 0.0 or not np.isfinite(G_right) or G_right < 0.0:
        raise ValueError("Eq 63 geometry factors must be finite and nonnegative")
    state = _prepare_fixed_boundary_eq63_charge_density_exact_pitch(
        invariant_speed_grid=invariant_speed_grid,
        local_speed_grid=local_speed_grid,
        pitch_grid=pitch_grid,
        left_H_U=left_H_U,
        right_H_U=right_H_U,
        left_parallel_temperature_J=left_parallel_temperature_J,
        right_parallel_temperature_J=right_parallel_temperature_J,
        mirror_ratio=mirror_ratio,
        particle_mass_kg=particle_mass_kg,
        charge_number=charge_number,
        speed_quadrature_order=speed_quadrature_order,
    )
   
    return _fixed_boundary_eq63_charge_density_exact_pitch_from_state(
        state=state,
        B_tilde=B,
        potential_drop_magnitude_J=P,
        left_geometry_G=G_left,
        right_geometry_G=G_right,
    )

def reconstruct_full_device_fixed_boundary_eq63(*, invariant_speed_grid: SpeedGrid, local_speed_grid: SpeedGrid, pitch_grid: PitchGrid, left_H_U: np.ndarray, right_H_U: np.ndarray, left_parallel_temperature_J: float, right_parallel_temperature_J: float, z_edges_m: np.ndarray, B_tilde_centers: np.ndarray, potential_drop_magnitude_centers_J: np.ndarray, left_throat_z_m: float, right_throat_z_m: float, left_wall_z_m: float, right_wall_z_m: float, mirror_ratio: float, half_length_m: float, geometry_factor_G: float, particle_mass_kg: float, quadrature_order: int = 3) -> tuple[np.ndarray, np.ndarray, np.ndarray, dict[str, object]]:
    """Reconstruct central left and right Eq 63 distributions after applying the axial potential
    
    Returned distribution arrays have shape (n_z, n_local_speed, n_pitch)
    The left branch occupies negative pitch and the right branch occupies positive pitch
    Expander reconstruction outside the central throats is intentionally deferred to the expander model
    """
    edges = np.asarray(z_edges_m, dtype=float)
    centers = 0.5 * (edges[:-1] + edges[1:])
    B = np.asarray(B_tilde_centers, dtype=float)
    P = np.asarray(potential_drop_magnitude_centers_J, dtype=float)
    invariant_faces = np.asarray(invariant_speed_grid.faces_m_s, dtype=float)
    invariant_v = np.asarray(invariant_speed_grid.centers_m_s, dtype=float)
    H_left = np.asarray(left_H_U, dtype=float)
    H_right = np.asarray(right_H_U, dtype=float)
    if edges.ndim != 1 or edges.size < 2 or np.any(~np.isfinite(edges)) or np.any(np.diff(edges) <= 0.0):
        raise ValueError("z_edges_m must be a finite strictly increasing grid")
    if B.shape != centers.shape or P.shape != centers.shape:
        raise ValueError("full device field and potential must match z_edges_m")
    if np.any(~np.isfinite(B)) or np.any(B <= 0.0):
        raise ValueError("B_tilde_centers must be positive and finite")
    if np.any(~np.isfinite(P)) or np.any(P < 0.0):
        raise ValueError("potential_drop_magnitude_centers_J must be finite and nonnegative")
    if H_left.shape != invariant_v.shape or H_right.shape != invariant_v.shape:
        raise ValueError("left_H_U and right_H_U must match the invariant speed grid")
    if np.any(~np.isfinite(H_left)) or np.any(H_left < 0.0) or np.any(~np.isfinite(H_right)) or np.any(H_right < 0.0):
        raise ValueError("Eq 63 H_U arrays must be finite and nonnegative")
    T_left = float(left_parallel_temperature_J)
    T_right = float(right_parallel_temperature_J)
    if not np.isfinite(T_left) or T_left <= 0.0 or not np.isfinite(T_right) or T_right <= 0.0:
        raise ValueError("Eq 63 parallel temperatures must be positive and finite")
    R = float(mirror_ratio)
    mass = float(particle_mass_kg)
    if not np.isfinite(R) or R <= 1.0:
        raise ValueError("mirror_ratio must be finite and greater than one")
    if not np.isfinite(mass) or mass <= 0.0:
        raise ValueError("particle_mass_kg must be positive and finite")
    left_G, right_G = _G_profiles(z_edges_m=edges, B_tilde_centers=B, left_throat_z_m=float(left_throat_z_m), right_throat_z_m=float(right_throat_z_m), half_length_m=float(half_length_m), mirror_ratio=R, target_G_at_throat=float(geometry_factor_G))
    speed_nodes, speed_weights = _cell_quadrature(np.asarray(local_speed_grid.faces_m_s, dtype=float), int(quadrature_order))
    local_speed_faces = np.asarray(local_speed_grid.faces_m_s, dtype=float)
    speed_measure = (local_speed_faces[1:] ** 3 - local_speed_faces[:-1] ** 3) / 3.0
    pitch_faces = np.asarray(pitch_grid.faces, dtype=float)
    pitch_widths = np.asarray(pitch_grid.widths, dtype=float)
    result_shape = (centers.size, local_speed_grid.centers_m_s.size, pitch_grid.centers.size)
    left_result = np.zeros(result_shape, dtype=float)
    right_result = np.zeros(result_shape, dtype=float)
    central_mask = (centers >= float(left_throat_z_m)) & (centers <= float(right_throat_z_m))
    for iz in np.flatnonzero(central_mask):
        v_nodes = np.asarray(speed_nodes, dtype=float)
        weights = np.asarray(speed_weights, dtype=float)
        K = 0.5 * mass * v_nodes**2
        # This sign convention stores the positive magnitude of the ion potential drop from the midplane
        U = K - float(P[iz])
        xi_min = _loss_cone_lower_abs_pitch(U, K, float(B[iz]), R)
        sampled_left_H = _sample_H_on_total_energy(total_energy_J=U, particle_mass_kg=mass, invariant_speed_faces_m_s=invariant_faces, H_U=H_left)
        sampled_right_H = _sample_H_on_total_energy(total_energy_J=U, particle_mass_kg=mass, invariant_speed_faces_m_s=invariant_faces, H_U=H_right)
        for ip in range(pitch_faces.size - 1):
            lo = float(pitch_faces[ip])
            hi = float(pitch_faces[ip + 1])
            if hi <= 0.0:
                abs_lo = max(-hi, 0.0)
                abs_hi = min(-lo, 1.0)
                if abs_hi > abs_lo:
                    lower = np.maximum(xi_min, abs_lo)
                    pitch_integral = _exact_eq63_pitch_integral(total_energy_J=U, local_kinetic_J=K, B_tilde=float(B[iz]), mirror_ratio=R, parallel_temperature_J=T_left, lower_abs_pitch=lower, upper_abs_pitch=abs_hi)
                    numerator = np.sum(weights * v_nodes**2 * sampled_left_H * float(left_G[iz]) * pitch_integral, axis=1)
                    left_result[iz, :, ip] = np.divide(numerator, speed_measure * pitch_widths[ip], out=np.zeros_like(numerator), where=(speed_measure * pitch_widths[ip]) > 0.0)
            elif lo >= 0.0:
                abs_lo = max(lo, 0.0)
                abs_hi = min(hi, 1.0)
                if abs_hi > abs_lo:
                    lower = np.maximum(xi_min, abs_lo)
                    pitch_integral = _exact_eq63_pitch_integral(total_energy_J=U, local_kinetic_J=K, B_tilde=float(B[iz]), mirror_ratio=R, parallel_temperature_J=T_right, lower_abs_pitch=lower, upper_abs_pitch=abs_hi)
                    numerator = np.sum(weights * v_nodes**2 * sampled_right_H * float(right_G[iz]) * pitch_integral, axis=1)
                    right_result[iz, :, ip] = np.divide(numerator, speed_measure * pitch_widths[ip], out=np.zeros_like(numerator), where=(speed_measure * pitch_widths[ip]) > 0.0)
            else:
                # Split a pitch cell that straddles zero so each directed branch keeps its physical sign
                negative_width = -lo
                positive_width = hi
                if negative_width > 0.0:
                    lower = np.maximum(xi_min, 0.0)
                    pitch_integral = _exact_eq63_pitch_integral(total_energy_J=U, local_kinetic_J=K, B_tilde=float(B[iz]), mirror_ratio=R, parallel_temperature_J=T_left, lower_abs_pitch=lower, upper_abs_pitch=min(negative_width, 1.0))
                    numerator = np.sum(weights * v_nodes**2 * sampled_left_H * float(left_G[iz]) * pitch_integral, axis=1)
                    left_result[iz, :, ip] = np.divide(numerator, speed_measure * pitch_widths[ip], out=np.zeros_like(numerator), where=(speed_measure * pitch_widths[ip]) > 0.0)
                if positive_width > 0.0:
                    lower = np.maximum(xi_min, 0.0)
                    pitch_integral = _exact_eq63_pitch_integral(total_energy_J=U, local_kinetic_J=K, B_tilde=float(B[iz]), mirror_ratio=R, parallel_temperature_J=T_right, lower_abs_pitch=lower, upper_abs_pitch=min(positive_width, 1.0))
                    numerator = np.sum(weights * v_nodes**2 * sampled_right_H * float(right_G[iz]) * pitch_integral, axis=1)
                    right_result[iz, :, ip] = np.divide(numerator, speed_measure * pitch_widths[ip], out=np.zeros_like(numerator), where=(speed_measure * pitch_widths[ip]) > 0.0)
    invariant_max = float(invariant_speed_grid.faces_m_s[-1])
    required_local_max = float(np.sqrt(invariant_max**2 + 2.0 * float(np.max(P[central_mask]) if np.any(central_mask) else 0.0) / mass))
    represented_local_max = float(local_speed_grid.faces_m_s[-1])
    local_speed_coverage_passed = bool(represented_local_max >= required_local_max * (1.0 - 128.0 * np.finfo(float).eps))
    diagnostics = {
        "model": "egedal_hot_fixed_magnetic_boundary_eq63_central_exact_pitch",
        "potential_model": "solved_central_Eq70_profile",
        "active_loss_boundary": "Lambda_M_equals_1_over_RM",
        "Eq71_72_loss_boundary_role": "confined_population_local_mapping_only",
        "speed_quadrature_order": int(quadrature_order),
        "pitch_integration_model": "analytic_erfcx_with_exact_magnetic_loss_cone_boundary",
        "expander_representation_role": "deferred_to_source_connected_expander_mapper",
        "expander_cells_zeroed_in_this_reconstruction": bool(np.any(~central_mask)),
        "invariant_speed_max_m_s": invariant_max,
        "required_local_physical_speed_max_m_s": required_local_max,
        "represented_local_physical_speed_max_m_s": represented_local_max,
        "local_speed_coverage_passed": local_speed_coverage_passed,
        "local_speed_tail_clipped": False,
        "forced_posterior_normalization_applied": False,
        "invariant_mapping": "U_and_mu_conserved_pointwise",
        "H_U_speed_sampling_model": "piecewise_constant_on_invariant_speed_cells",
        "left_wrong_sign_population_integral": 0.0,
        "right_wrong_sign_population_integral": 0.0,
    }
    left_result = np.maximum(left_result, 0.0)
    right_result = np.maximum(right_result, 0.0)
  
    return left_result + right_result, left_result, right_result, diagnostics
