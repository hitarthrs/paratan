"""
Cell integrated beam source placement on speed Lambda grids
"""
from __future__ import annotations
from collections.abc import Iterable
import numpy as np
from numpy.typing import ArrayLike
from source_model_revamp.fbis.beam_source_definition import BeamSourceDefinition, beam_speed_from_energy_m_s, particle_birth_rate_from_current_s, particle_birth_rate_from_power_s, source_rate_density_from_birth_rate_m3_s
from source_model_revamp.fbis.mirror_invariant_grid import LambdaGrid
from source_model_revamp.fbis.velocity_space_grid import SpeedGrid

def _require_1d_finite(name: str, values: ArrayLike) -> np.ndarray:
    arr = np.asarray(values, dtype=float)
    if arr.ndim != 1:
        raise ValueError(f"{name} must be 1D")
    if arr.size == 0:
        raise ValueError(f"{name} must be non empty")
    if np.any(~np.isfinite(arr)):
        raise ValueError(f"{name} must be finite")
 
    return arr

def _require_positive_scalar(name: str, value: float) -> float:
    val = float(value)
    if not np.isfinite(val) or val <= 0.0:
        raise ValueError(f"{name} must be positive and finite")
   
    return val

def _require_increasing_faces(name: str, faces: ArrayLike) -> np.ndarray:
    face_values = _require_1d_finite(name, faces)
    if face_values.size < 2:
        raise ValueError(f"{name} must contain at least two faces")
    if np.any(np.diff(face_values) <= 0.0):
        raise ValueError(f"{name} must be strictly increasing")
   
    return face_values

def cell_index_from_faces(faces: ArrayLike, value: float, *, allow_clamp: bool = False):
    """Return the index of the cell containing a value"""
    face_values = _require_increasing_faces("faces", faces)
    val = float(value)
    if not np.isfinite(val):
        raise ValueError("value must be finite")
    lo = float(face_values[0])
    hi = float(face_values[-1])
    if val < lo or val > hi:
        if not allow_clamp:
            raise ValueError(f"value {val} is outside cell-centered grid range [{lo}, {hi}]")
        val = min(max(val, lo), hi)
    if np.isclose(val, hi, rtol=0.0, atol=1.0e-14 * max(1.0, abs(hi))):
        return face_values.size - 2
    
    return int(np.searchsorted(face_values, val, side="right") - 1)

def speed_cell_index_for_speed(speed_grid: SpeedGrid, speed_m_s: float, *, allow_clamp: bool = False):
    """Speed cell index for a monoenergetic source"""
    return cell_index_from_faces(speed_grid.faces_m_s, speed_m_s, allow_clamp=allow_clamp)

def lambda_cell_index_for_lambda(lambda_grid: LambdaGrid, lambda_value: float, *, allow_clamp: bool = False):
    """Λ cell index for a mono pitch source"""
    return cell_index_from_faces(lambda_grid.faces, lambda_value, allow_clamp=allow_clamp)

def lambda_pitch_cell_measures(lambda_grid: LambdaGrid, B_tilde_birth: float):
    """
    Return two branch pitch widths corresponding to Λ cells at B_tilde_birth
    From ξ^2 = 1 - Λ * B
        Δ_ξ = integral B/sqrt(1 - Λ * B) dΛ = 2 * [sqrt(1 - B * Λ_lo) - sqrt(1 - B * Λ_hi)]
    """
    B = _require_positive_scalar("B_tilde_birth", B_tilde_birth)
    lo = np.asarray(lambda_grid.faces[:-1], dtype=float)
    hi = np.asarray(lambda_grid.faces[1:], dtype=float)
    accessible_hi = np.minimum(hi, 1.0 / B)
    accessible_lo = np.minimum(lo, 1.0 / B)
    width = 2.0 * (np.sqrt(np.maximum(1.0 - B * accessible_lo, 0.0)) - np.sqrt(np.maximum(1.0 - B * accessible_hi, 0.0)))
    width[hi <= 0.0] = 0.0
    width[lo >= 1.0 / B] = 0.0
    
    return width

def monoenergetic_lambda_source_density(speed_grid: SpeedGrid, lambda_grid: LambdaGrid, birth_rate_density_m3_s: float, birth_speed_m_s: float, lambda_birth: float, B_tilde_birth: float):
    """Cell integrated monoenergetic source on F(v, Lambda)"""
    source = np.zeros((speed_grid.centers_m_s.size, lambda_grid.centers.size), dtype=float)
    i_speed = speed_cell_index_for_speed(speed_grid, birth_speed_m_s)
    i_lambda = lambda_cell_index_for_lambda(lambda_grid, lambda_birth)
    pitch_measures = lambda_pitch_cell_measures(lambda_grid, B_tilde_birth=B_tilde_birth)
    cell_measure = 0.5 * speed_grid.shell_volumes_m3_s3[i_speed] * pitch_measures[i_lambda]
    if not np.isfinite(cell_measure) or cell_measure <= 0.0:
        raise ValueError("birth source cell has zero or invalid phase-space measure")
    source[i_speed, i_lambda] = float(birth_rate_density_m3_s) / cell_measure
    
    return source

def integrate_lambda_source_density(speed_grid: SpeedGrid, lambda_grid: LambdaGrid, source_v_lambda_per_s: ArrayLike, B_tilde_birth: float):
    """Integrate S(v, Lambda) over the birth Lambda pitch measure"""
    source = np.asarray(source_v_lambda_per_s, dtype=float)
    if source.shape != (speed_grid.centers_m_s.size, lambda_grid.centers.size):
        raise ValueError("source_v_lambda_per_s has inconsistent shape")
    pitch_measures = lambda_pitch_cell_measures(lambda_grid, B_tilde_birth=B_tilde_birth)
    cell_measures = 0.5 * speed_grid.shell_volumes_m3_s3[:, None] * pitch_measures[None, :]
    
    return float(np.sum(source * cell_measures))

def beam_source_definition_to_lambda_source(speed_grid: SpeedGrid, lambda_grid: LambdaGrid, source_definition: BeamSourceDefinition):
    """Convert a source definition to a cell integrated speed-Lambda source"""
    return monoenergetic_lambda_source_density(speed_grid=speed_grid, lambda_grid=lambda_grid, birth_rate_density_m3_s=source_definition.source_rate_density_m3_s, birth_speed_m_s=source_definition.speed_m_s, lambda_birth=source_definition.Lambda_birth, B_tilde_birth=source_definition.B_tilde_birth)

def beam_component_to_lambda_source(speed_grid: SpeedGrid, lambda_grid: LambdaGrid, beam_component, B_tilde_birth: float | None = None):
    """Create a Lambda grid source from a beam component like object"""
    if B_tilde_birth is None:
        if not hasattr(beam_component, "B_tilde_birth"):
            raise ValueError("B_tilde_birth must be supplied or present on beam_component")
        B_birth = float(beam_component.B_tilde_birth)
    else:
        B_birth = float(B_tilde_birth)
    
    return monoenergetic_lambda_source_density(speed_grid=speed_grid, lambda_grid=lambda_grid, birth_rate_density_m3_s=float(beam_component.source_rate_density_m3_s), birth_speed_m_s=float(beam_component.speed_m_s), lambda_birth=float(beam_component.Lambda_birth), B_tilde_birth=B_birth)

def sum_beam_components_to_lambda_source(speed_grid: SpeedGrid, lambda_grid: LambdaGrid, components: Iterable):
    """Sum cell integrated Lambda grid sources from multiple energy components"""
    total = np.zeros((speed_grid.centers_m_s.size, lambda_grid.centers.size), dtype=float)
    for component in components:
        total += beam_component_to_lambda_source(speed_grid=speed_grid, lambda_grid=lambda_grid, beam_component=component)
    return total

def lambda_source_from_power(speed_grid: SpeedGrid, lambda_grid: LambdaGrid, beam_power_W: float, beam_energy_J: float, pitch_angle_rad: float, B_tilde_birth: float, source_volume_m3: float, fast_ion_mass_kg: float):
    """Convenience builder for a monoenergetic Lambda grid source from beam power"""
    speed = float(beam_speed_from_energy_m_s(beam_energy_J=beam_energy_J, fast_ion_mass_kg=fast_ion_mass_kg))
    birth_rate = float(particle_birth_rate_from_power_s(beam_power_W=beam_power_W, beam_energy_J=beam_energy_J))
    birth_rate_density = float(source_rate_density_from_birth_rate_m3_s(birth_rate_s=birth_rate, volume_m3=source_volume_m3))
    from source_model_revamp.fbis.beam_source_definition import beam_lambda_birth
    Lambda_b = float(beam_lambda_birth(pitch_angle_rad=pitch_angle_rad, B_tilde_birth=B_tilde_birth))
    
    return monoenergetic_lambda_source_density(speed_grid=speed_grid, lambda_grid=lambda_grid, birth_rate_density_m3_s=birth_rate_density, birth_speed_m_s=speed, lambda_birth=Lambda_b, B_tilde_birth=B_tilde_birth)

def lambda_source_from_current(speed_grid: SpeedGrid, lambda_grid: LambdaGrid, beam_current_A: float, beam_energy_J: float, charge_state: float, pitch_angle_rad: float, B_tilde_birth: float, source_volume_m3: float, fast_ion_mass_kg: float):
    """Convenience builder for a monoenergetic Lambda grid source from equivalent ion current"""
    speed = float(beam_speed_from_energy_m_s(beam_energy_J=beam_energy_J, fast_ion_mass_kg=fast_ion_mass_kg))
    birth_rate = float(particle_birth_rate_from_current_s(beam_current_A=beam_current_A, charge_state=charge_state))
    birth_rate_density = float(source_rate_density_from_birth_rate_m3_s(birth_rate_s=birth_rate, volume_m3=source_volume_m3))
    from source_model_revamp.fbis.beam_source_definition import beam_lambda_birth
    Lambda_b = float(beam_lambda_birth(pitch_angle_rad=pitch_angle_rad, B_tilde_birth=B_tilde_birth))
    
    return monoenergetic_lambda_source_density(speed_grid=speed_grid, lambda_grid=lambda_grid, birth_rate_density_m3_s=birth_rate_density, birth_speed_m_s=speed, lambda_birth=Lambda_b, B_tilde_birth=B_tilde_birth)