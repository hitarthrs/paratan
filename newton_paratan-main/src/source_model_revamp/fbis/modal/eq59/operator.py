"""Assemble and solve the finite volume speed operator for Egedal Eq 59"""
from __future__ import annotations
from collections import OrderedDict
from dataclasses import dataclass
import numpy as np
from numpy.typing import ArrayLike
from scipy.linalg import solve_banded
from source_model_revamp.fbis.modal.types import ModalRosenbluthCoefficients
from source_model_revamp.fbis.modal.eq59.collision_operator_state import Eq59CollisionOperatorState
from source_model_revamp.fbis.modal.utils import _EPS
from source_model_revamp.fbis.modal.velocity_solve_cold import _source_cell_weights
from source_model_revamp.fbis.rosenbluth_potentials import isotropic_rosenbluth_potentials
from source_model_revamp.fbis.velocity_space_grid import SpeedGrid

@dataclass(frozen=True)
class _Eq59OperatorTemplate:
    """Cache the mode independent pieces of the discrete Eq 59 operator
    
    lower base_diagonal and upper describe speed flux divergence
    pitch and charge exchange coefficients are added to each mode diagonal
    """
    collision_operator_state: Eq59CollisionOperatorState
    source_speeds_m_s: np.ndarray
    charge_exchange_loss_frequency_s: np.ndarray
    lower: np.ndarray
    base_diagonal: np.ndarray
    upper: np.ndarray
    pitch_diagonal_coefficient: np.ndarray
    charge_exchange_diagonal_coefficient: np.ndarray
    source_cell_weights: np.ndarray
    slowing_time_s: float

_OPERATOR_TEMPLATE_CACHE: OrderedDict[tuple[object, ...], _Eq59OperatorTemplate] = OrderedDict()
_SOURCE_WEIGHT_CACHE: OrderedDict[tuple[object, ...], tuple[SpeedGrid, np.ndarray]] = OrderedDict()
_TEMPLATE_CACHE_SIZE = 24

def _cache_insert(cache: OrderedDict, key: tuple[object, ...], value: object) -> object:
    """Insert one cached value and enforce the bounded least recently used cache size"""
    cache[key] = value
    cache.move_to_end(key)
    while len(cache) > _TEMPLATE_CACHE_SIZE:
        cache.popitem(last=False)
  
    return value

def _roundoff_nonnegative(values: ArrayLike, name: str) -> np.ndarray:
    """Validate one coefficient array and clip negative roundoff to zero"""
    array = np.asarray(values, dtype=float)
    if np.any(~np.isfinite(array)):
        raise ValueError(f"{name} must be finite")
    scale = max(float(np.max(np.abs(array))) if array.size else 0.0, 1.0)
    tolerance = 256.0 * np.finfo(float).eps * scale
    if float(np.min(array)) < -tolerance:
        raise ValueError(f"{name} contains material negative values")
   
    return np.where(array < 0.0, 0.0, array)

def _density_normalized_rosenbluth_coefficients(speed_grid: SpeedGrid, f1_v: ArrayLike) -> ModalRosenbluthCoefficients:
    """Return density normalized Rosenbluth coefficients from the first modal speed amplitude
    
    h_tilde = −v² ∂h/∂v
    g_tilde_1 = ∂g/∂v
    g_tilde_2 = v² ∂²g/∂v²
    """
    v = np.asarray(speed_grid.centers_m_s, dtype=float)
    f1 = np.asarray(f1_v, dtype=float)
    if f1.shape != v.shape or np.any(~np.isfinite(f1)):
        raise ValueError("f1_v must have one finite value per speed cell")
    f1_scale = max(float(np.max(np.abs(f1))) if f1.size else 0.0, 1.0)
    negative_tolerance = 256.0 * np.finfo(float).eps * f1_scale
    if float(np.min(f1)) < -negative_tolerance:
        raise ValueError("first mode distribution contains material negative values")
    f1 = np.where(f1 < 0.0, 0.0, f1)
    potentials = isotropic_rosenbluth_potentials(v, f1)
    density = float(potentials.density_m3)
    if not np.isfinite(density) or density <= 0.0:
        raise ValueError("first mode distribution must have positive density for Rosenbluth normalization")
    dH = potentials.dH_dv / density
    dG = potentials.dG_dv / density
    d2G = potentials.d2G_dv2 / density
    h_tilde = _roundoff_nonnegative(-v**2 * dH, "h_tilde")
    g1 = _roundoff_nonnegative(dG, "g_tilde_1")
    g2 = _roundoff_nonnegative(v**2 * d2G, "g_tilde_2")
  
    return ModalRosenbluthCoefficients(h_tilde=h_tilde, g_tilde_1=g1, g_tilde_2=g2, density_normalization_m3=density)

def _face_average(values: np.ndarray) -> np.ndarray:
    """Map one cell centered coefficient array to arithmetic face values"""
    faces = np.empty(values.size + 1, dtype=float)
    faces[1:-1] = 0.5 * (values[:-1] + values[1:])
    faces[0] = values[0]
    faces[-1] = values[-1]
 
    return faces

def _scharfetter_gummel_bernoulli(argument: ArrayLike) -> np.ndarray:
    """Stable Bernoulli function B(x)=x/(exp(x)-1) for SG fluxes"""
    x = np.asarray(argument, dtype=float)
    result = np.empty_like(x)
    small = np.abs(x) < 1.0e-4
    large_positive = x > 50.0
    large_negative = x < -50.0
    ordinary = ~(small | large_positive | large_negative)
    xs = x[small]
    result[small] = 1.0 - 0.5 * xs + xs**2 / 12.0 - xs**4 / 720.0
    result[large_positive] = x[large_positive] * np.exp(-x[large_positive])
    result[large_negative] = -x[large_negative]
    result[ordinary] = x[ordinary] / np.expm1(x[ordinary])
 
    return result

def _scharfetter_gummel_face_coefficients(drift: float, diffusion: float, delta_v: float) -> tuple[float, float]:
    """Return left and right coefficients for F=b f + D df/dv"""
    b = float(drift)
    D = float(diffusion)
    dv = float(delta_v)
    if not np.isfinite(b) or not np.isfinite(D) or not np.isfinite(dv) or D < 0.0 or dv <= 0.0:
        raise ValueError("Scharfetter Gummel drift, diffusion, and spacing are invalid")
    if D == 0.0:
        if b > 0.0:
            return 0.0, b
        if b < 0.0:
            return b, 0.0
        return 0.0, 0.0
    b_dv_over_D = b * dv / D
    left = -(D / dv) * float(_scharfetter_gummel_bernoulli(np.asarray([b_dv_over_D]))[0])
    right = (D / dv) * float(_scharfetter_gummel_bernoulli(np.asarray([-b_dv_over_D]))[0])
   
    return left, right

def _validated_operator_arrays(speed_grid: SpeedGrid, collision_operator_state: Eq59CollisionOperatorState) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Validate that collision coefficient arrays exactly match the active speed grid"""
    speed = np.asarray(speed_grid.centers_m_s, dtype=float)
    state_speed = np.asarray(collision_operator_state.speed_grid.centers_m_s, dtype=float)
    state_faces = np.asarray(collision_operator_state.speed_grid.faces_m_s, dtype=float)
    if speed.shape != state_speed.shape or not np.allclose(speed, state_speed, rtol=0.0, atol=0.0) or not np.array_equal(np.asarray(speed_grid.faces_m_s, dtype=float), state_faces):
        raise ValueError("Eq 59 collision operator speed grid must match the solve grid exactly")
    shape = speed.shape
    arrays: list[np.ndarray] = []
    for name, values in (("ion_drag_velocity_cubed_m3_s3", collision_operator_state.ion_drag_velocity_cubed_m3_s3), ("ion_energy_diffusion_velocity_fourth_m4_s4", collision_operator_state.ion_energy_diffusion_velocity_fourth_m4_s4), ("ion_pitch_scattering_velocity_cubed_m3_s3", collision_operator_state.ion_pitch_scattering_velocity_cubed_m3_s3)):
        array = np.asarray(values, dtype=float)
        if array.shape != shape or np.any(~np.isfinite(array)) or np.any(array < 0.0):
            raise ValueError(f"{name} must contain one finite nonnegative value per speed cell")
        arrays.append(array)
   
    return arrays[0], arrays[1], arrays[2]

def _source_weights_for_components(speed_grid: SpeedGrid, source_speeds_m_s: np.ndarray) -> np.ndarray:
    """Return cached spherical shell source weights for each monoenergetic component"""
    speeds = np.asarray(source_speeds_m_s, dtype=float)
    key = (id(speed_grid), speeds.shape, speeds.tobytes())
    cached = _SOURCE_WEIGHT_CACHE.get(key)
    if cached is not None and cached[0] is speed_grid:
        _SOURCE_WEIGHT_CACHE.move_to_end(key)
        return cached[1]
    weights = np.asarray([_source_cell_weights(speed_grid, float(value)) for value in speeds], dtype=float)
 
    return _cache_insert(_SOURCE_WEIGHT_CACHE, key, (speed_grid, weights))[1]

def _prepare_eq59_operator_template(*, speed_grid: SpeedGrid, source_speeds_m_s: np.ndarray, collision_operator_state: Eq59CollisionOperatorState, charge_exchange_loss_frequency_s: ArrayLike | None) -> _Eq59OperatorTemplate:
    """Build the mode independent finite volume speed operator template
    
    Internal face fluxes use Scharfetter Gummel fitting and boundary face fluxes remain zero
    The template also stores pitch scattering charge exchange and source cell factors
    """
    v = np.asarray(speed_grid.centers_m_s, dtype=float)
    widths = np.asarray(speed_grid.widths_m_s, dtype=float)
    source_speeds = np.asarray(source_speeds_m_s, dtype=float)
    n = v.size
    tau_s = float(collision_operator_state.spitzer_slowing_down_time_s)
    if n < 2:
        raise ValueError("Eq 59 requires at least two speed cells")
    if not np.isfinite(tau_s) or tau_s <= 0.0:
        raise ValueError("Eq 59 slowing time must be positive and finite")
    if source_speeds.ndim != 1 or not np.all(np.isfinite(source_speeds)) or np.any(source_speeds <= 0.0):
        raise ValueError("source speeds must be one dimensional positive and finite")
    ion_drag, diffusion_center, pitch = _validated_operator_arrays(speed_grid, collision_operator_state)
    if charge_exchange_loss_frequency_s is None:
        charge_exchange_frequency = np.zeros(n, dtype=float)
    else:
        charge_exchange_frequency = np.asarray(charge_exchange_loss_frequency_s, dtype=float)
        if charge_exchange_frequency.shape != (n,) or np.any(~np.isfinite(charge_exchange_frequency)) or np.any(charge_exchange_frequency < 0.0):
            raise ValueError("charge exchange loss frequency must be finite and nonnegative on the Eq 59 speed grid")
    key = (id(collision_operator_state), source_speeds.shape, source_speeds.tobytes(), charge_exchange_frequency.tobytes())
    cached = _OPERATOR_TEMPLATE_CACHE.get(key)
    if cached is not None and cached.collision_operator_state is collision_operator_state:
        _OPERATOR_TEMPLATE_CACHE.move_to_end(key)
        return cached
    b_face = _face_average(v**3 + ion_drag)
    diffusion_face = _face_average(diffusion_center)
    lower = np.zeros(n - 1, dtype=float)
    diagonal = np.zeros(n, dtype=float)
    upper = np.zeros(n - 1, dtype=float)
    face_left_coefficients = np.zeros(n + 1, dtype=float)
    face_right_coefficients = np.zeros(n + 1, dtype=float)
    # Boundary face flux is zero because only internal faces receive fitted coefficients
    for face_index in range(1, n):
        delta_v = v[face_index] - v[face_index - 1]
        left_coefficient, right_coefficient = _scharfetter_gummel_face_coefficients(b_face[face_index], diffusion_face[face_index], delta_v)
        face_left_coefficients[face_index] = left_coefficient
        face_right_coefficients[face_index] = right_coefficient
    for index in range(n):
        if index > 0:
            lower[index - 1] -= face_left_coefficients[index]
            diagonal[index] -= face_right_coefficients[index]
        if index < n - 1:
            diagonal[index] += face_left_coefficients[index + 1]
            upper[index] += face_right_coefficients[index + 1]
    template = _Eq59OperatorTemplate(
        collision_operator_state=collision_operator_state,
        source_speeds_m_s=source_speeds,
        charge_exchange_loss_frequency_s=charge_exchange_frequency,
        lower=lower,
        base_diagonal=diagonal,
        upper=upper,
        pitch_diagonal_coefficient=pitch / np.maximum(v, _EPS) * widths,
        charge_exchange_diagonal_coefficient=tau_s * charge_exchange_frequency * v**2 * widths,
        source_cell_weights=_source_weights_for_components(speed_grid, source_speeds),
        slowing_time_s=tau_s,
    )
   
    return _cache_insert(_OPERATOR_TEMPLATE_CACHE, key, template)

def _assemble_eq59_mode_tridiagonal(*, speed_grid: SpeedGrid, source_coefficients_by_component: np.ndarray, source_speeds_m_s: np.ndarray, eigenvalue: ArrayLike, collision_operator_state: Eq59CollisionOperatorState, charge_exchange_loss_frequency_s: ArrayLike | None = None) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Assemble one modal tridiagonal Eq 59 system from the shared speed template
    
    A scalar λ_j or speed dependent λ_j(v) multiplies the pitch scattering sink for the selected mode
    """
    source_coefficients = np.asarray(source_coefficients_by_component, dtype=float)
    source_speeds = np.asarray(source_speeds_m_s, dtype=float)
    n = speed_grid.centers_m_s.size
    lam = np.asarray(eigenvalue, dtype=float)
    if lam.ndim == 0:
        lam = np.full(n, float(lam), dtype=float)
    if lam.shape != (n,) or np.any(~np.isfinite(lam)) or np.any(lam <= 0.0):
        raise ValueError("eigenvalue must be positive and finite with one value per speed cell")
    if source_coefficients.ndim != 1 or source_speeds.shape != source_coefficients.shape:
        raise ValueError("source coefficients and speeds must be matching 1D arrays")
    if not np.all(np.isfinite(source_coefficients)):
        raise ValueError("source coefficients must be finite")
    template = _prepare_eq59_operator_template(speed_grid=speed_grid, source_speeds_m_s=source_speeds, collision_operator_state=collision_operator_state, charge_exchange_loss_frequency_s=charge_exchange_loss_frequency_s)
    diagonal = template.base_diagonal.copy()
    for index in range(n):
        diagonal[index] -= lam[index] * template.pitch_diagonal_coefficient[index]
        diagonal[index] -= template.charge_exchange_diagonal_coefficient[index]
    rhs = np.zeros(n, dtype=float)
    for component_index, coefficient in enumerate(source_coefficients):
        rhs -= template.slowing_time_s * float(coefficient) * template.source_cell_weights[component_index]
   
    return template.lower, diagonal, template.upper, rhs

def _tridiagonal_banded_matrix(lower: np.ndarray, diagonal: np.ndarray, upper: np.ndarray) -> np.ndarray:
    """Pack three diagonals in the format required by scipy solve_banded"""
    banded = np.zeros((3, diagonal.size), dtype=float)
    banded[0, 1:] = upper
    banded[1] = diagonal
    banded[2, :-1] = lower
   
    return banded

def _solve_eq59_mode_finite_difference(*, speed_grid: SpeedGrid, source_coefficients_by_component: np.ndarray, source_speeds_m_s: np.ndarray, eigenvalue: ArrayLike, collision_operator_state: Eq59CollisionOperatorState, charge_exchange_loss_frequency_s: ArrayLike | None = None) -> np.ndarray:
    """Solve one Eq 59 modal speed amplitude with the tridiagonal finite volume operator"""
    lower, diagonal, upper, rhs = _assemble_eq59_mode_tridiagonal(speed_grid=speed_grid, source_coefficients_by_component=source_coefficients_by_component, source_speeds_m_s=source_speeds_m_s, eigenvalue=eigenvalue, collision_operator_state=collision_operator_state, charge_exchange_loss_frequency_s=charge_exchange_loss_frequency_s)
    banded = _tridiagonal_banded_matrix(lower, diagonal, upper)
    try:
        solution = solve_banded((1, 1), banded, rhs, overwrite_ab=True, overwrite_b=False, check_finite=False)
    except np.linalg.LinAlgError as exc:
        raise RuntimeError("Eq 59 tridiagonal operator is singular") from exc
    if not np.all(np.isfinite(solution)):
        raise RuntimeError("Eq 59 tridiagonal solve produced nonfinite values")
  
    return np.asarray(solution, dtype=float)

def _apply_tridiagonal(lower: np.ndarray, diagonal: np.ndarray, upper: np.ndarray, values: np.ndarray) -> np.ndarray:
    """Apply one tridiagonal operator without constructing a dense matrix"""
    result = diagonal * values
    result[1:] += lower * values[:-1]
    result[:-1] += upper * values[1:]
  
    return result

def _eq59_nonlinear_residual(*, speed_grid: SpeedGrid, modal_distribution: np.ndarray, source_coefficients_by_component_j: np.ndarray, source_speeds_m_s: np.ndarray, eigenvalues: np.ndarray, collision_operator_state: Eq59CollisionOperatorState, charge_exchange_loss_frequency_s: ArrayLike | None = None) -> tuple[float, float]:
    """Return maximum absolute and normwise backward residuals over all retained modes"""
    maximum_absolute = 0.0
    maximum_relative = 0.0
    tiny = np.finfo(float).tiny
    eigenvalue_array = np.asarray(eigenvalues, dtype=float)
    if eigenvalue_array.ndim == 1:
        eigenvalue_array = np.broadcast_to(eigenvalue_array[:, None], (eigenvalue_array.size, speed_grid.centers_m_s.size))
    expected_shape = (modal_distribution.shape[0], speed_grid.centers_m_s.size)
    if eigenvalue_array.shape != expected_shape:
        raise ValueError("eigenvalues must have shape (n_mode,) or (n_mode, n_speed_cell)")
    for mode_index, eigenvalue in enumerate(eigenvalue_array):
        lower, diagonal, upper, rhs = _assemble_eq59_mode_tridiagonal(speed_grid=speed_grid, source_coefficients_by_component=source_coefficients_by_component_j[:, mode_index], source_speeds_m_s=source_speeds_m_s, eigenvalue=eigenvalue, collision_operator_state=collision_operator_state, charge_exchange_loss_frequency_s=charge_exchange_loss_frequency_s)
        values = modal_distribution[mode_index]
        residual = _apply_tridiagonal(lower, diagonal, upper, values) - rhs
        absolute = float(np.max(np.abs(residual)))
        row_sum = np.abs(diagonal)
        row_sum[1:] += np.abs(lower)
        row_sum[:-1] += np.abs(upper)
        scale = float(np.max(row_sum)) * float(np.max(np.abs(values))) + float(np.max(np.abs(rhs)))
        relative = absolute / max(scale, tiny)
        maximum_absolute = max(maximum_absolute, absolute)
        maximum_relative = max(maximum_relative, relative)
   
    return maximum_absolute, maximum_relative
