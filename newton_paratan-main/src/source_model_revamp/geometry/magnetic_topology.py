"""Magnetic throat discovery and central well validation

Refine field maxima on each side and check the field between them
Also record sampled extrema in the two outer search intervals
"""
from __future__ import annotations
from collections.abc import Callable
from dataclasses import dataclass
import numpy as np
from scipy.optimize import minimize_scalar

@dataclass(frozen=True)
class MagneticThroatPair:
    """Independent left and right throats with central well checks

    Mirror ratios use the field at the supplied midplane
    Expander extrema are sampled locations rather than refined roots
    """
    left_z_m: float
    right_z_m: float
    left_B_T: float
    right_B_T: float
    left_mirror_ratio: float
    right_mirror_ratio: float
    monotonic_left: bool
    monotonic_right: bool
    left_expander_extrema_m: tuple[float, ...]
    right_expander_extrema_m: tuple[float, ...]

def _side_maximum(B_T: Callable[[float], float], start_m: float, stop_m: float, samples: int) -> tuple[float, float]:
    """Refine sampled interior maxima and select the strongest one

    Return values are ordered as position in meters, then field in tesla
    """
    grid = np.linspace(start_m, stop_m, max(int(samples), 101))
    values = np.asarray([B_T(float(z)) for z in grid], dtype=float)
    if np.any(~np.isfinite(values)) or np.any(values <= 0.0):
        raise ValueError("coil field must be finite and positive throughout throat search")
    candidates: list[tuple[float, float]] = []
    # Search interior peaks so an interval endpoint cannot become a throat
    for index in range(1, grid.size - 1):
        if values[index] >= values[index - 1] and values[index] >= values[index + 1]:
            lower, upper = sorted((float(grid[index - 1]), float(grid[index + 1])))
            result = minimize_scalar(lambda z: -B_T(float(z)), bounds=(lower, upper), method="bounded", options={"xatol": 1.0e-12})
            candidates.append((float(result.x), float(-result.fun)))
    if not candidates:
        raise ValueError("no interior magnetic throat was found on one side of the midplane")

    return max(candidates, key=lambda item: item[1])

def _extrema_locations(B_T: Callable[[float], float], z_min_m: float, z_max_m: float, samples: int) -> tuple[float, ...]:
    """Return sampled positions where adjacent field slopes change sign"""
    if z_max_m <= z_min_m:
        return ()
    grid = np.linspace(z_min_m, z_max_m, max(int(samples), 51))
    values = np.asarray([B_T(float(z)) for z in grid], dtype=float)
    slope = np.diff(values)
    locations: list[float] = []
    for i in range(1, slope.size):
        if slope[i - 1] * slope[i] < 0.0:
            locations.append(float(grid[i]))

    return tuple(locations)

def discover_magnetic_throats(B_T: Callable[[float], float], *, midplane_z_m: float, search_z_min_m: float, search_z_max_m: float, samples: int = 2001, requested_mirror_ratio: float | None = None, mirror_ratio_relative_tolerance: float = 5.0e-3) -> MagneticThroatPair:
    """Locate both throats and require a monotonic field rise away from the midplane

    Check each throat ratio against the requested value when one is supplied
    """
    if not search_z_min_m < midplane_z_m < search_z_max_m:
        raise ValueError("magnetic throat search domain must bracket the midplane")
    left_z, left_B = _side_maximum(B_T, midplane_z_m, search_z_min_m, samples // 2)
    right_z, right_B = _side_maximum(B_T, midplane_z_m, search_z_max_m, samples // 2)
    if left_z >= midplane_z_m or right_z <= midplane_z_m:
        raise ValueError("independent throat discovery did not bracket the midplane")
    B0 = float(B_T(midplane_z_m))
    if not np.isfinite(B0) or B0 <= 0.0:
        raise ValueError("midplane field must be finite and positive")
    left_ratio = left_B / B0
    right_ratio = right_B / B0
    # Both arrays increase in z, so the left and right field slopes have opposite signs
    left_grid = np.linspace(left_z, midplane_z_m, max(samples // 2, 101))
    right_grid = np.linspace(midplane_z_m, right_z, max(samples // 2, 101))
    left_values = np.asarray([B_T(float(z)) for z in left_grid])
    right_values = np.asarray([B_T(float(z)) for z in right_grid])
    field_scale = max(left_B, right_B, B0)
    tolerance = 1.0e-9 * field_scale
    monotonic_left = bool(np.all(np.diff(left_values) <= tolerance))
    monotonic_right = bool(np.all(np.diff(right_values) >= -tolerance))
    if not monotonic_left or not monotonic_right:
        raise ValueError("fitted field is not a single monotonic magnetic well between midplane and throats")
    if requested_mirror_ratio is not None:
        errors = (abs(left_ratio / requested_mirror_ratio - 1.0), abs(right_ratio / requested_mirror_ratio - 1.0))
        if max(errors) > float(mirror_ratio_relative_tolerance):
            raise ValueError("independent fitted throat ratio does not match requested mirror ratio:" f"left={left_ratio:.12g}, right={right_ratio:.12g}, " f"requested={requested_mirror_ratio:.12g}")

    return MagneticThroatPair(
        left_z_m=left_z,
        right_z_m=right_z,
        left_B_T=left_B,
        right_B_T=right_B,
        left_mirror_ratio=left_ratio,
        right_mirror_ratio=right_ratio,
        monotonic_left=monotonic_left,
        monotonic_right=monotonic_right,
        left_expander_extrema_m=_extrema_locations(B_T, search_z_min_m, left_z, samples // 2),
        right_expander_extrema_m=_extrema_locations(B_T, right_z, search_z_max_m, samples // 2),
    )
