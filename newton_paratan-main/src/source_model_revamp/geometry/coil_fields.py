"""On axis fields and field gradients for effective circular coils

Evaluate supplied coil strengths and sample their mirror metrics
Coil placement and strength fitting are handled by integration/config/geometry.py
"""
from __future__ import annotations
from dataclasses import dataclass
from collections.abc import Sequence
import math
import numpy as np
from numpy.typing import ArrayLike
from scipy.constants import mu_0

@dataclass(frozen=True)
class CircularCoilOnAxis:
    """Effective thin circular coil for on axis Bz calculations

    Positions and radii are in meters
    amp_turns is the signed product of current and turn count
    group labels the coil family
    """
    name: str
    z_center_m: float
    radius_m: float
    amp_turns: float
    group: str = "coil"

def circular_coil_Bz_on_axis_T(z_m: ArrayLike, coil: CircularCoilOnAxis):
    """Return the on axis circular loop field in tesla

    Use the thin loop Biot Savart expression with amp_turns = N * I
    """
    z = np.asarray(z_m, dtype=float)
    R = float(coil.radius_m)
    dz = z - float(coil.z_center_m)

    return mu_0 * float(coil.amp_turns) * R**2 / (2.0 * (R**2 + dz**2) ** 1.5)

def circular_coil_dBz_dz_on_axis_T_per_m(z_m: ArrayLike, coil: CircularCoilOnAxis):
    """Axial derivative of the on axis circular coil field"""
    z = np.asarray(z_m, dtype=float)
    R = float(coil.radius_m)
    dz = z - float(coil.z_center_m)

    return -3.0 * mu_0 * float(coil.amp_turns) * R**2 * dz / (2.0 * (R**2 + dz**2) ** 2.5)

def coil_set_Bz_on_axis_T(z_m: ArrayLike, coils: Sequence[CircularCoilOnAxis]):
    """Superpose on axis Bz from a set of effective circular coils"""
    z = np.asarray(z_m, dtype=float)
    terms = [circular_coil_Bz_on_axis_T(z, coil) for coil in coils]
    if not terms:
        return np.zeros_like(z, dtype=float)

    # Accumulate coil contributions at extended precision before returning float values
    return np.asarray(np.sum(np.asarray(terms, dtype=np.longdouble), axis=0, dtype=np.longdouble), dtype=float)

def coil_set_dBz_dz_on_axis_T_per_m(z_m: ArrayLike, coils: Sequence[CircularCoilOnAxis]):
    """Superpose dBz/dz from a set of effective circular coils"""
    z = np.asarray(z_m, dtype=float)
    total = np.zeros_like(z, dtype=float)
    for coil in coils:
        total = total + circular_coil_dBz_dz_on_axis_T_per_m(z, coil)

    return total

def mirror_metrics_from_coil_set(coils: Sequence[CircularCoilOnAxis], z_min_m: float, z_max_m: float, reference_z_m: float = 0.0, grid_points: int = 1001):
    """Return the sampled peak field, its location, and its ratio to the reference field

    The peak is selected from the supplied interval grid without further refinement
    """
    z = np.linspace(float(z_min_m), float(z_max_m), int(grid_points))
    B = coil_set_Bz_on_axis_T(z, coils)
    B_ref = float(coil_set_Bz_on_axis_T(np.asarray([reference_z_m]), coils)[0])
    imax = int(np.argmax(B)) if B.size else 0
    B_max = float(B[imax]) if B.size else 0.0
    mirror_ratio = float(B_max / B_ref) if B_ref != 0.0 else math.inf

    return {
        "reference_z_m": float(reference_z_m),
        "B_reference_T": B_ref,
        "B_max_T": B_max,
        "mirror_ratio": mirror_ratio,
        "max_location_z_m": float(z[imax]) if B.size else None,
        "z_min_m": float(z_min_m),
        "z_max_m": float(z_max_m),
        "grid_points": int(grid_points),
    }
