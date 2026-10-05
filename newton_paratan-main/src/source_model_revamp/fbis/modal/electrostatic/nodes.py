"""Exact midplane and throat nodes for the axial electrostatic closure"""
from __future__ import annotations
import numpy as np
from source_model_revamp.constants import ELECTRON_CHARGE_C
from source_model_revamp.fbis.modal.electrostatic.eq70 import _electron_density_fraction_eq70

def _exact_electrostatic_node_profiles(*, cell_zeta: np.ndarray, cell_B_tilde: np.ndarray, cell_potential_energy_J: np.ndarray, cell_electron_density_m3: np.ndarray, cell_ion_density_m3: np.ndarray, cell_fast_ion_density_m3: np.ndarray, cell_background_density_m3: np.ndarray, mirror_ratio: float, left_throat_potential_energy_J: float, right_throat_potential_energy_J: float, electron_parent_maxwellian_n0_m3: float, wall_barrier_energy_J: float, electron_temperature_J: float, exact_midplane_fast_density_m3: float, exact_midplane_potential_energy_J: float, exact_left_throat_fast_density_m3: float, exact_right_throat_fast_density_m3: float, exact_midplane_background_density_m3: float, exact_left_throat_background_density_m3: float, exact_right_throat_background_density_m3: float, exact_midplane_ion_density_m3: float | None = None, exact_left_throat_ion_density_m3: float | None = None, exact_right_throat_ion_density_m3: float | None = None) -> dict[str, object]:
    """Insert exact ζ = −1, 0, 1 anchors into cell centered electrostatic profiles

    Cell arrays have shape (n_z,) and density values use m⁻³
    Anchor electron densities are recomputed directly from Eq 70 rather than interpolated
    """
    z = np.asarray(cell_zeta, dtype=float)
    anchors = np.array((-1.0, 0.0, 1.0), dtype=float)
    nodes = np.unique(np.concatenate((z, anchors)))
    source_order = np.argsort(z)
    source_z = z[source_order]

    def interpolate(values: np.ndarray) -> np.ndarray:
        """Interpolate one cell centered axial vector onto the augmented node grid"""
        array = np.asarray(values, dtype=float)[source_order]
        return np.interp(nodes, source_z, array)

    node_B = interpolate(cell_B_tilde)
    node_phi = interpolate(cell_potential_energy_J)
    node_ne = interpolate(cell_electron_density_m3)
    node_ni = interpolate(cell_ion_density_m3)
    node_fast = interpolate(cell_fast_ion_density_m3)
    node_background = interpolate(cell_background_density_m3)
    anchor_indices = [int(np.flatnonzero(nodes == value)[0]) for value in anchors]
    left_index, midplane_index, right_index = anchor_indices
    node_B[[left_index, midplane_index, right_index]] = [float(mirror_ratio), 1.0, float(mirror_ratio)]
    node_phi[[left_index, midplane_index, right_index]] = [float(left_throat_potential_energy_J), float(exact_midplane_potential_energy_J), float(right_throat_potential_energy_J)]
    node_fast[[left_index, midplane_index, right_index]] = [float(exact_left_throat_fast_density_m3), float(exact_midplane_fast_density_m3), float(exact_right_throat_fast_density_m3)]
    node_background[[left_index, midplane_index, right_index]] = [float(exact_left_throat_background_density_m3), float(exact_midplane_background_density_m3), float(exact_right_throat_background_density_m3)]
    default_anchor_ion = node_fast[[left_index, midplane_index, right_index]] + node_background[[left_index, midplane_index, right_index]]
    node_ni[[left_index, midplane_index, right_index]] = [
        default_anchor_ion[0] if exact_left_throat_ion_density_m3 is None else float(exact_left_throat_ion_density_m3),
        default_anchor_ion[1] if exact_midplane_ion_density_m3 is None else float(exact_midplane_ion_density_m3),
        default_anchor_ion[2] if exact_right_throat_ion_density_m3 is None else float(exact_right_throat_ion_density_m3),
    ]
    for index in (left_index, midplane_index, right_index):
        node_ne[index] = float(electron_parent_maxwellian_n0_m3) * _electron_density_fraction_eq70(float(node_phi[index]), float(wall_barrier_energy_J), float(electron_temperature_J))
    node_absolute = node_ne - node_ni
    node_density_scale = np.maximum(node_ne, node_ni)
    node_relative = np.divide(node_absolute, node_density_scale, out=np.zeros_like(node_absolute), where=node_density_scale > 0.0)

    return {
        "zeta": nodes,
        "B_tilde": node_B,
        "potential_energy_J": node_phi,
        "potential_relative_to_midplane_V": -(node_phi - float(exact_midplane_potential_energy_J)) / ELECTRON_CHARGE_C,
        "electron_density_m3": node_ne,
        "ion_density_m3": node_ni,
        "fast_ion_density_m3": node_fast,
        "background_positive_charge_density_m3": node_background,
        "absolute_residual_profile_m3": node_absolute,
        "relative_residual_profile": node_relative,
        "left_throat_index": left_index,
        "midplane_index": midplane_index,
        "right_throat_index": right_index,
        "left_throat_node_exact": True,
        "midplane_node_exact": True,
        "right_throat_node_exact": True,
        "midplane_gauge_exact": bool(node_phi[midplane_index] == 0.0),
        "node_model": "cell_average_profile_with_direct_quasineutral_exact_midplane_and_throat_nodes",
    }

def _symmetric_input(zeta: np.ndarray, B_tilde: np.ndarray, volumes: np.ndarray) -> bool:
    """Return whether ζ is antisymmetric and B_tilde and cell volumes are symmetric"""
    scale_z = max(float(np.max(np.abs(zeta))), 1.0)
    scale_B = max(float(np.max(np.abs(B_tilde))), 1.0)
    scale_volume = max(float(np.max(np.abs(volumes))), 1.0)

    return bool(np.allclose(zeta, -zeta[::-1], rtol=0.0, atol=1.0e-11 * scale_z) and np.allclose(B_tilde, B_tilde[::-1], rtol=0.0, atol=1.0e-11 * scale_B) and np.allclose(volumes, volumes[::-1], rtol=0.0, atol=1.0e-11 * scale_volume))

def _symmetric_average(values: np.ndarray) -> np.ndarray:
    """Average a one dimensional axial profile with its reflected ordering"""
    array = np.asarray(values, dtype=float)
    return 0.5 * (array + array[::-1])

def _midplane_density_value(zeta: np.ndarray, density: np.ndarray) -> float:
    """Return density at an exact ζ = 0 node or average the two nearest cells"""
    z = np.asarray(zeta, dtype=float)
    n = np.asarray(density, dtype=float)
    exact = np.flatnonzero(np.isclose(z, 0.0, rtol=0.0, atol=1.0e-14))
    if exact.size:
        return float(np.mean(n[exact]))
    order = np.argsort(np.abs(z))
    count = min(2, z.size)

    return float(np.mean(n[order[:count]]))

def _enforce_midplane_reference(zeta: np.ndarray, phi: np.ndarray) -> np.ndarray:
    """Set exact ζ = 0 entries of the potential energy coordinate to zero"""
    result = np.asarray(phi, dtype=float).copy()
    exact = np.flatnonzero(np.isclose(np.asarray(zeta, dtype=float), 0.0, rtol=0.0, atol=1.0e-14))
    result[exact] = 0.0

    return result
