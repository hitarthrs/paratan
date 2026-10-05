"""
Straight neutral beam path geometry through axisymmetric mirror plasma

Coordinates
The mirror axis is the z axis, a beamline is represented by a finite straight segment
    x(s) = x_start + u * s
    0 <= s <= L
with "u" being the unit vector from the start point to the end point of the beam
"""
from __future__ import annotations
from dataclasses import dataclass
import math
import numpy as np
from source_model_revamp.numerical_quadrature import gauss_legendre_rule
from numpy.typing import ArrayLike
from source_model_revamp.geometry.axial_grid_geometry import cell_centers_from_faces, flux_tube_cell_volumes_m3
_EPS = 1.0e-14

@dataclass(frozen=True)
class BeamPathOnAxisymmetricGrid:
    """
    Geometry factors for one straight beam on an axial mirror grid
    
    Axial arrays have one value per magnetic grid cell
    physical_path_lengths_m stores centerline distance within each axial cell
    radial_overlap_fractions stores the path averaged beam footprint fraction inside the local plasma cross section
    effective_path_lengths_m is their product and controls attenuation
    cell_indices_along_beam contains only cells with positive effective length ordered by increasing beam path coordinate
    """
    coordinate_centers_m: np.ndarray
    coordinate_faces_m: np.ndarray
    B_tilde_centers: np.ndarray
    plasma_radius_m: np.ndarray
    cell_volumes_m3: np.ndarray
    reference_area_m2: float
    start_point_m: np.ndarray
    end_point_m: np.ndarray
    direction_unit: np.ndarray
    total_length_m: float
    cell_indices_along_beam: np.ndarray
    physical_path_lengths_m: np.ndarray
    radial_overlap_fractions: np.ndarray
    effective_path_lengths_m: np.ndarray
    beam_center_points_m: np.ndarray
    beam_radius_profile_m: np.ndarray
    beam_radius_m: float
    divergence_half_angle_rad: float
    radial_overlap_model: str

def unit_vector(vector: ArrayLike) -> np.ndarray:
    """Return the normalized direction of a finite vector"""
    vec = np.asarray(vector, dtype=float)
    norm = float(np.linalg.norm(vec))
    if norm <= 0.0:
        raise ValueError("Beam direction vector must have non zero length")
    
    return vec / norm

def circle_intersection_area(radius_a: float, radius_b: float, separation: float) -> float:
    """
    Return the geometric intersection area of two circles
    
    Radii and center separation use the same length unit and the returned area uses that unit squared
    """
    ra = max(float(radius_a), 0.0)
    rb = max(float(radius_b), 0.0)
    d = abs(float(separation))
    if ra <= 0.0 or rb <= 0.0:
        return 0.0
    if d >= ra + rb:
        return 0.0
    if d <= abs(ra - rb):
        return math.pi * min(ra, rb) ** 2
    cos_a = (d * d + ra * ra - rb * rb) / (2.0 * d * ra)
    cos_b = (d * d + rb * rb - ra * ra) / (2.0 * d * rb)
    cos_a = min(1.0, max(-1.0, cos_a))
    cos_b = min(1.0, max(-1.0, cos_b))
    term_a = ra * ra * math.acos(cos_a)
    term_b = rb * rb * math.acos(cos_b)
    term_c = 0.5 * math.sqrt(max(0.0, (-d + ra + rb) * (d + ra - rb) * (d - ra + rb) * (d + ra + rb)))

    return term_a + term_b - term_c

def beam_overlap_fraction(beam_center_radius_m: float, beam_radius_m: float, plasma_radius_m: float, *, model: str = "hard_edge_overlap") -> float:
    """
    Return the fraction of a circular beam footprint inside a circular plasma cross section
    
    The no clip aliases return one
    The centerline aliases test only the beam axis
    The hard edge overlap aliases divide the circle intersection area by the beam footprint area
    """
    model_norm = str(model).strip().lower()
    rho = abs(float(beam_center_radius_m))
    rb = max(float(beam_radius_m), 0.0)
    rp = max(float(plasma_radius_m), 0.0)
    if model_norm in {"none", "ignore", "no_radial_clip"}:
        return 1.0
    if model_norm in {"centerline", "ray", "beam_axis"} or rb <= 0.0:
        return 1.0 if rho <= rp else 0.0
    if model_norm in {"hard_edge", "hard_edge_overlap", "circle_overlap", "area_overlap"}:
        beam_area = math.pi * rb * rb
        if beam_area <= 0.0:
            return 1.0 if rho <= rp else 0.0
        
        return float(circle_intersection_area(rb, rp, rho) / beam_area)

def _s_interval_for_z_cell(start_z: float, direction_z: float, length_m: float, z0: float, z1: float):
    """
    Return the beam path coordinate interval lying inside one axial cell
    
    The result is clipped to the finite beam segment from zero to length_m
    None indicates that the segment does not cross the cell
    """
    z_lower = min(float(z0), float(z1))
    z_upper = max(float(z0), float(z1))
    if abs(float(direction_z)) <= _EPS:
        if z_lower <= float(start_z) < z_upper:
            return 0.0, float(length_m)
        return None
    s_a = (z_lower - float(start_z)) / float(direction_z)
    s_b = (z_upper - float(start_z)) / float(direction_z)
    s0 = max(0.0, min(s_a, s_b))
    s1 = min(float(length_m), max(s_a, s_b))
    if s1 <= s0:
        return None
    
    return s0, s1

def _cell_average_overlap_fraction(*, start_point_m: np.ndarray, direction_unit: np.ndarray, s0: float, s1: float, plasma_radius_m: float, beam_radius_m: float, tan_divergence: float, radial_overlap_model: str, quadrature_order: int = 16) -> float:
    """
    Return the path averaged beam and plasma overlap fraction in one axial cell
    
    Gauss Legendre quadrature integrates the footprint overlap along the beam path
    The beam radius grows linearly with path coordinate through tan_divergence
    """
    ds = float(s1) - float(s0)
    if ds <= 0.0:
        return 0.0
    nodes, weights = gauss_legendre_rule(int(max(2, quadrature_order)))
    s_values = 0.5 * (float(s0) + float(s1)) + 0.5 * ds * nodes
    total = 0.0
    for s, weight in zip(s_values, weights):
        point = start_point_m + direction_unit * float(s)
        rho = float(np.hypot(point[0], point[1]))
        beam_radius_here = max(float(beam_radius_m) + float(tan_divergence) * float(s), 0.0)
        total += float(weight) * beam_overlap_fraction(beam_center_radius_m=rho, beam_radius_m=beam_radius_here, plasma_radius_m=float(plasma_radius_m), model=radial_overlap_model)
  
    return min(1.0, max(0.0, 0.5 * total))

def beam_path_on_axisymmetric_grid(coordinate_faces_m: ArrayLike, B_tilde_centers: ArrayLike, *, start_point_m: ArrayLike, end_point_m: ArrayLike, reference_area_m2: float, plasma_radius_m: float, beam_radius_m: float = 0.0, divergence_half_angle_rad: float = 0.0, radial_overlap_model: str = "hard_edge_overlap", overlap_quadrature_order: int = 16) -> BeamPathOnAxisymmetricGrid:
    """
    Map a finite straight beam onto the axisymmetric magnetic grid
    
    coordinate_faces_m defines the axial cells and B_tilde_centers gives normalized field magnitude at their centers
    The local plasma radius scales as the supplied midplane radius divided by sqrt B_tilde
    Flux tube cell volumes use reference_area_m2 and the same magnetic profile
    The returned effective path length in each cell is the physical path integral of the radial overlap fraction
    """
    faces = np.asarray(coordinate_faces_m, dtype=float)
    if faces.ndim != 1 or faces.size < 2:
        raise ValueError("coordinate_faces_m must be a 1D array with at least two entries")
    B_tilde = np.asarray(B_tilde_centers, dtype=float)
    centers = cell_centers_from_faces(faces)
    if B_tilde.size != centers.size:
        raise ValueError("B_tilde_centers must have one value per axial cell")
    start = np.asarray(start_point_m, dtype=float)
    end = np.asarray(end_point_m, dtype=float)
    if start.shape != (3,) or end.shape != (3,):
        raise ValueError("Beam start_point_m and end_point_m must be length 3 vectors")
    delta = end - start
    length = float(np.linalg.norm(delta))
    if length <= 0.0:
        raise ValueError("Beam start and end points must be different")
    direction = delta / length
    volumes = flux_tube_cell_volumes_m3(faces, B_tilde, reference_area_m2=reference_area_m2)
    # Flux conservation expands the local radius as B_tilde decreases
    plasma_radius_profile = float(plasma_radius_m) / np.sqrt(np.maximum(B_tilde, np.finfo(float).tiny))
    physical_lengths = np.zeros(centers.size, dtype=float)
    overlap = np.zeros(centers.size, dtype=float)
    effective_lengths = np.zeros(centers.size, dtype=float)
    beam_centers = np.full((centers.size, 3), np.nan, dtype=float)
    beam_radii = np.zeros(centers.size, dtype=float)
    s_mid_values = np.full(centers.size, np.nan, dtype=float)
    tan_div = math.tan(max(float(divergence_half_angle_rad), 0.0))
    for i in range(centers.size):
        interval = _s_interval_for_z_cell(start[2], direction[2], length, faces[i], faces[i + 1])
        if interval is None:
            continue
        s0, s1 = interval
        ds = max(0.0, s1 - s0)
        if ds <= 0.0:
            continue
        s_mid = 0.5 * (s0 + s1)
        center_point = start + direction * s_mid
        beam_radius_here = max(float(beam_radius_m) + tan_div * s_mid, 0.0)
        frac = _cell_average_overlap_fraction(start_point_m=start, direction_unit=direction, s0=s0, s1=s1, plasma_radius_m=float(plasma_radius_profile[i]), beam_radius_m=float(beam_radius_m), tan_divergence=tan_div, radial_overlap_model=radial_overlap_model, quadrature_order=overlap_quadrature_order)
        physical_lengths[i] = ds
        overlap[i] = min(1.0, max(0.0, frac))
        effective_lengths[i] = physical_lengths[i] * overlap[i]
        beam_centers[i] = center_point
        beam_radii[i] = beam_radius_here
        s_mid_values[i] = s_mid
    # Attenuation must traverse intersected cells in beam path order rather than axial index order
    positive = np.flatnonzero(effective_lengths > 0.0)
    if positive.size:
        order = positive[np.argsort(s_mid_values[positive])]
    else:
        order = np.asarray([], dtype=int)

    return BeamPathOnAxisymmetricGrid(
        coordinate_centers_m=centers,
        coordinate_faces_m=faces,
        B_tilde_centers=B_tilde,
        plasma_radius_m=plasma_radius_profile,
        cell_volumes_m3=volumes,
        reference_area_m2=float(reference_area_m2),
        start_point_m=start,
        end_point_m=end,
        direction_unit=direction,
        total_length_m=length,
        cell_indices_along_beam=order.astype(int),
        physical_path_lengths_m=physical_lengths,
        radial_overlap_fractions=overlap,
        effective_path_lengths_m=effective_lengths,
        beam_center_points_m=beam_centers,
        beam_radius_profile_m=beam_radii,
        beam_radius_m=float(beam_radius_m),
        divergence_half_angle_rad=float(divergence_half_angle_rad),
        radial_overlap_model=str(radial_overlap_model),
    )
