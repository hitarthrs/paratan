"""Plasma facing surfaces and geometric exhaust path limits

Build candidate interfaces from the ParaTAN layout and apply role overrides
Find first material hits and cumulative radial access along supplied paths
"""
from __future__ import annotations
from collections.abc import Sequence
from dataclasses import dataclass, replace
from fnmatch import fnmatchcase
from typing import Any
import numpy as np
from source_model_revamp.geometry.device_domains import ParaTANLayout
from source_model_revamp.radial_profiles import canonical_radial_profile_model, radial_probability_cdf

PARTICLE_ROLES = frozenset({"absorbing", "collector", "end_ring", "reflecting", "excluded"})
ELECTRICAL_MODELS = frozenset({"floating", "grounded", "prescribed_bias"})
TERMINAL_PARTICLE_ROLES = frozenset({"absorbing", "collector", "end_ring"})

@dataclass(frozen=True)
class BoundaryOverride:
    """Particle or electrical settings selected by surface ID or component pattern

    The selector uses case sensitive shell wildcard matching
    A prescribed bias is in volts
    """
    selector: str
    particle_role: str | None = None
    electrical_model: str | None = None
    prescribed_bias_V: float | None = None

    def __post_init__(self) -> None:
        """Check the selector, supported roles, and required bias value"""
        if not self.selector.strip():
            raise ValueError("boundary override selector must not be empty")
        if self.particle_role is not None and self.particle_role not in PARTICLE_ROLES:
            raise ValueError(f"unsupported particle_role {self.particle_role!r}")
        if self.electrical_model is not None and self.electrical_model not in ELECTRICAL_MODELS:
            raise ValueError(f"unsupported electrical_model {self.electrical_model!r}")
        if self.electrical_model == "prescribed_bias":
            if self.prescribed_bias_V is None or not np.isfinite(self.prescribed_bias_V):
                raise ValueError("prescribed_bias electrical boundaries require finite prescribed_bias_V")
        elif self.prescribed_bias_V is not None:
            raise ValueError("prescribed_bias_V is valid only with electrical_model: prescribed_bias")

@dataclass(frozen=True)
class PlasmaFacingSurface:
    """Axisymmetric material interface bounding a connected vacuum volume

    Planes use an axial position and radial bounds
    Cylinders use a radius and axial bounds
    Cones use axial bounds and a radius at each end
    Geometry dimensions are in meters and unused fields remain None
    """
    surface_id: str
    component: str
    side: str
    geometry: str
    particle_role: str
    electrical_model: str
    prescribed_bias_V: float | None = None
    z_m: float | None = None
    radius_m: float | None = None
    radius_min_m: float | None = None
    radius_max_m: float | None = None
    z_min_m: float | None = None
    z_max_m: float | None = None
    radius_at_z_min_m: float | None = None
    radius_at_z_max_m: float | None = None
    connected_vacuum_region: str = ""
    automatically_detected: bool = True

    def __post_init__(self) -> None:
        """Validate particle and electrical roles and their bias settings"""
        if self.particle_role not in PARTICLE_ROLES:
            raise ValueError(f"unsupported particle role {self.particle_role!r}")
        if self.electrical_model not in ELECTRICAL_MODELS:
            raise ValueError(f"unsupported electrical model {self.electrical_model!r}")
        if self.electrical_model == "prescribed_bias":
            if self.prescribed_bias_V is None or not np.isfinite(self.prescribed_bias_V):
                raise ValueError("prescribed bias surface requires a finite prescribed_bias_V")
        elif self.prescribed_bias_V is not None:
            raise ValueError("prescribed_bias_V is valid only for a prescribed bias surface")

    @property
    def is_terminal(self) -> bool:
        """Return whether the particle role labels an absorbing terminal surface"""
        return self.particle_role in TERMINAL_PARTICLE_ROLES

    def as_dict(self) -> dict[str, Any]:
        """Return surface geometry and boundary settings as metadata"""
        return {
            "surface_id": self.surface_id,
            "component": self.component,
            "side": self.side,
            "geometry": self.geometry,
            "particle_role": self.particle_role,
            "electrical_model": self.electrical_model,
            "prescribed_bias_V": self.prescribed_bias_V,
            "z_m": self.z_m,
            "radius_m": self.radius_m,
            "radius_min_m": self.radius_min_m,
            "radius_max_m": self.radius_max_m,
            "z_min_m": self.z_min_m,
            "z_max_m": self.z_max_m,
            "radius_at_z_min_m": self.radius_at_z_min_m,
            "radius_at_z_max_m": self.radius_at_z_max_m,
            "connected_vacuum_region": self.connected_vacuum_region,
            "automatically_detected": self.automatically_detected,
        }

@dataclass(frozen=True)
class SurfaceIntersection:
    """One geometric material hit along an exhaust path

    path_parameter_m is distance from the path start
    z_m and radius_m locate the hit in machine coordinates
    """
    surface: PlasmaFacingSurface
    path_parameter_m: float
    z_m: float
    radius_m: float

@dataclass(frozen=True)
class SourceConnectedRadialSupport:
    """Radial support that remains connected to the throat along one exhaust path

    Edge arrays share the outward ordered z grid
    rho_limit_edges is the surviving radius divided by the nominal flux tube radius
    Survival probabilities use the selected radial profile
    Each segment records the material surface responsible for any new loss
    """
    side: str
    z_edges_outward_m: np.ndarray
    nominal_outer_radius_edges_m: np.ndarray
    rho_limit_edges: np.ndarray
    survival_probability_edges: np.ndarray
    segment_material_surface_ids: tuple[str | None, ...]
    radial_profile_model: str

    def __post_init__(self) -> None:
        """Validate outward order, cumulative radial limits, and survival probabilities"""
        branch = str(self.side).strip().lower()
        if branch not in {"left", "right"}:
            raise ValueError("radial support side must be left or right")
        z = np.asarray(self.z_edges_outward_m, dtype=float)
        radius = np.asarray(self.nominal_outer_radius_edges_m, dtype=float)
        rho = np.asarray(self.rho_limit_edges, dtype=float)
        survival = np.asarray(self.survival_probability_edges, dtype=float)
        if z.ndim != 1 or z.size < 2 or radius.shape != z.shape or rho.shape != z.shape or survival.shape != z.shape:
            raise ValueError("radial support edge arrays must have matching nonempty shapes")
        direction = -1.0 if branch == "left" else 1.0
        if np.any(~np.isfinite(z)) or np.any(direction * np.diff(z) <= 0.0):
            raise ValueError("radial support coordinates must be finite and ordered outward")
        if np.any(~np.isfinite(radius)) or np.any(radius <= 0.0):
            raise ValueError("radial support nominal radii must be positive and finite")
        if np.any(~np.isfinite(rho)) or np.any(rho < 0.0) or np.any(rho > 1.0):
            raise ValueError("radial support rho limits must lie inside zero to one")
        if np.any(np.diff(rho) > 256.0 * np.finfo(float).eps):
            raise ValueError("radial support rho limit must not increase outward")
        expected_survival = radial_probability_cdf(rho, self.radial_profile_model)
        if not np.allclose(survival, expected_survival, rtol=1.0e-13, atol=1.0e-15):
            raise ValueError("radial support survival probabilities do not match the radial profile")
        if len(self.segment_material_surface_ids) != z.size - 1:
            raise ValueError("radial support segment material surface IDs must match path segments")
        # Store normalized values during construction of this frozen dataclass
        object.__setattr__(self, "side", branch)
        object.__setattr__(self, "z_edges_outward_m", z)
        object.__setattr__(self, "nominal_outer_radius_edges_m", radius)
        object.__setattr__(self, "rho_limit_edges", rho)
        object.__setattr__(self, "survival_probability_edges", survival)
        object.__setattr__(self, "radial_profile_model", canonical_radial_profile_model(self.radial_profile_model))

    @property
    def effective_outer_radius_edges_m(self) -> np.ndarray:
        """Return nominal radii multiplied by the surviving radial limits"""
        return np.asarray(self.nominal_outer_radius_edges_m) * np.asarray(self.rho_limit_edges)

def _surface_radial_continuation_limit( surface: PlasmaFacingSurface, z_m: float, nominal_outer_radius_m: float, *, terminal_surface_id: str, tolerance_m: float) -> float | None:
    """Return the local aperture radius divided by the nominal flux tube radius

    Return None for an irrelevant, excluded, or selected terminal surface
    """
    # The selected terminal closes the path rather than removing its radial support early
    if surface.surface_id == terminal_surface_id or surface.particle_role == "excluded":
        return None
    radius = float(nominal_outer_radius_m)
    if radius <= 0.0:
        raise ValueError("nominal flux tube radius must be positive")
    if surface.geometry == "radial_cylinder" and surface.radius_m is not None:
        if surface.z_min_m is not None and z_m < float(surface.z_min_m) - tolerance_m:
            return None
        if surface.z_max_m is not None and z_m > float(surface.z_max_m) + tolerance_m:
            return None
        return float(surface.radius_m) / radius
    if (
        surface.geometry == "conical_frustum"
        and surface.z_min_m is not None
        and surface.z_max_m is not None
        and surface.radius_at_z_min_m is not None
        and surface.radius_at_z_max_m is not None
    ):
        if z_m < float(surface.z_min_m) - tolerance_m or z_m > float(surface.z_max_m) + tolerance_m:
            return None
        wall_radius = float(np.interp(z_m, [float(surface.z_min_m), float(surface.z_max_m)], [float(surface.radius_at_z_min_m), float(surface.radius_at_z_max_m)]))
        return wall_radius / radius
    if surface.geometry == "axial_plane" and surface.z_m is not None:
        if abs(z_m - float(surface.z_m)) > tolerance_m:
            return None
        # The inner edge of an annular face sets the open aperture
        aperture_radius = 0.0 if surface.radius_min_m is None else float(surface.radius_min_m)
        return aperture_radius / radius

    return None

def source_connected_radial_support( *, path_z_m: Sequence[float], nominal_outer_radius_m: Sequence[float], surfaces: Sequence[PlasmaFacingSurface], side: str, terminal_surface_id: str, radial_profile_model: str, tolerance_m: float = 1.0e-12) -> SourceConnectedRadialSupport:
    """Track the cumulative radial support from the throat toward one device end

    Evaluate material limits at supplied path nodes
    The surviving normalized radius cannot increase after a restriction
    """
    branch = str(side).strip().lower()
    if branch not in {"left", "right"}:
        raise ValueError("radial support side must be left or right")
    z = np.asarray(path_z_m, dtype=float)
    radius = np.asarray(nominal_outer_radius_m, dtype=float)
    if z.ndim != 1 or z.size < 2 or radius.shape != z.shape:
        raise ValueError("radial support path arrays must be matching vectors")
    direction = -1.0 if branch == "left" else 1.0
    if np.any(~np.isfinite(z)) or np.any(direction * np.diff(z) <= 0.0):
        raise ValueError("radial support path must be finite and ordered outward")
    if np.any(~np.isfinite(radius)) or np.any(radius <= 0.0):
        raise ValueError("radial support path radius must be positive and finite")
    model = canonical_radial_profile_model(radial_profile_model)
    scale = max(float(np.max(np.abs(z))), float(np.max(radius)), 1.0)
    tolerance = max(float(tolerance_m), 256.0 * np.finfo(float).eps * scale)
    rho = np.ones(z.size, dtype=float)
    edge_new_surface: list[str | None] = [None] * z.size
    current_limit = 1.0
    for index, (z_value, radius_value) in enumerate(zip(z, radius, strict=True)):
        local_limit = 1.0
        local_surface: str | None = None
        for surface in surfaces:
            if surface.side not in {branch, "both"}:
                continue
            candidate = _surface_radial_continuation_limit( surface, float(z_value), float(radius_value), terminal_surface_id=str(terminal_surface_id), tolerance_m=tolerance)
            if candidate is None:
                continue
            candidate_limit = float(np.clip(candidate, 0.0, 1.0))
            if candidate_limit < local_limit - 64.0 * np.finfo(float).eps:
                local_limit = candidate_limit
                local_surface = surface.surface_id
        # Radial support removed upstream cannot return when the vessel widens
        new_limit = min(current_limit, local_limit)
        if new_limit < current_limit - 64.0 * np.finfo(float).eps:
            edge_new_surface[index] = local_surface
        current_limit = new_limit
        rho[index] = current_limit
    if rho[0] < 1.0 - 1.0e-10:
        raise ValueError(f"{branch} radial support is already material limited at the magnetic throat")
    # Convert the radial cutoff to a surviving fraction using the prescribed profile
    survival = radial_probability_cdf(rho, model)
    segment_surface_ids: list[str | None] = []
    for index in range(z.size - 1):
        loss = float(survival[index] - survival[index + 1])
        surface_id = edge_new_surface[index + 1] if loss > 256.0 * np.finfo(float).eps else None
        if loss > 256.0 * np.finfo(float).eps and surface_id is None:
            raise ValueError("radial support loss could not be assigned to a material surface")
        segment_surface_ids.append(surface_id)

    return SourceConnectedRadialSupport(
        side=branch,
        z_edges_outward_m=z,
        nominal_outer_radius_edges_m=radius,
        rho_limit_edges=rho,
        survival_probability_edges=survival,
        segment_material_surface_ids=tuple(segment_surface_ids),
        radial_profile_model=model,
    )

@dataclass(frozen=True)
class PlasmaBoundarySelection:
    """Candidate surfaces, applied overrides, and the selected hit on each side

    Selection records geometry and boundary labels for downstream loss calculations
    """
    detection_model: str
    accessible_vacuum_regions: tuple[str, ...]
    candidate_surfaces: tuple[PlasmaFacingSurface, ...]
    explicit_overrides: tuple[BoundaryOverride, ...]
    left_selected_loss_surface: PlasmaFacingSurface
    right_selected_loss_surface: PlasmaFacingSurface
    selection_method: str
    valid: bool = True
    failure_reason: str = ""

    def as_metadata(self) -> dict[str, Any]:
        """Return candidate surfaces, overrides, and selected boundary settings"""
        return {
            "plasma_boundary_detection_model": self.detection_model,
            "plasma_accessible_vacuum_regions": list(self.accessible_vacuum_regions),
            "candidate_plasma_facing_surface_count": len(self.candidate_surfaces),
            "candidate_plasma_facing_surfaces": [item.as_dict() for item in self.candidate_surfaces],
            "explicit_boundary_overrides": [override.__dict__ for override in self.explicit_overrides],
            "left_selected_loss_surface": self.left_selected_loss_surface.as_dict(),
            "right_selected_loss_surface": self.right_selected_loss_surface.as_dict(),
            "loss_surface_selection_method": self.selection_method,
            "loss_surface_selection_valid": self.valid,
            "loss_surface_selection_failure_reason": self.failure_reason,
            "left_wall_electrical_model": self.left_selected_loss_surface.electrical_model,
            "right_wall_electrical_model": self.right_selected_loss_surface.electrical_model,
        }

def automatic_candidate_surfaces(layout: ParaTANLayout, *, default_particle_role: str = "absorbing", default_electrical_model: str = "floating") -> tuple[PlasmaFacingSurface, ...]:
    """Construct candidate material interfaces from the supported ParaTAN layout

    Use analytic planes, cylinders, and conical walls from the component dimensions
    """
    if default_particle_role not in PARTICLE_ROLES:
        raise ValueError(f"unsupported default particle role {default_particle_role!r}")
    if default_electrical_model not in ELECTRICAL_MODELS:
        raise ValueError(f"unsupported default electrical model {default_electrical_model!r}")
    common = dict(particle_role=default_particle_role, electrical_model=default_electrical_model)
    surfaces = [PlasmaFacingSurface("vacuum_vessel.central.radial_wall", "vacuum_vessel", "both", "radial_cylinder", radius_m=layout.central_cell_vacuum_radius_m, z_min_m=layout.central_cell_plasma_domain.z_min_m, z_max_m=layout.central_cell_plasma_domain.z_max_m, connected_vacuum_region="central_cell_vacuum", **common),]
    conical_domains = (layout.left_conical_vacuum_domain, layout.right_conical_vacuum_domain)
    if all(domain is not None for domain in conical_domains):
        left_conical, right_conical = conical_domains
        surfaces.extend([
            PlasmaFacingSurface("vacuum_vessel.left.conical_wall", "vacuum_vessel", "left", "conical_frustum", z_min_m=left_conical.z_min_m, z_max_m=left_conical.z_max_m, radius_at_z_min_m=layout.bottleneck_vacuum_radius_m, radius_at_z_max_m=layout.central_cell_vacuum_radius_m, connected_vacuum_region="left_conical_vacuum", **common),
            PlasmaFacingSurface("vacuum_vessel.right.conical_wall", "vacuum_vessel", "right", "conical_frustum", z_min_m=right_conical.z_min_m, z_max_m=right_conical.z_max_m, radius_at_z_min_m=layout.central_cell_vacuum_radius_m, radius_at_z_max_m=layout.bottleneck_vacuum_radius_m, connected_vacuum_region="right_conical_vacuum", **common),
        ])
    elif any(domain is not None for domain in conical_domains):
        raise ValueError("ParaTAN layout must provide both conical vacuum domains or neither")
    else:
        surfaces.extend([
            PlasmaFacingSurface("vacuum_vessel.left.perpendicular_transition_face", "vacuum_vessel", "left", "axial_plane", z_m=layout.central_cell_plasma_domain.z_min_m, radius_min_m=layout.bottleneck_vacuum_radius_m, radius_max_m=layout.central_cell_vacuum_radius_m, connected_vacuum_region="central_cell_vacuum", **common),
            PlasmaFacingSurface("vacuum_vessel.right.perpendicular_transition_face", "vacuum_vessel", "right", "axial_plane", z_m=layout.central_cell_plasma_domain.z_max_m, radius_min_m=layout.bottleneck_vacuum_radius_m, radius_max_m=layout.central_cell_vacuum_radius_m, connected_vacuum_region="central_cell_vacuum", **common),
        ])
    surfaces.extend([
        PlasmaFacingSurface("vacuum_vessel.left.bottleneck_radial_wall", "vacuum_vessel", "left", "radial_cylinder", radius_m=layout.bottleneck_vacuum_radius_m, z_min_m=layout.left_bottleneck_vacuum_domain.z_min_m, z_max_m=layout.left_bottleneck_vacuum_domain.z_max_m, connected_vacuum_region="left_bottleneck_vacuum", **common),
        PlasmaFacingSurface("vacuum_vessel.right.bottleneck_radial_wall", "vacuum_vessel", "right", "radial_cylinder", radius_m=layout.bottleneck_vacuum_radius_m, z_min_m=layout.right_bottleneck_vacuum_domain.z_min_m, z_max_m=layout.right_bottleneck_vacuum_domain.z_max_m, connected_vacuum_region="right_bottleneck_vacuum", **common),
        PlasmaFacingSurface("end_cell.left.upstream_annular_face", "end_cell", "left", "axial_plane", z_m=layout.left_end_cell_vacuum_domain.z_max_m, radius_min_m=layout.bottleneck_vacuum_radius_m, radius_max_m=layout.end_cell_vacuum_radius_m, connected_vacuum_region="left_end_cell_vacuum", **common),
        PlasmaFacingSurface("end_cell.right.upstream_annular_face", "end_cell", "right", "axial_plane", z_m=layout.right_end_cell_vacuum_domain.z_min_m, radius_min_m=layout.bottleneck_vacuum_radius_m, radius_max_m=layout.end_cell_vacuum_radius_m, connected_vacuum_region="right_end_cell_vacuum", **common),
        PlasmaFacingSurface("end_cell.left.axial_cap", "end_cell", "left", "axial_plane", z_m=layout.left_end_cell_vacuum_domain.z_min_m, radius_min_m=0.0, radius_max_m=layout.end_cell_vacuum_radius_m, connected_vacuum_region="left_end_cell_vacuum", **common),
        PlasmaFacingSurface("end_cell.right.axial_cap", "end_cell", "right", "axial_plane", z_m=layout.right_end_cell_vacuum_domain.z_max_m, radius_min_m=0.0, radius_max_m=layout.end_cell_vacuum_radius_m, connected_vacuum_region="right_end_cell_vacuum", **common),
        PlasmaFacingSurface("end_cell.left.radial_wall", "end_cell", "left", "radial_cylinder", radius_m=layout.end_cell_vacuum_radius_m, z_min_m=layout.left_end_cell_vacuum_domain.z_min_m, z_max_m=layout.left_end_cell_vacuum_domain.z_max_m, connected_vacuum_region="left_end_cell_vacuum", **common),
        PlasmaFacingSurface("end_cell.right.radial_wall", "end_cell", "right", "radial_cylinder", radius_m=layout.end_cell_vacuum_radius_m, z_min_m=layout.right_end_cell_vacuum_domain.z_min_m, z_max_m=layout.right_end_cell_vacuum_domain.z_max_m, connected_vacuum_region="right_end_cell_vacuum", **common),
    ])

    return tuple(surfaces)

def apply_boundary_overrides(surfaces: Sequence[PlasmaFacingSurface], overrides: Sequence[BoundaryOverride]) -> tuple[PlasmaFacingSurface, ...]:
    """Apply role overrides matched to surface IDs or component names

    Reject unmatched selectors and any surface claimed by multiple overrides
    """
    claimed: dict[str, str] = {}
    result = list(surfaces)
    for override in overrides:
        matched = [i for i, surface in enumerate(result) if fnmatchcase(surface.surface_id, override.selector) or fnmatchcase(surface.component, override.selector)]
        if not matched:
            raise ValueError(f"boundary override selector {override.selector!r} matched no candidate surface")
        for i in matched:
            surface = result[i]
            if surface.surface_id in claimed:
                raise ValueError(f"contradictory boundary overrides {claimed[surface.surface_id]!r} and" f"{override.selector!r} both match {surface.surface_id!r}")
            claimed[surface.surface_id] = override.selector
            electrical = override.electrical_model or surface.electrical_model
            bias = override.prescribed_bias_V if electrical == "prescribed_bias" else None
            result[i] = replace(surface, particle_role=override.particle_role or surface.particle_role, electrical_model=electrical, prescribed_bias_V=bias, automatically_detected=False)

    return tuple(result)

def first_trajectory_surface(*, start_z_m: float, start_radius_m: float, direction_z: float, direction_radius: float, surfaces: Sequence[PlasmaFacingSurface], side: str, tolerance_m: float = 1.0e-12) -> SurfaceIntersection:
    """Find the first nonexcluded material interface along a straight radial axial ray

    Normalize the direction so the returned path parameter is distance in meters
    Particle roles label the hit but do not change the geometric search
    """
    norm = float(np.hypot(direction_z, direction_radius))
    if not np.isfinite(norm) or norm <= 0.0:
        raise ValueError("trajectory direction must be finite and nonzero")
    dz = float(direction_z) / norm
    dr = float(direction_radius) / norm
    intersections: list[SurfaceIntersection] = []
    for surface in surfaces:
        if surface.side not in {side, "both"} or surface.particle_role == "excluded":
            continue
        distance: float | None = None
        if surface.geometry == "axial_plane" and surface.z_m is not None and abs(dz) > tolerance_m:
            distance = (float(surface.z_m) - start_z_m) / dz
            r_hit = start_radius_m + distance * dr
            if r_hit < -tolerance_m:
                distance = None
            if distance is not None and surface.radius_min_m is not None and r_hit < surface.radius_min_m - tolerance_m:
                distance = None
            if distance is not None and surface.radius_max_m is not None and r_hit > surface.radius_max_m + tolerance_m:
                distance = None
        elif surface.geometry == "radial_cylinder" and surface.radius_m is not None and abs(dr) > tolerance_m:
            distance = (float(surface.radius_m) - start_radius_m) / dr
            if distance is not None:
                z_hit = start_z_m + distance * dz
                if surface.z_min_m is not None and z_hit < surface.z_min_m - tolerance_m:
                    distance = None
                if surface.z_max_m is not None and z_hit > surface.z_max_m + tolerance_m:
                    distance = None
        elif (surface.geometry == "conical_frustum" and surface.z_min_m is not None and surface.z_max_m is not None and surface.radius_at_z_min_m is not None and surface.radius_at_z_max_m is not None):
            wall_slope = (surface.radius_at_z_max_m - surface.radius_at_z_min_m) / (surface.z_max_m - surface.z_min_m)
            # Equal ray and wall slopes have no isolated intersection
            denominator = dr - wall_slope * dz
            if abs(denominator) > tolerance_m:
                wall_at_start = surface.radius_at_z_min_m + wall_slope * (start_z_m - surface.z_min_m)
                distance = (wall_at_start - start_radius_m) / denominator
                z_hit = start_z_m + distance * dz
                if z_hit < surface.z_min_m - tolerance_m or z_hit > surface.z_max_m + tolerance_m:
                    distance = None
        if distance is None or distance <= tolerance_m:
            continue
        intersections.append(
        SurfaceIntersection(surface=surface, path_parameter_m=float(distance), z_m=float(start_z_m + distance * dz), radius_m=float(start_radius_m + distance * dr)))
    if not intersections:
        raise ValueError(f"{side} exhaust trajectory exits the modeled accessible geometry without reaching a terminal material surface")

    return min(intersections, key=lambda item: item.path_parameter_m)

def first_flux_tube_surface(*, path_z_m: Sequence[float], path_radius_m: Sequence[float], surfaces: Sequence[PlasmaFacingSurface], side: str, tolerance_m: float = 1.0e-12) -> SurfaceIntersection:
    """Find the first nonexcluded material hit along a supplied flux tube path

    Treat each pair of path nodes as a straight segment
    Return the hit with its accumulated path length
    """
    z = np.asarray(path_z_m, dtype=float)
    radius = np.asarray(path_radius_m, dtype=float)
    if z.ndim != 1 or radius.shape != z.shape or z.size < 2:
        raise ValueError("flux tube path arrays must be matching vectors with at least two points")
    if np.any(~np.isfinite(z)) or np.any(~np.isfinite(radius)) or np.any(radius < 0.0):
        raise ValueError("flux tube path arrays must be finite with nonnegative radius")
    direction = -1.0 if str(side).strip().lower() == "left" else 1.0 if str(side).strip().lower() == "right" else 0.0
    if direction == 0.0 or np.any(direction * np.diff(z) <= 0.0):
        raise ValueError("flux tube path must be strictly ordered toward the selected side")
    segment_lengths = np.hypot(np.diff(z), np.diff(radius))
    cumulative = np.concatenate(([0.0], np.cumsum(segment_lengths)))
    intersections: list[SurfaceIntersection] = []
    for index in range(z.size - 1):
        z0 = float(z[index])
        z1 = float(z[index + 1])
        r0 = float(radius[index])
        r1 = float(radius[index + 1])
        dz = z1 - z0
        dr = r1 - r0
        segment_length = float(segment_lengths[index])
        fraction_tolerance = tolerance_m / max(segment_length, tolerance_m)
        for surface in surfaces:
            if surface.side not in {side, "both"} or surface.particle_role == "excluded":
                continue
            fraction: float | None = None
            if surface.geometry == "axial_plane" and surface.z_m is not None:
                fraction = (float(surface.z_m) - z0) / dz
                if fraction < -fraction_tolerance or fraction > 1.0 + fraction_tolerance:
                    fraction = None
                if fraction is not None:
                    r_hit = r0 + fraction * dr
                    if surface.radius_min_m is not None and r_hit < float(surface.radius_min_m) - tolerance_m:
                        fraction = None
                    if surface.radius_max_m is not None and r_hit > float(surface.radius_max_m) + tolerance_m:
                        fraction = None
            elif surface.geometry == "radial_cylinder" and surface.radius_m is not None and abs(dr) > tolerance_m:
                fraction = (float(surface.radius_m) - r0) / dr
                if fraction < -fraction_tolerance or fraction > 1.0 + fraction_tolerance:
                    fraction = None
                if fraction is not None:
                    z_hit = z0 + fraction * dz
                    if surface.z_min_m is not None and z_hit < float(surface.z_min_m) - tolerance_m:
                        fraction = None
                    if surface.z_max_m is not None and z_hit > float(surface.z_max_m) + tolerance_m:
                        fraction = None
            elif surface.geometry == "conical_frustum" and surface.z_min_m is not None and surface.z_max_m is not None and surface.radius_at_z_min_m is not None and surface.radius_at_z_max_m is not None:
                wall_slope = (float(surface.radius_at_z_max_m) - float(surface.radius_at_z_min_m)) / (float(surface.z_max_m) - float(surface.z_min_m))
                wall_at_z0 = float(surface.radius_at_z_min_m) + wall_slope * (z0 - float(surface.z_min_m))
                denominator = dr - wall_slope * dz
                if abs(denominator) > tolerance_m:
                    fraction = (wall_at_z0 - r0) / denominator
                    if fraction < -fraction_tolerance or fraction > 1.0 + fraction_tolerance:
                        fraction = None
                    if fraction is not None:
                        z_hit = z0 + fraction * dz
                        if z_hit < float(surface.z_min_m) - tolerance_m or z_hit > float(surface.z_max_m) + tolerance_m:
                            fraction = None
            if fraction is None:
                continue
            clipped = min(max(float(fraction), 0.0), 1.0)
            # Ignore contact at the initial point but retain later segment boundary hits
            if index == 0 and clipped * segment_length <= tolerance_m:
                continue
            z_hit = z0 + clipped * dz
            r_hit = r0 + clipped * dr
            path_parameter = float(cumulative[index] + clipped * segment_lengths[index])
            intersections.append(SurfaceIntersection(surface=surface, path_parameter_m=path_parameter, z_m=z_hit, radius_m=r_hit))
    if not intersections:
        raise ValueError(f"{side} flux tube exits the modeled accessible geometry without reaching a material surface")

    return min(intersections, key=lambda item: item.path_parameter_m)

def build_boundary_selection(layout: ParaTANLayout, *, left_throat_z_m: float, right_throat_z_m: float, detection_mode: str, default_particle_role: str, default_electrical_model: str, overrides: Sequence[BoundaryOverride] = ()) -> PlasmaBoundarySelection:
    """Build candidate surfaces, apply overrides, and select axial exhaust hits

    Start an outward ray on the axis at each magnetic throat
    Selection records the first nonexcluded hit regardless of its particle role
    """
    if detection_mode != "automatic_with_overrides":
        raise ValueError("only plasma boundary detection_mode: automatic_with_overrides is supported")
    candidates = apply_boundary_overrides(automatic_candidate_surfaces(layout, default_particle_role=default_particle_role, default_electrical_model=default_electrical_model), overrides,)
    # Axis rays select the end interfaces separately from the outer flux tube wall hits
    left = first_trajectory_surface(start_z_m=left_throat_z_m, start_radius_m=0.0, direction_z=-1.0, direction_radius=0.0, surfaces=candidates, side="left").surface
    right = first_trajectory_surface(start_z_m=right_throat_z_m, start_radius_m=0.0, direction_z=1.0, direction_radius=0.0, surfaces=candidates, side="right").surface
    accessible_regions = ["central_cell_vacuum"]
    if layout.left_conical_vacuum_domain is not None and layout.right_conical_vacuum_domain is not None:
        accessible_regions.extend(["left_conical_vacuum", "right_conical_vacuum"])
    accessible_regions.extend(["left_bottleneck_vacuum", "right_bottleneck_vacuum", "left_end_cell_vacuum", "right_end_cell_vacuum"])

    return PlasmaBoundarySelection(
        detection_model=detection_mode,
        accessible_vacuum_regions=tuple(accessible_regions),
        candidate_surfaces=candidates,
        explicit_overrides=tuple(overrides),
        left_selected_loss_surface=left,
        right_selected_loss_surface=right,
        selection_method="first_terminal_interface_on_geometry-linked_axial_exhaust_ray",
    )
