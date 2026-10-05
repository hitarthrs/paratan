"""
Magnetic geometry, plasma boundary, and axial grid construction for the pipeline

The stage builds either the analytic Egedal mirror or the ParaTAN linked coil field and returns confined and full device flux tube geometry
"""
from __future__ import annotations
from dataclasses import replace
from typing import Any
import numpy as np
from numpy.typing import ArrayLike
from scipy.optimize import brentq
from source_model_revamp.geometry.axial_grid_geometry import cell_centers_from_faces, flux_tube_cell_volumes_m3
from source_model_revamp.geometry.background_profiles import BackgroundProfileSpecification, MAGNETICALLY_MAPPED_STARTUP_PROFILE_MODEL, build_ion_seed_reference_profiles
from source_model_revamp.geometry.coil_fields import CircularCoilOnAxis, coil_set_Bz_on_axis_T
from source_model_revamp.geometry.device_domains import AxialDomain, attach_confined_domain
from source_model_revamp.geometry.egedal_mirror_profile import egedal_B_tilde
from source_model_revamp.geometry.magnetic_topology import discover_magnetic_throats
from source_model_revamp.geometry.plasma_boundaries import PlasmaFacingSurface, build_boundary_selection
from source_model_revamp.integration.pipeline_types import GeometryStageResult
from source_model_revamp.integration.plasma_closure import uses_startup_seed_only
from source_model_revamp.integration.config.root import SourceModelRunConfig

def _build_effective_coils(config: SourceModelRunConfig) -> tuple[CircularCoilOnAxis, ...]:
    """Convert configured effective coil records into on axis circular coil field objects"""
    return tuple(CircularCoilOnAxis( name=coil.name, z_center_m=coil.z_center_m, radius_m=coil.radius_m, amp_turns=coil.amp_turns, group=coil.group) for coil in config.geometry.effective_coils)

def _scalar_or_array(values: ArrayLike, original: ArrayLike) -> float | np.ndarray:
    """Return a scalar for scalar input and otherwise preserve the original array shape"""
    result = np.asarray(values, dtype=float)
  
    return float(result) if np.ndim(original) == 0 else result

def _event_aligned_edges(z_min_m: float, z_max_m: float, target_spacing_m: float, events: list[float]) -> np.ndarray:
    """
    Build an axial edge grid that includes every physical event coordinate exactly
    
    Intervals between events are subdivided so their spacing does not exceed the requested target spacing
    """
    if not np.isfinite(target_spacing_m) or target_spacing_m <= 0.0:
        raise ValueError("event aligned grid spacing must be finite and positive")
 
    bounded_events = [float(value) for value in events if z_min_m <= value <= z_max_m and np.isfinite(value)]
    anchors = np.unique(np.asarray([z_min_m, *bounded_events, z_max_m], dtype=float))
    pieces: list[np.ndarray] = []
   
    for lower, upper in zip(anchors[:-1], anchors[1:], strict=True):
        count = max(1, int(np.ceil((upper - lower) / target_spacing_m)))
        segment = np.linspace(lower, upper, count + 1)
        pieces.append(segment if not pieces else segment[1:])
    edges = np.concatenate(pieces)

    return edges

def _flux_tube_radius_m(z_m: float, *, plasma_radius_m: float, B0_T: float, B_T_function) -> float:
    """Return the flux tube radius `a(z) = a0 sqrt(B0 / B(z))` in m"""
    field = float(B_T_function(float(z_m)))
  
    if not np.isfinite(field) or field <= 0.0:
        raise ValueError("flux tube material intersection requires positive finite magnetic field")
   
    return float(plasma_radius_m) * np.sqrt(float(B0_T) / field)

def _surface_radius_m(surface: PlasmaFacingSurface, z_m: float) -> float:
    """Return the selected plasma facing surface radius at one axial coordinate in m"""
    if surface.geometry == "radial_cylinder" and surface.radius_m is not None:
        return float(surface.radius_m)
  
    if surface.geometry == "conical_frustum" and surface.z_min_m is not None and surface.z_max_m is not None and surface.radius_at_z_min_m is not None and surface.radius_at_z_max_m is not None:
        fraction = (float(z_m) - float(surface.z_min_m)) / (float(surface.z_max_m) - float(surface.z_min_m))
        return float(surface.radius_at_z_min_m) + fraction * (float(surface.radius_at_z_max_m) - float(surface.radius_at_z_min_m))
  
    raise ValueError("surface does not define a radial material boundary")

def _first_flux_tube_material_event(*, throat_z_m: float, side: str, z_min_m: float, z_max_m: float, plasma_radius_m: float, B0_T: float, B_T_function, surfaces: tuple[PlasmaFacingSurface, ...]) -> tuple[float, float, PlasmaFacingSurface]:
    """
    Find the first axial intersection between the expanding flux tube and a selected material surface on one side of the device
    
    The search considers axial planes and cylindrical or frustum surfaces and returns the event coordinate, local flux tube radius, and selected surface
    """
    branch = str(side).strip().lower()
    direction = -1.0 if branch == "left" else 1.0 if branch == "right" else 0.0
   
    if direction == 0.0:
        raise ValueError("flux tube material event side must be left or right")
    domain_limit = float(z_min_m if branch == "left" else z_max_m)
    z_tolerance = 256.0 * np.finfo(float).eps * max(abs(float(throat_z_m)), abs(domain_limit), 1.0)
    candidates: list[tuple[float, float, float, PlasmaFacingSurface]] = []
  
    for surface in surfaces:
        if surface.side not in {branch, "both"} or surface.particle_role == "excluded":
            continue
       
        if surface.geometry == "axial_plane" and surface.z_m is not None:
            z_hit = float(surface.z_m)
            distance = direction * (z_hit - float(throat_z_m))
         
            if distance <= z_tolerance or direction * (domain_limit - z_hit) < -z_tolerance:
                continue
            radius = _flux_tube_radius_m(z_hit, plasma_radius_m=plasma_radius_m, B0_T=B0_T, B_T_function=B_T_function)
            if surface.radius_min_m is not None and radius < float(surface.radius_min_m) - z_tolerance:
                continue
            if surface.radius_max_m is not None and radius > float(surface.radius_max_m) + z_tolerance:
                continue
            candidates.append((distance, z_hit, radius, surface))
            continue
      
        if surface.geometry not in {"radial_cylinder", "conical_frustum"} or surface.z_min_m is None or surface.z_max_m is None:
            continue
      
        interval_low = max(float(z_min_m), float(surface.z_min_m), min(float(throat_z_m), domain_limit))
        interval_high = min(float(z_max_m), float(surface.z_max_m), max(float(throat_z_m), domain_limit))
      
        if interval_high <= interval_low:
            continue
        outward_start = interval_high if branch == "left" else interval_low
        outward_end = interval_low if branch == "left" else interval_high
      
        if direction * (outward_start - float(throat_z_m)) < -z_tolerance:
            outward_start = float(throat_z_m)
       
        scan = np.linspace(outward_start, outward_end, 513)
        residual = np.asarray([_flux_tube_radius_m(value, plasma_radius_m=plasma_radius_m, B0_T=B0_T, B_T_function=B_T_function) - _surface_radius_m(surface, value) for value in scan], dtype=float)
        root: float | None = None
        radius_scale = max(float(plasma_radius_m), float(np.max(np.abs(residual))), 1.0)
        radius_tolerance = 256.0 * np.finfo(float).eps * radius_scale
       
        for index in range(scan.size - 1):
            left_value = float(residual[index])
            right_value = float(residual[index + 1])
            if abs(left_value) <= radius_tolerance:
                root = float(scan[index])
                break
            if left_value * right_value < 0.0:
                a = float(min(scan[index], scan[index + 1]))
                b = float(max(scan[index], scan[index + 1]))
                root = float(brentq(lambda value: _flux_tube_radius_m(value, plasma_radius_m=plasma_radius_m, B0_T=B0_T, B_T_function=B_T_function) - _surface_radius_m(surface, value), a, b, xtol=8.0 * np.finfo(float).eps * max(abs(a), abs(b), 1.0), rtol=8.0 * np.finfo(float).eps))
                break
    
        if root is None and abs(float(residual[-1])) <= radius_tolerance:
            root = float(scan[-1])
        if root is None:
            continue
        distance = direction * (root - float(throat_z_m))
        if distance <= z_tolerance:
            continue
        candidates.append((distance, root, _flux_tube_radius_m(root, plasma_radius_m=plasma_radius_m, B0_T=B0_T, B_T_function=B_T_function), surface))
   
    if not candidates:
        raise ValueError(f"{branch} magnetic flux tube does not reach a modeled material surface")
    _, z_hit, radius_hit, surface = min(candidates, key=lambda item: item[0])
  
    return float(z_hit), float(radius_hit), surface

def build_geometry_stage(config: SourceModelRunConfig) -> GeometryStageResult:
    """
    Build the magnetic field, confined domain, full device domain, flux tube volumes, and startup or reference background geometry
    
    The analytic path uses the supplied Egedal field shape
    The ParaTAN path discovers magnetic throats from the fitted coil field, resolves plasma facing boundaries, finds the first material intersections beyond each throat, and builds a full device grid aligned to those physical events
    """
    g = config.geometry
    linked = False
    extra: dict[str, Any] = {}
    domains = None
    boundaries = None
    background = None
    full_edges: np.ndarray
    full_centers: np.ndarray
    full_B_T: np.ndarray
    full_B_tilde: np.ndarray
    full_cell_volumes: np.ndarray

    if g.model == "egedal_analytic":
        mirror_ratio = g.field_strength.mirror_ratio
        B0_T = g.field_strength.midplane_B_T
        half_length = float(g.half_length_m)
        orbit_midpoint = 0.0
        z_edges = np.linspace(-half_length, half_length, g.axial_bins + 1)
        z_centers = cell_centers_from_faces(z_edges)
        zeta_edges = z_edges / half_length
        zeta_centers = z_centers / half_length
        B_tilde = np.asarray( egedal_B_tilde(zeta_centers, mirror_ratio, g.transition_width, g.shape_exponent), dtype=float)

        def B_tilde_function(zeta: ArrayLike) -> float | np.ndarray:
            """Evaluate analytic normalized magnetic field on scalar or array ζ coordinates"""
            return _scalar_or_array(egedal_B_tilde(zeta, mirror_ratio, g.transition_width, g.shape_exponent), zeta)
       
        def B_T_function(z_m: ArrayLike) -> float | np.ndarray:
            """Evaluate analytic magnetic field magnitude in T on scalar or array axial coordinates"""
            return _scalar_or_array(B0_T * egedal_B_tilde(np.asarray(z_m, dtype=float) / half_length, mirror_ratio, g.transition_width, g.shape_exponent), z_m)

        full_edges = z_edges.copy()
        full_centers = z_centers.copy()
        full_B_tilde = B_tilde.copy()
        full_B_T = B0_T * full_B_tilde
        magnetic_model = "egedal_analytic_manual_benchmark_domain"
        extra.update({
                "requested_midplane_B_T": g.field_strength.midplane_B_T,
                "requested_mirror_throat_B_T": g.field_strength.mirror_throat_B_T,
                "requested_mirror_ratio": g.field_strength.mirror_ratio,
                "fitted_midplane_B_T": g.field_strength.midplane_B_T,
                "fitted_mirror_throat_B_T": g.field_strength.mirror_throat_B_T,
                "fitted_mirror_ratio": g.field_strength.mirror_ratio,
                "left_fitted_mirror_throat_z_m": -half_length,
                "right_fitted_mirror_throat_z_m": half_length,
                "coil_current_fit_model": "not_applicable_egedal_analytic",
                "coil_effective_amp_turns_are_internal_scaling": False,
                "transition_width": g.transition_width,
                "shape_exponent": g.shape_exponent,
                "domain_mode": "manual_benchmark",
                "explicit_device_domains_available": False,
            })
   
    elif g.model == "paratan_coil_field":
        if g.paratan_layout is None:
            raise ValueError(" linked ParaTAN field requires a derived ParaTAN layout")
        layout = g.paratan_layout
        coils = _build_effective_coils(config)

        def B_T_function(z_m: ArrayLike) -> float | np.ndarray:
            """Evaluate fitted ParaTAN linked on axis magnetic field magnitude in T"""
            value = np.abs(coil_set_Bz_on_axis_T(z_m, coils))
            return _scalar_or_array(value, z_m)

        B0_T = B_T_function(layout.midplane_z_m)
       
        if B0_T <= 0.0:
            raise ValueError("paratan_coil_field geometry produced nonpositive midplane B0")
       
        throats = discover_magnetic_throats( B_T_function, midplane_z_m=layout.midplane_z_m, search_z_min_m=layout.full_device_domain.z_min_m, search_z_max_m=layout.full_device_domain.z_max_m, samples=max(g.coil_field_grid_points, 401), requested_mirror_ratio=g.field_strength.mirror_ratio,)
        boundaries = build_boundary_selection( layout, left_throat_z_m=throats.left_z_m, right_throat_z_m=throats.right_z_m, detection_mode=config.plasma_boundaries.detection_mode, default_particle_role=config.plasma_boundaries.default_material_interface_role, default_electrical_model=config.plasma_boundaries.default_electrical_model, overrides=config.plasma_boundaries.overrides,)
        domains = attach_confined_domain( layout, throats.left_z_m, throats.right_z_m, left_boundary_z_m=boundaries.left_selected_loss_surface.z_m, right_boundary_z_m=boundaries.right_selected_loss_surface.z_m,)
        z_edges = np.linspace(throats.left_z_m, throats.right_z_m, g.axial_bins + 1)
        z_centers = cell_centers_from_faces(z_edges)
        half_length = 0.5 * (throats.right_z_m - throats.left_z_m)
        orbit_midpoint = 0.5 * (throats.right_z_m + throats.left_z_m)
        zeta_edges = (z_edges - orbit_midpoint) / half_length
        zeta_centers = (z_centers - orbit_midpoint) / half_length
        B_tilde = np.asarray(B_T_function(z_centers), dtype=float) / B0_T

        def B_tilde_function(zeta: ArrayLike) -> float | np.ndarray:
            """Evaluate fitted normalized magnetic field on scalar or array ζ coordinates"""
            z_m = orbit_midpoint + np.asarray(zeta, dtype=float) * half_length
            return _scalar_or_array(np.asarray(B_T_function(z_m), dtype=float) / B0_T, zeta)

        mirror_ratio = float(0.5 * (throats.left_mirror_ratio + throats.right_mirror_ratio))
        magnetic_model = "paratan_lf_hf_coil_field_direct_on_axis"
        linked = True
        left_flux_hit_z, left_flux_hit_radius, left_flux_hit_surface = _first_flux_tube_material_event(throat_z_m=throats.left_z_m, side="left", z_min_m=domains.full_device_domain.z_min_m, z_max_m=domains.full_device_domain.z_max_m, plasma_radius_m=g.plasma_radius_m, B0_T=B0_T, B_T_function=B_T_function, surfaces=boundaries.candidate_surfaces)
        right_flux_hit_z, right_flux_hit_radius, right_flux_hit_surface = _first_flux_tube_material_event(throat_z_m=throats.right_z_m, side="right", z_min_m=domains.full_device_domain.z_min_m, z_max_m=domains.full_device_domain.z_max_m, plasma_radius_m=g.plasma_radius_m, B0_T=B0_T, B_T_function=B_T_function, surfaces=boundaries.candidate_surfaces)
        material_events: list[float] = []
      
        for surface in boundaries.candidate_surfaces:
            for value in (surface.z_m, surface.z_min_m, surface.z_max_m):
                if value is not None:
                    material_events.append(float(value))
      
        events = [
            layout.midplane_z_m,
            layout.central_cell_plasma_domain.z_min_m,
            layout.central_cell_plasma_domain.z_max_m,
            throats.left_z_m,
            throats.right_z_m,
            float(boundaries.left_selected_loss_surface.z_m),
            float(boundaries.right_selected_loss_surface.z_m),
            left_flux_hit_z,
            right_flux_hit_z,
            *material_events,
            *throats.left_expander_extrema_m,
            *throats.right_expander_extrema_m,
        ]
        full_edges = _event_aligned_edges( domains.full_device_domain.z_min_m, domains.full_device_domain.z_max_m, float(np.min(np.diff(z_edges))), events)
        full_centers = cell_centers_from_faces(full_edges)
        full_B_T = np.asarray(B_T_function(full_centers), dtype=float)
        full_B_tilde = full_B_T / B0_T
        extra.update(dict(g.coil_fit_metadata))
        extra.update({
                "coil_count": len(coils),
                "coil_z_center_m": [coil.z_center_m for coil in coils],
                "coil_radius_m": [coil.radius_m for coil in coils],
                "coil_group": [coil.group for coil in coils],
                "coil_effective_amp_turns": [coil.amp_turns for coil in coils],
                "field_evaluation_model": "direct_fitted_coil_sum_in_physical_z",
                "field_fit_search_domain_is_separate_from_evaluation_domain": True,
                "left_fitted_mirror_throat_z_m": throats.left_z_m,
                "right_fitted_mirror_throat_z_m": throats.right_z_m,
                "left_fitted_mirror_throat_B_T": throats.left_B_T,
                "right_fitted_mirror_throat_B_T": throats.right_B_T,
                "left_fitted_mirror_ratio": throats.left_mirror_ratio,
                "right_fitted_mirror_ratio": throats.right_mirror_ratio,
                "fitted_mirror_throat_z_m": throats.right_z_m,
                "fitted_mirror_throat_inside_kinetic_domain": True,
                "field_monotonic_midplane_to_left_throat": throats.monotonic_left,
                "field_monotonic_midplane_to_right_throat": throats.monotonic_right,
                "left_expander_additional_field_extrema_m": throats.left_expander_extrema_m,
                "right_expander_additional_field_extrema_m": throats.right_expander_extrema_m,
                "expander_effective_potential_reflection_check_required": True,
                "left_flux_tube_first_material_hit_z_m": left_flux_hit_z,
                "right_flux_tube_first_material_hit_z_m": right_flux_hit_z,
                "left_flux_tube_first_material_hit_radius_m": left_flux_hit_radius,
                "right_flux_tube_first_material_hit_radius_m": right_flux_hit_radius,
                "left_flux_tube_first_material_surface_id": left_flux_hit_surface.surface_id,
                "right_flux_tube_first_material_surface_id": right_flux_hit_surface.surface_id,
                "full_device_population_grid_includes_flux_tube_material_intersections": True,
                "full_device_population_grid_includes_all_candidate_material_axial_events": True,
                "outer_flux_tube_first_material_hit_role": "geometry_diagnostic_not_whole_population_terminal",
                "kinetic_domain_mirror_ratio": mirror_ratio,
                "kinetic_domain_to_requested_mirror_ratio": mirror_ratio / g.field_strength.mirror_ratio,
                "kinetic_domain_mirror_ratio_relative_error": (mirror_ratio - g.field_strength.mirror_ratio) / g.field_strength.mirror_ratio,
                "domain_mode": "geometry_linked",
                "explicit_device_domains_available": True,
                "device_domains": domains.as_dict(),
                "paratan_layout": layout.as_dict(),
                **boundaries.as_metadata(),
            })
    else:
        raise ValueError(f"Unsupported geometry model {g.model!r}")
    if np.any(~np.isfinite(B_tilde)) or np.any(B_tilde <= 0.0):
        raise ValueError("geometry produced invalid B_tilde values")
    if np.any(~np.isfinite(full_B_tilde)) or np.any(full_B_tilde <= 0.0):
        raise ValueError("full device geometry produced invalid magnetic field values")

    midplane_area = g.midplane_area_m2
    cell_volumes = flux_tube_cell_volumes_m3(z_edges, B_tilde, midplane_area)
    full_cell_volumes = flux_tube_cell_volumes_m3(full_edges, full_B_tilde, midplane_area)
    radial_outer = g.plasma_radius_m / np.sqrt(B_tilde)
    radial_inner = np.zeros_like(radial_outer)
    full_flux_radius = g.plasma_radius_m / np.sqrt(full_B_tilde)
    full_B_T_edges = np.asarray(B_T_function(full_edges), dtype=float)
  
    if full_B_T_edges.shape != full_edges.shape or np.any(~np.isfinite(full_B_T_edges)) or np.any(full_B_T_edges <= 0.0):
        raise ValueError("full device edge magnetic field must be positive and finite")
  
    full_flux_radius_edges = g.plasma_radius_m * np.sqrt(float(B0_T) / full_B_T_edges)
    p = config.plasma_closure
   
    if uses_startup_seed_only(p):
        if p.background_deuterium_midplane_density_m3 + p.background_tritium_midplane_density_m3 <= 0.0:
            raise ValueError("NBI supported startup requires positive D plus T midplane density")
      
        support = domains.confined_orbit_domain if domains is not None else AxialDomain("confined_orbit_domain", float(z_edges[0]), float(z_edges[-1]), "magnetic mirror throat boundaries")
        midpoint = domains.midplane_z_m if domains is not None else orbit_midpoint
        sample_z = np.unique(np.concatenate((full_edges, full_centers, z_edges, z_centers, np.asarray([support.z_min_m, midpoint, support.z_max_m]))))
        sample_B = np.asarray(B_T_function(sample_z), dtype=float)
        reference_B = np.asarray(B_T_function(np.asarray([midpoint, support.z_min_m, support.z_max_m])), dtype=float)
        confining_throat_B_T = float(min(reference_B[1], reference_B[2]))
        denominator = 1.0 - float(reference_B[0]) / confining_throat_B_T
       
        if denominator <= 0.0:
            raise ValueError("magnetically mapped startup seed requires throat fields above the midplane field")
      
        seed_shape = np.where(support.contains(sample_z), np.sqrt(np.clip((1.0 - sample_B / confining_throat_B_T) / denominator, 0.0, None)), 0.0)
        seed_specification = BackgroundProfileSpecification(model="tabulated_axial_profile", tabulated_z_m=tuple(sample_z), tabulated_scale=tuple(seed_shape))
        background = build_ion_seed_reference_profiles(full_centers, full_cell_volumes, support, seed_specification, magnetic_midplane_z_m=midpoint, deuterium_midplane_density_m3=p.background_deuterium_midplane_density_m3, tritium_midplane_density_m3=p.background_tritium_midplane_density_m3, ion_temperature_keV=p.background_ion_temperature_keV, symmetry_tolerance=config.kinetic_electrostatic.modal_eq42_density_symmetry_tolerance)
        background = replace(background, model=MAGNETICALLY_MAPPED_STARTUP_PROFILE_MODEL)
        extra.update(background.as_metadata())
        extra.update({"startup_seed_profile_role": "initialization_only", "startup_seed_velocity_distribution_model": "midplane_maxwellian_with_empty_magnetic_loss_cone", "startup_seed_axial_mapping_model": "collisionless_energy_and_magnetic_moment_invariants_with_zero_startup_potential", "startup_seed_profile_user_selectable": False, "startup_seed_midplane_B_T": float(reference_B[0]), "startup_seed_left_throat_B_T": float(reference_B[1]), "startup_seed_right_throat_B_T": float(reference_B[2]), "startup_seed_confining_throat_B_T": confining_throat_B_T, "startup_seed_density_shape_definition": "sqrt((1-B(z)/B_confining_throat)/(1-B_midplane/B_confining_throat))", "startup_seed_electron_density_model": "quasineutral_D_plus_T", "startup_seed_electron_density_initial_guess_input_used": False})
    elif domains is not None and config.plasma_closure.model == "egedal_beam_plasma_quasineutral":
        background = build_ion_seed_reference_profiles(full_centers, full_cell_volumes, domains.central_cell_plasma_domain, p.background_profile, magnetic_midplane_z_m=domains.midplane_z_m, deuterium_midplane_density_m3=p.background_deuterium_midplane_density_m3, tritium_midplane_density_m3=p.background_tritium_midplane_density_m3, ion_temperature_keV=p.background_ion_temperature_keV, symmetry_tolerance=config.kinetic_electrostatic.modal_eq42_density_symmetry_tolerance)
        extra.update(background.as_metadata())

    if background is not None:
        extra.update({"plasma_closure_model": p.model, "electron_density_mode": p.electron_density_mode, "electron_density_initial_guess_m3": p.electron_density_initial_guess_m3, "configured_seed_reference_midpoint_positive_charge_m3": float(p.background_deuterium_midplane_density_m3 + p.background_tritium_midplane_density_m3), "electron_density_initial_guess_valid": p.electron_density_initial_guess_valid, "electron_density_initial_guess_source": p.electron_density_initial_guess_source, "electron_density_initial_guess_is_authoritative_physical_density": False, "electron_temperature_initial_guess_keV": p.electron_temperature_initial_guess_keV})

    metadata = {
        "geometry_model": g.model,
        "magnetic_field_model": magnetic_model,
        "linked_to_paratan_geometry": linked,
        "B0_T": B0_T,
        "mirror_ratio": mirror_ratio,
        "plasma_radius_m": g.plasma_radius_m,
        "half_length_m": half_length,
        "confined_orbit_midpoint_m": orbit_midpoint,
        "axial_bins": g.axial_bins,
        "source_volume_m3": float(np.sum(cell_volumes)),
        "coordinate_edges_m": z_edges,
        "coordinate_centers_m": z_centers,
        "z_edges_cm": 100.0 * z_edges,
        "z_centers_cm": 100.0 * z_centers,
        "zeta_edges": zeta_edges,
        "zeta_centers": zeta_centers,
        "B_tilde_centers": B_tilde,
        "cell_volumes_m3": cell_volumes,
        "magnetic_field_visual_coordinate_m": full_centers,
        "magnetic_field_visual_B_tilde": full_B_tilde,
        "magnetic_field_visual_B_T": full_B_T,
        "magnetic_field_visual_flux_tube_radius_m": full_flux_radius,
        "magnetic_field_visual_flux_tube_radius_edges_m": full_flux_radius_edges,
        "source_geometry_linkage_source_z_min_m": float(z_edges[0]),
        "source_geometry_linkage_source_z_max_m": float(z_edges[-1]),
        "full_device_z_edges_m": full_edges,
        "full_device_z_centers_m": full_centers,
        "full_device_cell_volumes_m3": full_cell_volumes,
        **extra,
    }

    return GeometryStageResult(
        z_edges_m=z_edges,
        z_centers_m=z_centers,
        zeta_edges=zeta_edges,
        zeta_centers=zeta_centers,
        B_tilde_centers=B_tilde,
        B_tilde_midpoints=B_tilde.copy(),
        B0_T=B0_T,
        mirror_ratio=mirror_ratio,
        half_length_m=half_length,
        plasma_radius_m=g.plasma_radius_m,
        radial_outer_radius_m_by_z=radial_outer,
        radial_inner_radius_m_by_z=radial_inner,
        midplane_area_m2=midplane_area,
        cell_volumes_m3=cell_volumes,
        volume_m3=float(np.sum(cell_volumes)),
        magnetic_field_model=magnetic_model,
        linked_to_paratan_geometry=linked,
        metadata=metadata,
        B_tilde_function=B_tilde_function,
        full_device_z_edges_m=full_edges,
        full_device_z_centers_m=full_centers,
        full_device_B_T_centers=full_B_T,
        full_device_B_tilde_centers=full_B_tilde,
        full_device_cell_volumes_m3=full_cell_volumes,
        domains=domains,
        plasma_boundaries=boundaries,
        background_profiles=background,
        B_T_function=B_T_function,
    )
