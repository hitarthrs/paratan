"""better quality source model figures"""
from __future__ import annotations
import argparse
from collections.abc import Mapping
from dataclasses import dataclass
from pathlib import Path
from typing import Any
import numpy as np
from source_model_revamp.neutrons.events import CorrelatedNeutronEventBank
from source_model_revamp.neutrons.events.spatial_sampling import sample_axisymmetric_cell_positions
from source_model_revamp.radial_profiles import KOTELNIKOV_PARABOLIC_FLUX_K1, radial_probability_cdf, radial_source_shape
from source_model_revamp.plotting.common import MEV_TO_J, MissingMetadata, _as_array, _beam_path_arrays, _beam_projected_polygon, _centers_from_edges, _configure_matplotlib, _first_scalar, _plasma_radius_profile, _source_inner_radius_profile, _source_rate_matrix, _volumes_for_length
from source_model_revamp.plotting.model_view import build_model_view

better_RC = {"font.family": "serif", "font.serif": ["DejaVu Serif"], "mathtext.fontset": "stix", "font.size": 9.5, "axes.titlesize": 11.0, "axes.labelsize": 9.5, "xtick.labelsize": 8.5, "ytick.labelsize": 8.5, "legend.fontsize": 8.0, "axes.linewidth": 0.8, "lines.linewidth": 1.5, "savefig.bbox": "tight", "savefig.facecolor": "white", "figure.facecolor": "white"}
BLACK = "#202020"
BLUE = "#0072B2"
SKY = "#56B4E9"
GREEN = "#009E73"
ORANGE = "#E69F00"
VERMILION = "#D55E00"
PURPLE = "#CC79A7"
GRAY = "#7A7A7A"
LIGHT_GRAY = "#E8E8E8"
PALE_BLUE = "#DCEEF8"
PALE_GREEN = "#D8EFE5"
PALE_ORANGE = "#F8E5D3"

def _save_better(fig: Any, output_path: Path) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.canvas.draw()
    fig.savefig(output_path, dpi=320)
    plt = _configure_matplotlib()
    plt.close(fig)

def _mapping(value: Any) -> Mapping[str, Any]:
    return value if isinstance(value, Mapping) else {}

def _domain_from_mapping(mapping: Mapping[str, Any], key: str) -> tuple[float, float] | None:
    value = _mapping(mapping.get(key))
    try:
        lower = float(value["z_min_m"])
        upper = float(value["z_max_m"])
    except Exception:
        return None
    if not np.isfinite(lower) or not np.isfinite(upper) or upper <= lower:
        return None
    return lower, upper

def _layout(metadata: Mapping[str, Any]) -> Mapping[str, Any]:
    value = _mapping(metadata.get("paratan_layout"))
    if value:
        return value
    return _mapping(metadata.get("device_domains"))

def _layout_domain(metadata: Mapping[str, Any], key: str) -> tuple[float, float] | None:
    return _domain_from_mapping(_layout(metadata), key)

def _within(z: np.ndarray, interval: tuple[float, float] | None) -> np.ndarray:
    if interval is None:
        return np.zeros(z.shape, dtype=bool)
    return (z >= interval[0]) & (z <= interval[1])

def _vessel_envelope(metadata: Mapping[str, Any]) -> tuple[np.ndarray, np.ndarray]:
    layout = _mapping(metadata.get("paratan_layout"))
    if not layout:
        return np.asarray([], dtype=float), np.asarray([], dtype=float)
    keys = ("central_cell_plasma_domain", "left_conical_vacuum_domain", "right_conical_vacuum_domain", "left_bottleneck_vacuum_domain", "right_bottleneck_vacuum_domain", "left_end_cell_vacuum_domain", "right_end_cell_vacuum_domain", "full_device_domain")
    domains = {key: _domain_from_mapping(layout, key) for key in keys}
    full = domains["full_device_domain"]
    if full is None:
        return np.asarray([], dtype=float), np.asarray([], dtype=float)
    radii = np.asarray([ layout.get("central_cell_vacuum_radius_m", np.nan), layout.get("bottleneck_vacuum_radius_m", np.nan), layout.get("end_cell_vacuum_radius_m", np.nan),], dtype=float,)
    if not np.all(np.isfinite(radii)) or np.any(radii <= 0.0):
        return np.asarray([], dtype=float), np.asarray([], dtype=float)
    central_radius, bottleneck_radius, end_radius = radii
    z = np.linspace(full[0], full[1], 1800)
    radius = np.full(z.shape, np.nan, dtype=float)
    central = domains["central_cell_plasma_domain"]
    left_conical = domains["left_conical_vacuum_domain"]
    right_conical = domains["right_conical_vacuum_domain"]
    left_bottle = domains["left_bottleneck_vacuum_domain"]
    right_bottle = domains["right_bottleneck_vacuum_domain"]
    left_end = domains["left_end_cell_vacuum_domain"]
    right_end = domains["right_end_cell_vacuum_domain"]
    radius[_within(z, left_bottle) | _within(z, right_bottle)] = bottleneck_radius
    radius[_within(z, central)] = central_radius
    if left_conical is not None:
        mask = _within(z, left_conical)
        radius[mask] = np.interp(z[mask], [left_conical[0], left_conical[1]], [bottleneck_radius, central_radius])
    if right_conical is not None:
        mask = _within(z, right_conical)
        radius[mask] = np.interp(z[mask], [right_conical[0], right_conical[1]], [central_radius, bottleneck_radius])
    radius[_within(z, left_end) | _within(z, right_end)] = end_radius
    return z, radius

def _throat_positions(metadata: Mapping[str, Any]) -> tuple[float, ...]:
    values: list[float] = []
    for keys in (("left_fitted_mirror_throat_z_m",), ("right_fitted_mirror_throat_z_m", "fitted_mirror_throat_z_m")):
        value = _first_scalar(metadata, keys)
        if value is not None and all(abs(value - other) > 1.0e-12 for other in values):
            values.append(value)
    return tuple(values)

def _plot_domain_shading(ax: Any, metadata: Mapping[str, Any]) -> None:
    central = _layout_domain(metadata, "central_cell_plasma_domain")
    if central is not None:
        ax.axvspan(central[0], central[1], color=PALE_BLUE, alpha=0.42, zorder=0)
    for value in _throat_positions(metadata):
        ax.axvline(value, color=GRAY, linestyle="--", linewidth=0.95, zorder=1)

def _coil_centerlines(ax: Any, metadata: Mapping[str, Any]) -> None:
    z = _as_array(metadata, "coil_z_center_m", ndim=1)
    radius = _as_array(metadata, "coil_radius_m", ndim=1)
    groups = metadata.get("coil_group", ())
    if z.size == 0 or radius.size != z.size:
        return
    seen: set[str] = set()
    for index, (z_value, radius_value) in enumerate(zip(z, radius, strict=True)):
        group = str(groups[index]) if index < len(groups) else "coil"
        is_hf = "hf" in group.lower()
        color = PURPLE if is_hf else ORANGE
        label = "HF coil centerline" if is_hf else "LF coil centerline"
        if label in seen:
            label = "_nolegend_"
        else:
            seen.add(label)
        ax.scatter([z_value, z_value], [radius_value, -radius_value], marker="s", s=24.0, color=color, edgecolor="white", linewidth=0.45, label=label, zorder=5)

def _significant_axial_source_mask(rate_by_z: np.ndarray) -> np.ndarray:
    rates = np.asarray(rate_by_z, dtype=float)
    if rates.ndim != 1 or rates.size == 0:
        return np.zeros(rates.shape, dtype=bool)
    rates = np.where(np.isfinite(rates) & (rates > 0.0), rates, 0.0)
    total = float(np.sum(rates))
    maximum = float(np.max(rates, initial=0.0))
    if total <= 0.0 or maximum <= 0.0:
        return np.zeros(rates.shape, dtype=bool)
    relative_mask = rates >= maximum * 1.0e-7
    cumulative = np.cumsum(rates) / total
    lower = int(np.searchsorted(cumulative, 5.0e-5, side="left"))
    upper = int(np.searchsorted(cumulative, 1.0 - 5.0e-5, side="left"))
    cumulative_mask = np.zeros(rates.shape, dtype=bool)
    cumulative_mask[max(0, lower) : min(rates.size, upper + 1)] = True
    active = relative_mask & cumulative_mask
    if not np.any(active):
        active = cumulative_mask
    return active

def _source_radial_profile_model(metadata: Mapping[str, Any]) -> str:
    value = metadata.get("neutron_radial_source_profile_model")
    if value is None:
        value = metadata.get("expander_radial_profile_model")
    return KOTELNIKOV_PARABOLIC_FLUX_K1 if value is None else str(value)

def _source_connected_rho_edges(metadata: Mapping[str, Any], edge_count: int) -> np.ndarray:
    for key in ("neutron_source_connected_radial_rho_limit_edges", "expander_source_connected_radial_rho_limit_edges_full_device"):
        value = _as_array(metadata, key, ndim=1)
        if value.size:
            if value.shape != (edge_count,):
                raise MissingMetadata(f"{key} must match the neutron axial edge grid")
            return np.asarray(value, dtype=float)
    return np.ones(edge_count, dtype=float)

def _source_envelope(metadata: Mapping[str, Any]) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    matrix, z_edges, _ = _source_rate_matrix(metadata)
    z_centers = _centers_from_edges(z_edges)
    outer = _plasma_radius_profile(metadata, z_centers)
    inner = _source_inner_radius_profile(z_centers)
    if outer.size != z_centers.size or inner.size != z_centers.size:
        raise MissingMetadata("missing source radial envelope")
    rho_edges = _source_connected_rho_edges(metadata, z_edges.size)
    rho_centers = 0.5 * (rho_edges[:-1] + rho_edges[1:])
    effective_outer = outer * rho_centers
    active = _significant_axial_source_mask(np.sum(matrix, axis=1))
    return z_centers, effective_outer, active

def _add_source_volume(ax: Any, metadata: Mapping[str, Any], label: str) -> None:
    try:
        z_centers, outer, active = _source_envelope(metadata)
    except MissingMetadata:
        return
    if not np.any(active):
        return
    ax.fill_between(z_centers, -outer, outer, where=active, step="mid", color=PALE_GREEN, alpha=0.78, label=label, zorder=2)
    matrix, _, _ = _source_rate_matrix(metadata)
    nonzero = np.sum(matrix, axis=1) > 0.0
    if np.any(nonzero):
        upper = np.where(nonzero, outer, np.nan)
        ax.plot(z_centers, upper, color=GREEN, linestyle=":", linewidth=1.0, label="Nonzero neutron source support", zorder=3)
        ax.plot(z_centers, -upper, color=GREEN, linestyle=":", linewidth=1.0, zorder=3)

@dataclass(frozen=True)
class _GeometryShellSurface:
    geometry: str
    surface_id: str
    component: str
    side: str
    z_min_m: float | None = None
    z_max_m: float | None = None
    radius_m: float | None = None
    radius_at_z_min_m: float | None = None
    radius_at_z_max_m: float | None = None
    z_m: float | None = None
    radius_min_m: float | None = None
    radius_max_m: float | None = None

def _finite_float(value: Any) -> float | None:
    try:
        number = float(value)
    except Exception:
        return None
    return number if np.isfinite(number) else None

def _geometry_shell_surfaces(metadata: Mapping[str, Any]) -> tuple[_GeometryShellSurface, ...]:
    raw = metadata.get("candidate_plasma_facing_surfaces")
    if not isinstance(raw, (list, tuple)):
        return ()
    surfaces: list[_GeometryShellSurface] = []
    seen: set[str] = set()
    for item in raw:
        if not isinstance(item, Mapping):
            continue
        geometry = str(item.get("geometry", "")).strip().lower()
        surface_id = str(item.get("surface_id", "")).strip()
        if not geometry:
            continue
        dedupe_key = surface_id or repr(sorted(item.items(), key=lambda pair: str(pair[0])))
        if dedupe_key in seen:
            continue
        seen.add(dedupe_key)
        component = str(item.get("component", "")).strip()
        side = str(item.get("side", "")).strip()
        if geometry == "radial_cylinder":
            z_min = _finite_float(item.get("z_min_m"))
            z_max = _finite_float(item.get("z_max_m"))
            radius = _finite_float(item.get("radius_m"))
            if None in (z_min, z_max, radius) or z_max <= z_min or radius <= 0.0:
                continue
            surfaces.append(_GeometryShellSurface(geometry=geometry, surface_id=surface_id, component=component, side=side, z_min_m=z_min, z_max_m=z_max, radius_m=radius))
            continue
        if geometry == "conical_frustum":
            z_min = _finite_float(item.get("z_min_m"))
            z_max = _finite_float(item.get("z_max_m"))
            radius_min = _finite_float(item.get("radius_at_z_min_m"))
            radius_max = _finite_float(item.get("radius_at_z_max_m"))
            if (
                None in (z_min, z_max, radius_min, radius_max)
                or z_max <= z_min
                or radius_min <= 0.0
                or radius_max <= 0.0
            ):
                continue
            surfaces.append(_GeometryShellSurface(geometry=geometry, surface_id=surface_id, component=component, side=side, z_min_m=z_min, z_max_m=z_max, radius_at_z_min_m=radius_min, radius_at_z_max_m=radius_max))
            continue
        if geometry == "axial_plane":
            z_value = _finite_float(item.get("z_m"))
            radius_min = _finite_float(item.get("radius_min_m"))
            radius_max = _finite_float(item.get("radius_max_m"))
            if z_value is None:
                continue
            radius_min = 0.0 if radius_min is None else max(0.0, radius_min)
            radius_max = radius_min if radius_max is None else max(radius_min, radius_max)
            if radius_max <= 0.0:
                continue
            surfaces.append(_GeometryShellSurface(geometry=geometry, surface_id=surface_id, component=component, side=side, z_m=z_value, radius_min_m=radius_min, radius_max_m=radius_max))
    return tuple(sorted(surfaces, key=lambda surface: (float(surface.z_m if surface.z_m is not None else surface.z_min_m or 0.0), surface.geometry, surface.surface_id),))

def _geometry_shell_radius_at_z(surfaces: tuple[_GeometryShellSurface, ...], z_m: float) -> float | None:
    matches: list[float] = []
    for surface in surfaces:
        if surface.geometry == "radial_cylinder":
            if (
                surface.z_min_m is not None
                and surface.z_max_m is not None
                and surface.radius_m is not None
                and surface.z_min_m - 1.0e-12 <= z_m <= surface.z_max_m + 1.0e-12
            ):
                matches.append(surface.radius_m)
            continue
        if surface.geometry == "conical_frustum":
            if (
                surface.z_min_m is not None
                and surface.z_max_m is not None
                and surface.radius_at_z_min_m is not None
                and surface.radius_at_z_max_m is not None
                and surface.z_min_m - 1.0e-12 <= z_m <= surface.z_max_m + 1.0e-12
            ):
                matches.append(float(np.interp(z_m, [surface.z_min_m, surface.z_max_m], [surface.radius_at_z_min_m, surface.radius_at_z_max_m],)))
            continue
        if (
            surface.geometry == "axial_plane"
            and surface.z_m is not None
            and surface.radius_max_m is not None
            and abs(surface.z_m - z_m) <= 1.0e-12
        ):
            matches.append(surface.radius_max_m)
    if not matches:
        return None
    return float(max(matches))

def _plot_geometry_shell_cross_section(ax: Any, metadata: Mapping[str, Any], label: str) -> bool:
    surfaces = _geometry_shell_surfaces(metadata)
    if not surfaces:
        return False
    labeled = False
    for surface in surfaces:
        plot_label = label if not labeled else "_nolegend_"
        if (
            surface.geometry == "radial_cylinder"
            and surface.z_min_m is not None
            and surface.z_max_m is not None
            and surface.radius_m is not None
        ):
            z = np.asarray([surface.z_min_m, surface.z_max_m], dtype=float)
            r = np.full(2, surface.radius_m, dtype=float)
            ax.plot(z, r, color=BLACK, linewidth=1.1, label=plot_label, zorder=0)
            ax.plot(z, -r, color=BLACK, linewidth=1.1, zorder=0)
            labeled = True
            continue
        if (
            surface.geometry == "conical_frustum"
            and surface.z_min_m is not None
            and surface.z_max_m is not None
            and surface.radius_at_z_min_m is not None
            and surface.radius_at_z_max_m is not None
        ):
            z = np.asarray([surface.z_min_m, surface.z_max_m], dtype=float)
            r = np.asarray([surface.radius_at_z_min_m, surface.radius_at_z_max_m], dtype=float)
            ax.plot(z, r, color=BLACK, linewidth=1.1, label=plot_label, zorder=0)
            ax.plot(z, -r, color=BLACK, linewidth=1.1, zorder=0)
            labeled = True
            continue
        if (
            surface.geometry == "axial_plane"
            and surface.z_m is not None
            and surface.radius_max_m is not None
        ):
            radius_min = 0.0 if surface.radius_min_m is None else surface.radius_min_m
            radius_max = max(radius_min, surface.radius_max_m)
            ax.plot([surface.z_m, surface.z_m], [radius_min, radius_max], color=BLACK, linewidth=1.0, label=plot_label, zorder=0)
            if radius_min > 0.0:
                ax.plot([surface.z_m, surface.z_m], [-radius_max, -radius_min], color=BLACK, linewidth=1.0, zorder=0)
            else:
                ax.plot([surface.z_m, surface.z_m], [-radius_max, radius_max], color=BLACK, linewidth=1.0, zorder=0)
            labeled = True
    return labeled

def _plot_ring_3d(ax: Any, z_cm: float, radius_cm: float, *, color: str, linewidth: float, alpha: float, linestyle: str = "-") -> None:
    theta = np.linspace(0.0, 2.0 * np.pi, 161)
    ax.plot(np.full(theta.size, z_cm), radius_cm * np.cos(theta), radius_cm * np.sin(theta), color=color, linewidth=linewidth, alpha=alpha, linestyle=linestyle)

def _plot_cylindrical_shell_3d(ax: Any, z_min_cm: float, z_max_cm: float, radius_cm: float) -> None:
    axial_samples = np.linspace(z_min_cm, z_max_cm, 6)
    for z_value in axial_samples:
        _plot_ring_3d(ax, float(z_value), radius_cm, color=GRAY, linewidth=0.85, alpha=0.46)
    z = np.linspace(z_min_cm, z_max_cm, 160)
    for angle in np.linspace(0.0, 2.0 * np.pi, 8, endpoint=False):
        ax.plot(z, np.full(z.size, radius_cm * np.cos(angle)), np.full(z.size, radius_cm * np.sin(angle)), color=GRAY, linewidth=0.85, alpha=0.42)

def _plot_conical_shell_3d(ax: Any, z_min_cm: float, z_max_cm: float, radius_min_cm: float, radius_max_cm: float) -> None:
    axial_samples = np.linspace(z_min_cm, z_max_cm, 6)
    radii = np.interp(axial_samples, [z_min_cm, z_max_cm], [radius_min_cm, radius_max_cm])
    for z_value, radius_value in zip(axial_samples, radii, strict=True):
        _plot_ring_3d(ax, float(z_value), float(radius_value), color=GRAY, linewidth=0.85, alpha=0.46)
    z = np.linspace(z_min_cm, z_max_cm, 180)
    radius = np.interp(z, [z_min_cm, z_max_cm], [radius_min_cm, radius_max_cm])
    for angle in np.linspace(0.0, 2.0 * np.pi, 8, endpoint=False):
        ax.plot(z, radius * np.cos(angle), radius * np.sin(angle), color=GRAY, linewidth=0.85, alpha=0.42)

def _plot_axial_face_3d(ax: Any, z_cm: float, radius_min_cm: float, radius_max_cm: float) -> None:
    _plot_ring_3d(ax, z_cm, radius_max_cm, color=GRAY, linewidth=0.85, alpha=0.46)
    if radius_min_cm > 0.0:
        _plot_ring_3d(ax, z_cm, radius_min_cm, color=GRAY, linewidth=0.75, alpha=0.42, linestyle="--")
    for angle in np.linspace(0.0, 2.0 * np.pi, 8, endpoint=False):
        ax.plot([z_cm, z_cm], [radius_min_cm * np.cos(angle), radius_max_cm * np.cos(angle)], [radius_min_cm * np.sin(angle), radius_max_cm * np.sin(angle)], color=GRAY, linewidth=0.75, alpha=0.38)

def _plot_geometry_shell_wireframe(ax: Any, metadata: Mapping[str, Any]) -> tuple[_GeometryShellSurface, ...]:
    surfaces = _geometry_shell_surfaces(metadata)
    for surface in surfaces:
        if (
            surface.geometry == "radial_cylinder"
            and surface.z_min_m is not None
            and surface.z_max_m is not None
            and surface.radius_m is not None
        ):
            _plot_cylindrical_shell_3d(
                ax,
                100.0 * surface.z_min_m,
                100.0 * surface.z_max_m,
                100.0 * surface.radius_m,
            )
            continue
        if (
            surface.geometry == "conical_frustum"
            and surface.z_min_m is not None
            and surface.z_max_m is not None
            and surface.radius_at_z_min_m is not None
            and surface.radius_at_z_max_m is not None
        ):
            _plot_conical_shell_3d(
                ax,
                100.0 * surface.z_min_m,
                100.0 * surface.z_max_m,
                100.0 * surface.radius_at_z_min_m,
                100.0 * surface.radius_at_z_max_m,
            )
            continue
        if (
            surface.geometry == "axial_plane"
            and surface.z_m is not None
            and surface.radius_max_m is not None
        ):
            _plot_axial_face_3d(
                ax,
                100.0 * surface.z_m,
                100.0 * (0.0 if surface.radius_min_m is None else surface.radius_min_m),
                100.0 * surface.radius_max_m,
            )
    return surfaces

def _beam_plot_styles(beams: tuple[Any, ...]) -> Mapping[str, tuple[str, str]]:
    colors = {"deuterium": BLUE, "tritium": VERMILION}
    linestyles = ("-", "--", "-.", ":")
    counts: dict[str, int] = {}
    styles: dict[str, tuple[str, str]] = {}
    for beam in beams:
        count = counts.get(beam.species_id, 0)
        counts[beam.species_id] = count + 1
        styles[beam.beam_id] = (colors.get(beam.species_id, BLACK), linestyles[count % len(linestyles)])
    return styles

def _beam_plot_label(beam: Any) -> str:
    species = "D" if beam.species_id == "deuterium" else "T" if beam.species_id == "tritium" else beam.species_id
    angle = beam.effective_injection_angle_deg
    if angle is None:
        angle = beam.configured_injection_angle_deg
    angle_text = "" if angle is None else f", {angle:.1f}°"
    return f"{beam.beam_id}: {species}, {beam.energy_keV:.0f} keV, {beam.power_W / 1.0e6:.3g} MW{angle_text}"

def _plot_beam_projection_panel(ax: Any, metadata: Mapping[str, Any], view: Any, styles: Mapping[str, tuple[str, str]], coordinate_index: int, full_device: bool) -> None:
    coordinate_name = "x" if coordinate_index == 0 else "y"
    radial_values: list[float] = []
    z_values: list[np.ndarray] = []
    if full_device:
        has_direct_shell = _plot_geometry_shell_cross_section(ax, metadata, "Plasma facing shell")
        if not has_direct_shell:
            vessel_z, vessel_radius = _vessel_envelope(metadata)
            if vessel_z.size:
                ax.fill_between(vessel_z, -vessel_radius, vessel_radius, where=np.isfinite(vessel_radius), color=LIGHT_GRAY, alpha=0.45, label="Inferred chamber envelope", zorder=0)
                ax.plot(vessel_z, vessel_radius, color=BLACK, linewidth=1.15)
                ax.plot(vessel_z, -vessel_radius, color=BLACK, linewidth=1.15)
                radial_values.append(float(np.nanmax(vessel_radius)))
                z_values.append(vessel_z)
        field_z = _as_array(metadata, "magnetic_field_visual_coordinate_m", ndim=1)
        flux_radius = _as_array(metadata, "magnetic_field_visual_flux_tube_radius_m", ndim=1)
        if field_z.size >= 2 and flux_radius.size == field_z.size:
            ax.fill_between(field_z, -flux_radius, flux_radius, color=PALE_BLUE, alpha=0.82, label="Modeled magnetic flux tube", zorder=1)
            ax.plot(field_z, flux_radius, color=BLUE, linewidth=1.3)
            ax.plot(field_z, -flux_radius, color=BLUE, linewidth=1.3)
            radial_values.append(float(np.nanmax(flux_radius)))
            z_values.append(field_z)
        _add_source_volume(ax, metadata, "Neutron source volume")
        _coil_centerlines(ax, metadata)
        for index, throat in enumerate(_throat_positions(metadata)):
            ax.axvline(throat, color=GRAY, linestyle="--", linewidth=1.0, label="Mirror throats" if index == 0 else "_nolegend_", zorder=3)
    else:
        confined_grid = view.grid("confined_kinetic")
        confined_radius = _plasma_radius_profile(metadata, confined_grid.centers_m)
        if confined_radius.size == confined_grid.centers_m.size:
            ax.fill_between(confined_grid.centers_m, -confined_radius, confined_radius, color=PALE_BLUE, alpha=0.86, label="Confined flux tube", zorder=0)
            ax.plot(confined_grid.centers_m, confined_radius, color=BLUE, linewidth=1.25)
            ax.plot(confined_grid.centers_m, -confined_radius, color=BLUE, linewidth=1.25)
            radial_values.append(float(np.nanmax(confined_radius)))
            z_values.append(confined_grid.centers_m)
        _plot_geometry_shell_cross_section(ax, metadata, "Plasma facing shell")
        _add_source_volume(ax, metadata, "Neutron source volume")
    for beam in view.beams:
        centers, radii = _beam_path_arrays(beam)
        polygon = _beam_projected_polygon(beam, coordinate_index)
        color, linestyle = styles[beam.beam_id]
        if polygon.size:
            ax.fill(polygon[:, 0], polygon[:, 1], color=color, alpha=0.13, zorder=4)
        ax.plot(centers[:, 2], centers[:, coordinate_index], color=color, linestyle=linestyle, linewidth=2.0, label=_beam_plot_label(beam), zorder=5)
        ax.scatter([centers[0, 2], centers[-1, 2]], [centers[0, coordinate_index], centers[-1, coordinate_index]], s=26, facecolor="white", edgecolor=color, zorder=6)
        radial_values.append(float(np.nanmax(np.abs(centers[:, coordinate_index])) + np.nanmax(radii)))
        z_values.append(centers[:, 2])
    if z_values:
        all_z = np.concatenate(z_values)
        lower = float(np.nanmin(all_z))
        upper = float(np.nanmax(all_z))
        span = max(upper - lower, 0.1)
        if full_device:
            ax.set_xlim(lower, upper)
        else:
            ax.set_xlim(lower - 0.12 * span, upper + 0.12 * span)
    radial = max(radial_values, default=0.1)
    ax.set_ylim(-1.12 * radial, 1.12 * radial)
    ax.axhline(0.0, color=BLACK, linewidth=0.6, alpha=0.5)
    ax.set_xlabel("Axial coordinate z [m]")
    ax.set_ylabel(f"Transverse coordinate {coordinate_name} [m]")
    ax.grid(True, alpha=0.16)

def _geometry_shell_axial_bounds_m(surfaces: tuple[_GeometryShellSurface, ...],) -> tuple[float, float] | None:
    bounds: list[float] = []
    for surface in surfaces:
        if surface.z_m is not None:
            bounds.append(float(surface.z_m))
        if surface.z_min_m is not None:
            bounds.append(float(surface.z_min_m))
        if surface.z_max_m is not None:
            bounds.append(float(surface.z_max_m))
    if not bounds:
        return None
    return float(min(bounds)), float(max(bounds))

def _geometry_shell_radius_max_m(surfaces: tuple[_GeometryShellSurface, ...]) -> float:
    radii: list[float] = []
    for surface in surfaces:
        if surface.radius_m is not None:
            radii.append(float(surface.radius_m))
        if surface.radius_at_z_min_m is not None:
            radii.append(float(surface.radius_at_z_min_m))
        if surface.radius_at_z_max_m is not None:
            radii.append(float(surface.radius_at_z_max_m))
        if surface.radius_max_m is not None:
            radii.append(float(surface.radius_max_m))
    return max(radii, default=0.0)

def _full_device_context_limits_m(metadata: Mapping[str, Any], view: Any | None = None) -> tuple[float, float, float]:
    z_lower: float | None = None
    z_upper: float | None = None
    radial_limit: float = 0.0

    def extend_z(values: np.ndarray) -> None:
        nonlocal z_lower, z_upper
        if values.size == 0:
            return
        finite = values[np.isfinite(values)]
        if finite.size == 0:
            return
        local_lower = float(np.min(finite))
        local_upper = float(np.max(finite))
        z_lower = local_lower if z_lower is None else min(z_lower, local_lower)
        z_upper = local_upper if z_upper is None else max(z_upper, local_upper)

    field_z = _as_array(metadata, "magnetic_field_visual_coordinate_m", ndim=1)
    flux_radius = _as_array(metadata, "magnetic_field_visual_flux_tube_radius_m", ndim=1)
    if field_z.size and flux_radius.size == field_z.size:
        extend_z(field_z)
        finite_radius = flux_radius[np.isfinite(flux_radius)]
        if finite_radius.size:
            radial_limit = max(radial_limit, float(np.max(finite_radius)))

    vessel_z, vessel_radius = _vessel_envelope(metadata)
    if vessel_z.size and vessel_radius.size == vessel_z.size:
        extend_z(vessel_z)
        finite_radius = vessel_radius[np.isfinite(vessel_radius)]
        if finite_radius.size:
            radial_limit = max(radial_limit, float(np.max(finite_radius)))

    surfaces = _geometry_shell_surfaces(metadata)
    shell_bounds = _geometry_shell_axial_bounds_m(surfaces)
    if shell_bounds is not None:
        extend_z(np.asarray(shell_bounds, dtype=float))
        radial_limit = max(radial_limit, _geometry_shell_radius_max_m(surfaces))

    layout_full = _layout_domain(metadata, "full_device_domain")
    if layout_full is not None:
        extend_z(np.asarray(layout_full, dtype=float))

    throats = np.asarray(_throat_positions(metadata), dtype=float)
    extend_z(throats)

    if view is not None:
        for beam in view.beams:
            centers, radii = _beam_path_arrays(beam)
            if centers.size:
                extend_z(centers[:, 2])
                radial_limit = max(radial_limit, float(np.nanmax(np.abs(centers[:, 0])) + np.nanmax(radii)), float(np.nanmax(np.abs(centers[:, 1])) + np.nanmax(radii)))

    if z_lower is None or z_upper is None:
        raise MissingMetadata("unable to determine full device plotting limits")

    span = max(z_upper - z_lower, 0.1)
    padding = 0.04 * span
    if radial_limit <= 0.0:
        radial_limit = 0.25 * span
    return z_lower - padding, z_upper + padding, radial_limit

def _active_source_interval_m(metadata: Mapping[str, Any]) -> tuple[float, float] | None:
    try:
        matrix, z_edges, _ = _source_rate_matrix(metadata)
    except MissingMetadata:
        return None
    active = np.flatnonzero(_significant_axial_source_mask(np.sum(matrix, axis=1)))
    if active.size == 0:
        return None
    return float(z_edges[active[0]]), float(z_edges[active[-1] + 1])

def plot_full_device_cross_section(metadata: Mapping[str, Any], output_path: Path, args: argparse.Namespace) -> None:
    plt = _configure_matplotlib()
    view = build_model_view(metadata)
    if not view.beams:
        raise MissingMetadata("no enabled beam records")
    styles = _beam_plot_styles(view.beams)
    z_lower, z_upper, radial_limit = _full_device_context_limits_m(metadata, view)
    source_interval = _active_source_interval_m(metadata)
    source_matrix, source_edges, _ = _source_rate_matrix(metadata)
    nonzero_source = np.flatnonzero(np.sum(source_matrix, axis=1) > 0.0)
    nonzero_source_span_text = "unavailable"
    if nonzero_source.size:
        nonzero_source_span_text = f"{100.0 * (source_edges[nonzero_source[-1] + 1] - source_edges[nonzero_source[0]]):.1f} cm"
    physical_domain = _layout_domain(metadata, "full_device_domain")
    physical_device_length_cm = 100.0 * (physical_domain[1] - physical_domain[0]) if physical_domain is not None else 100.0 * (z_upper - z_lower)
    throat_positions = _throat_positions(metadata)
    throat_spacing_text = "unavailable"
    if len(throat_positions) >= 2:
        throat_spacing_text = f"{100.0 * (max(throat_positions) - min(throat_positions)):.1f} cm"
    source_span_text = "unavailable"
    if source_interval is not None:
        source_span_text = f"{100.0 * (source_interval[1] - source_interval[0]):.1f} cm"

    with plt.rc_context(better_RC):
        fig, axes = plt.subplots(1, 2, figsize=(11.6, 4.9), sharey=True, constrained_layout=True)
        _plot_beam_projection_panel(axes[0], metadata, view, styles, 0, True)
        _plot_beam_projection_panel(axes[1], metadata, view, styles, 1, True)
        for ax in axes:
            ax.set_xlim(z_lower, z_upper)
            ax.set_ylim(-1.08 * radial_limit, 1.08 * radial_limit)
        axes[0].set_title("(a) Full device x-z cross section", loc="left")
        axes[1].set_title("(b) Full device y-z cross section", loc="left")

        handles: list[Any] = []
        labels: list[str] = []
        for ax in axes:
            local_handles, local_labels = ax.get_legend_handles_labels()
            for handle, label in zip(local_handles, local_labels, strict=True):
                if label and label != "_nolegend_" and label not in labels:
                    handles.append(handle)
                    labels.append(label)
        fig.legend(handles, labels, loc="outside upper center", ncol=3, frameon=True)

        axes[1].text(
            0.985,
            0.03,
            f"Full device length: {physical_device_length_cm:.1f} cm\n"
            f"Throat spacing: {throat_spacing_text}\n"
            f"Dominant source span: {source_span_text}\n"
            f"Nonzero source span: {nonzero_source_span_text}",
            transform=axes[1].transAxes,
            ha="right",
            va="bottom",
            fontsize=8.0,
            bbox={"facecolor": "white", "alpha": 0.88, "edgecolor": "none"},
        )
        fig.suptitle("Full device geometry cross section with source and beam context")
        _save_better(fig, output_path)

def _step(ax: Any, edges: np.ndarray, values: np.ndarray, **kwargs: Any) -> None:
    if edges.size == values.size + 1:
        ax.stairs(values, edges, **kwargs)
        return
    centers = _centers_from_edges(edges)
    if centers.size == values.size:
        ax.plot(centers, values, **kwargs)

def _positive_range(*arrays: np.ndarray) -> tuple[float, float] | None:
    positive_parts = []
    for array in arrays:
        values = np.asarray(array, dtype=float)
        positive_parts.append(values[np.isfinite(values) & (values > 0.0)])
    available = [part for part in positive_parts if part.size]
    if not available:
        return None
    positive = np.concatenate(available)
    return float(np.min(positive)), float(np.max(positive))

def _log_profile(values: np.ndarray) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    return np.where(np.isfinite(values) & (values > 0.0), values, np.nan)


def _weighted_quantile(values: np.ndarray, weights: np.ndarray, quantiles: tuple[float, ...]) -> np.ndarray:
    x = np.asarray(values, dtype=float)
    w = np.asarray(weights, dtype=float)
    valid = np.isfinite(x) & np.isfinite(w) & (w > 0.0)
    if not np.any(valid):
        return np.full(len(quantiles), np.nan, dtype=float)
    x = x[valid]
    w = w[valid]
    order = np.argsort(x)
    x = x[order]
    w = w[order]
    cumulative = np.cumsum(w)
    cumulative /= cumulative[-1]
    return np.interp(np.asarray(quantiles, dtype=float), cumulative, x)

def _field_crossings(z_m: np.ndarray, b_over_b0: np.ndarray, target: float, lower_m: float, upper_m: float) -> tuple[float, ...]:
    z = np.asarray(z_m, dtype=float)
    b = np.asarray(b_over_b0, dtype=float)
    valid = np.isfinite(z) & np.isfinite(b)
    z = z[valid]
    b = b[valid]
    if z.size < 2 or not np.isfinite(target):
        return ()
    order = np.argsort(z)
    z = z[order]
    b = b[order]
    residual = b - float(target)
    crossings: list[float] = []
    for index in range(z.size - 1):
        z0 = float(z[index])
        z1 = float(z[index + 1])
        y0 = float(residual[index])
        y1 = float(residual[index + 1])
        if z1 < lower_m or z0 > upper_m:
            continue
        if y0 == 0.0:
            crossings.append(z0)
            continue
        if y1 == 0.0:
            crossings.append(z1)
            continue
        if y0 * y1 > 0.0:
            continue
        fraction = -y0 / (y1 - y0)
        crossing = z0 + fraction * (z1 - z0)
        if lower_m <= crossing <= upper_m:
            crossings.append(float(crossing))
    crossings.sort()
    deduplicated: list[float] = []
    for value in crossings:
        if not deduplicated or abs(value - deduplicated[-1]) > 1.0e-8:
            deduplicated.append(value)
    return tuple(deduplicated)

def _beam_magnetic_turning_markers(metadata: Mapping[str, Any], view: Any) -> tuple[tuple[str, tuple[float, ...]], ...]:
    confined = view.grid("confined_kinetic")
    field_z = _as_array(metadata, "magnetic_field_visual_coordinate_m", ndim=1)
    field_b = _as_array(metadata, "magnetic_field_visual_B_tilde", ndim=1)
    if field_z.size == 0 or field_b.size != field_z.size or confined.edges_m.size < 2:
        return ()
    lower_m = float(confined.edges_m[0])
    upper_m = float(confined.edges_m[-1])
    markers: list[tuple[str, tuple[float, ...]]] = []
    for beam in view.beams:
        rate = np.asarray(beam.axial_birth_rate_s, dtype=float)
        lam = np.asarray(beam.birth_lambda_profile, dtype=float)
        if rate.shape != lam.shape:
            continue
        support = np.isfinite(rate) & (rate > 0.0) & np.isfinite(lam) & (lam > 0.0) & (lam <= 1.0)
        if not np.any(support):
            continue
        median_lambda = float(_weighted_quantile(lam[support], rate[support], (0.50,))[0])
        if not np.isfinite(median_lambda) or median_lambda <= 0.0:
            continue
        crossings = _field_crossings(field_z, field_b, 1.0 / median_lambda, lower_m, upper_m)
        if crossings:
            markers.append((beam.species_id, crossings))
    return tuple(markers)

def _add_beam_magnetic_turning_markers(ax: Any, metadata: Mapping[str, Any], view: Any) -> None:
    species_colors = {"deuterium": BLUE, "tritium": ORANGE}
    for species_id, crossings in _beam_magnetic_turning_markers(metadata, view):
        color = species_colors.get(species_id, GRAY)
        for crossing in crossings:
            ax.axvline(crossing, color=color, linestyle=":", linewidth=1.05, alpha=0.82)

def _full_device_flux_radius_at_z(metadata: Mapping[str, Any], z_m: float, full_edges_m: np.ndarray) -> float:
    z_visual = _as_array(metadata, "magnetic_field_visual_coordinate_m", ndim=1)
    radius_visual = _as_array(metadata, "magnetic_field_visual_flux_tube_radius_m", ndim=1)
    radius_edges = _as_array(metadata, "magnetic_field_visual_flux_tube_radius_edges_m", ndim=1)
    if radius_edges.shape == full_edges_m.shape and np.all(np.isfinite(radius_edges)) and np.all(radius_edges > 0.0):
        return float(np.interp(float(z_m), full_edges_m, radius_edges))
    if z_visual.size >= 2 and radius_visual.shape == z_visual.shape and np.all(np.isfinite(z_visual)) and np.all(np.isfinite(radius_visual)) and np.all(radius_visual > 0.0):
        return float(np.interp(float(z_m), z_visual, radius_visual))
    plasma_radius = _first_scalar(metadata, ("plasma_radius_m",))
    if plasma_radius is None or plasma_radius <= 0.0:
        raise MissingMetadata("missing full device magnetic flux tube radius")
    return float(plasma_radius)

def _radial_source_location(metadata: Mapping[str, Any], view: Any, location: str) -> tuple[int, float]:
    full = view.grid("full_device_population")
    if location == "midplane":
        target = _first_scalar(metadata, ("magnetic_midplane_z_m",))
        target_z = 0.0 if target is None else float(target)
    elif location == "throat":
        target = _first_scalar(metadata, ("right_fitted_mirror_throat_z_m", "fitted_mirror_throat_z_m"))
        if target is None:
            raise MissingMetadata("missing right mirror throat position")
        target_z = float(target)
    elif location == "end_cell":
        layout = _layout(metadata)
        target = _finite_float(layout.get("right_end_cell_center_m"))
        if target is None:
            raise MissingMetadata("missing right end cell center")
        target_z = target
    else:
        raise ValueError(f"unknown radial source location {location!r}")
    index = int(np.argmin(np.abs(full.centers_m - target_z)))
    return index, float(full.centers_m[index])

def _better_radial_neutron_source_figure(metadata: Mapping[str, Any], location: str) -> Any:
    plt = _configure_matplotlib()
    view = build_model_view(metadata)
    full = view.grid("full_device_population")
    total_rate = np.asarray(view.fusion_total_axial_neutron_rate_s, dtype=float)
    if total_rate.shape != full.centers_m.shape:
        raise MissingMetadata("full device axial neutron rate is unavailable")
    index, z_center = _radial_source_location(metadata, view, location)
    z_lo = float(full.edges_m[index])
    z_hi = float(full.edges_m[index + 1])
    dz = z_hi - z_lo
    if dz <= 0.0:
        raise MissingMetadata("invalid full device axial cell width")
    nominal_radius = _full_device_flux_radius_at_z(metadata, z_center, full.edges_m)
    rho_edges = _source_connected_rho_edges(metadata, full.edges_m.size)
    rho_limit = float(np.interp(z_center, full.edges_m, rho_edges))
    rho_limit = float(np.clip(rho_limit, 0.0, 1.0))
    if rho_limit <= 0.0:
        raise MissingMetadata(f"radial neutron source support is empty at {location}")
    model = _source_radial_profile_model(metadata)
    support = float(radial_probability_cdf(rho_limit, model))
    if not np.isfinite(support) or support <= 0.0:
        raise MissingMetadata(f"radial neutron source normalization is unavailable at {location}")
    local_rate = float(total_rate[index])
    mean_density = local_rate / (np.pi * nominal_radius**2 * dz)
    rho = np.linspace(0.0, rho_limit, 600)
    radius_m = nominal_radius * rho
    density = mean_density * radial_source_shape(rho, model) / support
    with plt.rc_context(better_RC):
        fig, ax = plt.subplots(figsize=(7.6, 5.2), constrained_layout=True)
        ax.plot(100.0 * radius_m, density, color=GREEN, linewidth=2.0, label="Radial neutron source density")
        ax.axhline(mean_density, color=GRAY, linestyle="--", linewidth=1.25, label="Axial cell mean density")
        ax.set_xlim(0.0, 100.0 * nominal_radius * rho_limit)
        ax.set_xlabel("Radius [cm]")
        ax.set_ylabel(r"Neutron source density [m$^{-3}$ s$^{-1}$]")
        ax.set_title(f"Volume normalized radial neutron source density at {location.replace('_', ' ')}", loc="left")
        ax.grid(True, alpha=0.16)
        ax.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
        ax.legend(loc="upper right")
        ax.text(0.02, 0.04, f"z = {z_center:.3f} m\nAxial cell rate = {local_rate:.3e} n/s", transform=ax.transAxes, ha="left", va="bottom", fontsize=8.0, bbox={"facecolor": "white", "alpha": 0.86, "edgecolor": "none"})
        return fig

def plot_better_radial_neutron_source_midplane(metadata: Mapping[str, Any], output_path: Path, args: argparse.Namespace) -> None:
    fig = _better_radial_neutron_source_figure(metadata, "midplane")
    _save_better(fig, output_path)

def plot_better_radial_neutron_source_throat(metadata: Mapping[str, Any], output_path: Path, args: argparse.Namespace) -> None:
    fig = _better_radial_neutron_source_figure(metadata, "throat")
    _save_better(fig, output_path)

def plot_better_radial_neutron_source_end_cell(metadata: Mapping[str, Any], output_path: Path, args: argparse.Namespace) -> None:
    fig = _better_radial_neutron_source_figure(metadata, "end_cell")
    _save_better(fig, output_path)

def plot_better_fast_ion_density_profiles(metadata: Mapping[str, Any], output_path: Path, args: argparse.Namespace) -> None:
    plt = _configure_matplotlib()
    view = build_model_view(metadata)
    confined = view.grid("confined_kinetic")
    full = view.grid("full_device_population")
    species_colors = {"deuterium": BLUE, "tritium": ORANGE}
    plotted: list[np.ndarray] = []
    with plt.rc_context(better_RC):
        fig, ax = plt.subplots(figsize=(9.6, 5.2), constrained_layout=True)
        for fast in view.fast_species:
            if fast.full_device_density_m3.size == full.centers_m.size:
                values = _log_profile(fast.full_device_density_m3)
                edges = full.edges_m
            elif fast.local_density_m3.size == confined.centers_m.size:
                values = _log_profile(fast.local_density_m3)
                edges = confined.edges_m
            else:
                continue
            if not np.any(np.isfinite(values)):
                continue
            plotted.append(values)
            _step(ax, edges, values, color=species_colors.get(fast.species_id, GRAY), linewidth=1.8, label=f"Fast {fast.symbol}")
        density_range = _positive_range(*plotted)
        if density_range is None:
            raise MissingMetadata("missing fast ion density profiles")
        ax.set_yscale("log")
        ax.set_ylim(density_range[0] / 1.8, density_range[1] * 1.8)
        ax.set_xlim(float(full.edges_m[0]), float(full.edges_m[-1]))
        ax.set_xlabel("Axial position z [m]")
        ax.set_ylabel(r"Density [m$^{-3}$]")
        ax.set_title("Fast ion density profiles", loc="left")
        ax.axvline(_first_scalar(metadata, ("magnetic_midplane_z_m",)) or 0.0, color=BLACK, linewidth=0.8, alpha=0.55)
        _add_beam_magnetic_turning_markers(ax, metadata, view)
        ax.grid(True, which="both", alpha=0.16)
        ax.legend(loc="center left", bbox_to_anchor=(1.01, 0.5))
        _save_better(fig, output_path)

def plot_better_axial_neutron_source_density(metadata: Mapping[str, Any], output_path: Path, args: argparse.Namespace) -> None:
    plt = _configure_matplotlib()
    view = build_model_view(metadata)
    full = view.grid("full_device_population")
    total = _log_profile(view.fusion_total_axial_neutron_rate_density_m3_s)
    if total.size != full.centers_m.size or not np.any(np.isfinite(total)):
        raise MissingMetadata("missing full device neutron source density")
    reaction_profiles: dict[str, np.ndarray] = {}
    for reaction in ("dd_n", "dt_n"):
        profiles = tuple(np.asarray(component.axial_neutron_rate_density_m3_s, dtype=float) for component in view.fusion_components if component.reaction == reaction)
        if profiles:
            reaction_profiles[reaction] = _log_profile(np.sum(profiles, axis=0))
    with plt.rc_context(better_RC):
        fig, ax = plt.subplots(figsize=(9.6, 5.2), constrained_layout=True)
        _step(ax, full.edges_m, total, color=BLUE, linewidth=2.0, label="Total neutron source density")
        if "dd_n" in reaction_profiles:
            _step(ax, full.edges_m, reaction_profiles["dd_n"], color=ORANGE, linewidth=1.5, linestyle=":", label="DD source density")
        if "dt_n" in reaction_profiles:
            _step(ax, full.edges_m, reaction_profiles["dt_n"], color=GREEN, linewidth=1.5, linestyle="--", label="DT source density")
        source_range = _positive_range(total, *reaction_profiles.values())
        if source_range is None:
            raise MissingMetadata("missing neutron source density profiles")
        ax.set_yscale("log")
        ax.set_ylim(source_range[0] / 1.8, source_range[1] * 1.8)
        ax.set_xlim(float(full.edges_m[0]), float(full.edges_m[-1]))
        ax.set_xlabel("Axial position z [m]")
        ax.set_ylabel(r"Neutron source density [m$^{-3}$ s$^{-1}$]")
        ax.set_title("Volume normalized axial neutron production", loc="left")
        ax.axvline(_first_scalar(metadata, ("magnetic_midplane_z_m",)) or 0.0, color=BLACK, linewidth=0.8, alpha=0.55)
        _add_beam_magnetic_turning_markers(ax, metadata, view)
        ax.grid(True, which="both", alpha=0.16)
        ax.legend(loc="center left", bbox_to_anchor=(1.01, 0.5))
        _save_better(fig, output_path)

def plot_better_axial_profiles(metadata: Mapping[str, Any], output_path: Path, args: argparse.Namespace) -> None:
    plt = _configure_matplotlib()
    view = build_model_view(metadata)
    confined_grid = view.grid("confined_kinetic")
    full_grid = view.grid("full_device_population")
    electrostatic = view.electrostatic
    if electrostatic is None:
        raise MissingMetadata("missing electrostatic plotting profiles")
    species_colors = {"deuterium": BLUE, "tritium": VERMILION}
    with plt.rc_context(better_RC):
        fig, axes = plt.subplots(4, 1, figsize=(9.6, 10.2), sharex=True, constrained_layout=True)
        ax_b, ax_phi, ax_density, ax_source = axes

        field_z = _as_array(metadata, "magnetic_field_visual_coordinate_m", ndim=1)
        field_b = _as_array(metadata, "magnetic_field_visual_B_tilde", ndim=1)
        if field_z.size == 0 or field_b.size != field_z.size:
            raise MissingMetadata("missing full device magnetic field profile")
        ax_b.plot(field_z, field_b, color=BLACK, linewidth=1.8)
        ax_b.set_ylabel("$B/B_0$")
        ax_b.set_title("(a) Magnetic field", loc="left")
        _plot_domain_shading(ax_b, metadata)

        _step(ax_phi, confined_grid.edges_m, electrostatic.potential_relative_to_midplane_V / 1.0e3, color=VERMILION, linewidth=1.6, label=r"$\Phi(z)$")
        ax_phi.axhline(0.0, color=BLACK, linewidth=0.65)
        ax_phi.set_ylabel("Potential [kV]")
        ax_phi.set_title("(b) Shared ambipolar potential on the confined domain", loc="left")
        if electrostatic.converged is not None:
            ax_phi.text(0.015, 0.94, "Eq 70 converged" if electrostatic.converged else "Eq 70 not converged", transform=ax_phi.transAxes, va="top", ha="left", fontsize=8.0)
        _plot_domain_shading(ax_phi, metadata)

        density_values: list[np.ndarray] = []
        startup_seed_only = str(metadata.get("background_population_scope", "")).strip().lower() == "startup_seed_initialization_only"
        for fast in view.fast_species:
            if fast.full_device_density_m3.size != full_grid.centers_m.size:
                continue
            color = species_colors.get(fast.species_id, GRAY)
            values = _log_profile(fast.full_device_density_m3)
            density_values.append(values)
            _step(ax_density, full_grid.edges_m, values, color=color, linewidth=1.55, linestyle="--", label=f"Fast {fast.symbol}")
        ion_values = _log_profile(electrostatic.ion_density_m3)
        electron_values = _log_profile(electrostatic.electron_density_m3)
        density_values.extend((ion_values, electron_values))
        _step(ax_density, confined_grid.edges_m, ion_values, color=GRAY, linewidth=1.25, linestyle="-.", label="Total positive ion charge")
        _step(ax_density, confined_grid.edges_m, electron_values, color=BLACK, linewidth=1.35, linestyle=":", label="$n_e$")
        density_range = _positive_range(*density_values)
        if density_range is None:
            raise MissingMetadata("missing axial density profiles")
        ax_density.set_yscale("log")
        ax_density.set_ylim(density_range[0] / 2.0, density_range[1] * 2.0)
        ax_density.set_ylabel("Density [m$^{-3}$]")
        ax_density.set_title("(c) Final represented D/T populations and confined quasineutrality", loc="left")
        if startup_seed_only:
            ax_density.text(0.015, 0.06, "Startup D/T seed omitted from final population curves\ninitialization only", transform=ax_density.transAxes, va="bottom", ha="left", fontsize=7.6)
        ax_density.legend(loc="upper right", ncol=2, fontsize=7.2)
        _plot_domain_shading(ax_density, metadata)

        source_values: list[np.ndarray] = []
        total_neutron = _log_profile(view.fusion_total_axial_neutron_rate_density_m3_s)
        if total_neutron.size == full_grid.centers_m.size:
            source_values.append(total_neutron)
            _step(ax_source, full_grid.edges_m, total_neutron, color=BLACK, linewidth=1.8, label="Total neutron production")
        for reaction, color, linestyle, label in (("dd_n", SKY, ":", "DD neutron production"), ("dt_n", GREEN, "-.", "DT neutron production")):
            profiles = tuple(component.axial_neutron_rate_density_m3_s for component in view.fusion_components if component.reaction == reaction)
            if not profiles:
                continue
            values = _log_profile(np.sum(profiles, axis=0))
            source_values.append(values)
            _step(ax_source, full_grid.edges_m, values, color=color, linewidth=1.4, linestyle=linestyle, label=label)
        for species_id, profile in view.beam_profiles_by_species.items():
            color = species_colors.get(species_id, GRAY)
            symbol = "D" if species_id == "deuterium" else "T" if species_id == "tritium" else species_id
            values = _log_profile(profile.axial_birth_rate_density_m3_s)
            source_values.append(values)
            _step(ax_source, confined_grid.edges_m, values, color=color, linewidth=1.5, linestyle="--", label=f"{symbol} beam births")
        source_range = _positive_range(*source_values)
        if source_range is None:
            raise MissingMetadata("missing beam birth and neutron source profiles")
        if source_range[1] / source_range[0] > 50.0:
            ax_source.set_yscale("log")
            ax_source.set_ylim(source_range[0] / 2.0, source_range[1] * 2.0)
        ax_source.set_ylabel("Rate density [m$^{-3}$ s$^{-1}$]")
        ax_source.set_xlabel("Axial coordinate z [m]")
        ax_source.set_title("(d) Species resolved beam births and neutron production", loc="left")
        ax_source.legend(loc="upper right", ncol=2, fontsize=7.4)
        _plot_domain_shading(ax_source, metadata)

        for ax in axes:
            ax.grid(True, which="both", alpha=0.16)
        ax_source.set_xlim(float(np.min(field_z)), float(np.max(field_z)))
        _save_better(fig, output_path)

def _current_component_label(component_id: str) -> str:
    labels = {'fast_deuterium': 'Fast D', 'fast_tritium': 'Fast T', 'prompt_deuterium': 'Prompt D', 'prompt_tritium': 'Prompt T'}
    return labels.get(component_id, component_id.replace("_", " ").title())

def plot_better_current_balance(metadata: Mapping[str, Any], output_path: Path, args: argparse.Namespace) -> None:
    plt = _configure_matplotlib()
    view = build_model_view(metadata)
    balance = view.current_balance
    if balance is None:
        raise MissingMetadata("Pass 11 current balance metadata is unavailable")
    labels = [_current_component_label(component.component_id) for component in balance.components]
    currents = np.asarray([component.current_A for component in balance.components], dtype=float)
    colors = [BLUE if component.species_id == "deuterium" else VERMILION if component.species_id == "tritium" else GRAY for component in balance.components]
    summary_labels: list[str] = []
    summary_values: list[float] = []
    for label, value in (("Ion target", balance.target_ion_current_A), ("Electron", balance.electron_current_A), ("Left ion", balance.left_ion_current_A), ("Right ion", balance.right_ion_current_A)):
        if value is not None:
            summary_labels.append(label)
            summary_values.append(value)
    with plt.rc_context(better_RC):
        fig, axes = plt.subplots(1, 2, figsize=(10.2, 4.8), constrained_layout=True)
        axes[0].barh(labels, currents, color=colors, alpha=0.82)
        axes[0].invert_yaxis()
        axes[0].set_xlabel("Current magnitude [A]")
        axes[0].set_title("(a) Represented ion loss current components", loc="left")
        axes[0].grid(True, axis="x", alpha=0.16)
        axes[1].bar(summary_labels, summary_values, color=(BLACK, GRAY, BLUE, VERMILION)[:len(summary_values)], alpha=0.82)
        axes[1].set_ylabel("Current magnitude [A]")
        axes[1].set_title("(b) Total and side resolved current balance", loc="left")
        axes[1].tick_params(axis="x", rotation=20)
        axes[1].grid(True, axis="y", alpha=0.16)
        residual = "unavailable" if balance.current_relative_residual is None else f"{balance.current_relative_residual:.3e}"
        barrier = "unavailable" if balance.wall_barrier_energy_keV is None else f"{balance.wall_barrier_energy_keV:.3g} keV"
        asymmetry = "not applicable" if not balance.end_asymmetry_applicable else "unavailable" if balance.end_asymmetry_relative_error is None else f"{balance.end_asymmetry_relative_error:.3e}"
        status = "converged" if balance.converged else "not converged"
        axes[1].text(0.02, 0.97, f"Status: {status}\nRelative residual: {residual}\nShared wall barrier: {barrier}\nEnd asymmetry: {asymmetry}", transform=axes[1].transAxes, va="top", ha="left", fontsize=8.2, bbox={"facecolor": "white", "alpha": 0.86, "edgecolor": "none"})
        fig.suptitle("Pass 11 total represented confined plasma current balance")
        _save_better(fig, output_path)

def _event_bank(metadata: Mapping[str, Any]) -> CorrelatedNeutronEventBank | None:
    value = metadata.get("_correlated_event_bank")
    return value if isinstance(value, CorrelatedNeutronEventBank) else None

def _sample_binned_events(metadata: Mapping[str, Any], count: int, seed: int) -> tuple[np.ndarray, np.ndarray, float]:
    matrix, z_edges, energy_edges = _source_rate_matrix(metadata)
    z_centers = _centers_from_edges(z_edges)
    outer = _plasma_radius_profile(metadata, z_centers)
    inner = _source_inner_radius_profile(z_centers)
    if outer.size != z_centers.size or inner.size != z_centers.size:
        raise MissingMetadata("missing source radial envelope")
    flat_rates = np.asarray(matrix, dtype=float).ravel()
    positive = np.isfinite(flat_rates) & (flat_rates > 0.0)
    if not np.any(positive):
        raise MissingMetadata("source has no positive energy and axial cells")
    probabilities = flat_rates[positive] / np.sum(flat_rates[positive])
    flat_indices = np.flatnonzero(positive)
    rng = np.random.default_rng(int(seed))
    chosen = rng.choice(flat_indices, size=int(count), replace=True, p=probabilities)
    axial_indices, energy_indices = np.unravel_index(chosen, matrix.shape)
    rho_edges = _source_connected_rho_edges(metadata, z_edges.size)
    outer_edges = _as_array(metadata, "magnetic_field_visual_flux_tube_radius_edges_m", ndim=1)
    if outer_edges.size != z_edges.size:
        outer_edges = None
    positions = sample_axisymmetric_cell_positions(axial_indices, z_edges, inner, outer, rng, radial_profile_model=_source_radial_profile_model(metadata), source_connected_radial_rho_limit_edges=rho_edges, radial_outer_radius_m_by_edge=outer_edges)
    energy = _centers_from_edges(energy_edges)[energy_indices] / MEV_TO_J
    total_rate = _first_scalar(metadata,("total_neutron_rate_s", "neutron_fusion_profile_total_rate_s"),)
    return positions, energy, 0.0 if total_rate is None else total_rate

def _weighted_resample_indices(weights: np.ndarray, count: int, seed: int) -> np.ndarray:
    weights = np.asarray(weights, dtype=float)
    positive = np.isfinite(weights) & (weights > 0.0)
    if not np.any(positive):
        raise MissingMetadata("event bank has no positive probability weights")
    probabilities = weights[positive] / np.sum(weights[positive])
    indices = np.flatnonzero(positive)
    rng = np.random.default_rng(int(seed))
    return rng.choice(indices, size=int(count), replace=True, p=probabilities)

def _wireframe(ax: Any, z: np.ndarray, radius: np.ndarray, *, color: str = GRAY, ring_alpha: float = 0.38, line_alpha: float = 0.34, linewidth: float = 0.75) -> None:
    valid = np.isfinite(z) & np.isfinite(radius) & (radius >= 0.0)
    z = z[valid]
    radius = radius[valid]
    if z.size < 2:
        return
    selected = np.unique(np.linspace(0, z.size - 1, min(9, z.size), dtype=int))
    theta = np.linspace(0.0, 2.0 * np.pi, 121)
    for index in selected:
        ax.plot(np.full(theta.size, z[index]), radius[index] * np.cos(theta), radius[index] * np.sin(theta), color=color, alpha=ring_alpha, linewidth=linewidth)
    for angle in np.linspace(0.0, 2.0 * np.pi, 6, endpoint=False):
        ax.plot(z, radius * np.cos(angle), radius * np.sin(angle), color=color, alpha=line_alpha, linewidth=linewidth)

def _active_source_limits(metadata: Mapping[str, Any]) -> tuple[float, float]:
    matrix, z_edges, _ = _source_rate_matrix(metadata)
    active = np.flatnonzero(_significant_axial_source_mask(np.sum(matrix, axis=1)))
    if active.size == 0:
        raise MissingMetadata("source has no significant axial cells")
    lower = float(z_edges[active[0]])
    upper = float(z_edges[active[-1] + 1])
    span = max(upper - lower, 0.1)
    return lower - 0.08 * span, upper + 0.08 * span

def _active_shell_radius_max_m(surfaces: tuple[_GeometryShellSurface, ...], z_lower_m: float, z_upper_m: float) -> float:
    radii: list[float] = []
    for surface in surfaces:
        if (
            surface.geometry == "radial_cylinder"
            and surface.z_min_m is not None
            and surface.z_max_m is not None
            and surface.radius_m is not None
            and surface.z_max_m >= z_lower_m
            and surface.z_min_m <= z_upper_m
        ):
            radii.append(float(surface.radius_m))
            continue
        if (
            surface.geometry == "conical_frustum"
            and surface.z_min_m is not None
            and surface.z_max_m is not None
            and surface.radius_at_z_min_m is not None
            and surface.radius_at_z_max_m is not None
            and surface.z_max_m >= z_lower_m
            and surface.z_min_m <= z_upper_m
        ):
            local_lower = max(z_lower_m, surface.z_min_m)
            local_upper = min(z_upper_m, surface.z_max_m)
            local_radii = np.interp([local_lower, local_upper], [surface.z_min_m, surface.z_max_m], [surface.radius_at_z_min_m, surface.radius_at_z_max_m])
            radii.extend(float(value) for value in local_radii)
            continue
        if (
            surface.geometry == "axial_plane"
            and surface.z_m is not None
            and surface.radius_max_m is not None
            and z_lower_m <= surface.z_m <= z_upper_m
        ):
            radii.append(float(surface.radius_max_m))
    return max(radii, default=0.0)

def plot_better_source_cloud(metadata: Mapping[str, Any], output_path: Path, args: argparse.Namespace) -> None:
    plt = _configure_matplotlib()
    from matplotlib.lines import Line2D
    from matplotlib.ticker import MaxNLocator

    bank = _event_bank(metadata)
    sample_count = max(1, int(getattr(args, "better_cloud_samples", 30000)))
    seed = int(getattr(args, "seed", 12345))
    if bank is not None:
        chosen = _weighted_resample_indices(bank.normalized_weights, sample_count, seed)
        positions = np.asarray(bank.positions_m, dtype=float)[chosen]
        energy_MeV = np.asarray(bank.energies_J, dtype=float)[chosen] / MEV_TO_J
        total_rate = float(bank.physical_total_rate_s)
    else:
        positions, energy_MeV, total_rate = _sample_binned_events(metadata, sample_count, seed)

    if positions.shape[0] == 0:
        raise MissingMetadata("source cloud has no events")
    support_positions = np.asarray([], dtype=float).reshape(0, 3)
    if bank is not None:
        bank_positions = np.asarray(bank.positions_m, dtype=float)
        bank_axial = np.asarray(bank.axial_cell_indices, dtype=int)
        source_matrix, _, _ = _source_rate_matrix(metadata)
        axial_rate = np.sum(source_matrix, axis=1)
        significant = _significant_axial_source_mask(axial_rate)
        weak_cells = np.flatnonzero((axial_rate > 0.0) & ~significant)
        if weak_cells.size:
            support_mask = np.isin(bank_axial, weak_cells)
            support_positions = bank_positions[support_mask]

    color_mode = str(getattr(args, "better_cloud_color", "energy")).lower()
    if color_mode == "rate_density":
        matrix, source_z_edges, _ = _source_rate_matrix(metadata)
        source_z_centers = _centers_from_edges(source_z_edges)
        cell_volumes = _volumes_for_length(metadata, source_z_centers.size)
        local_density = np.divide(np.sum(matrix, axis=1), cell_volumes, out=np.zeros(source_z_centers.size, dtype=float), where=cell_volumes > 0.0)
        maximum_density = float(np.max(local_density, initial=0.0))
        if maximum_density <= 0.0:
            raise MissingMetadata("source rate density is unavailable")
        point_density = np.interp(positions[:, 2], source_z_centers, local_density, left=0.0, right=0.0)
        nominal_outer = _plasma_radius_profile(metadata, source_z_centers)
        source_outer_edges = _as_array(metadata, "magnetic_field_visual_flux_tube_radius_edges_m", ndim=1)
        if source_outer_edges.size == source_z_edges.size:
            point_outer = np.interp(positions[:, 2], source_z_edges, source_outer_edges)
        else:
            point_outer = np.interp(positions[:, 2], source_z_centers, nominal_outer, left=nominal_outer[0], right=nominal_outer[-1])
        point_radius = np.sqrt(positions[:, 0] ** 2 + positions[:, 1] ** 2)
        point_rho = np.divide(point_radius, point_outer, out=np.zeros_like(point_radius), where=point_outer > 0.0)
        rho_edges = _source_connected_rho_edges(metadata, source_z_edges.size)
        point_rho_max = np.interp(positions[:, 2], source_z_edges, rho_edges)
        radial_normalization = radial_probability_cdf(point_rho_max, _source_radial_profile_model(metadata))
        point_density *= np.divide(radial_source_shape(point_rho, _source_radial_profile_model(metadata)), radial_normalization, out=np.zeros_like(point_density), where=radial_normalization > 0.0)
        point_color = np.log10(np.maximum(point_density / maximum_density, 1.0e-12))
        color_label = r"$\log_{10}$ local neutron rate density / maximum"
    else:
        point_color = energy_MeV
        color_label = "Neutron energy [MeV]"

    view = build_model_view(metadata)
    field_z = _as_array(metadata, "magnetic_field_visual_coordinate_m", ndim=1)
    flux_radius = _as_array(metadata, "magnetic_field_visual_flux_tube_radius_m", ndim=1)
    if field_z.size < 2 or flux_radius.size != field_z.size:
        raise MissingMetadata("missing modeled magnetic flux tube profile")

    context_lower, context_upper, device_radial_limit_m = _full_device_context_limits_m(metadata, view)
    throats = np.asarray(_throat_positions(metadata), dtype=float)
    surfaces = _geometry_shell_surfaces(metadata)
    vessel_z, vessel_radius = _vessel_envelope(metadata)

    flux_keep = (np.isfinite(field_z) & np.isfinite(flux_radius) & (flux_radius >= 0.0) & (field_z >= context_lower) & (field_z <= context_upper))
    flux_z = field_z[flux_keep]
    flux_r = flux_radius[flux_keep]
    if flux_z.size < 2:
        raise MissingMetadata("modeled magnetic flux tube does not cover the full device interval")

    positions_cm = positions * 100.0
    flux_z_cm = 100.0 * flux_z
    flux_r_cm = 100.0 * flux_r
    point_radius_cm = np.sqrt(positions_cm[:, 0] ** 2 + positions_cm[:, 1] ** 2)
    radial_limit_cm = 1.08 * max(100.0 * device_radial_limit_m, float(np.max(point_radius_cm, initial=0.0)), float(np.max(flux_r_cm, initial=0.0)), 1.0e-6)
    exaggeration = max(1.0, float(getattr(args, "better_cloud_radial_exaggeration", 4.0)))

    with plt.rc_context(better_RC):
        fig = plt.figure(figsize=(11.8, 6.8), constrained_layout=False)
        ax = fig.add_subplot(111, projection="3d")

        scatter = ax.scatter(positions_cm[:, 2], positions_cm[:, 0], positions_cm[:, 1], c=point_color, cmap="viridis", s=3.0, alpha=0.50, depthshade=False, rasterized=True)
        if support_positions.size:
            support_cm = 100.0 * support_positions
            ax.scatter(support_cm[:, 2], support_cm[:, 0], support_cm[:, 1], s=18.0, facecolors="none", edgecolors=BLACK, linewidths=0.75, alpha=0.80, depthshade=False)

        if surfaces:
            _plot_geometry_shell_wireframe(ax, metadata)
            shell_legend_label = "Plasma facing shell"
        elif vessel_z.size >= 2 and vessel_radius.size == vessel_z.size:
            _wireframe(ax, 100.0 * vessel_z, 100.0 * vessel_radius, color=BLACK, ring_alpha=0.22, line_alpha=0.24, linewidth=0.72)
            shell_legend_label = "Inferred chamber envelope"
        else:
            shell_legend_label = None

        _wireframe(ax, flux_z_cm, flux_r_cm, color=GRAY, ring_alpha=0.38, line_alpha=0.34, linewidth=0.75)

        for throat in throats:
            throat_radius_cm = 100.0 * float(np.interp(throat, field_z, flux_radius))
            _plot_ring_3d(ax, 100.0 * float(throat), throat_radius_cm, color=VERMILION, linewidth=1.1, alpha=0.90, linestyle=":")

        axis_z_lower_cm = 100.0 * context_lower
        axis_z_upper_cm = 100.0 * context_upper
        ax.set_xlim(axis_z_lower_cm, axis_z_upper_cm)
        ax.set_ylim(-radial_limit_cm, radial_limit_cm)
        ax.set_zlim(-radial_limit_cm, radial_limit_cm)

        box_aspect = (float(axis_z_upper_cm - axis_z_lower_cm), 2.0 * radial_limit_cm * exaggeration, 2.0 * radial_limit_cm * exaggeration)
        try:
            ax.set_box_aspect(box_aspect, zoom=1.22)
        except TypeError:
            ax.set_box_aspect(box_aspect)

        ax.set_proj_type("ortho")
        ax.view_init(elev=float(getattr(args, "better_cloud_elevation", 18.0)), azim=float(getattr(args, "better_cloud_azimuth", -62.0)))
        ax.set_xlabel("Axial position z [cm]", labelpad=7)
        ax.set_ylabel("x [cm]", labelpad=5)
        ax.set_zlabel("y [cm]", labelpad=5)
        ax.set_title("Neutron source cloud inside the full modeled device", pad=12)
        ax.grid(True, alpha=0.14)
        ax.xaxis.set_major_locator(MaxNLocator(nbins=6))
        ax.yaxis.set_major_locator(MaxNLocator(nbins=5))
        ax.zaxis.set_major_locator(MaxNLocator(nbins=5))
        ax.tick_params(axis="x", labelsize=7.4, pad=0)
        ax.tick_params(axis="y", labelsize=7.4, pad=0)
        ax.tick_params(axis="z", labelsize=7.4, pad=0)
        for axis in (ax.xaxis, ax.yaxis, ax.zaxis):
            axis.pane.set_alpha(0.025)

        ax.text2D(0.025, 0.95, f"{positions.shape[0]:,} weighted samples\n" f"Physical rate {total_rate:.3e} n/s\n", transform=ax.transAxes, va="top", ha="left", fontsize=8.2, bbox={"facecolor": "white", "alpha": 0.88, "edgecolor": "none"})

        legend_handles = []
        if shell_legend_label is not None:
            legend_handles.append(Line2D([0], [0], color=BLACK, linewidth=1.0, label=shell_legend_label))
       
        ax.legend(handles=legend_handles, loc="lower left", bbox_to_anchor=(0.02, 0.035), fontsize=7.4)
        ax.set_position([0.03, 0.18, 0.93, 0.73])
        colorbar = fig.colorbar(scatter, ax=ax, orientation="horizontal", shrink=0.82, pad=0.025, aspect=42)
        colorbar.set_label(color_label)
        _save_better(fig, output_path)

def _histogram_density(values: np.ndarray, weights: np.ndarray, edges: np.ndarray) -> np.ndarray:
    histogram, _ = np.histogram(values, bins=edges, weights=weights)
    return np.divide(histogram, np.diff(edges), out=np.zeros_like(histogram, dtype=float), where=np.diff(edges) > 0.0)

def _beam_axis(metadata: Mapping[str, Any], beam_id: str | None) -> tuple[np.ndarray, str]:
    view = build_model_view(metadata)
    if beam_id:
        beam = view.beam(beam_id)
    elif len(view.beams) == 1:
        beam = view.beams[0]
    elif not view.beams:
        raise MissingMetadata("no enabled beam is available as a direction reference")
    else:
        raise MissingMetadata("--direction-beam-id is required when multiple beams are enabled")
    vector = np.asarray(beam.direction_unit, dtype=float)
    norm = float(np.linalg.norm(vector))
    if norm <= 0.0:
        raise MissingMetadata(f"beam direction is not defined for {beam.beam_id}")
    return vector / norm, beam.beam_id

def _direction_axis(metadata: Mapping[str, Any], mode: str, beam_id: str | None = None) -> tuple[np.ndarray, str, str]:
    if mode == "beam":
        axis, selected_beam_id = _beam_axis(metadata, beam_id)
        return axis, r"$\mu_{\mathrm{beam}}$", f"{selected_beam_id} axis"
    return np.asarray([0.0, 0.0, 1.0]), r"$\mu_z$", "machine axis"

def _weighted_mean(values: np.ndarray, weights: np.ndarray) -> float:
    weight_sum = float(np.sum(weights))
    if weight_sum <= 0.0:
        return float("nan")
    return float(np.sum(values * weights) / weight_sum)

def _reaction_spectrum(ax: Any, energy: np.ndarray, physical_weight: np.ndarray, reaction: np.ndarray, key: str, label: str, color: str, bin_count: int) -> None:
    mask = reaction == key
    if not np.any(mask):
        ax.text(0.5, 0.5, f"No {label} events", transform=ax.transAxes, ha="center", va="center")
        ax.set_title(f"{label} neutron spectrum")
        ax.set_xlabel("Neutron energy [MeV]")
        ax.set_ylabel("$dR/dE$ [s$^{-1}$ MeV$^{-1}$]")
        return
    local_energy = energy[mask]
    local_weight = physical_weight[mask]
    lower = float(np.min(local_energy))
    upper = float(np.max(local_energy))
    span = max(upper - lower, 0.02)
    edges = np.linspace(max(0.0, lower - 0.06 * span), upper + 0.06 * span, bin_count + 1)
    density = _histogram_density(local_energy, local_weight, edges)
    ax.stairs(density, edges, color=color, linewidth=1.7, fill=True, alpha=0.17)
    mean_energy = _weighted_mean(local_energy, local_weight)
    rate = float(np.sum(local_weight))
    nominal_energy = 2.45 if key == "dd_n" else 14.10
    ax.axvline(nominal_energy, color=BLACK, linestyle=":", linewidth=1.0, label="_nolegend_")
    ax.axvline(mean_energy, color=GRAY, linestyle="--", linewidth=1.0, label="Weighted mean")
    ax.text(0.03, 0.95, f"Rate {rate:.3e} s$^{{-1}}$\nWeighted mean {mean_energy:.3f} MeV", transform=ax.transAxes, va="top", ha="left", fontsize=8.2)
    ax.set_xlabel("Neutron energy [MeV]")
    ax.set_ylabel("$dR/dE$ [s$^{-1}$ MeV$^{-1}$]")
    ax.set_title(f"{label} neutron spectrum")
    ax.grid(True, alpha=0.16)
    ax.legend(loc="upper right", frameon=False, fontsize=7.2)

def plot_better_neutron_energy_direction(metadata: Mapping[str, Any], output_path: Path, args: argparse.Namespace) -> None:
    plt = _configure_matplotlib()
    from matplotlib.colors import LogNorm

    bank = _event_bank(metadata)
    if bank is None:
        load_error = metadata.get("_correlated_event_bank_load_error")
        suffix = "" if load_error is None else f": {load_error}"
        raise MissingMetadata("correlated_neutron_events.npz is required for the joint energy and direction figure" + suffix)
    weights = np.asarray(bank.normalized_weights, dtype=float)
    positive = np.isfinite(weights) & (weights > 0.0)
    if not np.any(positive):
        raise MissingMetadata("correlated event bank has no positive probability weights")
    energy = np.asarray(bank.energies_J, dtype=float)[positive] / MEV_TO_J
    direction = np.asarray(bank.directions, dtype=float)[positive]
    probability = weights[positive]
    physical_weight = probability * float(bank.physical_total_rate_s)
    reaction = np.asarray(bank.reaction_keys)[positive]
    axis_mode = str(getattr(args, "direction_reference", "machine")).lower()
    beam_id = getattr(args, "direction_beam_id", None)
    axis_vector, mu_symbol, axis_label = _direction_axis(metadata, axis_mode, beam_id)
    mu = direction @ axis_vector
    energy_bins = max(40, int(getattr(args, "better_energy_bins", 240)))
    reaction_bins = max(50, min(140, energy_bins // 2))
    dd_mask = reaction == "dd_n"
    dd_energy = energy
    dd_weight = physical_weight
    dd_reaction = reaction
    if np.any(dd_mask):
        dd_rate = float(np.sum(physical_weight[dd_mask]))
        dd_count = max(int(getattr(args, "better_neutron_dd_min_samples", 4000)), min(int(getattr(args, "better_neutron_display_samples", 40000)), int(np.count_nonzero(dd_mask))))
        local_weights = physical_weight[dd_mask]
        local_indices = np.flatnonzero(dd_mask)
        rng = np.random.default_rng(int(getattr(args, "seed", 12345)) + 2718)
        chosen = rng.choice(local_indices, size=dd_count, replace=True, p=local_weights / np.sum(local_weights))
        dd_energy = energy[chosen]
        dd_weight = np.full(dd_count, dd_rate / dd_count, dtype=float)
        dd_reaction = reaction[chosen]
    energy_lower = max(0.0, float(np.min(energy)) - 0.20)
    energy_upper = float(np.max(energy)) + 0.20
    energy_edges = np.linspace(energy_lower, energy_upper, energy_bins + 1)
    mu_edges = np.linspace(-1.0, 1.0, 101)

    with plt.rc_context(better_RC):
        fig = plt.figure(figsize=(9.8, 8.1), constrained_layout=True)
        grid = fig.add_gridspec(2, 2, height_ratios=(0.72, 1.30))
        ax_dd = fig.add_subplot(grid[0, 0])
        ax_dt = fig.add_subplot(grid[0, 1])
        ax_joint = fig.add_subplot(grid[1, :])

        _reaction_spectrum(ax_dd, dd_energy, dd_weight, dd_reaction, "dd_n", "DD", BLUE, reaction_bins)
        _reaction_spectrum(ax_dt, energy, physical_weight, reaction, "dt_n", "DT", VERMILION, reaction_bins)

        joint, e_joint_edges, mu_joint_edges = np.histogram2d(energy, mu, bins=(energy_edges, mu_edges), weights=physical_weight)
        cell_area = np.diff(energy_edges)[:, None] * np.diff(mu_edges)[None, :]
        joint = np.divide(joint, cell_area, out=np.zeros_like(joint), where=cell_area > 0.0)
        positive_joint = joint[joint > 0.0]
        if positive_joint.size == 0:
            raise MissingMetadata("joint neutron energy and direction distribution is empty")
        vmax = float(np.max(positive_joint))
        vmin = max(vmax * 1.0e-6, float(np.min(positive_joint)))
        mesh = ax_joint.pcolormesh(e_joint_edges, mu_joint_edges, joint.T, shading="auto", cmap="magma", norm=LogNorm(vmin=vmin, vmax=vmax), rasterized=True)
        ax_joint.set_xlabel("Neutron energy [MeV]")
        ax_joint.set_ylabel(f"Direction cosine {mu_symbol}")
        ax_joint.set_ylim(-1.0, 1.0)
        ax_joint.set_title("Joint laboratory energy and direction distribution " f"relative to the {axis_label}", loc="left")
        colorbar = fig.colorbar(mesh, ax=ax_joint, pad=0.012)
        colorbar.set_label(r"$d^2R/(dE\,d\mu)$ [s$^{-1}$ MeV$^{-1}$]")
        ax_joint.text(0.01, 0.02, f"{bank.event_count:,} correlated events, physical rate " f"{bank.physical_total_rate_s:.3e} n/s", transform=ax_joint.transAxes, va="bottom", ha="left", color="white", fontsize=8.0, bbox={"facecolor": "black", "alpha": 0.45, "edgecolor": "none"})

        _save_better(fig, output_path)