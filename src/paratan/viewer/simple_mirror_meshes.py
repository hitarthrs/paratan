"""Build analytic component meshes for the Paratan simple mirror from YAML."""
from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np
import pyvista as pv
import yaml

from src.paratan.viewer.components import ComponentMesh, color_for_material, component_label
from src.paratan.viewer.revolve import Profile

# Keep in sync with src.paratan.geometry.core.SIMPLE_MIRROR_VV_STYLES (avoid OpenMC import).
_SIMPLE_MIRROR_VV_STYLES = frozenset({"conical", "perpendicular"})


def load_simple_mirror_yaml(path: str | Path) -> dict[str, Any]:
    with Path(path).open() as f:
        return yaml.safe_load(f)


def _vv_style(vv: dict[str, Any]) -> str:
    style = str(vv.get("geometry_style", "conical") or "conical").lower()
    if style not in _SIMPLE_MIRROR_VV_STYLES:
        raise ValueError(
            "vacuum_vessel.geometry_style must be one of "
            f"{sorted(_SIMPLE_MIRROR_VV_STYLES)}, got {vv.get('geometry_style')!r}"
        )
    return style


def _vv_z_extents(vv: dict[str, Any]) -> dict[str, float]:
    mid = float(vv.get("axial_midplane", 0.0))
    first = float(vv["outer_axial_length"]) / 2.0
    second = float(vv["central_axial_length"]) / 2.0
    left_bn = float(vv["left_bottleneck_length"])
    right_bn = float(vv["right_bottleneck_length"])
    return {
        "mid": mid,
        "central_left": mid - first,
        "central_right": mid + first,
        "cone_left_outer": mid - second,  # bottle/cone junction (left)
        "cone_right_outer": mid + second,
        "left_end": mid - second - left_bn,
        "right_end": mid + second + right_bn,
    }


def _perpendicular_z_extents(vv: dict[str, Any], thickness: float = 0.0) -> dict[str, float]:
    """Stepped VV extents; structural thickness grows the central cylinder axially.

    Matches :func:`src.paratan.models.simple_mirror_volumes._profile` and the
    ports builder's perpendicular structural nesting (outermost ends fixed).
    """
    mid = float(vv.get("axial_midplane", 0.0))
    half = float(vv["central_axial_length"]) / 2.0
    left_bn = float(vv["left_bottleneck_length"])
    right_bn = float(vv["right_bottleneck_length"])
    h = half + float(thickness)
    left_end = mid - half - left_bn
    right_end = mid + half + right_bn
    if mid - h <= left_end or mid + h >= right_end:
        raise ValueError("Structural thickness exhausts a vessel bottleneck.")
    return {
        "mid": mid,
        "central_left": mid - h,
        "central_right": mid + h,
        "left_end": left_end,
        "right_end": right_end,
    }


def _vessel_outline(r_central: float, r_bottle: float, z: dict[str, float]) -> list[tuple[float, float]]:
    """Left-to-right (r, z) outline of the conical vessel wall, bottle ends included."""
    return [
        (r_bottle, z["left_end"]),
        (r_bottle, z["cone_left_outer"]),
        (r_central, z["central_left"]),
        (r_central, z["central_right"]),
        (r_bottle, z["cone_right_outer"]),
        (r_bottle, z["right_end"]),
    ]


def _perpendicular_vessel_outline(
    r_central: float, r_bottle: float, z: dict[str, float]
) -> list[tuple[float, float]]:
    """Left-to-right (r, z) outline of the stepped (perpendicular) vessel wall."""
    return [
        (r_bottle, z["left_end"]),
        (r_bottle, z["central_left"]),
        (r_central, z["central_left"]),
        (r_central, z["central_right"]),
        (r_bottle, z["central_right"]),
        (r_bottle, z["right_end"]),
    ]


def _conical_vessel_solid(r_central: float, r_bottle: float, z: dict[str, float]) -> Profile:
    """Filled conical VV (plasma / vessel interior) matching single_vacuum_vessel_region."""
    wall = _vessel_outline(r_central, r_bottle, z)
    return Profile.polygon([(0.0, wall[0][1]), *wall, (0.0, wall[-1][1])])


def _perpendicular_vessel_solid(
    r_central: float, r_bottle: float, z: dict[str, float]
) -> Profile:
    """Filled stepped VV matching perpendicular_vacuum_vessel_region."""
    wall = _perpendicular_vessel_outline(r_central, r_bottle, z)
    return Profile.polygon([(0.0, wall[0][1]), *wall, (0.0, wall[-1][1])])


def _conical_vessel_shell(
    r_c_in: float, r_b_in: float, r_c_out: float, r_b_out: float, z: dict[str, float]
) -> Profile:
    """Structural shell between two conical VV surfaces."""
    outer = _vessel_outline(r_c_out, r_b_out, z)
    inner = _vessel_outline(r_c_in, r_b_in, z)
    return Profile.polygon([*outer, *inner[::-1]])


def _perpendicular_vessel_shell(
    r_c_in: float,
    r_b_in: float,
    z_in: dict[str, float],
    r_c_out: float,
    r_b_out: float,
    z_out: dict[str, float],
) -> Profile:
    """Structural shell between two stepped VV surfaces (possibly different step z)."""
    outer = _perpendicular_vessel_outline(r_c_out, r_b_out, z_out)
    inner = _perpendicular_vessel_outline(r_c_in, r_b_in, z_in)
    return Profile.polygon([*outer, *inner[::-1]])


def _surround_vessel_profile(r_out: float, z0: float, z1: float, wall) -> Profile:
    """Cylindrical blanket outside the vessel, including its tapered ends.

    ``wall`` must be a left-to-right outline with strictly increasing ``z``
    (conical). Stepped / perpendicular vessels use
    :func:`_surround_stepped_vessel_profile` instead.
    """
    radii, heights = np.asarray(wall).T
    if np.any(np.diff(heights) <= 0):
        raise ValueError("conical surround wall must have strictly increasing z")
    inner = []
    if z0 < heights[0]:
        inner.extend([(0.0, z0), (0.0, min(z1, heights[0]))])
    lo, hi = max(z0, heights[0]), min(z1, heights[-1])
    if lo < hi:
        zs = sorted({lo, hi, *[z for z in heights if lo < z < hi]})
        inner.extend((float(np.interp(z, heights, radii)), z) for z in zs)
    if z1 > heights[-1]:
        inner.extend([(0.0, max(z0, heights[-1])), (0.0, z1)])
    return Profile.polygon([(r_out, z0), (r_out, z1), *inner[::-1]])


def _surround_stepped_vessel_profile(
    r_out: float,
    z0: float,
    z1: float,
    wall: list[tuple[float, float]],
) -> Profile:
    """Cylindrical blanket outside a stepped (perpendicular) vessel wall.

    Walks the outline so vertical step faces become radial bore edges; regions
    beyond the vessel ends are solid out to the axis (same convention as the
    conical surround helper).
    """
    if z1 <= z0:
        raise ValueError("surround needs z1 > z0")
    w0, w1 = float(wall[0][1]), float(wall[-1][1])

    def clip_seg(
        r_a: float, z_a: float, r_b: float, z_b: float
    ) -> list[tuple[float, float]]:
        """Clip one outline segment to ``[z0, z1]`` (inclusive)."""
        if abs(z_b - z_a) <= 1e-12:
            # Vertical (radial) step: keep if the step plane is in-range.
            if z0 - 1e-12 <= z_a <= z1 + 1e-12:
                return [(r_a, z_a), (r_b, z_b)]
            return []
        # Parametric clip in z.
        za, zb, ra, rb = z_a, z_b, r_a, r_b
        if za > zb:
            za, zb, ra, rb = zb, za, rb, ra
            flipped = True
        else:
            flipped = False
        if zb < z0 - 1e-12 or za > z1 + 1e-12:
            return []
        t0 = 0.0 if za >= z0 else (z0 - za) / (zb - za)
        t1 = 1.0 if zb <= z1 else (z1 - za) / (zb - za)
        p0 = (ra + t0 * (rb - ra), za + t0 * (zb - za))
        p1 = (ra + t1 * (rb - ra), za + t1 * (zb - za))
        return [p1, p0] if flipped else [p0, p1]

    inner: list[tuple[float, float]] = []
    if z0 < w0 - 1e-12:
        inner.extend([(0.0, z0), (0.0, min(z1, w0))])
    for (r_a, z_a), (r_b, z_b) in zip(wall[:-1], wall[1:]):
        for pt in clip_seg(float(r_a), float(z_a), float(r_b), float(z_b)):
            if not inner or abs(inner[-1][0] - pt[0]) > 1e-9 or abs(inner[-1][1] - pt[1]) > 1e-9:
                inner.append(pt)
    if z1 > w1 + 1e-12:
        start = (0.0, max(z0, w1))
        if not inner or abs(inner[-1][0] - start[0]) > 1e-9 or abs(inner[-1][1] - start[1]) > 1e-9:
            inner.append(start)
        inner.append((0.0, z1))
    if len(inner) < 2:
        return Profile.rect(0.0, r_out, z0, z1)
    return Profile.polygon([(r_out, z0), (r_out, z1), *inner[::-1]])


def _add(
    components: list[ComponentMesh],
    *,
    name: str,
    group: str,
    material: str,
    meta: dict[str, Any],
    n_theta: int = 96,
    profile: Profile | None = None,
    mesh: pv.DataSet | None = None,
    opacity: float = 1.0,
) -> None:
    if mesh is None:
        mesh = profile.revolve(n_theta=n_theta)
    components.append(
        ComponentMesh(
            name=name,
            group=group,
            material=material,
            mesh=mesh,
            meta=meta,
            color=color_for_material(material),
            opacity=opacity,
            profile=profile,
            label=component_label(name, group, material, meta),
        )
    )


def build_simple_mirror_components(
    input_data: dict[str, Any],
    *,
    n_theta: int = 96,
) -> list[ComponentMesh]:
    """YAML dict → list of ComponentMesh (no OpenMC build required)."""
    components: list[ComponentMesh] = []
    vv = input_data["vacuum_vessel"]
    style = _vv_style(vv)
    r_c0 = float(vv["central_radius"])
    r_b0 = float(vv["bottleneck_radius"])
    conical = style == "conical"
    z = _vv_z_extents(vv) if conical else _perpendicular_z_extents(vv)

    if conical:
        profile = _conical_vessel_solid(r_c0, r_b0, z)
    else:
        profile = _perpendicular_vessel_solid(r_c0, r_b0, z)

    _add(
        components,
        name="vacuum_vessel_plasma",
        group="VV",
        material="vacuum",
        profile=profile,
        n_theta=n_theta,
        meta={
            "central_radius_cm": r_c0,
            "bottleneck_radius_cm": r_b0,
            "geometry_style": style,
        },
        opacity=0.25,
    )

    structure = vv.get("structure") or {}
    structural_thicknesses = [float(p["thickness"]) for p in structure.values()]
    structural_names = list(structure.keys())
    structural_materials = [str(p["material"]) for p in structure.values()]
    central_radii = np.cumsum([r_c0] + structural_thicknesses)
    bottle_radii = np.cumsum([r_b0] + structural_thicknesses)
    # Cumulative structural thickness for perpendicular axial growth of the step.
    cum_thickness = np.cumsum([0.0] + structural_thicknesses)

    for i, layer_name in enumerate(structural_names):
        group = "FW" if "first_wall" in layer_name.lower() or layer_name == "first_wall" else "VV"
        if conical:
            shell = _conical_vessel_shell(
                float(central_radii[i]),
                float(bottle_radii[i]),
                float(central_radii[i + 1]),
                float(bottle_radii[i + 1]),
                z,
            )
        else:
            shell = _perpendicular_vessel_shell(
                float(central_radii[i]),
                float(bottle_radii[i]),
                _perpendicular_z_extents(vv, float(cum_thickness[i])),
                float(central_radii[i + 1]),
                float(bottle_radii[i + 1]),
                _perpendicular_z_extents(vv, float(cum_thickness[i + 1])),
            )
        _add(
            components,
            name=f"vv_{layer_name}",
            group=group,
            material=structural_materials[i],
            profile=shell,
            n_theta=n_theta,
            meta={"thickness_cm": structural_thicknesses[i], "layer_index": i},
        )

    vv_outer_r = float(central_radii[-1]) if len(central_radii) else r_c0
    vv_outer_bottle = float(bottle_radii[-1]) if len(bottle_radii) else r_b0
    outer_wall = (
        _vessel_outline(vv_outer_r, vv_outer_bottle, z)
        if conical
        else _perpendicular_vessel_outline(
            vv_outer_r,
            vv_outer_bottle,
            _perpendicular_z_extents(vv, float(cum_thickness[-1])),
        )
    )
    surround = _surround_vessel_profile if conical else _surround_stepped_vessel_profile

    # --- Central cell radial layers (outside VV) ---
    cc = input_data["central_cell"]
    cc_half = float(cc["axial_length"]) / 2.0
    r_in = vv_outer_r
    for i, layer in enumerate(cc["layers"]):
        thick = float(layer["thickness"])
        r_out = r_in + thick
        mat = str(layer["material"])
        profile = (
            surround(r_out, -cc_half, cc_half, outer_wall)
            if i == 0
            else Profile.rect(r_in, r_out, -cc_half, cc_half)
        )
        _add(
            components,
            name=f"cc_layer_{i}_{mat}",
            group="CC",
            material=mat,
            profile=profile,
            n_theta=n_theta,
            meta={
                "thickness_cm": thick,
                "layer_index": i,
                "r_inner_cm": r_in,
                "r_outer_cm": r_out,
                "axial_length_cm": float(cc["axial_length"]),
            },
        )
        r_in = r_out
    cc_outer = r_in

    # --- LF coils ---
    lf = input_data["lf_coil"]
    shell_t = lf["shell_thicknesses"]
    inner = lf["inner_dimensions"]
    for i, z0 in enumerate(lf["positions"]):
        z0 = float(z0)
        r0 = cc_outer
        # magnet pack (inner annular region)
        r_front = r0 + float(shell_t["front"])
        r_mag_out = r_front + float(inner["radial_thickness"])
        z_half = float(inner["axial_length"]) / 2.0
        mag = Profile.rect(r_front, r_mag_out, z0 - z_half, z0 + z_half)
        _add(
            components,
            name=f"lf_coil_{i}_magnet",
            group="LF",
            material=str(lf["materials"]["magnet"]),
            profile=mag,
            n_theta=n_theta,
            meta={"z0_cm": z0, "radial_thickness_cm": float(inner["radial_thickness"])},
        )
        # outer shield box (approximate as larger annulus covering shell)
        r_out = r_mag_out + float(shell_t["back"])
        z_shell = z_half + float(shell_t["axial"])
        shield = Profile.rect(r0, r_out, z0 - z_shell, z0 + z_shell)
        # subtract magnet volume visually by using shell-only: outer annulus + ends
        # Keep full outer for simplicity; opacity distinguishes.
        _add(
            components,
            name=f"lf_coil_{i}_shield",
            group="LF",
            material=str(lf["materials"]["shield"]),
            profile=shield,
            n_theta=n_theta,
            meta={"z0_cm": z0},
            opacity=0.55,
        )

    # --- HF coils ---
    hf = input_data["hf_coil"]
    casing_thicks = np.array([float(L["thickness"]) for L in hf["casing_layers"]])
    cc_half_for_hf = float(cc["axial_length"]) / 2.0
    hf_z0_offset = 0.0
    if float(vv["central_axial_length"]) >= float(cc["axial_length"]):
        hf_z0_offset = float(vv["central_axial_length"]) / 2.0 - float(cc["axial_length"]) / 2.0
    hf_coil_center_z0 = (
        cc_half_for_hf
        + float(hf["shield"]["shield_central_cell_gap"])
        + float(hf["shield"]["axial_thickness"][0])
        + float(np.sum(casing_thicks))
        + float(hf["magnet"]["axial_thickness"]) / 2.0
        + hf_z0_offset
    )
    r_hf0 = vv_outer_bottle + float(hf["shield"]["radial_gap_before_casing"])
    mag_dr = float(hf["magnet"]["radial_thickness"])
    mag_dz = float(hf["magnet"]["axial_thickness"])
    # Match nested_cylindrical_shells: the magnet starts outside all inward
    # casing/shield layers; each surrounding layer expands radially and axially.
    front = list(casing_thicks) + [float(hf["shield"]["radial_thickness"][0])]
    back = list(casing_thicks) + [float(hf["shield"]["radial_thickness"][1])]
    axial = list(casing_thicks) + [float(hf["shield"]["axial_thickness"][0])]
    for side, sign in (("right", 1.0), ("left", -1.0)):
        zc = sign * hf_coil_center_z0
        ri = r_hf0 + sum(front)
        ro = ri + mag_dr
        zlo, zhi = zc - mag_dz / 2.0, zc + mag_dz / 2.0
        _add(
            components, name=f"hf_coil_{side}_magnet", group="HF",
            material=str(hf["magnet"]["material"]),
            profile=Profile.rect(ri, ro, zlo, zhi), n_theta=n_theta,
            meta={"z0_cm": zc, "bore_start_r_cm": ri},
        )
        layers = [(f"casing_{j}", str(layer["material"]))
                  for j, layer in enumerate(hf["casing_layers"])]
        layers.append(("shield", str(hf["shield"]["material"])))
        for j, (label, material) in enumerate(layers):
            new_ri, new_ro = ri - front[j], ro + back[j]
            new_zlo, new_zhi = zlo - axial[j], zhi + axial[j]
            # Expanded annulus minus the previous one: a ring in the (r, z) plane.
            ring = Profile.shell(Profile.rect(new_ri, new_ro, new_zlo, new_zhi),
                                 Profile.rect(ri, ro, zlo, zhi))
            _add(
                components, name=f"hf_coil_{side}_{label}", group="HF",
                material=material, profile=ring, n_theta=n_theta,
                meta={"inner_radius_cm": new_ri, "outer_radius_cm": new_ro,
                      "z_min_cm": new_zlo, "z_max_cm": new_zhi},
                opacity=0.7 if label.startswith("casing") else 1.0,
            )
            ri, ro, zlo, zhi = new_ri, new_ro, new_zlo, new_zhi

    # --- End cells ---
    end = input_data.get("end_cell") or {}
    if end:
        casing_thicks = np.array(
            [float(L["thickness"]) for L in hf["casing_layers"]]
        )
        hf_z0_offset_end = hf_z0_offset + 10.0
        hf_coil_center_z0_end = (
            cc_half_for_hf
            + float(hf["shield"]["shield_central_cell_gap"])
            + float(hf["shield"]["axial_thickness"][0])
            + float(np.sum(casing_thicks))
            + float(hf["magnet"]["axial_thickness"]) / 2.0
            + hf_z0_offset_end
        )
        end_z0 = (
            hf_coil_center_z0_end
            + float(hf["shield"]["axial_thickness"][0])
            + float(np.sum(casing_thicks))
            + float(hf["magnet"]["axial_thickness"]) / 2.0
            + float(end.get("shell_thickness", 2))
            + float(end.get("axial_length", 50)) / 2.0
        )
        r_end = float(end.get("diameter", 100)) / 2.0
        L = float(end.get("axial_length", 50))
        t = float(end.get("shell_thickness", 2))
        for side, sign in (("right", 1.0), ("left", -1.0)):
            zc = sign * end_z0
            z_in0, z_in1 = zc - L / 2.0, zc + L / 2.0
            # OpenMC cylinder_with_shell: closed can (barrel + full end disks).
            # Previous annular-only mesh left a see-through hole at each circular face.
            ro, zo0, zo1 = r_end + t, z_in0 - t, z_in1 + t
            shell = Profile.polygon([
                (0.0, zo0), (ro, zo0), (ro, zo1), (0.0, zo1),
                (0.0, z_in1), (r_end, z_in1), (r_end, z_in0), (0.0, z_in0),
            ])
            inner = Profile.rect(0.0, r_end, z_in0, z_in1)
            _add(
                components,
                name=f"end_cell_{side}_shell",
                group="ends",
                material=str(end.get("shell_material", "stainless")),
                profile=shell,
                n_theta=n_theta,
                meta={"shell_thickness_cm": t, "z0_cm": zc, "axial_length_cm": L},
            )
            _add(
                components,
                name=f"end_cell_{side}_inner",
                group="ends",
                material=str(end.get("inner_material", "vacuum")),
                profile=inner,
                n_theta=n_theta,
                meta={"z0_cm": zc, "axial_length_cm": L},
                opacity=0.15,
            )

    # --- Ports: two exterior arms meeting the plasma boundary ---
    ports = input_data.get('ports') or {}
    if ports.get('ports_present', False):
        from src.paratan.viewer.port_openings import PortOpenings
        openings = PortOpenings(input_data)
        inner_radius = float(ports.get('inner_radius',25))
        x_limit = float(ports.get('port_x_limits',250))
        radii = [inner_radius, *(inner_radius + np.cumsum(ports.get('port_layers_thicknesses',[])))]
        materials = ['vacuum', *ports.get('port_layers_materials',[])]
        if len(radii) != len(materials):
            raise ValueError('Each port wall thickness must have a material')
        for side, angle in (('left',45.),('right',-45.)):
            for i,(radius,material) in enumerate(zip(radii,materials)):
                ri = 0 if i==0 else radii[i-1]
                tube = openings.arms(ri,radius,angle,max(64,n_theta//2))
                _add(components,name=f'port_{side}_core' if i==0 else f'port_{side}_wall_{i}',
                     group='ports',material=str(material),mesh=tube,
                     meta={'tilt_deg':angle,'r_inner_cm':ri,'r_outer_cm':radius,'x_limits_cm':[-x_limit,x_limit], 'plasma_clipped': True},
                     # Vacuum is an inspectable simulation cell, not a visible
                     # plug: its end disks otherwise cover the hollow bore.
                     opacity=0. if i==0 else 1.)

    return components


def components_to_multiblock(
    components: list[ComponentMesh],
) -> pv.MultiBlock:
    mb = pv.MultiBlock()
    for c in components:
        # Store pick metadata as field data
        mesh = c.mesh.copy(deep=True)
        mesh.field_data["component_name"] = np.array([c.name])
        mesh.field_data["component_group"] = np.array([c.group])
        mesh.field_data["component_material"] = np.array([c.material])
        mb.append(mesh, c.name)
    return mb


def build_from_yaml(
    path: str | Path,
    *,
    n_theta: int = 96,
) -> list[ComponentMesh]:
    return build_simple_mirror_components(load_simple_mirror_yaml(path), n_theta=n_theta)
