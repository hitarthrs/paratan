"""Analytic PyVista primitives for revolution / annular solids (cm)."""
from __future__ import annotations

import numpy as np
import pyvista as pv


def _with_normals(mesh: pv.PolyData, *, feature_angle: float = 30.0) -> pv.PolyData:
    """Smooth curved walls; keep sharp splits at flat caps / cut faces."""
    return mesh.compute_normals(
        cell_normals=False,
        point_normals=True,
        split_vertices=True,
        feature_angle=feature_angle,
        auto_orient_normals=True,
        inplace=False,
    )


def annular_cylinder(
    r_inner: float,
    r_outer: float,
    z_min: float,
    z_max: float,
    *,
    n_theta: int = 64,
    n_radial: int = 2,
) -> pv.PolyData:
    """Hollow (or solid if r_inner<=0) cylinder aligned with +z."""
    if r_outer <= 0:
        raise ValueError("r_outer must be positive")
    if z_max <= z_min:
        raise ValueError("z_max must exceed z_min")
    r0 = max(r_inner, 0.0)
    if r0 >= r_outer:
        raise ValueError("r_inner must be < r_outer")
    # Constant-radius frustum shell == annular cylinder
    return frustum_shell(r0, r_outer, r0, r_outer, z_min, z_max, n_theta=n_theta)


def solid_disk_slab(
    radius: float,
    z_min: float,
    z_max: float,
    *,
    n_theta: int = 64,
    n_radial: int = 16,
) -> pv.PolyData:
    """Thick disk with radial×theta quads — no triangle-fan ridge shading on faces."""
    if radius <= 0:
        raise ValueError("radius must be positive")
    if z_max <= z_min:
        raise ValueError("z_max must exceed z_min")
    disc = pv.Disc(
        center=(0.0, 0.0, float(z_min)),
        inner=0.0,
        outer=float(radius),
        r_res=max(4, int(n_radial)),
        c_res=max(16, int(n_theta)),
        normal=(0.0, 0.0, 1.0),
    )
    slab = disc.extrude((0.0, 0.0, float(z_max - z_min)), capping=True)
    return _with_normals(slab.triangulate().clean())


def solid_cylinder(
    radius: float,
    z_min: float,
    z_max: float,
    *,
    n_theta: int = 64,
) -> pv.PolyData:
    """Solid cylinder: smooth side wall + ridgeless end caps."""
    height = z_max - z_min
    # Side only (no caps) — PyVista capped Cylinder uses a triangle fan that ridges.
    side = pv.Cylinder(
        center=(0.0, 0.0, 0.5 * (z_min + z_max)),
        direction=(0.0, 0.0, 1.0),
        radius=radius,
        height=height,
        resolution=n_theta,
        capping=False,
    )
    # Thin numerical eps so caps sit flush without z-fight.
    eps = max(1e-4, 1e-6 * height)
    bot = solid_disk_slab(radius, z_min - eps, z_min + eps, n_theta=n_theta)
    top = solid_disk_slab(radius, z_max - eps, z_max + eps, n_theta=n_theta)
    return merge_polydata([side.triangulate(), bot, top])


def frustum_shell(
    r_inner_bot: float,
    r_outer_bot: float,
    r_inner_top: float,
    r_outer_top: float,
    z_bot: float,
    z_top: float,
    *,
    n_theta: int = 64,
) -> pv.PolyData:
    """Annular truncated cone (shell) between two z planes."""
    theta = np.linspace(0.0, 2.0 * np.pi, n_theta, endpoint=False)
    ct, st = np.cos(theta), np.sin(theta)

    def ring(r: float, z: float) -> np.ndarray:
        return np.column_stack([r * ct, r * st, np.full_like(ct, z)])

    # Outer surface quads + inner surface quads + end annuli
    pts: list[np.ndarray] = []
    faces: list[list[int]] = []

    def add_quad_strip(a: np.ndarray, b: np.ndarray) -> None:
        base = len(pts)
        pts.extend(a)
        pts.extend(b)
        n = len(a)
        for i in range(n):
            j = (i + 1) % n
            faces.append([4, base + i, base + j, base + n + j, base + n + i])

    add_quad_strip(ring(r_outer_bot, z_bot), ring(r_outer_top, z_top))
    if r_inner_bot > 0 or r_inner_top > 0:
        add_quad_strip(ring(r_inner_top, z_top), ring(r_inner_bot, z_bot))  # inward winding
        # bottom annulus
        add_quad_strip(ring(r_inner_bot, z_bot), ring(r_outer_bot, z_bot))
        # top annulus
        add_quad_strip(ring(r_outer_top, z_top), ring(r_inner_top, z_top))
    else:
        # Solid frustum: outer wall + ridgeless disc caps (no triangle-fan).
        wall = _with_normals(
            pv.PolyData(np.vstack(pts), np.hstack(faces)).clean().triangulate()
        )
        dz = max(1e-3, 1e-6 * abs(z_top - z_bot))
        bot = solid_disk_slab(
            float(r_outer_bot), z_bot - dz, z_bot, n_theta=n_theta
        )
        top = solid_disk_slab(
            float(r_outer_top), z_top, z_top + dz, n_theta=n_theta
        )
        return merge_polydata([wall, bot, top])

    points = np.vstack(pts)
    faces_arr = np.hstack(faces)
    return _with_normals(pv.PolyData(points, faces_arr).clean().triangulate())


def solid_frustum(
    r_bot: float,
    r_top: float,
    z_bot: float,
    z_top: float,
    *,
    n_theta: int = 64,
) -> pv.PolyData:
    return frustum_shell(0.0, r_bot, 0.0, r_top, z_bot, z_top, n_theta=n_theta)


def closed_cylindrical_shell(
    r_inner: float,
    r_outer: float,
    z_inner_min: float,
    z_inner_max: float,
    shell_thickness: float,
    *,
    n_theta: int = 64,
) -> pv.PolyData:
    """Hollow cylinder with full end disks — matches ``cylinder_with_shell``.

    Barrel spans ``[z_inner_min, z_inner_max]``; end caps of axial thickness
    ``shell_thickness`` close both ends out to ``r_outer`` (no see-through hole).
    """
    if r_outer <= r_inner:
        raise ValueError("r_outer must exceed r_inner")
    if z_inner_max <= z_inner_min:
        raise ValueError("z_inner_max must exceed z_inner_min")
    if shell_thickness <= 0:
        raise ValueError("shell_thickness must be positive")

    z_outer_min = z_inner_min - shell_thickness
    z_outer_max = z_inner_max + shell_thickness
    parts = [
        # Radial wall over the inner cavity height
        annular_cylinder(r_inner, r_outer, z_inner_min, z_inner_max, n_theta=n_theta),
        # Full end plates (r = 0 → r_outer) — disc slabs, not fan-capped cylinders
        solid_disk_slab(r_outer, z_outer_min, z_inner_min, n_theta=n_theta),
        solid_disk_slab(r_outer, z_inner_max, z_outer_max, n_theta=n_theta),
    ]
    return merge_polydata(parts)


def finish_clipped_mesh(mesh: pv.PolyData) -> pv.PolyData:
    """Clean clip output so the planar cut face doesn't shade as ridges."""
    if mesh.n_cells == 0:
        return mesh
    cleaned = mesh.triangulate().clean()
    return _with_normals(cleaned, feature_angle=20.0)


def merge_polydata(parts: list[pv.PolyData]) -> pv.PolyData:
    if not parts:
        raise ValueError("no parts to merge")
    out = parts[0]
    for p in parts[1:]:
        out = out.merge(p)
    return _with_normals(out.clean().triangulate())
