"""Exact sections of solids of revolution in cylindrical (r, phi, z) coordinates.

Every axisymmetric component (solid / hollow cylinder, cone frustum, nested
shell, closed can ...) is fully described by its planar (r, z) cross-section,
a :class:`Profile`.  A :class:`Section` restricts a solid to a sub-box of
(phi, z).  Meshes are produced analytically from that pair, so a cut is

* capped: the planar faces at the phi / z limits are real faces, never a see-through shell,
* crisp: normals are analytic (smooth around phi, sharp across corners), no
  vertex-normal averaging and no coincident duplicate surfaces,
* generic: nothing here knows about any particular component.
"""
from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass

import numpy as np
import pyvista as pv

TWO_PI = 2.0 * np.pi
_EPS = 1e-9


@dataclass(frozen=True)
class Section:
    """Kept sub-box of (phi, z).  ``phi`` is the kept arc in radians (span <= 2*pi)."""

    phi: tuple[float, float] = (0.0, TWO_PI)
    z: tuple[float, float] = (-np.inf, np.inf)

    @property
    def is_full(self) -> bool:
        return self.phi[1] - self.phi[0] >= TWO_PI - _EPS and not np.isfinite(self.z[0]) and not np.isfinite(self.z[1])

    @classmethod
    def without_wedge(cls, start: float, width: float, z=(-np.inf, np.inf)) -> "Section":
        """Remove the wedge ``[start, start + width]`` (radians); keep the rest."""
        return cls(phi=(start + width, start + TWO_PI), z=z)

    def key(self) -> tuple:
        return tuple(round(float(v), 6) for v in (*self.phi, *self.z))


FULL = Section()


def polylines(lines: list[np.ndarray]) -> pv.PolyData:
    """Several open polylines as one ``PolyData``."""
    if not lines:
        return pv.PolyData()
    cells, start = [], 0
    for line in lines:
        cells.extend([len(line), *range(start, start + len(line))])
        start += len(line)
    return pv.PolyData(np.vstack(lines), lines=np.asarray(cells))


def _merge_collinear(edges):
    """Join boundary edges that continue in a straight line (drops slab-split vertices)."""
    key = lambda r, z: (round(r, 7), round(z, 7))  # noqa: E731
    starts = defaultdict(list)
    for e in edges:
        starts[key(e[0], e[1])].append(e)
    used, loops = set(), []
    for e in edges:
        if id(e) in used:
            continue
        loop, cur = [], e
        while cur is not None and id(cur) not in used:
            used.add(id(cur))
            loop.append(cur)
            cur = next((n for n in starts[key(cur[2], cur[3])] if id(n) not in used), None)
        loops.append(loop)
    merged = []
    for loop in loops:
        def bends(a, b):
            ax, az, bx, bz = a[2] - a[0], a[3] - a[1], b[2] - b[0], b[3] - b[1]
            return abs(ax * bz - az * bx) > 1e-9 * max(1.0, np.hypot(ax, az) * np.hypot(bx, bz))
        k = next((i for i in range(len(loop)) if bends(loop[i - 1], loop[i])), 0)
        loop = loop[k:] + loop[:k]
        run = list(loop[0])
        for e in loop[1:]:
            if bends(tuple(run), e):
                merged.append(tuple(run))
                run = list(e)
            else:
                run[2], run[3] = e[2], e[3]
        merged.append(tuple(run))
    return merged


def _subtract(a: list[tuple[float, float]], b: list[tuple[float, float]]) -> list[tuple[float, float]]:
    """Interval set difference ``a - b``."""
    out: list[tuple[float, float]] = []
    for lo, hi in a:
        pieces = [(lo, hi)]
        for blo, bhi in b:
            nxt = []
            for plo, phi in pieces:
                if bhi <= plo + _EPS or blo >= phi - _EPS:
                    nxt.append((plo, phi))
                    continue
                if blo > plo + _EPS:
                    nxt.append((plo, blo))
                if bhi < phi - _EPS:
                    nxt.append((bhi, phi))
            pieces = nxt
        out.extend(pieces)
    return out


class Profile:
    """Planar region in the (r, z) half-plane; loops are closed polylines (outer + holes)."""

    def __init__(self, loops):
        self.loops = [np.asarray(loop, dtype=float).reshape(-1, 2) for loop in loops]
        if not self.loops or any(len(loop) < 3 for loop in self.loops):
            raise ValueError("a profile needs closed loops of at least three points")
        self._a = np.vstack(self.loops)
        self._b = np.vstack([np.roll(loop, -1, axis=0) for loop in self.loops])
        if (self._a[:, 0] < -_EPS).any():
            raise ValueError("profile radii must be non-negative")
        self.bounds = (
            float(self._a[:, 0].min()), float(self._a[:, 0].max()),
            float(self._a[:, 1].min()), float(self._a[:, 1].max()),
        )

    # ---- constructors ---------------------------------------------------
    @classmethod
    def rect(cls, r_in: float, r_out: float, z0: float, z1: float) -> "Profile":
        """Hollow (``r_in > 0``) or solid (``r_in <= 0``) cylinder."""
        r0 = max(float(r_in), 0.0)
        if r_out <= r0 or z1 <= z0:
            raise ValueError("rect needs r_out > r_in and z1 > z0")
        return cls([[(r0, z0), (r_out, z0), (r_out, z1), (r0, z1)]])

    @classmethod
    def polygon(cls, points) -> "Profile":
        return cls([points])

    @classmethod
    def shell(cls, outer: "Profile", cavity: "Profile") -> "Profile":
        """``outer`` with the region of ``cavity`` removed (cavity must lie inside)."""
        return cls([*outer.loops, *cavity.loops])

    # ---- decomposition into trapezoidal slabs ---------------------------
    def slabs(self, z_range=(-np.inf, np.inf)) -> list[tuple[float, ...]]:
        """Trapezoids ``(z0, z1, rl0, rr0, rl1, rr1)`` covering the region in ``z_range``."""
        zlo, zhi = max(z_range[0], self.bounds[2]), min(z_range[1], self.bounds[3])
        if zhi - zlo <= _EPS:
            return []
        a, b = self._a, self._b
        az, bz = a[:, 1], b[:, 1]
        lo, hi = np.minimum(az, bz), np.maximum(az, bz)
        sloped = hi - lo > _EPS
        zs = np.unique(np.clip(np.concatenate([az, [zlo, zhi]]), zlo, zhi))
        zs = zs[np.concatenate([[True], np.diff(zs) > _EPS])]

        def r_at(i: int, z: float) -> float:
            return a[i, 0] + (z - az[i]) * (b[i, 0] - a[i, 0]) / (bz[i] - az[i])

        out = []
        for zl, zh in zip(zs[:-1], zs[1:]):
            zm = 0.5 * (zl + zh)
            idx = np.nonzero(sloped & (lo <= zm) & (zm < hi))[0]
            if len(idx) % 2:
                raise ValueError("profile loops are not closed")
            order = idx[np.argsort([r_at(i, zm) for i in idx])]
            for left, right in zip(order[::2], order[1::2]):
                slab = (zl, zh, r_at(left, zl), r_at(right, zl), r_at(left, zh), r_at(right, zh))
                if slab[3] - slab[2] > _EPS or slab[5] - slab[4] > _EPS:
                    out.append(slab)
        return out

    @staticmethod
    def boundary_edges(slabs) -> list[tuple[float, float, float, float]]:
        """Boundary of the slab union as ``(r0, z0, r1, z1)`` with the region on the left."""
        edges = []
        above, below = defaultdict(list), defaultdict(list)
        for zl, zh, rl0, rr0, rl1, rr1 in slabs:
            edges.append((rr0, zl, rr1, zh))
            edges.append((rl1, zh, rl0, zl))
            above[round(zl, 9)].append((rl0, rr0))
            below[round(zh, 9)].append((rl1, rr1))
        for z in set(above) | set(below):
            for lo, hi in _subtract(above[z], below[z]):  # faces looking down
                edges.append((lo, z, hi, z))
            for lo, hi in _subtract(below[z], above[z]):  # faces looking up
                edges.append((hi, z, lo, z))
        return edges

    def outline(self, section: Section = FULL, n_theta: int = 96) -> pv.PolyData | None:
        """Edge lines (corner rings, cap outlines) to draw over :meth:`revolve`."""
        return _outline(self, section, n_theta)

    # ---- meshing ----------------------------------------------------------
    def revolve(self, section: Section = FULL, n_theta: int = 96) -> pv.PolyData | None:
        """Mesh of the solid restricted to ``section``; ``None`` if nothing remains."""
        phi0, phi1 = section.phi
        span = phi1 - phi0
        if span <= _EPS:
            return None
        full = span >= TWO_PI - _EPS
        n = max(3, int(np.ceil(n_theta * min(span, TWO_PI) / TWO_PI)))
        slabs = self.slabs(section.z)
        if not slabs:
            return None
        ang = np.linspace(phi0, phi0 + min(span, TWO_PI), n + 1)
        c, s = np.cos(ang), np.sin(ang)
        points: list[np.ndarray] = []
        normals: list[np.ndarray] = []
        tris: list[np.ndarray] = []
        count = 0

        def emit(p, nrm, t):
            nonlocal count
            points.append(p)
            normals.append(nrm)
            tris.append(t + count)
            count += len(p)

        i = np.arange(n)
        for r0, z0, r1, z1 in self.boundary_edges(slabs):
            dr, dz = r1 - r0, z1 - z0
            length = np.hypot(dr, dz)
            if length < _EPS or (r0 < _EPS and r1 < _EPS):
                continue
            nr, nz = dz / length, -dr / length  # outward normal in the (r, z) plane
            p = np.vstack([
                np.column_stack([r0 * c, r0 * s, np.full(n + 1, z0)]),
                np.column_stack([r1 * c, r1 * s, np.full(n + 1, z1)]),
            ])
            nrm = np.tile(np.column_stack([nr * c, nr * s, np.full(n + 1, nz)]), (2, 1))
            quad = []
            if r0 >= _EPS:
                quad.append(np.column_stack([i, i + 1, n + 2 + i]))
            if r1 >= _EPS:
                quad.append(np.column_stack([i, n + 2 + i, n + 1 + i]))
            emit(p, nrm, np.vstack(quad))

        if not full:
            for phic, sign in ((ang[0], -1.0), (ang[-1], 1.0)):
                cc, ss = np.cos(phic), np.sin(phic)
                tangent = sign * np.array([-ss, cc, 0.0])
                for zl, zh, rl0, rr0, rl1, rr1 in slabs:
                    rz = np.array([(rl0, zl), (rr0, zl), (rr1, zh), (rl1, zh)])
                    p = np.column_stack([rz[:, 0] * cc, rz[:, 0] * ss, rz[:, 1]])
                    order = [(0, 1, 2), (0, 2, 3)] if sign < 0 else [(0, 2, 1), (0, 3, 2)]
                    keep = []
                    for t in order:
                        u, v = rz[t[1]] - rz[t[0]], rz[t[2]] - rz[t[0]]
                        if abs(u[0] * v[1] - u[1] * v[0]) > 1e-12:
                            keep.append(t)
                    if keep:
                        emit(p, np.tile(tangent, (4, 1)), np.array(keep))
        if not tris:
            return None
        tri = np.vstack(tris)
        faces = np.column_stack([np.full(len(tri), 3), tri]).ravel()
        mesh = pv.PolyData(np.vstack(points).astype(np.float32), faces)
        # Register as the normals attribute only: ``point_data[...] = `` would also make
        # them the active scalars, and WebGL renderers then colour the surface by them.
        normals_array = pv.convert_array(np.ascontiguousarray(np.vstack(normals), dtype=np.float32), name="Normals")
        mesh.GetPointData().SetNormals(normals_array)
        return mesh


def _outline(profile: Profile, section: Section, n_theta: int) -> pv.PolyData | None:
    """Crisp edge lines of a section: corner rings plus the planar cap outlines."""
    phi0, phi1 = section.phi
    span = min(phi1 - phi0, TWO_PI)
    slabs = profile.slabs(section.z)
    if span <= _EPS or not slabs:
        return None
    edges = _merge_collinear(profile.boundary_edges(slabs))
    n = max(3, int(np.ceil(n_theta * span / TWO_PI)))
    ang = np.linspace(phi0, phi0 + span, n + 1)
    corners = {(round(e[0], 7), round(e[1], 7)) for e in edges}
    lines = [np.column_stack([r * np.cos(ang), r * np.sin(ang), np.full(n + 1, z)])
             for r, z in corners if r > _EPS]
    if span < TWO_PI - _EPS:
        for phic in (ang[0], ang[-1]):
            c, s = np.cos(phic), np.sin(phic)
            lines += [np.array([[r0 * c, r0 * s, z0], [r1 * c, r1 * s, z1]]) for r0, z0, r1, z1 in edges]
    return polylines(lines)


def section_closed_mesh(mesh: pv.PolyData, section: Section) -> pv.PolyData | None:
    """Fallback for closed, non-axisymmetric meshes (tilted ports): capped plane clips."""
    out = mesh
    z0, z1 = section.z
    if np.isfinite(z0):
        out = _clip_side(out, (0, 0, 1), (0, 0, z0), +1)
    if np.isfinite(z1) and out is not None:
        out = _clip_side(out, (0, 0, 1), (0, 0, z1), -1)
    span = section.phi[1] - section.phi[0]
    if out is not None and span < TWO_PI - _EPS:
        a, b = section.phi[1], section.phi[0] + TWO_PI  # removed wedge is (a, b)
        na = (-np.sin(a), np.cos(a), 0.0)  # positive side: angles just past a
        nb = (-np.sin(b), np.cos(b), 0.0)  # positive side: angles just past b
        origin = (0.0, 0.0, 0.0)
        if TWO_PI - span <= np.pi:
            # wedge = {past a} and {before b}; keep the disjoint complement.
            first = _clip_side(out, na, origin, -1)
            rest = _clip_side(out, na, origin, +1)
            rest = _clip_side(rest, nb, origin, +1)
            parts = [m for m in (first, rest) if m is not None]
            out = parts[0].merge(parts[1]) if len(parts) == 2 else (parts[0] if parts else None)
        else:
            out = _clip_side(out, na, origin, -1)
            out = _clip_side(out, nb, origin, +1)
    return out


def _clip_side(mesh, normal, origin, side: int):
    """Capped clip keeping the ``side`` (+1 / -1) of the plane along ``normal``."""
    if mesh is None:
        return None
    n = np.asarray(normal, float) * float(side)
    clipped = mesh.clip_closed_surface(normal=tuple(n), origin=origin)
    return clipped if clipped.n_cells else None
