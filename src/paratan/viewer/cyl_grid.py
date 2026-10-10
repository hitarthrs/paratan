"""Cylindrical (r, phi, z) bin grids and their azimuthal / axial slices.

A cylindrical mesh tally is a box in (r, phi, z).  :class:`CylGrid` describes
that box for *any* annular or solid cylinder (its ``Profile`` bounds), and cuts
planar slices out of it at an arbitrary azimuth or height:

* ``plane_phi``  - an (r, z) map at a fixed azimuth,
* ``plane_z``    - an (r, phi) map at a fixed height,
* ``volume``     - the whole grid restricted to a :class:`~revolve.Section`.

Bin values come from ``field``: either an ``(nr, nphi, nz)`` array (for example
a statepoint reshaped to the mesh) or a callable ``f(r_frac, phi, z_frac)``
evaluated at bin centres.  Planes only evaluate the bins they touch, so large
configured meshes stay cheap.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pyvista as pv

from src.paratan.viewer.revolve import FULL, TWO_PI, Profile, Section, polylines

_EPS = 1e-9


def demo_heating_field(r_frac, phi, z_frac):
    """Invented attenuation + two axial hot spots (arbitrary units); no OpenMC data."""
    return 5 + 95 * np.exp(-2.4 * r_frac) * (
        .2 + .8 * np.exp(-((z_frac - .3) / .18) ** 2) + .45 * np.exp(-((z_frac - .78) / .12) ** 2)
    ) * (.85 + .15 * np.cos(phi))


def validate_dims(dimensions) -> tuple[int, int, int]:
    if not isinstance(dimensions, (list, tuple)) or len(dimensions) != 3:
        raise ValueError('Mesh dimensions must contain three positive integers')
    if any(isinstance(n, bool) or not isinstance(n, int) or n < 1 for n in dimensions):
        raise ValueError('Mesh dimensions must contain three positive integers')
    return tuple(int(n) for n in dimensions)


_polylines = polylines


def _thin(values: np.ndarray, limit: int) -> np.ndarray:
    if len(values) <= limit:
        return values
    return values[np.unique(np.linspace(0, len(values) - 1, limit).round().astype(int))]


@dataclass(frozen=True)
class CylGrid:
    r: np.ndarray      # radial bin edges [cm]
    phi: np.ndarray    # azimuthal bin edges [rad], 0 .. 2*pi
    z: np.ndarray      # axial bin edges [cm]
    logical: tuple[int, int, int]  # configured (r, phi, z) bin counts

    @classmethod
    def from_profile(cls, profile: Profile, dimensions, *, max_bins=(256, 360, 256)) -> "CylGrid":
        """Grid spanning the cylinder that bounds ``profile``.

        Configured counts above ``max_bins`` are displayed coarser; the
        configured counts stay in ``logical``.
        """
        logical = validate_dims(dimensions)
        nr, nphi, nz = (min(n, m) for n, m in zip(logical, max_bins))
        rmin, rmax, zmin, zmax = profile.bounds
        return cls(np.linspace(rmin, rmax, nr + 1), np.linspace(0.0, TWO_PI, nphi + 1),
                   np.linspace(zmin, zmax, nz + 1), logical)

    @property
    def shape(self) -> tuple[int, int, int]:
        return len(self.r) - 1, len(self.phi) - 1, len(self.z) - 1

    # ---- bin lookup -------------------------------------------------------
    def locate_phi(self, angle):
        a = np.mod(angle, TWO_PI)
        return np.clip(np.searchsorted(self.phi, a, side='right') - 1, 0, self.shape[1] - 1)

    def locate_z(self, z):
        return np.clip(np.searchsorted(self.z, z, side='right') - 1, 0, self.shape[2] - 1)

    def phi_centre(self, j: int) -> float:
        return float(.5 * (self.phi[j] + self.phi[j + 1]))

    def z_centre(self, k: int) -> float:
        return float(.5 * (self.z[k] + self.z[k + 1]))

    def bin_text_phi(self, angle: float) -> str:
        j = int(self.locate_phi(angle))
        lo, hi = np.degrees(self.phi[j]), np.degrees(self.phi[j + 1])
        return f'azimuthal bin {j + 1}/{self.shape[1]} ({lo:.1f}°–{hi:.1f}°)'

    def bin_text_z(self, z: float) -> str:
        k = int(self.locate_z(z))
        return f'axial bin {k + 1}/{self.shape[2]} ({self.z[k]:.1f}–{self.z[k + 1]:.1f} cm)'

    # ---- sampling helpers ---------------------------------------------------
    def _values(self, field, ri, pj, zk):
        """Bin values for broadcastable index arrays."""
        nr, nphi, nz = self.shape
        if field is None:
            return None
        if callable(field):
            return field((ri + .5) / nr, (pj + .5) * TWO_PI / nphi, (zk + .5) / nz)
        return np.asarray(field)[ri, pj, zk]

    def _z_samples(self, z_range) -> np.ndarray:
        lo, hi = max(z_range[0], self.z[0]), min(z_range[1], self.z[-1])
        if hi - lo <= _EPS:
            return np.empty(0)
        zs = np.unique(np.clip(np.concatenate([self.z, [lo, hi]]), lo, hi))
        return zs[np.concatenate([[True], np.diff(zs) > _EPS])]

    def _phi_samples(self, arc, step: float = np.radians(4.0)) -> np.ndarray:
        """Angles over the kept arc: every bin edge in it plus <= ``step`` subdivisions."""
        p0 = arc[0]
        p1 = p0 + min(arc[1] - arc[0], TWO_PI)
        shifts = np.arange(np.floor(p0 / TWO_PI) - 1, np.floor(p1 / TWO_PI) + 2) * TWO_PI
        edges = (self.phi[None, :] + shifts[:, None]).ravel()
        base = np.unique(np.concatenate([[p0, p1], edges[(edges > p0 + _EPS) & (edges < p1 - _EPS)]]))
        out = [base[:1]]
        for lo, hi in zip(base[:-1], base[1:]):
            out.append(np.linspace(lo, hi, max(1, int(np.ceil((hi - lo) / step - 1e-10))) + 1)[1:])
        return np.concatenate(out)

    # ---- slices -------------------------------------------------------------
    def plane_phi(self, angle: float, field=None, section: Section = FULL, lift: float = 0.0,
                  name: str = 'value') -> pv.PolyData | None:
        """(r, z) map at azimuth ``angle`` [rad], clipped to the section's z range."""
        zs = self._z_samples(section.z)
        if len(zs) < 2:
            return None
        nr = self.shape[0]
        c, s = np.cos(angle), np.sin(angle)
        R, Z = np.meshgrid(self.r, zs, indexing='ij')
        points = np.column_stack([(R * c).ravel(order='F'), (R * s).ravel(order='F'), Z.ravel(order='F')])
        mesh = pv.StructuredGrid()
        mesh.points = points
        mesh.dimensions = (nr + 1, len(zs), 1)
        j = int(self.locate_phi(angle))
        values = self._values(field, np.arange(nr)[:, None], j, self.locate_z(.5 * (zs[:-1] + zs[1:]))[None, :])
        if values is not None:
            mesh.cell_data[name] = np.asarray(values, float).ravel(order='F')
        out = mesh.extract_surface()
        out.points = out.points + lift * np.array([-s, c, 0.0])
        return out

    def plane_z(self, z: float, field=None, section: Section = FULL, lift: float = 0.0,
                name: str = 'value') -> pv.PolyData | None:
        """(r, phi) map at height ``z`` [cm], clipped to the section's kept arc."""
        if z < self.z[0] - _EPS or z > self.z[-1] + _EPS or z < section.z[0] or z > section.z[1] + _EPS:
            return None
        phis = self._phi_samples(section.phi)
        if len(phis) < 2:
            return None
        nr = self.shape[0]
        R, P = np.meshgrid(self.r, phis, indexing='ij')
        points = np.column_stack([(R * np.cos(P)).ravel(order='F'), (R * np.sin(P)).ravel(order='F'),
                                  np.full(R.size, z + lift)])
        mesh = pv.StructuredGrid()
        mesh.points = points
        mesh.dimensions = (nr + 1, len(phis), 1)
        k = int(self.locate_z(z))
        values = self._values(field, np.arange(nr)[:, None], self.locate_phi(.5 * (phis[:-1] + phis[1:]))[None, :], k)
        if values is not None:
            mesh.cell_data[name] = np.asarray(values, float).ravel(order='F')
        out = mesh.extract_surface()
        return out

    def volume(self, field=None, section: Section = FULL, name: str = 'value') -> pv.PolyData | None:
        """Outer surface of the grid restricted to ``section`` (cut faces show the bins)."""
        zs = self._z_samples(section.z)
        phis = self._phi_samples(section.phi)
        if len(zs) < 2 or len(phis) < 2:
            return None
        # Build only the boundary, never the nr*nphi*nz interior lattice.
        # Bin values remain cell-associated; curved subdivisions share a bin.
        surfaces = []
        P, Z = np.meshgrid(phis, zs, indexing='ij')
        pj = self.locate_phi(.5 * (phis[:-1] + phis[1:]))[:, None]
        zk = self.locate_z(.5 * (zs[:-1] + zs[1:]))[None, :]
        for radius, ri in ((self.r[0], 0), (self.r[-1], self.shape[0]-1)):
            if radius <= _EPS:
                continue
            wall = pv.StructuredGrid(radius*np.cos(P), radius*np.sin(P), Z)
            values = self._values(field, ri, pj, zk)
            if values is not None:
                wall.cell_data[name] = np.asarray(values, float).ravel(order='F')
            surfaces.append(wall.extract_surface())
        for z in (zs[0], zs[-1]):
            surfaces.append(self.plane_z(z, field, section, name=name))
        if section.phi[1]-section.phi[0] < TWO_PI-_EPS:
            for phi in section.phi:
                surfaces.append(self.plane_phi(phi, field, section, name=name))
        surfaces = [m for m in surfaces if m is not None and m.n_cells]
        return pv.merge(surfaces, merge_points=False) if surfaces else None

    # ---- bin boundary lines --------------------------------------------------
    def lines_phi(self, angle: float, section: Section = FULL, lift: float = 0.0, limit=(64, 64)) -> pv.PolyData:
        zs = self._z_samples(section.z)
        if len(zs) < 2:
            return pv.PolyData()
        c, s = np.cos(angle), np.sin(angle)
        n = np.array([-s, c, 0.0]) * lift

        def pt(r, z):
            return np.column_stack([r * c, r * s, z]) + n

        lines = [pt(np.full(2, r), zs[[0, -1]]) for r in _thin(self.r, limit[0])]
        lines += [pt(self.r[[0, -1]], np.full(2, z)) for z in _thin(self.z[(self.z >= zs[0] - _EPS) & (self.z <= zs[-1] + _EPS)], limit[1])]
        return _polylines(lines)

    def lines_z(self, z: float, section: Section = FULL, lift: float = 0.0, limit=(64, 64)) -> pv.PolyData:
        if z < self.z[0] - _EPS or z > self.z[-1] + _EPS or z < section.z[0] or z > section.z[1] + _EPS:
            return pv.PolyData()
        phis = self._phi_samples(section.phi)
        edges = phis[np.isin(np.round(phis, 9), np.round(np.mod(self.phi, TWO_PI), 9)) |
                     np.isin(np.round(np.mod(phis, TWO_PI), 9), np.round(self.phi, 9))]
        height = z + lift
        lines = [np.column_stack([r * np.cos(phis), r * np.sin(phis), np.full(len(phis), height)])
                 for r in _thin(self.r, limit[0])]
        lines += [np.column_stack([self.r[[0, -1]] * np.cos(p), self.r[[0, -1]] * np.sin(p), np.full(2, height)])
                  for p in _thin(edges, limit[1])]
        return _polylines(lines)

    def cage(self, section: Section = FULL, limit=(16, 32, 24)) -> pv.PolyData:
        """Bin boundaries drawn on the cylinder's outer / inner walls and end planes."""
        zs = self._z_samples(section.z)
        if len(zs) < 2:
            return pv.PolyData()
        phis = self._phi_samples(section.phi)
        z_lines = _thin(zs, limit[2])
        lines = []
        for r in (self.r[0], self.r[-1]):
            if r < _EPS:
                continue
            for z in z_lines:
                lines.append(np.column_stack([r * np.cos(phis), r * np.sin(phis), np.full(len(phis), z)]))
        edge = _thin(phis[np.isin(np.round(np.mod(phis, TWO_PI), 9), np.round(self.phi, 9))], limit[1])
        for p in edge:
            for r in (self.r[0], self.r[-1]):
                if r >= _EPS:
                    lines.append(np.array([[r * np.cos(p), r * np.sin(p), zs[0]], [r * np.cos(p), r * np.sin(p), zs[-1]]]))
        for z in (zs[0], zs[-1]):
            for r in _thin(self.r, limit[0])[1:-1]:
                lines.append(np.column_stack([r * np.cos(phis), r * np.sin(phis), np.full(len(phis), z)]))
            for p in edge:
                lines.append(np.array([[self.r[0] * np.cos(p), self.r[0] * np.sin(p), z],
                                       [self.r[-1] * np.cos(p), self.r[-1] * np.sin(p), z]]))
        if section.phi[1] - section.phi[0] < TWO_PI - _EPS:
            cage = _polylines(lines)
            for p in (section.phi[0], section.phi[1]):
                cage = cage.merge(self.lines_phi(p, section, limit=(limit[0], limit[2])))
            return cage
        return _polylines(lines)
