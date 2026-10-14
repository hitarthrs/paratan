"""Cylindrical tally grids, overlays, inspection helpers, and result datasets.

Merged from former cyl_grid / tally_overlay / inspection / results modules.
Arrays use (r, phi, z) order. Native heating is eV/source.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from uuid import uuid4
import json

import numpy as np
import pyvista as pv
import yaml

from src.paratan.viewer.revolve import FULL, TWO_PI, Profile, Section, polylines


# --- cylindrical grids ---

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


# --- overlay presentation ---


def lift_boundary(mesh, grid, section, distance=.12):
    """Lift a display surface into free space without changing bins or values.

    Radial/end boundaries move outwards; cut faces move into the removed wedge.
    The small offset avoids depth fighting with the underlying component.
    """
    out=mesh.copy(deep=True)
    points=out.points.copy()
    radius=np.linalg.norm(points[:,:2],axis=1)
    angle=np.mod(np.arctan2(points[:,1],points[:,0]),2*np.pi)
    delta=np.where(np.isclose(radius,grid.r[-1]),distance,
                   np.where(np.isclose(radius,grid.r[0]),-distance,0.))
    points[:,:2]*=(1+delta/np.maximum(radius,1e-9))[:,None]
    points[:,2]+=np.where(np.isclose(points[:,2],grid.z[-1]),distance,
                         np.where(np.isclose(points[:,2],grid.z[0]),-distance,0.))
    if section.phi[1]-section.phi[0] < 2*np.pi-1e-9:
        for phi,sign in ((section.phi[0],-1),(section.phi[1],1)):
            difference=np.arctan2(np.sin(angle-phi),np.cos(angle-phi))
            face=np.isclose(difference,0,atol=1e-7)&(radius>1e-9)
            points[face,:2]+=sign*distance*np.array([-np.sin(phi),np.cos(phi)])
    out.points=points
    return out


def describe_overlay(grid, dataset=None, quantity='mean'):
    summary = {
        'active': True,
        'source': 'Mesh configuration',
        'bins': ' × '.join(map(str, grid.logical)) + ' bins (r · φ · z)',
        'extent': f'r {grid.r[0]:g}–{grid.r[-1]:g} cm · z {grid.z[0]:g}–{grid.z[-1]:g} cm',
        'range': '',
    }
    if dataset is not None:
        values = dataset.values(quantity)
        finite = values[np.isfinite(values)]
        summary['source'] = 'DEMO · invented values' if dataset.source == 'synthetic' else 'OpenMC results'
        summary['quantity'] = dataset.scalar_name(quantity)
        if finite.size:
            summary['range'] = f'{finite.min():.3g} – {finite.max():.3g}'
    return summary


# --- inspection / previews ---
def tally_configuration(data, component):
    """Return only entries the simple-model builder attaches to this component."""
    if '_tally_entries' in component.meta:
        return component.meta.get('_tally_path', ''), component.meta['_tally_entries']
    name = component.name
    blocks = []
    path = "No tally binding for this component in the simple-model builder."
    if component.group == 'CC':
        index = int(component.meta['layer_index'])
        config = data['central_cell'].get('tallies') or {}
        path = 'central_cell.tallies'
        if index == 0 and config.get('breeder'):
            blocks.append(config['breeder'])
        blocks.extend(e for e in config.get('layer_tallies', []) if e.get('position') == index)
    elif component.group in ('HF', 'LF') and name.endswith('_magnet'):
        section = 'hf_coil' if component.group == 'HF' else 'lf_coil'
        key = section + '_tallies'
        path = f'{section}.{key}'
        blocks.append(data[section].get(key) or {})
    elif component.group == 'ends' and name.endswith('_shell'):
        path = 'end_cell.end_cell_tallies'
        blocks.append(data['end_cell'].get('end_cell_tallies') or {})
    entries = []
    for block in blocks:
        for kind in ('cell_tallies', 'mesh_tallies'):
            for entry in block.get(kind) or []:
                entries.append(dict(entry, kind=kind))
    return path, entries


def tally_text(path, entries):
    lines = ['CONFIGURATION ONLY · no simulation results loaded', f'YAML: {path}']
    if not entries:
        lines.append('No active tallies attached to this component.')
    for i, e in enumerate(entries):
        lines.extend([f"\n{i+1}. {'Cell total' if e['kind'] == 'cell_tallies' else 'Spatial mesh'}",
                      'Scores: ' + ', '.join(e.get('scores', [])),
                      'Filters: ' + (', '.join(e.get('filters', [])) or 'All particles / energies'),
                      'Nuclides: ' + (', '.join(e.get('nuclides', [])) or 'All nuclides')])
        if e['kind'] == 'mesh_tallies':
            lines.append(f"Bins (r, φ, z): {e.get('dimensions', 'Missing dimensions')}")
    return '\n'.join(lines)


def setup_example(component):
    entry = {'cell_tallies': [{'scores': ['flux'], 'filters': ['neutron_filter']}],
             'mesh_tallies': [{'scores': ['heating'], 'dimensions': [20, 1, 30]}]}
    if component.name.startswith('blanket_'):
        key = 'central_cell' if component.group == 'CC' else 'end_plug'
        section = 'blanket' if key == 'central_cell' else 'central_cylinder'
        i = int(component.meta['layer_index'])
        config = {'breeder': entry} if i == 0 else {'layer_tallies': [dict(position=i, **entry)]}
        return 'Merge into the existing YAML section; then rebuild the OpenMC model.\n' + yaml.safe_dump({key: {section: {'tallies': config}}}, sort_keys=False)
    if component.name.startswith(('lf_coil_central_', 'lf_coil_left_', 'lf_coil_right_', 'hf_coil_left_inward_', 'hf_coil_left_outward_', 'hf_coil_right_inward_', 'hf_coil_right_outward_')):
        key = 'central_cell' if component.name.startswith('lf_coil_central_') else 'end_plug'
        section = 'hf_coil' if component.group == 'HF' else 'lf_coil'
        return yaml.safe_dump({key: {section: {section + '_tallies': entry}}}, sort_keys=False)
    if component.group == 'CC':
        i = int(component.meta['layer_index'])
        example = {'central_cell': {'tallies': {'layer_tallies': [dict(position=i, **entry)]}}}
    elif component.group in ('HF', 'LF'):
        section = 'hf_coil' if component.group == 'HF' else 'lf_coil'
        example = {section: {section + '_tallies': entry}}
    elif component.group == 'ends':
        example = {'end_cell': {'end_cell_tallies': entry}}
    else:
        return 'This builder does not attach tallies to this component. Select a CC layer, coil magnet, or end-cell shell.'
    return 'Merge into the existing YAML section; then rebuild the OpenMC model.\n' + yaml.safe_dump(example, sort_keys=False)


def mesh_preview(component, dimensions):
    """Bounded preview; sample large grids without changing configured bin counts."""
    if not isinstance(dimensions, (list, tuple)) or len(dimensions) != 3:
        raise ValueError('Mesh dimensions must contain three positive integers')
    if any(isinstance(n, bool) or not isinstance(n, int) or n < 1 for n in dimensions):
        raise ValueError('Mesh dimensions must contain three positive integers')
    bounds = component.mesh.bounds
    radii = np.hypot(component.mesh.points[:, 0], component.mesh.points[:, 1])
    r = np.linspace(radii.min(), radii.max(), min(dimensions[0], 16) + 1)
    phi = np.linspace(0, 2*np.pi, min(max(dimensions[1], 12), 32) + 1)
    z = np.linspace(bounds[4], bounds[5], min(dimensions[2], 24) + 1)
    # Draw the cylindrical boundary cage directly: smooth rings, axial
    # generators, and radial spokes on the end planes. No triangle edges.
    points, lines = [], []

    def add_line(coords):
        start = len(points)
        points.extend(coords)
        lines.extend([len(coords), *range(start, start + len(coords))])

    smooth_phi = np.linspace(0, 2*np.pi, 145)
    for radius in (r[0], r[-1]):
        for height in z:
            add_line(np.column_stack((radius*np.cos(smooth_phi), radius*np.sin(smooth_phi),
                                      np.full_like(smooth_phi, height))))
        for angle in phi[:-1]:
            add_line([[radius*np.cos(angle), radius*np.sin(angle), z[0]],
                      [radius*np.cos(angle), radius*np.sin(angle), z[-1]]])
    for height in (z[0], z[-1]):
        for radius in r[1:-1]:
            add_line(np.column_stack((radius*np.cos(smooth_phi), radius*np.sin(smooth_phi),
                                      np.full_like(smooth_phi, height))))
        for angle in phi[:-1]:
            add_line([[r[0]*np.cos(angle), r[0]*np.sin(angle), height],
                      [r[-1]*np.cos(angle), r[-1]*np.sin(angle), height]])
    return pv.PolyData(np.asarray(points), lines=np.asarray(lines))


def demo_heating_mesh(component, dimensions):
    """Synthetic cell-wise heating on cylindrical tally bins (arbitrary units).

    Angular subdivisions only tessellate curved cells; every subdivision of a
    logical tally bin receives the same value. No statepoint data are loaded.
    """
    mesh_preview(component, dimensions)  # shared input validation
    nr, nphi, nz = dimensions
    if nr * nphi * nz > 200000:
        raise ValueError('Demo supports at most 200,000 logical tally bins')
    bounds = component.mesh.bounds
    radii = np.hypot(component.mesh.points[:, 0], component.mesh.points[:, 1])
    r = np.linspace(radii.min(), radii.max(), nr + 1)
    subdivisions = max(1, int(np.ceil(72 / nphi)))
    angles = np.linspace(0, 2*np.pi, nphi * subdivisions + 1)
    z = np.linspace(bounds[4], bounds[5], nz + 1)
    rr, pp, zz = np.meshgrid(r, angles, z, indexing='ij')
    grid = pv.StructuredGrid(rr*np.cos(pp), rr*np.sin(pp), zz)
    radial = (np.arange(nr) + .5) / nr
    axial = (np.arange(nz) + .5) / nz
    angular = (np.arange(nphi) + .5) * 2*np.pi / nphi
    values = demo_heating_field(radial[:, None, None], angular[None, :, None], axial[None, None, :])
    values = np.repeat(values, subdivisions, axis=1)
    grid.cell_data['Demo heating [a.u.]'] = values.ravel(order='F')
    return grid


# --- result datasets ---



@dataclass(frozen=True)
class CylindricalMesh:
    r: np.ndarray
    phi: np.ndarray
    z: np.ndarray
    origin: tuple[float, float, float] = (0., 0., 0.)

    def __post_init__(self):
        for name in ('r', 'phi', 'z'):
            a = np.array(getattr(self, name), dtype=float, copy=True)
            if a.ndim != 1 or len(a) < 2 or not np.isfinite(a).all() or not (np.diff(a) > 0).all():
                raise ValueError(f'{name} edges must be finite and strictly increasing')
            a.setflags(write=False)
            object.__setattr__(self, name, a)
        if self.r[0] < 0 or self.phi[0] < -1e-9 or self.phi[-1] > 2*np.pi + 1e-9:
            raise ValueError('Invalid radial or azimuthal extent')
        if len(self.origin) != 3 or not np.isfinite(self.origin).all():
            raise ValueError('Invalid cylindrical mesh origin')

    @property
    def shape(self):
        return tuple(len(a)-1 for a in (self.r, self.phi, self.z))


@dataclass(frozen=True)
class TallyDataset:
    component_id: str
    mesh: CylindricalMesh
    mean: np.ndarray
    std_dev: np.ndarray | None = None
    score: str = 'heating'
    nuclide: str = 'total'
    units: str = 'eV/source'
    source: str = 'openmc'
    provenance: dict = field(default_factory=dict)
    id: str = field(default_factory=lambda: str(uuid4()))

    def __post_init__(self):
        for key in ('mean', 'std_dev'):
            value = getattr(self, key)
            if value is None:
                continue
            arr = np.array(value, dtype=float, copy=True)
            if arr.shape != self.mesh.shape or not np.isfinite(arr).all():
                raise ValueError(f'{key} must be finite with shape {self.mesh.shape}')
            if key == 'std_dev' and np.any(arr < 0):
                raise ValueError('Standard deviations cannot be negative')
            arr.setflags(write=False)
            object.__setattr__(self, key, arr)
        if self.source not in ('synthetic', 'openmc'):
            raise ValueError('Unknown results source')

    def values(self, quantity='mean'):
        if quantity == 'mean':
            return self.mean
        if self.std_dev is None:
            raise ValueError('This dataset has no uncertainty estimates')
        if quantity == 'std_dev':
            return self.std_dev
        if quantity == 'relative_error':
            return np.divide(self.std_dev, np.abs(self.mean), out=np.full_like(self.mean, np.nan), where=self.mean != 0)
        raise ValueError('Unknown result quantity')

    def scalar_name(self, quantity):
        units = 'fraction' if quantity == 'relative_error' else self.units
        prefix = 'DEMO ' if self.source == 'synthetic' else ''
        return f'{prefix}{self.score} {quantity.replace("_", " ")} [{units}]'

    def save(self, path):
        metadata = {'version': 1, 'component_id': self.component_id, 'score': self.score,
                    'nuclide': self.nuclide, 'units': self.units, 'source': self.source,
                    'provenance': self.provenance, 'id': self.id, 'has_std_dev': self.std_dev is not None}
        np.savez_compressed(path, r=self.mesh.r, phi=self.mesh.phi, z=self.mesh.z,
                            origin=self.mesh.origin, mean=self.mean,
                            std_dev=self.std_dev if self.std_dev is not None else np.empty(0),
                            metadata=json.dumps(metadata))

    @classmethod
    def load(cls, path):
        with np.load(path, allow_pickle=False) as f:
            meta = json.loads(str(f['metadata']))
            if meta.pop('version') != 1:
                raise ValueError('Unsupported results version')
            has_std = meta.pop('has_std_dev')
            mesh = CylindricalMesh(f['r'], f['phi'], f['z'], tuple(f['origin']))
            return cls(mesh=mesh, mean=f['mean'], std_dev=f['std_dev'] if has_std else None, **meta)


def synthetic_dataset(component_id, grid):
    nr, nphi, nz = grid.shape
    r, p, z = np.meshgrid((np.arange(nr)+.5)/nr, (grid.phi[:-1]+grid.phi[1:])/2,
                          (np.arange(nz)+.5)/nz, indexing='ij')
    values = demo_heating_field(r, p, z)
    return TallyDataset(component_id, CylindricalMesh(grid.r, grid.phi, grid.z), values,
                        units='a.u.', source='synthetic', provenance={'generator': 'demo-heating-v1'})


def dataset_from_tally(tally, component_id, *, score='heating', nuclide='total', filter_bins=None, provenance=None):
    """Read one explicit score/nuclide/filter selection; never silently sum bins."""
    import openmc
    filter_bins = dict(filter_bins or {})
    meshes = [(i, f) for i, f in enumerate(tally.filters) if isinstance(f, openmc.MeshFilter)]
    if len(meshes) != 1 or not isinstance(meshes[0][1].mesh, openmc.CylindricalMesh):
        raise ValueError('Select a tally with exactly one cylindrical MeshFilter')
    mi, mf = meshes[0]
    if score not in tally.scores or nuclide not in tally.nuclides:
        raise ValueError('Score or nuclide is not available in this tally')
    si, ni = list(tally.scores).index(score), list(tally.nuclides).index(nuclide)
    selectors = []
    for i, f in enumerate(tally.filters):
        if i == mi:
            selectors.append(slice(None))
        else:
            if f.num_bins > 1 and i not in filter_bins:
                raise ValueError(f'Filter {i} ({type(f).__name__}) has {f.num_bins} bins; select a bin index explicitly')
            j = int(filter_bins.get(i, 0))
            if not 0 <= j < f.num_bins:
                raise ValueError(f'Invalid bin for filter {i}')
            selectors.append(j)
    if any(i == mi or not 0 <= i < len(tally.filters) for i in filter_bins):
        raise ValueError('Invalid extra-filter selection')
    m = mf.mesh
    origin = np.asarray(m.origin, float)
    translation = getattr(mf, 'translation', None)
    if translation is not None:
        origin = origin + np.asarray(translation, float)
    mesh = CylindricalMesh(m.r_grid, m.phi_grid, m.z_grid, tuple(origin))
    arrays = []
    for quantity in ('mean', 'std_dev'):
        shaped = tally.get_reshaped_data(value=quantity, expand_dims=False)
        flat = shaped[tuple(selectors)+(ni, si)]
        arr = np.empty(mesh.shape, dtype=float)
        # Actual OpenMC bins encode ordering, including its radial-fast convention.
        for value, index in zip(flat, mf.bins):
            arr[tuple(int(i)-1 for i in index)] = value
        arrays.append(arr)
    units = 'eV/source' if score in ('heating', 'heating-local') else 'native OpenMC units/source'
    info = dict(provenance or {}, tally_id=int(tally.id), tally_name=tally.name,
                mesh_id=int(m.id), filter_bins={str(k): v for k,v in filter_bins.items()},
                normalization='per source particle', array_order='r,phi,z')
    return TallyDataset(component_id, mesh, *arrays, score=score, nuclide=nuclide, units=units, provenance=info)


def load_statepoint(path, component_id, tally_id, **selection):
    import openmc
    with openmc.StatePoint(str(path), autolink=False) as sp:
        if tally_id not in sp.tallies:
            raise ValueError(f'Tally {tally_id} not found')
        return dataset_from_tally(sp.tallies[tally_id], component_id,
                                  provenance={'statepoint': str(Path(path).resolve())}, **selection)
