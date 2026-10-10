"""Validated cylindrical results, independent of Trame and rendering.

Arrays always use (r, phi, z) order. OpenMC filter ordering is handled only by
its adapter. Native heating is eV/source; no power/volume normalization is implicit.
"""
from __future__ import annotations
from dataclasses import dataclass, field
from pathlib import Path
from uuid import uuid4
import json
import numpy as np


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
    from src.paratan.viewer.cyl_grid import demo_heating_field
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
