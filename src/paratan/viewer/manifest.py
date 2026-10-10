"""Versioned, JSON-safe device/model contract shared with simulation builders.

No Trame, VTK, or OpenMC imports. Cell IDs come from a built simulation, never
from guessed numeric ranges. Analytic display geometry is explicitly identified.
"""
from __future__ import annotations
from dataclasses import dataclass, asdict, field
from hashlib import sha256
import json
from pathlib import Path

SCHEMA_VERSION = 1


def input_fingerprint(data):
    return sha256(json.dumps(data, sort_keys=True, separators=(',', ':'), default=str).encode()).hexdigest()


def device_kind(data):
    vv = data.get('vacuum_vessel') or {}
    return 'tandem' if 'central_cell' in vv and 'end_plug' in vv else 'simple_mirror'


def component_id(device, name):
    return f'{device}:{name}'


@dataclass
class ComponentRecord:
    id: str
    name: str
    label: str
    group: str
    material: str
    cell_ids: list[int] = field(default_factory=list)
    geometry: dict = field(default_factory=dict)
    tallies: list[dict] = field(default_factory=list)


@dataclass
class ModelManifest:
    device: str
    input_hash: str
    components: list[ComponentRecord]
    version: int = SCHEMA_VERSION
    provenance: dict = field(default_factory=dict)

    def to_dict(self):
        return asdict(self)

    def write(self, path):
        Path(path).write_text(json.dumps(self.to_dict(), indent=2, allow_nan=False))

    @classmethod
    def read(cls, path):
        data = json.loads(Path(path).read_text())
        if data.get('version') != SCHEMA_VERSION:
            raise ValueError('Unsupported model manifest version')
        records = [ComponentRecord(**r) for r in data['components']]
        if len({r.id for r in records}) != len(records):
            raise ValueError('Duplicate component IDs in manifest')
        return cls(data['device'], data['input_hash'], records, provenance=data.get('provenance', {}))

    def find(self, name):
        return next((r for r in self.components if r.name == name or r.id == name), None)


def display_manifest(data, components):
    device = device_kind(data)
    records = []
    for c in components:
        records.append(ComponentRecord(component_id(device, c.name), c.name, c.label, c.group, c.material,
            geometry={'representation': 'analytic_preview', 'bounds_cm': [float(v) for v in c.mesh.bounds],
                      'axisymmetric': c.profile is not None},
            tallies=list(c.meta.get('_tally_entries', []))))
    return ModelManifest(device, input_fingerprint(data), records,
                         provenance={'geometry': 'Display geometry; simulation boolean exclusions may differ.'})


def export_simulation_manifest(data, cells, tally_bindings, path):
    """Export actual cell/tally IDs without requiring the visualization stack."""
    device = device_kind(data)
    records = []
    for cell in cells:
        name = cell.name or f'cell_{cell.id}'
        material = getattr(cell.fill, 'name', '') or str(getattr(cell.fill, 'id', 'void'))
        record = ComponentRecord(component_id(device, name), name, name, 'simulation', material,
                                 cell_ids=[int(cell.id)], geometry={'representation': 'openmc_csg'})
        record.tallies = [dict(b) for b in tally_bindings if b['cell_id'] == cell.id]
        records.append(record)
    manifest = ModelManifest(device, input_fingerprint(data), records,
                             provenance={'source': 'OpenMC model builder', 'length_unit': 'cm'})
    manifest.write(path)
    return manifest
