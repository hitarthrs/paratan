"""Tally configuration inspection and lightweight cylindrical previews."""
from __future__ import annotations

import numpy as np
import pyvista as pv
import yaml

from src.paratan.viewer.cyl_grid import demo_heating_field


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
