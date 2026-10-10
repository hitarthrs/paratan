"""Presentation metadata shared by all cylindrical tally overlays."""
import numpy as np


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
