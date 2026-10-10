"""Validate saved views before applying any changes to the live scene."""
from copy import deepcopy
import numpy as np


def validate_view(view, components):
    if not isinstance(view, dict):
        raise ValueError('A saved view must be an object')
    result = deepcopy(view)
    camera = np.asarray(view.get('camera'), dtype=float)
    if camera.shape != (3,3) or not np.isfinite(camera).all():
        raise ValueError('Invalid camera')
    if np.linalg.norm(camera[0]-camera[1]) < 1e-9 or np.linalg.norm(camera[2]) < 1e-9:
        raise ValueError('Invalid camera direction')
    result['camera'] = camera.tolist()
    for key, default, lo, hi in [('explode',0,0,1),('demo_opacity',1,0,1),
                                 ('parallel_scale',1,1e-12,np.inf),('view_angle',30,1,179)]:
        value = float(view.get(key,default))
        if not np.isfinite(value) or not lo <= value <= hi:
            raise ValueError(f'Invalid {key}')
        result[key] = value
    if view.get('section','none') not in ('none','half','quarter'):
        raise ValueError('Unknown section')
    if view.get('result_quantity','mean') not in ('mean','std_dev','relative_error'):
        raise ValueError('Unknown result quantity')
    index = int(view.get('tally_mesh_index',0))
    if index < 0:
        raise ValueError('Invalid tally mesh index')
    result['tally_mesh_index'] = index
    selection = view.get('selected','')
    if selection and selection not in components:
        raise ValueError('Saved component does not exist in this model')
    for key in ('groups','opacity','slice'):
        if not isinstance(view.get(key,{}),dict):
            raise ValueError(f'Invalid {key}')
    for name,value in view.get('opacity',{}).items():
        if name not in components or not np.isfinite(float(value)) or not 0 <= float(value) <= 1:
            raise ValueError('Invalid component opacity')
    for key in ('phi','z'):
        if not np.isfinite(float(view.get('slice',{}).get(key,0))):
            raise ValueError('Invalid slice position')
    if not isinstance(view.get('hidden',[]),list) or any(n not in components for n in view.get('hidden',[])):
        raise ValueError('Invalid hidden component list')
    return result
