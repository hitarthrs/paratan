"""Device adapters: the viewer core consumes components, never a device YAML schema."""
from __future__ import annotations
from typing import Protocol
import numpy as np
from src.paratan.viewer.manifest import device_kind
from src.paratan.viewer.simple_mirror_meshes import (
    build_simple_mirror_components, _add, _vv_z_extents, _vessel_outline,
    _conical_vessel_solid, _conical_vessel_shell,
)
from src.paratan.viewer.revolve import Profile
from src.paratan.viewer.inspection import tally_configuration


class DeviceAdapter(Protocol):
    kind: str
    def build_components(self, data: dict, n_theta: int): ...


def _entries(config):
    return [dict(e, kind=kind) for kind in ('cell_tallies', 'mesh_tallies')
            for e in (config or {}).get(kind) or []]


class SimpleMirrorAdapter:
    kind = 'simple_mirror'
    def build_components(self, data, n_theta):
        parts = build_simple_mirror_components(data, n_theta=n_theta)
        for c in parts:
            path, entries = tally_configuration(data, c)
            c.meta['_tally_path'], c.meta['_tally_entries'] = path, entries
        return parts



class TandemAdapter:
    """Axisymmetric display profiles following TandemMachineBuilder dimensions.

    Boolean exclusions involving neighboring parts remain an explicitly marked
    display approximation. OpenMC IDs/bindings are supplied by its manifest.
    """
    kind = 'tandem'

    def build_components(self, data, n_theta):
        parts = []
        vv = data['vacuum_vessel']
        sep = float(vv['central_cell_end_plug_separation_distance'])
        mid = (vv['central_cell']['central_axis_length']/2 + sep + vv['end_plug']['central_axis_length']/2)
        midplanes = {'central': 0., 'left': -mid, 'right': mid}
        outlines, outer_r, outer_bottle = {}, {}, {}
        hf_centers = {}

        def add(name, group, material, profile, meta=None, config=None, path='', opacity=1.):
            _add(parts, name=name, group=group, material=material, profile=profile,
                 meta=dict(meta or {}, _tally_entries=_entries(config), _tally_path=path),
                 n_theta=n_theta, opacity=opacity)
            parts[-1].label = name.replace('_', ' ').capitalize()

        for key in ('central','left','right'):
            section = vv['central_cell'] if key == 'central' else vv['end_plug']
            fw = data['first_wall']['central_cell' if key=='central' else 'end_plug']['layers']
            left_len, right_len = ((sep/10,sep/10) if key=='central' else
                                  (vv['end_axial_distance'], .9*sep) if key=='left' else
                                  (.9*sep,vv['end_axial_distance']))
            z = _vv_z_extents({'axial_midplane':midplanes[key], 'outer_axial_length':section['outer_axial_length'],
                'central_axial_length':section['central_axis_length'], 'left_bottleneck_length':left_len,
                'right_bottleneck_length':right_len})
            rc, rb = float(section['central_radius']), float(vv['bottleneck radius'])
            add(f'{key}_vv_cell', 'VV', 'vacuum', _conical_vessel_solid(rc,rb,z), opacity=.15)
            for j, layer in enumerate(fw):
                t = float(layer['thickness'])
                add(f'fw_{key}_layer_{j}', 'FW', layer['material'], _conical_vessel_shell(rc,rb,rc+t,rb+t,z),
                    {'layer_index':j,'thickness_cm':t})
                rc, rb = rc+t, rb+t
            outlines[key] = _vessel_outline(rc,rb,z)
            outer_r[key],outer_bottle[key] = rc,rb
            if key == 'central':
                section_data = data['central_cell'].get('blanket')
                if section_data is None:
                    envelope = data['central_cell']['test_region']
                    section_data = {'axial_length': envelope['axial_length'], 'layers': [
                        {'thickness': envelope['radial_thickness'], 'material': 'vacuum'}]}
            else:
                section_data = data['end_plug']['central_cylinder']
            length = float(section_data['axial_length'])
            z0,z1 = midplanes[key]-length/2,midplanes[key]+length/2
            ri=rc
            for j, layer in enumerate(section_data['layers']):
                ro = ri+float(layer['thickness'])
                if j == 0:
                    contour = outlines[key]
                    zs = sorted(set([z0,z1]+[z for r,z in contour if z0<z<z1]))
                    inner = [(float(np.interp(z,[p[1] for p in contour],[p[0] for p in contour])),z) for z in zs]
                    profile=Profile.polygon([(ro,z0),(ro,z1),*inner[::-1]])
                else:
                    profile=Profile.rect(ri,ro,z0,z1)
                tallies=section_data.get('tallies') or {}
                config=tallies.get('breeder',{}) if j==0 else next((e for e in tallies.get('layer_tallies',[]) if e.get('position')==j),{})
                path=('central_cell.blanket' if key=='central' else 'end_plug.central_cylinder')+'.tallies'
                add(f'blanket_{key}_layer_{j}', 'CC' if key=='central' else 'EP',layer['material'],profile,
                    {'layer_index':j,'z0_cm':midplanes[key],'r_inner_cm':ri,'r_outer_cm':ro},config,path)
                if key=='central' and 'blanket' not in data['central_cell']:
                    parts[-1].label = 'Test region envelope (modules not modeled)'
                    parts[-1].meta['representation'] = 'Test region envelope; module geometry is not included'
                ri=ro
            coil=data['central_cell' if key=='central' else 'end_plug'].get('lf_coil') or {}
            dims=coil.get('inner_dimensions',{})
            shell=coil.get('shell_thicknesses',{})
            for j,pos in enumerate(coil.get('positions',[])):
                zc=midplanes[key]+pos
                front,back,axial=(float(shell.get(k,0)) for k in ('front','back','axial'))
                r0=ri+front;r1=r0+dims['radial_thickness'];lo=zc-dims['axial_length']/2;hi=zc+dims['axial_length']/2
                cavity=Profile.rect(r0,r1,lo,hi)
                add(f'lf_coil_{key}_{j}_magnet','LF',coil['materials']['magnet'],cavity,{'z0_cm':zc},
                    coil.get('lf_coil_tallies'),f"{'central_cell' if key=='central' else 'end_plug'}.lf_coil.lf_coil_tallies")
                add(f'lf_coil_{key}_{j}_shield','LF',coil['materials']['shield'],
                    Profile.shell(Profile.rect(ri,r1+back,lo-axial,hi+axial),cavity),opacity=.55)
        hf_data=data['end_plug'].get('hf_coil') or {}
        for key in ('left','right'):
            hf=hf_data.get(key)
            if not hf:
                continue
            mag,shield=hf['magnet'],hf['shield'];th=[float(l['thickness']) for l in hf.get('casing_layers',[])]
            for direction in ('inward','outward'):
                sign=(1 if key=='left' else -1) * (1 if direction=='inward' else -1)
                zc=midplanes[key]+sign*(data['end_plug']['central_cylinder']['axial_length']/2+
                    shield['shield_central_cell_gap']+shield['axial_thickness'][0]+sum(th)+mag['axial_thickness']/2)
                if direction=='outward':hf_centers[key]=zc
                ri=outer_bottle[key]+shield['radial_thickness'][0]+shield['radial_gap_before_casing']+sum(th)
                ro=ri+mag['radial_thickness'];lo=zc-mag['axial_thickness']/2;hi=zc+mag['axial_thickness']/2
                add(f'hf_coil_{key}_{direction}_magnet','HF',mag['material'],Profile.rect(ri,ro,lo,hi),{'z0_cm':zc},
                    hf_data.get('hf_coil_tallies'),'end_plug.hf_coil.hf_coil_tallies')
                for j,layer in enumerate(hf.get('casing_layers',[])):
                    t=th[j];cavity=Profile.rect(ri,ro,lo,hi)
                    ri,ro,lo,hi=ri-t,ro+t,lo-t,hi+t
                    add(f'hf_coil_{key}_{direction}_casing_{j}','HF',layer['material'],Profile.shell(Profile.rect(ri,ro,lo,hi),cavity),opacity=.7)
                base=outer_bottle[key]+shield['radial_gap_before_casing']
                cavity=Profile.rect(base+shield['radial_thickness'][0],base+shield['radial_thickness'][0]+mag['radial_thickness']+2*sum(th),lo,hi)
                outer=Profile.rect(base,cavity.bounds[1]+shield['radial_thickness'][1],lo-shield['axial_thickness'][0],hi+shield['axial_thickness'][0])
                add(f'hf_coil_{key}_{direction}_shield','HF',shield['material'],Profile.shell(outer,cavity))
        end=data.get('end_cell') or {}
        for key,center in hf_centers.items():
            if not end:continue
            hf=hf_data[key];sign=-1 if key=='left' else 1;t=float(end['shell_thickness']);L=float(end['axial_length']);ro=float(end['diameter'])/2
            zc=center+sign*(hf['shield']['axial_thickness'][0]+sum(l['thickness'] for l in hf['casing_layers'])+hf['magnet']['axial_thickness']/2+t+L/2+5)
            lo,hi=zc-L/2,zc+L/2
            profile=Profile.polygon([(0,lo-t),(ro,lo-t),(ro,hi+t),(0,hi+t),(0,hi),(ro-t,hi),(ro-t,lo),(0,lo)])
            add(f'end_cell_{key}_shell','ends',end['shell_material'],profile,{'z0_cm':zc,'axial_length_cm':L})
            add(f'end_cell_{key}_inner','ends',end['inner_material'],Profile.rect(0,ro-t,lo,hi),opacity=.15)
        return parts


ADAPTERS = {'simple_mirror': SimpleMirrorAdapter(), 'tandem': TandemAdapter()}


def adapter_for(data):
    return ADAPTERS[device_kind(data)]
