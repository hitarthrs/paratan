"""Port-only visual subtraction: plasma exclusion, bores and bounded mesh size."""
from pathlib import Path
import unittest
import numpy as np
from src.paratan.viewer.model import build_model
from src.paratan.viewer.port_openings import PortOpenings


class PortOpeningTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        path=Path(__file__).resolve().parents[1]/'input_files/simple_parametric_input_new.yaml'
        cls.model=build_model(path.read_text(),path.name,48)
        cls.openings=PortOpenings(cls.model.data)

    def test_ports_stop_at_plasma_in_every_cut(self):
        for mode, meshes in self.model.meshes.items():
            for name, mesh in meshes.items():
                if name.startswith('port_'):
                    self.assertGreaterEqual(self.openings.outside_plasma(mesh.points).min(),-.05,(mode,name))
                    self.assertLess(mesh.n_cells,5000,(mode,name))

    def test_central_layers_and_cut_faces_have_openings(self):
        for mode, meshes in self.model.meshes.items():
            for name,mesh in meshes.items():
                component=self.model.by_name[name]
                if component.group in {'CC','FW','VV','LF'} and name!='vacuum_vessel_plasma':
                    self.assertGreaterEqual(self.openings.outside_ports(mesh.points).min(),-.75,(mode,name))
            self.assertLess(sum(mesh.n_cells for mesh in meshes.values()),100000,mode)

    def test_bore_parameters_are_in_geometry_cache_key(self):
        import yaml
        data=yaml.safe_load(yaml.safe_dump(self.model.data))
        data['ports']['inner_radius']+=4
        other=build_model(yaml.safe_dump(data),'changed bore',48)
        name='cc_layer_0_eutectic_breeding_material'
        self.assertFalse(np.array_equal(self.model.meshes['half'][name].points,other.meshes['half'][name].points))

    def test_visible_port_walls_leave_the_bore_open_in_every_cut(self):
        for side,sign in [('left',1),('right',-1)]:
            self.assertEqual(self.model.by_name[f'port_{side}_core'].opacity,0)
            # Look along the bore, slightly below the y=0 section plane.
            # A filled end disk or incorrect cut cap would intercept this ray.
            start=np.array([-400.,-5.,-sign*400.])
            end=np.array([400.,-5.,sign*400.])
            for mode,meshes in self.model.meshes.items():
                for layer in (1,2):
                    name=f'port_{side}_wall_{layer}'
                    hits,_=meshes[name].ray_trace(start,end)
                    self.assertEqual(len(hits),0,(mode,name,hits))


if __name__=='__main__':
    unittest.main()
