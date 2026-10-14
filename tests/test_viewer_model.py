"""Model building (all cutaways preloaded), display names and upload parsing."""
import unittest
import math
from pathlib import Path

from src.paratan.viewer.app import MODES, build_model
from src.paratan.viewer.components import component_label, material_label, meta_row

INPUT = Path(__file__).resolve().parents[1] / 'input_files/simple_parametric_input_new.yaml'


class ModelTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.model = build_model(INPUT.read_text(), INPUT.name, 48)

    def test_every_cutaway_is_precomputed(self):
        for mode in MODES:
            self.assertEqual(set(self.model.meshes[mode]), set(self.model.by_name), mode)
        for name in self.model.by_name:
            full = self.model.meshes['none'][name].volume
            for mode in ('half', 'quarter'):
                self.assertLess(self.model.meshes[mode][name].volume, full * 1.001, (mode, name))

    def test_profiled_parts_get_edge_lines(self):
        self.assertTrue(all(name in self.model.outlines['quarter']
                            for name, c in self.model.by_name.items() if c.profile is not None))

    def test_breeder_fills_space_around_outer_vessel(self):
        model = self.model
        vv, cc = model.data['vacuum_vessel'], model.data['central_cell']
        wall = sum(layer['thickness'] for layer in vv['structure'].values())
        rc = vv['central_radius'] + wall
        end = cc['axial_length'] / 2
        ro = rc + cc['layers'][0]['thickness']
        style = str(vv.get('geometry_style', 'conical') or 'conical').lower()
        if style == 'perpendicular':
            # CC cylinder minus outer stepped vessel (central cylinder radius).
            self.assertLessEqual(end, vv['central_axial_length'] / 2 + wall)
            expected = math.pi * (ro**2 - rc**2) * 2 * end
        else:
            rb = vv['bottleneck_radius'] + wall
            straight_half, cone_end = vv['outer_axial_length']/2, vv['central_axial_length']/2
            re = rc + (rb-rc) * (end-straight_half)/(cone_end-straight_half)
            cavity = math.pi * (2*straight_half*rc**2 +
                               2*(end-straight_half)*(rc**2 + rc*re + re**2)/3)
            expected = math.pi * ro**2 * 2*end - cavity
        comp = model.by_name['cc_layer_0_eutectic_breeding_material']
        self.assertAlmostEqual(comp.mesh.volume / expected, 1, delta=.01)
        for mode, fraction in [('none',1), ('half',.5), ('quarter',.75)]:
            # Port penetrations are subtracted separately by model preparation;
            # verify the underlying blanket solid and all its analytic cuts.
            from src.paratan.viewer.model import MODE_SECTIONS
            mesh = comp.profile.revolve(MODE_SECTIONS[mode], n_theta=48)
            self.assertEqual(mesh.clean(tolerance=1e-7, absolute=True).n_open_edges, 0)
            self.assertAlmostEqual(mesh.volume / expected, fraction, delta=.01)

    def test_malformed_input_is_rejected(self):
        for text in ('- 1\n- 2', 'just text'):
            with self.assertRaises(ValueError):
                build_model(text, 'bad.yaml', 24)
        with self.assertRaises(KeyError):
            build_model('vacuum_vessel: {}\n', 'bad.yaml', 24)


class LabelTests(unittest.TestCase):
    def test_names_are_readable(self):
        parts = {c.name: c for c in build_model(INPUT.read_text(), INPUT.name, 24).components}
        self.assertEqual(parts['cc_layer_1_BW_FNSF'].label, 'Layer 2 · Back wall (FNSF)')
        self.assertEqual(parts['hf_coil_right_casing_1'].label, 'Right · Casing 2')
        self.assertEqual(parts['lf_coil_0_magnet'].label, 'Coil 1 · Magnet')
        self.assertEqual(parts['vv_first_wall'].label, 'First wall')
        for comp in parts.values():
            self.assertNotIn('_', comp.label, comp.name)
            self.assertTrue(comp.label[0].isupper(), comp.label)

    def test_materials_and_rows(self):
        self.assertEqual(material_label('water_cooled_wc'), 'Water-cooled WC')
        self.assertEqual(material_label('some_new_material'), 'Some new material')
        self.assertEqual(meta_row('r_inner_cm', 78.0), {'label': 'Inner radius', 'value': '78', 'unit': 'cm'})
        self.assertEqual(meta_row('layer_index', 0)['value'], '1')
        self.assertEqual(component_label('unknown_thing', 'CC', 'x', {}), 'Unknown thing')


if __name__ == '__main__':
    unittest.main()
