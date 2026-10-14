"""Check tally bindings against the simple-model builder's target rules."""
import unittest
import numpy as np
from pathlib import Path

from src.paratan.viewer.results import demo_heating_mesh, mesh_preview, setup_example, tally_configuration
from src.paratan.viewer.simple_mirror_meshes import build_simple_mirror_components, load_simple_mirror_yaml


class InspectionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.data = load_simple_mirror_yaml(Path(__file__).resolve().parents[1] / 'input_files/simple_parametric_input_new.yaml')
        cls.parts = {c.name: c for c in build_simple_mirror_components(cls.data, n_theta=24)}

    def test_hf_configuration_only_targets_magnets(self):
        _, entries = tally_configuration(self.data, self.parts['hf_coil_right_magnet'])
        self.assertEqual([e['scores'] for e in entries], [['heating'], ['flux']])
        self.assertEqual(entries[1]['filters'], ['neutron_filter', 'fast_energies_filter'])
        _, entries = tally_configuration(self.data, self.parts['hf_coil_right_shield'])
        self.assertEqual(entries, [])

    def test_cc_breeder_and_specific_layer_bindings(self):
        data = dict(self.data)
        data['central_cell'] = dict(data['central_cell'], tallies={
            'breeder': {'cell_tallies': [{'scores': ['H3-production']}]},
            'layer_tallies': [{'position': 1, 'cell_tallies': [{'scores': ['flux']}]}],
        })
        _, breeder = tally_configuration(data, self.parts['cc_layer_0_eutectic_breeding_material'])
        _, wall = tally_configuration(data, self.parts['cc_layer_1_BW_FNSF'])
        self.assertEqual(breeder[0]['scores'], ['H3-production'])
        self.assertEqual(wall[0]['scores'], ['flux'])
        self.assertIn('position: 1', setup_example(self.parts['cc_layer_1_BW_FNSF']))

    def test_preview_is_bounded_and_validates_dimensions(self):
        component = self.parts['hf_coil_right_magnet']
        mesh = mesh_preview(component, [1000, 1000, 1000])
        self.assertLess(mesh.n_points, 20000)
        self.assertAlmostEqual(mesh.bounds[4], component.mesh.bounds[4])
        self.assertAlmostEqual(mesh.bounds[5], component.mesh.bounds[5])
        for bad in ([0, 1, 2], [1, 2], [1.5, 1, 2], [True, 1, 2]):
            with self.assertRaises(ValueError):
                mesh_preview(component, bad)

    def test_demo_heating_preserves_logical_bins(self):
        component = self.parts['hf_coil_right_magnet']
        mesh = demo_heating_mesh(component, [3, 1, 4])
        values = mesh.cell_data['Demo heating [a.u.]'].reshape((3, 72, 4), order='F')
        self.assertTrue(np.isfinite(values).all())
        self.assertTrue((values > 0).all())
        self.assertTrue(np.allclose(values, values[:, :1, :]))
        self.assertFalse(np.allclose(values[0], values[-1]))
        self.assertAlmostEqual(mesh.bounds[4], component.mesh.bounds[4])
        self.assertAlmostEqual(mesh.bounds[5], component.mesh.bounds[5])
