"""FFH viewer: blanket inlet/outlet angles, meshes, and model smoke."""
from __future__ import annotations

import math
import unittest
from pathlib import Path

import yaml

from src.paratan.viewer.port_openings import (
    angular_distance_deg,
    build_blanket_inlet_outlet_meshes,
    resolve_port_angles,
)
from src.paratan.viewer.manifest import device_kind
from src.paratan.viewer.model import build_model


FFH_YAML = Path(__file__).resolve().parents[1] / "input_files/ffh_parametric_input.yaml"


class AngleResolverTests(unittest.TestCase):
    def test_auto_four_ports_clear_nbi_exclusion(self):
        angles = resolve_port_angles(4, None, (0.0, 180.0), 25.0)
        self.assertEqual(len(angles), 4)
        for a in angles:
            for forbidden in (0.0, 180.0):
                self.assertGreaterEqual(angular_distance_deg(a, forbidden), 25.0)

    def test_manual_conflict_raises(self):
        with self.assertRaises(ValueError):
            resolve_port_angles(2, [0.0, 90.0], (0.0, 180.0), 25.0)

    def test_manual_safe_angles(self):
        angles = resolve_port_angles(4, [45.0, 135.0, 225.0, 315.0], (0.0, 180.0), 25.0)
        self.assertEqual(angles, [45.0, 135.0, 225.0, 315.0])


class PipeMeshTests(unittest.TestCase):
    def test_four_ports_yield_sixteen_components(self):
        data = yaml.safe_load(FFH_YAML.read_text())
        parts = build_blanket_inlet_outlet_meshes(data, n_theta=32)
        # 4 azimuths × (inlet+outlet) × (wall+flow) = 16
        self.assertEqual(len(parts), 16)
        self.assertEqual(sum(1 for p in parts if p.group == "IO"), 16)
        self.assertEqual(sum(1 for p in parts if p.name.endswith("_wall")), 8)
        self.assertEqual(sum(1 for p in parts if p.name.endswith("_flow")), 8)
        for p in parts:
            self.assertGreater(p.mesh.n_cells, 0)
            self.assertEqual(p.mesh.n_open_edges, 0)


class FfhModelTests(unittest.TestCase):
    def test_device_kind_and_build(self):
        text = FFH_YAML.read_text()
        data = yaml.safe_load(text)
        self.assertEqual(device_kind(data), "ffh")
        model = build_model(text, FFH_YAML.name, 32)
        self.assertEqual(model.device_label, "FFH")
        io = [c for c in model.components if c.group == "IO"]
        self.assertEqual(len(io), 16)
        self.assertTrue(any(c.name.startswith("inlet_") for c in io))
        self.assertTrue(any(c.name.startswith("outlet_") for c in io))
        # Cutaways include at least one IO mesh
        self.assertTrue(any(name.startswith("inlet_") for name in model.meshes["none"]))


if __name__ == "__main__":
    unittest.main()
