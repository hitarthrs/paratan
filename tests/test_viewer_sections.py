"""Exact sections of revolved solids and cylindrical tally-grid slices."""
import unittest
from pathlib import Path

import numpy as np

from src.paratan.viewer.results import CylGrid, demo_heating_field
from src.paratan.viewer.revolve import FULL, Profile, Section
from src.paratan.viewer.simple_mirror_meshes import build_from_yaml


def welded_open_edges(mesh):
    """Open edges after merging the duplicated seam / strip-corner points."""
    return mesh.clean(tolerance=1e-7, absolute=False).n_open_edges


class SectionTests(unittest.TestCase):
    def test_full_revolve_volume_and_watertight(self):
        profile = Profile.rect(10, 20, 0, 30)
        mesh = profile.revolve(FULL, n_theta=180)
        self.assertAlmostEqual(mesh.volume, np.pi * 300 * 30, delta=0.002 * np.pi * 300 * 30)
        self.assertEqual(welded_open_edges(mesh), 0)

    def test_wedge_and_axial_cuts_are_capped_and_exact(self):
        profile = Profile.rect(10, 20, 0, 30)
        section = Section.without_wedge(0.3, np.pi / 2, z=(5, 12))
        mesh = profile.revolve(section, n_theta=180)
        expected = 0.75 * np.pi * 300 * 7
        self.assertAlmostEqual(mesh.volume, expected, delta=0.002 * expected)
        self.assertEqual(welded_open_edges(mesh), 0)  # real cap faces, not an open shell
        self.assertAlmostEqual(mesh.bounds[4], 5)
        self.assertAlmostEqual(mesh.bounds[5], 12)

    def test_solid_on_axis_and_cone_sections_close(self):
        profile = Profile.polygon([(0, 0), (5, 0), (5, 10), (8, 15), (8, 25), (5, 30), (5, 40), (0, 40)])
        for section in (FULL, Section.without_wedge(0, np.pi / 2, z=(3, 35)), Section(z=(12, 18))):
            self.assertEqual(welded_open_edges(profile.revolve(section)), 0)

    def test_ring_profile_matches_outer_minus_cavity(self):
        ring = Profile.shell(Profile.rect(5, 30, 0, 40), Profile.rect(10, 20, 10, 30))
        expected = np.pi * (900 - 25) * 40 - np.pi * 300 * 20
        self.assertAlmostEqual(ring.revolve(FULL, n_theta=180).volume, expected, delta=0.002 * expected)
        half = ring.revolve(Section.without_wedge(1.0, np.pi), n_theta=180)
        self.assertAlmostEqual(half.volume, 0.5 * expected, delta=0.002 * expected)

    def test_empty_section_returns_none(self):
        profile = Profile.rect(10, 20, 0, 30)
        self.assertIsNone(profile.revolve(Section(z=(40, 50))))

    def test_analytic_normals_are_unit_and_outward(self):
        mesh = Profile.rect(10, 20, 0, 30).revolve(FULL)
        normals = mesh.point_data['Normals']
        self.assertTrue(np.allclose(np.linalg.norm(normals, axis=1), 1.0))
        outer = np.isclose(np.hypot(mesh.points[:, 0], mesh.points[:, 1]), 20) & (normals[:, 2] == 0)
        radial = mesh.points[outer, :2] / 20
        self.assertTrue(np.allclose(normals[outer, :2], radial))

    def test_normals_are_never_colour_scalars(self):
        # WebGL renderers colour a mesh by its active scalars: normals there paint a rainbow.
        mesh = Profile.rect(10, 20, 0, 30).revolve(Section.without_wedge(0, np.pi / 2))
        self.assertIsNone(mesh.active_scalars_name)
        self.assertIsNotNone(mesh.GetPointData().GetNormals())
        self.assertEqual(mesh.points.dtype, np.float32)

    def test_every_viewer_component_with_a_profile_sections_cleanly(self):
        path = Path(__file__).resolve().parents[1] / 'input_files/simple_parametric_input_new.yaml'
        section = Section.without_wedge(0.0, np.pi / 2, z=(-np.inf, 50))
        for comp in build_from_yaml(path, n_theta=48):
            if comp.profile is None:
                continue
            mesh = comp.profile.revolve(section, n_theta=48)
            if mesh is None:
                continue
            self.assertLessEqual(mesh.volume, comp.mesh.volume * 1.001, comp.name)
            self.assertGreater(mesh.volume, 0, comp.name)


class CylGridTests(unittest.TestCase):
    def setUp(self):
        self.grid = CylGrid.from_profile(Profile.rect(10, 20, 0, 30), [5, 8, 6])

    def test_bin_lookup_wraps_and_clamps(self):
        g = self.grid
        self.assertEqual(int(g.locate_phi(np.radians(50))), 1)
        self.assertEqual(int(g.locate_phi(np.radians(50) + 2 * np.pi)), 1)
        self.assertEqual(int(g.locate_z(0.0)), 0)
        self.assertEqual(int(g.locate_z(30.0)), 5)
        self.assertEqual(int(g.locate_z(7.0)), 1)

    def test_phi_plane_reads_the_bins_at_that_azimuth(self):
        values = np.arange(5 * 8 * 6, dtype=float).reshape(5, 8, 6)
        angle = np.radians(50)  # azimuthal bin 1
        plane = self.grid.plane_phi(angle, values, name='v')
        self.assertEqual(plane.n_cells, 5 * 6)
        # Extraction can reorder cells; check the value at each physical bin.
        x,y,z = plane.cell_centers().points.T
        ri = np.searchsorted(self.grid.r, np.hypot(x,y), side='right')-1
        zk = self.grid.locate_z(z)
        np.testing.assert_array_equal(plane.cell_data['v'], values[ri,1,zk])
        radii = np.hypot(plane.points[:, 0], plane.points[:, 1])
        self.assertAlmostEqual(radii.min(), 10)
        self.assertAlmostEqual(radii.max(), 20)
        angles = np.arctan2(plane.points[:, 1], plane.points[:, 0])
        self.assertTrue(np.allclose(angles, angle))

    def test_z_plane_reads_the_bins_at_that_height(self):
        values = np.arange(5 * 8 * 6, dtype=float).reshape(5, 8, 6)
        plane = self.grid.plane_z(12.0, values, name='v')  # axial bin 2
        self.assertTrue(np.allclose(plane.points[:, 2], 12.0))
        cells = plane.cell_data['v']
        self.assertTrue(set(np.unique(cells)) <= set(values[:, :, 2].ravel()))
        self.assertEqual(len(np.unique(cells)), 5 * 8)

    def test_slices_respect_the_section(self):
        section = Section.without_wedge(0.0, np.pi / 2, z=(-np.inf, 12.0))
        self.assertIsNone(self.grid.plane_z(20.0, demo_heating_field, section))
        plane = self.grid.plane_z(12.0, demo_heating_field, section)
        angles = np.mod(np.arctan2(plane.points[:, 1], plane.points[:, 0]), 2 * np.pi)
        self.assertGreaterEqual(angles.min() + 1e-9, 0)  # kept arc is 90..360 degrees
        self.assertFalse(((angles > 0.1) & (angles < np.pi / 2 - 0.1)).any())
        side = self.grid.plane_phi(1.0, demo_heating_field, section)
        self.assertAlmostEqual(side.bounds[5], 12.0)

    def test_volume_exposes_cut_faces_with_bin_values(self):
        section = Section.without_wedge(0.0, np.pi / 2)
        volume = self.grid.volume(demo_heating_field, section, name='v')
        full = self.grid.volume(demo_heating_field, FULL, name='v')
        self.assertGreater(volume.n_cells, 0)
        self.assertTrue(np.isfinite(volume.cell_data['v']).all())
        self.assertLess(volume.area, full.area + 2 * 10 * 30 * 2)  # wedge removed, two cut faces added

    def test_dimension_validation_and_display_cap(self):
        for bad in ([0, 1, 2], [1, 2], [1.5, 1, 2], [True, 1, 2]):
            with self.assertRaises(ValueError):
                CylGrid.from_profile(Profile.rect(1, 2, 0, 1), bad)
        big = CylGrid.from_profile(Profile.rect(1, 2, 0, 1), [1000, 1000, 1000])
        self.assertEqual(big.logical, (1000, 1000, 1000))
        self.assertLessEqual(max(big.shape), 360)


if __name__ == '__main__':
    unittest.main()
