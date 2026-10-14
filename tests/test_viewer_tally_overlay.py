"""Boundary-only cylindrical rendering preserves bin identity and provenance."""
import unittest
import numpy as np
from src.paratan.viewer.results import CylGrid
from src.paratan.viewer.results import synthetic_dataset
from src.paratan.viewer.results import describe_overlay, lift_boundary
from src.paratan.viewer.revolve import Section


class TallyOverlayTests(unittest.TestCase):
    def test_display_lift_preserves_bin_values_and_source_geometry(self):
        grid=CylGrid(np.array([1.,2.,4.]),np.array([0.,np.pi,2*np.pi]),
                     np.array([0.,3.,8.]),(2,2,2))
        section=Section.without_wedge(0,np.pi/2)
        surface=grid.volume(np.arange(8).reshape(2,2,2),section)
        original=surface.points.copy()
        lifted=lift_boundary(surface,grid,section)
        np.testing.assert_array_equal(surface.points,original)
        np.testing.assert_array_equal(surface.cell_data['value'],lifted.cell_data['value'])
        self.assertGreater(np.max(np.linalg.norm(lifted.points-original,axis=1)),.1)

    def test_surface_values_follow_nonuniform_cylindrical_bins(self):
        grid=CylGrid(np.array([1.,2.,4.]),np.array([0.,np.pi,2*np.pi]),
                     np.array([0.,3.,8.]),(2,2,2))
        field=np.arange(8).reshape(2,2,2)
        surface=grid.volume(field)
        x,y,z=surface.cell_centers().points.T
        ri=np.clip(np.searchsorted(grid.r,np.hypot(x,y),side='right')-1,0,1)
        pj=grid.locate_phi(np.arctan2(y,x))
        zk=grid.locate_z(z)
        np.testing.assert_array_equal(surface.cell_data['value'],field[ri,pj,zk])

    def test_large_mesh_draws_boundaries_without_an_interior_lattice(self):
        grid=CylGrid(np.linspace(1,10,101),np.linspace(0,2*np.pi,91),
                     np.linspace(0,50,101),(100,90,100))
        surface=grid.volume()
        self.assertLess(surface.n_points,50000)

    def test_summary_explicitly_labels_synthetic_values(self):
        grid=CylGrid(np.array([1.,2.]),np.array([0.,2*np.pi]),np.array([0.,4.]),(1,1,1))
        summary=describe_overlay(grid,synthetic_dataset('simple:coil',grid))
        self.assertIn('invented',summary['source'])
        self.assertIn('DEMO',summary['quantity'])
        self.assertTrue(summary['range'])


if __name__=='__main__': unittest.main()
