"""Behavioral checks for independent actions, result ordering and device adapters."""
import contextlib
import io
import json
import tempfile
import unittest
from pathlib import Path
import numpy as np
from src.paratan.viewer.model import SceneState, SceneActions
from src.paratan.viewer.manifest import ModelManifest, export_simulation_manifest
from src.paratan.viewer.results import CylindricalMesh, TallyDataset, dataset_from_tally
from src.paratan.viewer.model import build_model

ROOT = Path(__file__).resolve().parents[1]


class StateTests(unittest.TestCase):
    def setUp(self):
        self.state = SceneState()
        self.actions = SceneActions(self.state)
        self.actions.install(['a','b'], ['CC','HF'])

    def test_actions_and_reset_do_not_require_a_renderer(self):
        self.actions.select('a')
        self.actions.hide_selected()
        self.actions.set('vis_CC', False)
        self.actions.set('explode', .5)
        self.state.solo = True
        self.state.demo_heating = True
        self.actions.reset()
        self.assertFalse(self.state.hidden)
        self.assertTrue(all(self.state.groups.values()))
        self.assertFalse(self.state.solo or self.state.demo_heating or self.state.has_selection)
        self.assertEqual(self.state.explode, 0)

    def test_invalid_actions_leave_scene_unchanged(self):
        old = self.state.snapshot()
        for key,value in [('explode',float('nan')),('section_mode','bad'),('demo_opacity',2)]:
            with self.assertRaises(ValueError):self.actions.set(key,value)
        with self.assertRaises(ValueError):self.actions.select('missing')
        self.assertEqual(self.state.snapshot(), old)


class ResultTests(unittest.TestCase):
    def test_validation_uncertainty_and_safe_roundtrip(self):
        mesh=CylindricalMesh([1,2,4],[0,2*np.pi],[0,5,9],(0,0,100))
        values=np.array([[[0,10]],[[20,30]]],float)
        dataset=TallyDataset('tandem:coil',mesh,values,np.ones_like(values))
        self.assertTrue(np.isnan(dataset.values('relative_error')[0,0,0]))
        self.assertAlmostEqual(dataset.values('relative_error')[1,0,0],.05)
        self.assertFalse(dataset.mean.flags.writeable)
        with tempfile.TemporaryDirectory() as tmp:
            path=Path(tmp)/'result.npz';dataset.save(path);copy=TallyDataset.load(path)
            np.testing.assert_array_equal(copy.mean, dataset.mean)
            self.assertEqual(copy.mesh.origin, dataset.mesh.origin)
            self.assertEqual(copy.id, dataset.id)
        with self.assertRaises(ValueError):TallyDataset('x',mesh,np.zeros((1,1,1)))
        with self.assertRaises(ValueError):CylindricalMesh([2,1],[0,1],[0,1])

    def test_actual_openmc_filter_order_and_explicit_bin_selection(self):
        import openmc
        mesh=openmc.CylindricalMesh(r_grid=[1,2,4],phi_grid=[0,1,2*np.pi],z_grid=[0,3,8],origin=[0,0,20])
        mf=openmc.MeshFilter(mesh)
        t=openmc.Tally(name='ordering_fixture')
        t.filters=[openmc.ParticleFilter(['neutron','photon']),mf,openmc.EnergyFilter([0,1,10])]
        t.nuclides=['total'];t.scores=['flux','heating']
        t._mean=np.arange(64.).reshape(32,1,2)
        t._std_dev=np.ones((32,1,2))
        with self.assertRaisesRegex(ValueError,'select a bin'):dataset_from_tally(t,'simple:coil')
        d=dataset_from_tally(t,'simple:coil',filter_bins={0:1,2:0})
        for i,index in enumerate(mf.bins):
            self.assertEqual(d.mean[tuple(j-1 for j in index)],t._mean[16+2*i,0,1])
        self.assertEqual(d.mesh.origin,(0.,0.,20.))
        self.assertEqual(d.units,'eV/source')
        with self.assertRaises(ValueError):dataset_from_tally(t,'simple:coil',filter_bins={0:5,2:0})


class DeviceTests(unittest.TestCase):
    def test_tandem_and_vns_use_shared_core(self):
        for filename in ('tandem_parametric_input.yaml','tandem_vns_parametric_input.yaml'):
            with self.subTest(filename=filename):
                p=ROOT/'input_files'/filename
                model=build_model(p.read_text(),p.name,24)
                self.assertEqual(model.manifest.device,'tandem')
                self.assertIn('EP', {c.group for c in model.components})
                self.assertEqual(set(model.meshes['none']),set(model.meshes['quarter']))
                self.assertEqual(len({r.id for r in model.manifest.components}),len(model.components))
                self.assertEqual(len([c for c in model.components if c.name.startswith('hf_coil') and c.name.endswith('_magnet')]),4)
                with tempfile.TemporaryDirectory() as tmp:
                    f=Path(tmp)/'manifest.json';model.manifest.write(f)
                    self.assertEqual(ModelManifest.read(f).input_hash,model.manifest.input_hash)

    def test_simulation_manifest_contains_actual_ids_and_mesh_coordinates(self):
        import openmc
        from src.paratan.tallies.tandem_tallies import TallyBuilder
        cell=openmc.Cell(name='hf_coil_right_magnet',region=+openmc.ZCylinder(r=10)&-openmc.ZCylinder(r=20)&+openmc.ZPlane(z0=30)&-openmc.ZPlane(z0=40))
        builder=TallyBuilder()
        builder.add_descriptors([{'type':'hf_coil','location':'right','description':'magnet','cell':cell,
                                  'mesh_tallies':[{'scores':['heating'],'dimensions':[3,1,4]}]}])
        with tempfile.TemporaryDirectory() as tmp:
            m=export_simulation_manifest({'vacuum_vessel':{}},[cell],builder.get_viewer_bindings(),Path(tmp)/'manifest.json')
            record=m.components[0]
            self.assertEqual(record.cell_ids,[cell.id])
            self.assertEqual(record.tallies[0]['tally_id'],builder.get_tallies()[0].id)
            self.assertEqual(len(record.tallies[0]['mesh']['r']),4)
            self.assertEqual(record.tallies[0]['mesh']['origin'][2],30)


class ApplicationTests(unittest.TestCase):
    def test_tally_button_overlays_geometry_and_selects_its_own_mesh(self):
        app=self.app
        app.reset_view()
        app.select_component('hf_coil_right_magnet')
        app.set_section_mode('half')
        app.view_cylindrical_tally(1)
        self.assertEqual(app.state.tally_mesh_index,1)
        self.assertTrue(app.state.preview_tally)
        self.assertFalse(app.state.solo or app.state.demo_heating)
        self.assertEqual(app.state.section_mode,'quarter')
        self.assertFalse(app.scene.hide_body)
        self.assertTrue(app.body['quarter']['hf_coil_right_magnet'].GetPickable())
        self.assertIn('tally-grid',app.pl.renderer.actors)
        context=app.body['quarter']['hf_coil_left_magnet']
        self.assertAlmostEqual(context.GetProperty().GetOpacity(),1.)
        app.view_cylindrical_tally(0)
        self.assertEqual(app.state.tally_mesh_index,0)
        self.assertTrue(app.state.demo_heating)
        self.assertIn('tally-map',app.pl.renderer.actors)
        self.assertIn('invented',app.scene_state.tally_summary['source'])
        self.assertFalse(app.scene.hide_body)
        app.reset_view()
        self.assertAlmostEqual(app.body['none']['hf_coil_left_magnet'].GetProperty().GetOpacity(),1.)

    """Headless integration: current controls share state and retain actors."""
    @classmethod
    def setUpClass(cls):
        from src.paratan.viewer.app import ViewerApp
        cls.tmp = tempfile.TemporaryDirectory()
        path=Path(cls.tmp.name)/'input.yaml'
        path.write_text((ROOT/'input_files/simple_parametric_input_new.yaml').read_text())
        cls.app=ViewerApp(path,n_theta=24,server_name='viewer-architecture-test',render_mode='client')

    @classmethod
    def tearDownClass(cls):
        cls.app.pl.close()
        cls.tmp.cleanup()

    def setUp(self):
        self.app.reset_view()

    def test_model_install_stages_actors_without_intermediate_renders(self):
        from unittest.mock import patch
        app=self.app
        view=app.view
        app.view=None  # Publication is a separate, final step.
        try:
            with patch.object(app.pl,'render') as render:
                app._install(app.model)
                # PyVista's camera-position setter may draw the final scene,
                # but adding actors must not draw once per part/cutaway.
                self.assertLessEqual(render.call_count,1)
            self.assertTrue(all(app.scene.body[mode] for mode in ('none','half','quarter')))
        finally:
            app.view=view

    def test_upload_preserves_browser_actor_and_property_identity(self):
        import yaml
        app=self.app
        original=app.model
        name='hf_coil_right_magnet'
        actor=app.scene.body['none'][name]
        ids=(actor.__this__,actor.GetProperty().__this__,actor.GetMapper().__this__)
        changed=yaml.safe_load(yaml.safe_dump(original.data))
        changed.pop('end_cell',None)
        changed['ports']['ports_present']=False
        try:
            result=app.load_upload({'name':'short.yaml','content':yaml.safe_dump(changed)})
            self.assertEqual(result['phase'],'success')
            current=app.scene.body['none'][name]
            self.assertEqual(ids,(current.__this__,current.GetProperty().__this__,current.GetMapper().__this__))
            self.assertFalse(any(n.startswith('port_') or n.startswith('end_cell') for n in app.model.by_name))
            self.assertIs(current.GetMapper().GetInput(),app.model.meshes['none'][name])
        finally:
            app._install(original)

    def test_failed_upload_returns_error_and_allows_retry(self):
        app=self.app
        previous=app.model
        result=app.load_upload({'name':'invalid.yaml','content':'not a model'})
        self.assertEqual(result['phase'],'error')
        self.assertFalse(app.state.upload_busy)
        self.assertIs(app.model,previous)
        result=app.load_upload({'name':previous.name,'content':(ROOT/'input_files/simple_parametric_input_new.yaml').read_text()})
        self.assertEqual(result['phase'],'success')

    def test_browser_batch_updates_do_not_overwrite_each_other(self):
        app=self.app
        app.server.state.ready()
        app.select_component('hf_coil_right_magnet')
        app.server.state.flush()
        app.server.state.update({'demo_heating':True,'slice_phi_on':True})
        app.server.state.flush()
        self.assertTrue(app.scene_state.demo_heating)
        self.assertTrue(app.scene_state.slice_phi_on)
        self.assertIn('tally-slice-phi',app.pl.renderer.actors)
        app.load_upload([{'name':'tandem_parametric_input.yaml',
            'content':(ROOT/'input_files/tandem_parametric_input.yaml').read_bytes()}])
        app.server.state.flush()
        self.assertEqual(app.model.manifest.device,'tandem')
        # Restore the class fixture for the following cases.
        from src.paratan.viewer.model import build_model
        app._install(build_model((ROOT/'input_files/simple_parametric_input_new.yaml').read_text(),'input.yaml',24))
        app.server.state.flush()

    def test_invalid_browser_batch_is_atomic(self):
        app=self.app
        app.server.state.ready();app.server.state.flush()
        previous=app.scene_state.explode
        app.server.state.update({'explode':.9,'demo_opacity':4})
        app.server.state.flush()
        self.assertEqual(app.scene_state.explode,previous)
        self.assertEqual(app.scene_state.demo_opacity,1)
        self.assertIn('must be between',app.scene_state.action_status)

    def test_example_devices_share_the_loader(self):
        app=self.app
        app.load_example('tandem')
        self.assertEqual(app.model.manifest.device,'tandem')
        self.assertEqual(len(app.model.components),71)  # includes 12 end-plug port layers
        self.assertFalse(app.state.device_under_construction)
        self.assertEqual(len([c for c in app.model.components if c.group=='ports']),12)
        app.load_example('tandem_vns')
        self.assertTrue(any('Test region envelope' in c.label for c in app.model.components))
        app.load_example('simple')
        self.assertEqual(app.model.manifest.device,'simple_mirror')

    def test_upload_notifications_cover_loading_success_and_failure(self):
        from unittest.mock import patch
        app = self.app
        client_only = app.server.state.trame__client_only
        for key in ('upload_busy', 'upload_phase', 'upload_notice', 'load_status'):
            self.assertIn(key, client_only)
        build = build_model
        phases = []
        def observe_build(*args):
            phases.append((app.state.upload_phase, app.state.upload_busy, app.state.load_status))
            return build(*args)
        with patch('src.paratan.viewer.app.build_model', side_effect=observe_build):
            app.load_upload({'name': 'uploaded.yaml', 'content':
                (ROOT/'input_files/simple_parametric_input_new.yaml').read_text() + '\n# notification test\n'})
        self.assertEqual(phases[0][:2], ('loading', True))
        self.assertIn('uploaded.yaml', phases[0][2])
        self.assertEqual(app.state.upload_phase, 'success')
        self.assertFalse(app.state.upload_busy)
        self.assertTrue(app.state.upload_notice)
        self.assertIn('Loaded uploaded.yaml', app.state.load_status)
        loaded_model = app.model
        for content in ['', 'invalid: [', 'x' * 2_000_001]:
            app.load_upload({'name': 'bad.yaml', 'content': content})
            self.assertEqual(app.state.upload_phase, 'error')
            self.assertFalse(app.state.upload_busy)
            self.assertTrue(app.state.upload_notice)
            self.assertIs(app.model, loaded_model)
        app.load_upload(None)
        self.assertFalse(app.state.upload_busy)
        self.assertEqual(app.state.upload_phase, 'error')

    def test_repeat_upload_reuses_prepared_model(self):
        from unittest.mock import patch
        app = self.app
        content = (ROOT/'input_files/simple_parametric_input_new.yaml').read_text()
        app.load_upload({'name':'original.yaml','content':content})
        meshes = app.model.meshes
        with patch('src.paratan.viewer.app.build_model', side_effect=AssertionError('rebuilt cached geometry')):
            app.load_upload({'name':'renamed.yaml','content':content})
        self.assertIs(app.model.meshes, meshes)
        self.assertEqual(app.model.name, 'renamed.yaml')

    def test_reload_skips_transfer_to_outgoing_browser(self):
        import asyncio
        from types import SimpleNamespace
        from unittest.mock import Mock
        from src.paratan.viewer.app import ViewerApp
        reload = Mock()
        target = SimpleNamespace(_page_reload=reload, view=object(),
            server=SimpleNamespace(protocol=object()), _pending_scene=True, _pending_camera=True)
        async def run():
            ViewerApp._reload_client(target)
            self.assertFalse(target._pending_scene)
            self.assertFalse(target._pending_camera)
            await asyncio.sleep(0)
            reload.exec.assert_called_once_with()
        asyncio.run(run())

    def test_cutaways_are_preloaded_and_retained(self):
        app = self.app
        app.load_example('simple')
        self.assertTrue(app.state.cuts_ready)
        self.assertTrue(all(app.body[m] for m in ('none','half','quarter')))
        half = app.body['half']['hf_coil_right_magnet']
        app.set_section_mode('half')
        app.set_section_mode('quarter')
        app.set_section_mode('half')
        self.assertIs(app.body['half']['hf_coil_right_magnet'], half)

    def test_display_changes_keep_dynamic_mesh_and_body_actors(self):
        app=self.app
        app.select_component('hf_coil_right_magnet')
        app.state.demo_heating=True;app._refresh()
        actor=app.pl.renderer.actors['tally-map']
        mapper=actor.GetMapper()
        data=mapper.GetInput()
        body=app.body['none']['hf_coil_right_magnet']
        dataset=app.scene.active_dataset
        app.state.demo_opacity=.4;app.state.explode=.3;app._refresh()
        self.assertIs(app.pl.renderer.actors['tally-map'],actor)
        self.assertIs(app.body['none']['hf_coil_right_magnet'],body)
        self.assertIs(actor.GetMapper().GetInput(),data)
        self.assertIs(app.scene.active_dataset,dataset)
        self.assertAlmostEqual(actor.GetProperty().GetOpacity(),.4)
        app.state.slice_phi_on=True;app.state.slice_phi=90;app._refresh()
        self.assertIn('tally-slice-phi',app.pl.renderer.actors)
        app.aim_at_cut()
        app.reset_view()
        self.assertFalse(app.scene_state.hidden or app.scene_state.solo)
        self.assertIsNone(app.scene.active_dataset)
        self.assertNotIn('tally-map',app.pl.renderer.actors)

    def test_saved_view_restore_and_invalid_import_are_atomic(self):
        app=self.app
        app.select_component('hf_coil_right_magnet');app.toggle_solo()
        app.state.section_mode='quarter';app.state.slice_z_on=True;app.state.demo_heating=True
        app.state.preset_name='Round trip';app.save_preset()
        app.reset_view();app.restore_preset()
        self.assertEqual(app.scene_state.selected_name,'hf_coil_right_magnet')
        self.assertTrue(app.scene_state.solo and app.scene_state.slice_z_on and app.scene_state.demo_heating)
        self.assertEqual(app.scene_state.section_mode,'quarter')
        before=app.scene_state.snapshot()
        payload={'version':1,'input_hash':app.model.fingerprint,'views':{'bad':{'camera':[[0,0,1],[0,0,0],[0,1,0]],'demo_opacity':2}}}
        app.state.preset_import=json.dumps(payload);app.import_presets()
        after=app.scene_state.snapshot()
        # The status is allowed to change, all scene settings must be preserved.
        before.pop('action_status');after.pop('action_status')
        self.assertEqual(before,after)
        self.assertNotIn('bad',app.presets)

    def test_real_dataset_uses_same_slices_with_nonuniform_bins(self):
        app=self.app;app.select_component('hf_coil_right_magnet')
        mesh=CylindricalMesh([61,100,211],[0,.3,2*np.pi],[0,20,100],(0,0,327))
        d=TallyDataset('simple_mirror:hf_coil_right_magnet',mesh,np.arange(8.).reshape(2,2,2)+1,np.ones((2,2,2)))
        app.scene.datasets[d.id]=d;app.state.result_id=d.id;app.state.slice_phi_on=True;app.state.slice_phi=45
        app._refresh()
        self.assertEqual(app.state.result_source,'openmc')
        surface=app.pl.renderer.actors['tally-slice-phi'].GetMapper().GetInput()
        self.assertEqual(surface.GetNumberOfCells(),4)
        grid,_=app.scene._grid_for_selection()
        self.assertEqual(int(grid.locate_phi(.4)),1)
        self.assertEqual(grid.z[0],327)
        app.state.result_quantity='std_dev';app._refresh()
        self.assertEqual(app.scene.active_dataset.id,d.id)


def write_statepoint_fixture(path, tally_id, mesh, means, sigma=.2):
    """Small documented OpenMC HDF5 layout fixture; not a transport simulation."""
    import h5py
    from openmc.statepoint import _VERSION_STATEPOINT
    flat=np.asarray(means).ravel(order='F')
    n=10
    with h5py.File(path,'w') as f:
        f.attrs['filetype']=np.bytes_('statepoint')
        f.attrs['version']=np.array([_VERSION_STATEPOINT,0])
        f.attrs['tallies_present']=1
        tg=f.create_group('tallies');tg.attrs['n_tallies']=1;tg.attrs['ids']=[tally_id]
        mg=tg.create_group(f'meshes/mesh {mesh.id}')
        mg['type']=np.bytes_('cylindrical')
        for name in ('r_grid','phi_grid','z_grid','origin'):mg[name]=getattr(mesh,name)
        fg=tg.create_group('filters/filter 707');fg['type']=np.bytes_('mesh');fg['bins']=mesh.id
        t=tg.create_group(f'tally {tally_id}')
        t['name']=np.bytes_('HDF5 integration fixture');t['n_realizations']=n
        t['estimator']=np.bytes_('tracklength');t['n_filters']=1;t['filters']=[707]
        t['nuclides']=np.array(['total'],dtype='S');t['score_bins']=np.array(['heating'],dtype='S')
        results=np.stack([flat*n,n*(flat**2+(n-1)*sigma**2)],axis=-1)
        t['results']=results[:,None,:]


class StatepointIntegrationTests(unittest.TestCase):
    def test_hdf5_loading_manifest_binding_and_rendering(self):
        import openmc
        from src.paratan.tallies.tandem_tallies import TallyBuilder
        from src.paratan.viewer.app import ViewerApp
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);input_path=root/'input.yaml'
            input_path.write_text((ROOT/'input_files/simple_parametric_input_new.yaml').read_text())
            app=ViewerApp(input_path,n_theta=24,server_name='statepoint-integration',render_mode='client')
            try:
                app.select_component('hf_coil_right_magnet')
                profile=app._comp().profile
                ri,ro,zlo,zhi=profile.bounds
                cell=openmc.Cell(name='hf_coil_right_magnet',region=+openmc.ZCylinder(r=ri)&-openmc.ZCylinder(r=ro)&+openmc.ZPlane(z0=zlo)&-openmc.ZPlane(z0=zhi))
                tb=TallyBuilder();tb.add_descriptors([{'cell':cell,'type':'hf_coil','location':'right','description':'magnet',
                    'mesh_tallies':[{'scores':['heating'],'dimensions':[3,1,4]}]}])
                manifest=export_simulation_manifest(app.model.data,[cell],tb.get_viewer_bindings(),root/'simulation-manifest.json')
                tally=tb.get_tallies()[0];mesh=tally.filters[0].mesh
                values=np.arange(12.).reshape(3,1,4)+1
                write_statepoint_fixture(root/'statepoint.h5',tally.id,mesh,values)
                s=app.state;s.statepoint_path=str(root/'statepoint.h5');s.simulation_manifest_path=str(root/'simulation-manifest.json')
                s.statepoint_tally_id=tally.id;app.load_results()
                self.assertIn('Loaded tally',s.result_status)
                np.testing.assert_allclose(app.scene.active_dataset.mean,values)
                np.testing.assert_allclose(app.scene.active_dataset.std_dev,.2)
                self.assertEqual(s.result_source,'openmc')
                self.assertIn('tally-map',app.pl.renderer.actors)
                s.slice_z_on=True;s.result_quantity='relative_error';app._refresh()
                self.assertIn('tally-slice-z',app.pl.renderer.actors)
                # Mismatched model identity is rejected without replacing a valid result.
                active=s.result_id
                bad=manifest.to_dict();bad['input_hash']='wrong';(root/'bad.json').write_text(json.dumps(bad))
                s.simulation_manifest_path=str(root/'bad.json');app.load_results()
                self.assertIn('does not match',s.result_status)
                self.assertEqual(s.result_id,active)
            finally:app.pl.close()


class FactoryLifetimeTests(unittest.TestCase):
    def test_server_retains_application_and_device_actions(self):
        import gc
        from src.paratan.viewer.app import create_app
        server=create_app(ROOT/'input_files/simple_parametric_input_new.yaml',n_theta=24,
                          server_name='viewer-lifetime-regression',render_mode='client')
        try:
            gc.collect()
            server.state.ready()
            server.controller.load_example('tandem')
            server.state.flush()
            self.assertEqual(server._paratan_viewer.model.manifest.device,'tandem')
            self.assertEqual(server.state.n_components,71)
            self.assertFalse(server.state.device_under_construction)
        finally:server._paratan_viewer.pl.close()
