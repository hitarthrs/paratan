"""Application services and Trame bridge for the Paratan viewer platform.

Performance model
-----------------
* Cutaways are cached and their actors are preloaded on model load.
  Switching cuts updates transforms without generating or downloading meshes.
* A slice isolates the selected component: one small mesh is rebuilt per slider
  move instead of the whole machine.
"""
from __future__ import annotations

import asyncio
import base64
import hashlib
import io
import json
from collections import OrderedDict
from time import perf_counter
from pathlib import Path
from typing import Any

import numpy as np
import pyvista as pv
import yaml
from PIL import Image
from trame.app import get_server

from src.paratan.viewer.components import ComponentMesh, color_for_material
from src.paratan.viewer.model import Model, build_model, GROUPS, MODES, MAX_UPLOAD_BYTES, BODY_STYLE, EDGE_COLOR, SLICE_WEDGE, is_shown
from src.paratan.viewer.model import _set_shown
from src.paratan.viewer.ui import ViewerUI
from src.paratan.viewer.scene import SceneController
from src.paratan.viewer.state import SceneState, SceneActions, TrameStateBridge
from src.paratan.viewer.manifest import ModelManifest, component_id
from src.paratan.viewer.results import load_statepoint, TallyDataset
from src.paratan.viewer.presets import validate_view
from src.paratan.viewer.inspection import setup_example, tally_configuration
from src.paratan.viewer.labels import GROUP_LABELS, material_label, meta_row
from src.paratan.viewer.materials_catalog import MaterialCatalog, default_catalog

def _rgb_css(rgb: tuple[float, float, float]) -> str:
    r, g, b = (int(255 * c) for c in rgb)
    return f"rgb({r},{g},{b})"

class ViewerApp(ViewerUI):
    def __init__(self, input_yaml: str | Path, *, n_theta: int, server_name: str, render_mode: str):
        self.n_theta = n_theta
        self.server = get_server(name=server_name, client_type="vue3")
        self.scene_state = SceneState()
        self.actions = SceneActions(self.scene_state)
        self.state = TrameStateBridge(self.server.state, self.scene_state, self.actions)
        self.ctrl = self.server.controller
        self.presets: dict[str, dict] = {}
        self.view = None
        self._flush_scheduled = False
        self._pending_camera = False
        self._pending_scene = False
        self.render_mode = render_mode
        self._cut_warming = False
        self._upload_started = None
        self._model_cache = OrderedDict()
        self.materials_catalog: MaterialCatalog = default_catalog()

        pl = self.pl = pv.Plotter(notebook=False, off_screen=True)
        pl.set_background("white")
        pl.add_axes()
        pl.remove_all_lights()
        pl.add_light(pv.Light(position=(1.0, -0.35, 0.85), light_type="scene light", intensity=0.9))
        pl.add_light(pv.Light(position=(-0.7, 0.6, 0.25), light_type="scene light", intensity=0.4))
        pl.enable_mesh_picking(callback=self._on_pick, left_clicking=True, use_actor=True,
                               show=False, show_message=False)

        self.scene = SceneController(self.scene_state, pl, n_theta, self._push)
        path = Path(input_yaml)
        # Saved views are persisted next to the file the viewer was started with only.
        self.preset_owner = path.name
        self.preset_path = path.resolve().with_suffix(".viewer.json")
        self._init_state()
        self._install(build_model(path.read_text(), path.name, n_theta))
        self._bind()
        self._build_ui(render_mode)
        # _install runs before plotter_ui exists; publish once the view is wired.
        self._refresh()
        self._push(geometry=True, camera=True)
        self.ctrl.on_server_bind.add(self._bind_upload_http)

    @property
    def model(self):
        return self.scene.model

    @model.setter
    def model(self, value):
        self.scene.model = value

    @property
    def hidden(self):
        return self.scene_state.hidden

    @hidden.setter
    def hidden(self, value):
        self.scene_state.hidden = value

    @property
    def opacity_override(self):
        return self.scene_state.opacity

    @opacity_override.setter
    def opacity_override(self, value):
        self.scene_state.opacity = value

    @property
    def body(self):
        return self.scene.body

    @property
    def edges(self):
        return self.scene.edges

    def _comp(self, *args, **kwargs):
        return self.scene._comp(*args, **kwargs)

    def _slice_comp(self, *args, **kwargs):
        return self.scene._slice_comp(*args, **kwargs)

    def _slicing(self, *args, **kwargs):
        return self.scene._slicing(*args, **kwargs)

    def _slice_section(self, *args, **kwargs):
        return self.scene._slice_section(*args, **kwargs)

    def _active_section(self, *args, **kwargs):
        return self.scene._active_section(*args, **kwargs)

    def _effective_opacity(self, *args, **kwargs):
        return self.scene._effective_opacity(*args, **kwargs)

    def _sync_scene(self, *args, **kwargs):
        result = self.scene._sync_scene(*args, **kwargs)
        self.state.publish_scene()
        return result

    def _refresh(self, *args, **kwargs):
        result = self.scene._refresh(*args, **kwargs)
        dataset = self.scene.active_dataset
        self.state.result_source = dataset.source if dataset is not None else ''
        self.state.result_summary = (f'{dataset.score} · {dataset.nuclide} · {dataset.units} · {dataset.mesh.shape}' if dataset else '')
        self.state.result_choices = [{'title': f'{d.score} · {d.nuclide} · tally {d.provenance.get("tally_id", "demo")}', 'value': d.id}
                                    for d in self.scene.datasets.values()
                                    if d.component_id == component_id(self.model.manifest.device, self.state.selected_name)]
        self.state.publish_scene()
        return result

    # ------------------------------------------------------------------ state
    def _init_state(self) -> None:
        s = self.state
        s.trame__title = "Paratan · Simple Mirror"
        for g in GROUPS:
            setattr(s, f"vis_{g}", True)
            setattr(s, f"items_{g}", [])
        s.update({
            "n_components": 0, "model_name": "", "load_status": "", "upload_file": None,
            "upload_busy": False, "upload_phase": "idle", "upload_notice": False,
            "upload_drag_over": False,
            "cuts_ready": False,
            "legend_items": [], "explode": 0.0, "panel": "geometry", "search": "",
            "selected_name": "", "selected_label": "", "has_selection": False, "selected_opacity": 1.0,
            "solo": False, "section_mode": "none",
            "has_slice": False, "slice_phi_on": False, "slice_phi": 0.0, "slice_z_on": False,
            "slice_z": 0.0, "slice_z_min": 0.0, "slice_z_max": 1.0, "slice_info": "",
            "tally_example": "", "inspector_title": "", "inspector_material": "", "inspector_group": "",
            "inspector_rows": [], "inspector_tallies": [], "inspector_path": "",
            "tally_meshes": [], "tally_mesh_index": 0, "preview_tally": False, "demo_heating": False,
            "heating_available": False, "demo_opacity": 1.0,
            "preset_name": "", "preset_selected": "", "preset_import": "", "preset_names": [],
            "action_status": "", "png_url": "", "presets_url": "",
            "device_label": "", "device_under_construction": False, "manifest_url": "", "result_source": "", "result_status": "",
            "result_summary": "", "result_choices": [], "statepoint_path": "", "simulation_manifest_path": "",
            "statepoint_tally_id": 0, "statepoint_score": "heating", "statepoint_nuclide": "total",
            "statepoint_filter_bins": "{}", "result_quantity": "mean", "result_id": "",
            "materials_panel_open": False, "materials_items": [], "materials_search": "",
            "materials_source": "", "materials_status": "", "materials_selected_key": "",
            "materials_detail_label": "", "materials_detail_density": "", "materials_detail_id": "—",
            "materials_composition": [],
        })
        # The browser may show immediate reading feedback, but must never send
        # stale loading text back over the server's final success/error state.
        s.client_only('upload_busy', 'upload_phase', 'upload_notice', 'load_status', 'upload_drag_over')

    # ------------------------------------------------------------------ model
    def _install(self, model: Model) -> None:
        """Swap in a model: rebuild the scene once and reset every view setting."""
        pl, s = self.pl, self.state
        self._model_cache[model.fingerprint] = model
        self._model_cache.move_to_end(model.fingerprint)
        while len(self._model_cache) > 2:
            self._model_cache.popitem(last=False)
        self.scene.install_model(model)
        self._cut_warming = False
        # The compact original meshes are cheap to preload together. All cut
        # actors are resident before interaction, so switching needs no meshes.
        for mode in MODES:
            self.scene._ensure_mode(mode)
        s.cuts_ready = True
        self.actions.install(model.by_name, GROUPS)
        for g in GROUPS:
            setattr(s, f"items_{g}", [{"name": c.name, "label": c.label, "material": material_label(c.material)}
                                      for c in model.components if c.group == g])
            setattr(s, f"vis_{g}", True)
        materials = sorted({c.material for c in model.components}, key=str.lower)
        s.legend_items = [{"name": material_label(m), "color": _rgb_css(color_for_material(m))} for m in materials]
        self._publish_model_materials(materials)
        self._clear_material_detail()
        s.n_components = len(model.components)
        s.model_name = model.name
        s.device_label = model.device_label
        s.device_under_construction = model.manifest.device == 'tandem'
        s.manifest_url = "data:application/json;base64," + base64.b64encode(json.dumps(model.manifest.to_dict()).encode()).decode()
        s.explode, s.section_mode = 0.0, "none"
        s.preview_tally = s.demo_heating = False
        self._clear_selection_state()
        self.presets = {}
        s.preset_selected = ""
        if model.name == self.preset_owner and self.preset_path.exists():
            try:
                saved = json.loads(self.preset_path.read_text())
                if saved.get("input_hash") == model.fingerprint:
                    self.presets = saved.get("views", {})
                else:
                    s.action_status = "Input changed; saved views were not loaded."
            except (ValueError, OSError):
                s.action_status = "Could not read saved views."
        self._publish_presets(write=False)
        # Show the new actors before fitting: the camera fits visible props only.
        self._sync_scene()
        pl.view_isometric(render=False)
        pl.reset_camera(render=False)
        if self.view is not None:
            self._refresh()
            self._push(geometry=True, camera=True)

    # ------------------------------------------------------------ view logic
    def _push(self, *, geometry: bool = False, camera: bool = False) -> None:
        """Request a scene sync; requests made while handling one event coalesce into one.

        Camera-only changes avoid scene publication. Scene changes use Trame
        serialization caching; geometry stays resident where supported.
        """
        if self.view is None:
            return
        self._pending_camera = self._pending_camera or camera
        self._pending_scene = self._pending_scene or geometry or not camera
        if self._flush_scheduled:
            return
        try:
            loop = asyncio.get_running_loop()
        except RuntimeError:  # no event loop (scripts / tests): sync immediately
            self._flush()
            return
        self._flush_scheduled = True
        loop.call_soon(self._flush)

    def _flush(self) -> None:
        self._flush_scheduled = False
        camera, self._pending_camera = self._pending_camera, False
        if camera and hasattr(self.view, "push_camera"):
            self.view.push_camera()
        scene, self._pending_scene = self._pending_scene, False
        if scene:
            self.view.update()

    # ----------------------------------------------------------------- camera
    def _aim(self, azimuth: float | None, elevation: float, bounds=None) -> None:
        """Look at the scene from ``azimuth`` / ``elevation`` [deg], keeping the current distance."""
        pl = self.pl
        if bounds is not None:
            pl.reset_camera(bounds=bounds)
        focal = np.array(pl.camera.focal_point, float)
        offset = np.array(pl.camera.position, float) - focal
        dist = float(np.linalg.norm(offset)) or 1.0
        az = np.radians(azimuth) if azimuth is not None else np.arctan2(offset[1], offset[0])
        el = np.radians(elevation)
        pl.camera.position = tuple(focal + dist * np.array([np.cos(az) * np.cos(el), np.sin(az) * np.cos(el), np.sin(el)]))
        pl.camera.up = (0.0, 0.0, 1.0)
        self._push(camera=True)

    def aim_at_cut(self) -> None:
        comp, s = self._slice_comp(), self.state
        if comp is not None and s.slice_phi_on:
            self._aim(float(s.slice_phi) + np.degrees(0.5 * SLICE_WEDGE), 30.0, bounds=comp.mesh.bounds)
        elif comp is not None and s.slice_z_on:
            self._aim(None, 55.0, bounds=comp.mesh.bounds)
        elif s.section_mode != "none":
            width = np.pi if s.section_mode == "half" else 0.5 * np.pi
            self._aim(float(np.degrees(0.5 * width)), 30.0)

    def cam_side(self) -> None:
        # Look squarely at the y=0 cut face, like the engineering reference.
        self._aim(90.0, 0.0)
        self.pl.enable_parallel_projection()
        self.pl.reset_camera()
        self._push(camera=True)

    def cam_end(self) -> None:
        self.pl.view_xy()
        self.pl.reset_camera()
        self._push(camera=True)

    def cam_iso(self) -> None:
        self.pl.view_isometric()
        self.pl.reset_camera()
        self._push(camera=True)

    def cam_fit_selected(self) -> None:
        comp = self._comp()
        if comp is None:
            self.reset_view()
            return
        self.pl.reset_camera(bounds=comp.mesh.bounds)
        self._push(camera=True)

    # -------------------------------------------------------------- selection
    def _clear_selection_state(self) -> None:
        self.actions.clear()
        self.state.publish_scene()
        self.state.update({
            "selected_name": "", "selected_label": "", "selected_opacity": 1.0, "has_selection": False,
            "solo": False, "tally_example": "", "tally_meshes": [], "inspector_rows": [],
            "inspector_tallies": [], "inspector_path": "", "has_slice": False,
            "slice_phi_on": False, "slice_z_on": False, "heating_available": False, "slice_info": ""})

    def clear_selection(self) -> None:
        self._clear_selection_state()
        self._refresh()

    def select_component(self, name: str) -> None:
        comp = self._comp(name)
        if comp is None:
            return
        s = self.state
        if name != s.selected_name:
            s.preview_tally = s.demo_heating = False
        self.actions.select(name)
        s.publish_scene()
        s.selected_name, s.selected_label = name, comp.label
        s.selected_opacity = self._effective_opacity(name)
        s.has_selection = True
        path, entries = tally_configuration(self.model.data, comp)
        s.inspector_title = comp.label
        s.inspector_material = material_label(comp.material)
        s.inspector_group = GROUP_LABELS[comp.group]
        s.inspector_path = path if entries or "." in path else ""
        s.inspector_rows = [meta_row(k, v) for k, v in comp.meta.items() if not k.startswith("_")]
        s.inspector_tallies = [
            {"index": i, "title": f"{'Spatial mesh' if e['kind'] == 'mesh_tallies' else 'Cell total'} {i + 1}",
             "scores": e.get("scores", []),
             "filters": [f.replace("_filter", "").replace("_", " ").capitalize() for f in e.get("filters", [])]
                        or ["All particles / energies"],
             "nuclides": ", ".join(e.get("nuclides", [])) or "All nuclides",
             "dimensions": e.get("dimensions", []) if e["kind"] == "mesh_tallies" else []}
            for i, e in enumerate(entries)]
        s.tally_example = setup_example(comp)
        s.tally_meshes = [{"title": f"Mesh {i + 1}: {' × '.join(str(d) for d in e.get('dimensions', []))}  (r × φ × z)",
                           "value": i} for i, e in enumerate(entries) if e["kind"] == "mesh_tallies"]
        s.tally_mesh_index = s.tally_meshes[0]["value"] if s.tally_meshes else 0
        s.has_slice = comp.profile is not None
        if comp.profile is not None:
            zmin, zmax = comp.profile.bounds[2:]
            s.slice_z_min, s.slice_z_max = round(zmin, 2), round(zmax, 2)
            if not zmin <= float(s.slice_z) <= zmax:
                s.slice_z = round(0.5 * (zmin + zmax), 2)
        else:
            s.slice_phi_on = s.slice_z_on = False
        self._refresh()

    def _on_pick(self, picked: Any = None, *args: Any, **_kwargs: Any) -> None:
        for candidate in (picked, *args):
            key = getattr(candidate, "name", None)
            if isinstance(key, str) and key.split("|", 1)[0] in MODES and "|" in key:
                self.select_component(key.split("|", 1)[1])
                return

    def _on_view_click(self, event: Any = None, **_kwargs: Any) -> None:
        """Browser-side pick (local rendering): map the clicked prop's scene id to a component."""
        picks = event if isinstance(event, list) else [event]
        ids = {str(p.get("remoteId")) for p in picks if isinstance(p, dict) and p.get("remoteId") is not None}
        if not ids or self.view is None:
            return
        for name, actor in self.body[self.state.section_mode].items():
            if str(self.view.get_scene_object_id(actor)) in ids:
                self.select_component(name)
                return

    def hide_selected(self) -> None:
        if self.state.selected_name:
            self.actions.hide_selected()
            self._refresh()

    def show_all(self) -> None:
        self.hidden.clear()
        for g in GROUPS:
            setattr(self.state, f"vis_{g}", True)
        self.state.solo = False
        self._refresh()

    def toggle_solo(self) -> None:
        if self.state.selected_name:
            self.state.solo = not bool(self.state.solo)
            self._refresh()

    def reset_view(self) -> None:
        s = self.state
        self.actions.reset()
        s.publish_scene()
        s.demo_heating = s.preview_tally = False
        s.explode, s.section_mode = 0.0, "none"
        self.hidden.clear()
        self.opacity_override.clear()
        for g in GROUPS:
            setattr(s, f"vis_{g}", True)
        self._clear_selection_state()
        self._refresh()
        self.pl.view_isometric()
        self.pl.reset_camera()
        self._push(camera=True)

    # -------------------------------------------------------------- cut / slice
    def view_cylindrical_tally(self, index=None) -> None:
        """Overlay the chosen mesh on a quarter cut of the geometry."""
        if not self.state.tally_meshes:
            return
        s = self.state
        chosen = int(s.tally_mesh_index if index is None else index)
        if chosen not in {item['value'] for item in s.tally_meshes}:
            return
        if chosen != s.tally_mesh_index:
            s.result_id = ''
        s.tally_mesh_index = chosen
        s.panel = 'inspect'
        s.preview_tally = True
        _, entries = tally_configuration(self.model.data, self._comp())
        s.demo_heating = not bool(s.result_id) and 'heating' in entries[chosen].get('scores', [])
        s.solo = False
        s.slice_phi_on = s.slice_z_on = False
        s.section_mode = 'quarter'
        self._refresh()
        self.aim_at_cut()

    def view_hf_tally(self) -> None:
        name = next((c.name for c in self.model.components
                     if c.group == 'HF' and 'magnet' in c.name and c.profile is not None), None)
        if name:
            self.select_component(name)
            self.view_cylindrical_tally()

    def set_section_mode(self, mode: str) -> None:
        self.state.section_mode = mode if mode in MODES else "none"
        self._refresh()
        self.aim_at_cut()

    def step_slice(self, axis: str, direction: int) -> None:
        """Move a slice to the next / previous bin centre of the selected component's grid."""
        grid = self.scene._grid_for_selection()[0]
        if grid is None:
            return
        s, forward = self.state, int(direction) > 0
        if axis == "phi":
            # Bin centres in degrees, repeated one turn either side so stepping wraps around.
            centres = np.degrees([grid.phi_centre(j) for j in range(grid.shape[1])])
            centres = np.concatenate([centres - 360.0, centres, centres + 360.0])
            current, tol = float(s.slice_phi), 1e-2 * 360.0 / grid.shape[1]
        else:
            centres = np.array([grid.z_centre(k) for k in range(grid.shape[2])])
            current, tol = float(s.slice_z), 1e-2 * (grid.z[-1] - grid.z[0]) / grid.shape[2]
        # Strictly the next centre in the pressed direction (the stored value is rounded).
        ahead = centres[centres > current + tol] if forward else centres[centres < current - tol]
        target = (ahead.min() if forward else ahead.max()) if len(ahead) else (centres.max() if forward else centres.min())
        if axis == "phi":
            s.slice_phi = round(float(target) % 360.0, 3)
        else:
            s.slice_z = round(float(target), 3)

    # --------------------------------------------------------- views & export
    def _publish_presets(self, *, write: bool = True) -> None:
        payload = json.dumps({"version": 1, "input_hash": self.model.fingerprint, "views": self.presets}, indent=2)
        if write and self.model.name == self.preset_owner:
            self.preset_path.write_text(payload)
        self.state.preset_names = sorted(self.presets)
        self.state.presets_url = ("data:application/json;base64," + base64.b64encode(payload.encode()).decode()
                                  if self.presets else "")

    def save_preset(self) -> None:
        s, pl = self.state, self.pl
        name = str(s.preset_name or "").strip()
        if not name:
            s.action_status = "Enter a name for this view."
            return
        self.presets[name] = {
            "camera": [list(v) for v in pl.camera_position],
            "parallel_projection": bool(pl.camera.parallel_projection), "parallel_scale": pl.camera.parallel_scale,
            "view_angle": pl.camera.view_angle, "section": s.section_mode, "explode": float(s.explode),
            "slice": {"phi_on": bool(s.slice_phi_on), "phi": float(s.slice_phi),
                      "z_on": bool(s.slice_z_on), "z": float(s.slice_z)},
            "groups": {g: bool(getattr(s, f"vis_{g}")) for g in GROUPS},
            "hidden": sorted(self.hidden), "opacity": dict(self.opacity_override),
            "selected": s.selected_name, "solo": bool(s.solo), "preview_tally": bool(s.preview_tally),
            "tally_mesh_index": int(s.tally_mesh_index), "demo_heating": bool(s.demo_heating),
            "demo_opacity": float(s.demo_opacity), "result_id": s.result_id, "result_quantity": s.result_quantity}
        try:
            self._publish_presets()
            s.preset_selected = name
            s.action_status = f"Saved “{name}”."
        except OSError as exc:
            s.action_status = f"Could not save view: {exc}"

    def restore_preset(self) -> None:
        s = self.state
        saved = self.presets.get(s.preset_selected)
        if not saved:
            return
        try:  # validate before touching the scene
            saved = validate_view(saved, self.model.by_name)
            camera = np.asarray(saved["camera"], dtype=float)
            if camera.shape != (3, 3) or not np.isfinite(camera).all():
                raise ValueError("Invalid camera")
            mode = saved.get("section", "none")  # views from the old half-space cuts open uncut
            explode = float(saved.get("explode", 0))
            sl = saved.get("slice") or {}
            phi, z = float(sl.get("phi", 0.0)), float(sl.get("z", 0.0))
            overrides = {n: float(v) for n, v in saved.get("opacity", {}).items() if n in self.model.by_name}
            scale = float(saved.get("parallel_scale", 1))
            if mode not in MODES or not 0 <= explode <= 1 or scale <= 0 or not np.isfinite([scale, phi, z]).all():
                raise ValueError("Invalid view settings")
            if any(not 0 <= v <= 1 for v in overrides.values()):
                raise ValueError("Invalid opacity")
            self.hidden = {n for n in saved.get("hidden", []) if n in self.model.by_name}
            self.opacity_override = overrides
            for g in GROUPS:
                setattr(s, f"vis_{g}", bool(saved.get("groups", {}).get(g, True)))
            s.section_mode, s.explode = mode, explode
            self._clear_selection_state()
            if saved.get("selected") in self.model.by_name:
                self.select_component(saved["selected"])
            if s.has_slice:
                s.slice_phi, s.slice_z = phi, z
                s.slice_phi_on, s.slice_z_on = bool(sl.get("phi_on")), bool(sl.get("z_on"))
            s.solo = bool(saved.get("solo", False))
            s.preview_tally = bool(saved.get("preview_tally", False))
            s.tally_mesh_index = int(saved.get("tally_mesh_index", 0))
            s.demo_heating = bool(saved.get('demo_heating', False))
            s.demo_opacity = float(saved.get('demo_opacity', 1))
            saved_result = saved.get('result_id', '')
            s.result_id = saved_result if saved_result in self.scene.datasets else ''
            s.result_quantity = saved.get('result_quantity', 'mean')
            self._refresh()
            self.pl.camera_position = camera.tolist()
            self.pl.camera.parallel_projection = bool(saved.get("parallel_projection", False))
            self.pl.camera.parallel_scale = scale
            self.pl.camera.view_angle = float(saved.get("view_angle", 30))
            self._push(camera=True)
            s.action_status = f"Restored “{s.preset_selected}”."
        except (ValueError, TypeError, KeyError) as exc:
            s.action_status = f"Could not restore view: {exc}"

    def delete_preset(self) -> None:
        self.presets.pop(self.state.preset_selected, None)
        try:
            self._publish_presets()
            self.state.preset_selected = ""
            self.state.action_status = "View deleted."
        except OSError as exc:
            self.state.action_status = f"Could not save views: {exc}"

    def import_presets(self) -> None:
        s = self.state
        try:
            payload = json.loads(s.preset_import)
            if payload.get("version") != 1 or payload.get("input_hash") != self.model.fingerprint:
                raise ValueError("Views must belong to this input configuration")
            views = payload["views"]
            if not isinstance(views, dict) or any(not isinstance(v, dict) for v in views.values()):
                raise ValueError("Invalid views")
            validated = {name: validate_view(view, self.model.by_name) for name,view in views.items()}
            self.presets.update(validated)
            self._publish_presets()
            s.action_status = "Views imported."
        except (ValueError, KeyError, TypeError, OSError) as exc:
            s.action_status = f"Import failed: {exc}"

    def export_png(self) -> None:
        try:
            self.pl.render()
            pixels = self.pl.screenshot(return_img=True, scale=2)
            out = io.BytesIO()
            Image.fromarray(pixels).save(out, format="PNG")
            self.state.png_url = "data:image/png;base64," + base64.b64encode(out.getvalue()).decode()
            self.state.action_status = "Image ready; use Download PNG."
        except Exception as exc:  # noqa: BLE001 - surfaced to the user
            self.state.action_status = f"Image export failed: {exc}"

    def _sync_camera(self, camera, **_kwargs) -> None:
        if isinstance(camera, dict):
            self.scene_state.camera = dict(camera)
            pl = self.pl
            pl.camera_position = [camera["position"], camera["focalPoint"], camera["viewUp"]]
            pl.camera.parallel_projection = camera.get("parallelProjection", False)
            pl.camera.parallel_scale = camera.get("parallelScale", pl.camera.parallel_scale)
            pl.camera.view_angle = camera.get("viewAngle", pl.camera.view_angle)

    # ------------------------------------------------------------------ upload
    def _bind_upload_http(self, web_server) -> None:
        """Transfer files outside the scene websocket, with a definite response."""
        from aiohttp import web

        async def upload_yaml(request):
            content = bytearray()
            while chunk := await request.content.read(64 * 1024):
                content.extend(chunk)
                if len(content) > MAX_UPLOAD_BYTES:
                    return web.json_response({'phase':'error','message':'Input file is larger than 2 MB.'}, status=413)
            with self.server.state:
                result = self.load_upload({'name':request.query.get('name','input.yaml'), 'content':content})
            return web.json_response(result)

        web_server.app.router.add_post('/paratan-upload', upload_yaml)

    def _upload_status(self, phase: str, message: str) -> dict:
        self.state.upload_phase = phase
        self.state.upload_busy = phase == 'loading'
        self.state.load_status = message
        self.state.upload_notice = phase in ('success', 'error')
        self.state.flush()
        return {'phase': phase, 'message': message}

    def load_example(self, kind: str) -> None:
        files = {'simple': 'simple_parametric_input_new.yaml', 'tandem': 'tandem_parametric_input.yaml',
                 'tandem_vns': 'tandem_vns_parametric_input.yaml'}
        if kind not in files:
            self._upload_status('error', 'Unknown example device.')
            return
        path = Path(__file__).resolve().parents[3] / 'input_files' / files[kind]
        try:
            self.load_upload({'name': path.name, 'content': path.read_bytes()})
        except OSError as exc:
            self._upload_status('error', f'Could not read example: {exc}')

    def load_upload(self, upload) -> dict:
        s = self.state
        if not upload:
            return self._upload_status('error', 'No input file was received. Please choose a YAML file.')
        item = upload[0] if isinstance(upload, (list, tuple)) else upload
        if not isinstance(item, dict):
            return self._upload_status('error', 'Could not read the selected input file.')
        name, content = str(item.get("name", "input.yaml")), item.get("content")
        if isinstance(content, str):
            content = content.encode()
        if not isinstance(content, (bytes, bytearray)) or not content:
            return self._upload_status('error', f"Could not read {name}: the file is empty or unreadable.")
        if len(content) > MAX_UPLOAD_BYTES:
            return self._upload_status('error', f"{name} is larger than {MAX_UPLOAD_BYTES // 1_000_000} MB.")
        self._upload_status('loading', f"Loading {name} · building geometry…")
        self._upload_started = perf_counter()
        try:
            text = bytes(content).decode('utf-8')
            fingerprint = hashlib.sha256(text.encode()).hexdigest()
            model = self._model_cache.get(fingerprint)
            if model is None:
                model = build_model(text, name, self.n_theta)
            else:
                from dataclasses import replace
                model = replace(model, name=name)
        except KeyError as exc:
            self._upload_started = None
            return self._upload_status('error', f"Could not load {name}: missing input entry {exc}.")
        except Exception as exc:  # noqa: BLE001 - any malformed input is reported, not raised
            self._upload_started = None
            return self._upload_status('error', f"Could not load {name}: {exc}")
        try:
            self._install(model)
        except Exception as exc:  # noqa: BLE001 - make renderer errors visible too
            self._upload_started = None
            return self._upload_status('error', f"Could not display {name}: {exc}")
        elapsed = perf_counter() - self._upload_started
        self._upload_started = None
        self._upload_status('success', f"Loaded {name} · {len(model.components)} parts · prepared in {elapsed:.1f}s")
        s.upload_file = None
        return {'phase': s.upload_phase, 'message': s.load_status}

    def _reload_client(self) -> None:
        """Have the browser resubscribe to the scene from scratch after a model swap.

        trame-vtk's browser cache does not share in-flight array requests, so streaming a
        whole new model into a live view re-requests each array many times (measured ~58x
        after an upload). Over a slow link parts stay missing until that storm drains, so
        they appear to vanish. A fresh page subscription fetches each array about once.
        """
        reload = getattr(self, "_page_reload", None)
        if reload is None or self.view is None or not getattr(self.server, "protocol", None):
            return
        try:
            loop = asyncio.get_running_loop()
        except RuntimeError:  # scripts / tests: nothing to reload
            return
        # Do not publish the new mesh to the outgoing subscription: that starts
        # thousands of array requests just before we disconnect it. The fresh
        # subscription reads the installed render window directly.
        self._pending_scene = False
        self._pending_camera = False
        loop.call_soon(reload.exec)

    def load_results(self) -> None:
        """Load an explicitly bound cylindrical tally; never infer IDs from labels."""
        s = self.state
        comp = self._slice_comp()
        if comp is None:
            s.result_status = 'Select an axisymmetric component first.'
            return
        try:
            manifest = ModelManifest.read(s.simulation_manifest_path)
            if manifest.device != self.model.manifest.device or manifest.input_hash != self.model.manifest.input_hash:
                raise ValueError('Simulation manifest does not match this device input')
            # HF casing unions are display pieces; results must map to a named simulation component.
            record = manifest.find(comp.name)
            tally_id = int(s.statepoint_tally_id)
            if record is None or not any(t.get('tally_id') == tally_id and t.get('kind') == 'mesh_tallies' for t in record.tallies):
                raise ValueError('That tally is not bound to the selected component in the simulation manifest')
            selections = json.loads(s.statepoint_filter_bins or '{}')
            if not isinstance(selections, dict):
                raise ValueError('Filter selections must be a JSON object of filter index to bin index')
            dataset = load_statepoint(s.statepoint_path, component_id(manifest.device, comp.name), tally_id,
                                     score=s.statepoint_score, nuclide=s.statepoint_nuclide,
                                     filter_bins={int(k): int(v) for k,v in selections.items()})
            if not np.allclose([dataset.mesh.phi[0], dataset.mesh.phi[-1]], [0, 2*np.pi]):
                raise ValueError('Only full-azimuth cylindrical meshes are supported in this first reader')
            if not np.allclose(dataset.mesh.origin[:2], [0, 0]):
                raise ValueError('Off-axis cylindrical meshes are not supported yet')
            expected_meshes = [t.get('mesh') for t in record.tallies if t.get('tally_id') == tally_id]
            if not expected_meshes or not expected_meshes[0]:
                raise ValueError('Simulation manifest is missing the tally mesh coordinates')
            expected = expected_meshes[0]
            for key in ('r', 'phi', 'z'):
                actual = getattr(dataset.mesh, key)
                if actual.shape != np.shape(expected[key]) or not np.allclose(actual, expected[key]):
                    raise ValueError('Statepoint mesh coordinates do not match the simulation manifest')
            if not np.allclose(dataset.mesh.origin, expected['origin']):
                raise ValueError('Statepoint mesh origin does not match the simulation manifest')
            self.scene.datasets[dataset.id] = dataset
            s.demo_heating = False
            s.result_id, s.result_quantity = dataset.id, 'mean'
            s.has_slice = True
            s.slice_z_min = float(dataset.mesh.z[0] + dataset.mesh.origin[2])
            s.slice_z_max = float(dataset.mesh.z[-1] + dataset.mesh.origin[2])
            s.slice_z = .5 * (s.slice_z_min + s.slice_z_max)
            self._refresh()
            s.result_status = f'Loaded tally {tally_id}: {dataset.mesh.shape} bins, {dataset.units}'
        except (ValueError, TypeError, KeyError, OSError, ImportError) as exc:
            s.result_status = f'Could not load results: {exc}'

    def clear_results(self) -> None:
        self.state.result_id = ''
        self.state.demo_heating = False
        self._refresh()
        self.state.result_status = ''

    # ---------------------------------------------------------------- bindings
    def _bind(self) -> None:
        state, ctrl = self.state, self.ctrl
        for name in ("reset_view", "clear_selection", "hide_selected", "toggle_solo", "cam_side", "cam_end",
                     "cam_iso", "cam_fit_selected", "set_section_mode", "step_slice", "save_preset",
                     "restore_preset", "delete_preset", "import_presets", "export_png", "select_component",
                     "load_results", "clear_results", "load_example", "view_hf_tally", "view_cylindrical_tally",
                     "toggle_materials_panel", "select_material"):
            setattr(ctrl, name, getattr(self, name))
        ctrl.show_all_components = self.show_all
        ctrl.face_slice = self.aim_at_cut

        @ctrl.add("on_server_ready")
        def _ready(**_kwargs: Any) -> None:
            self._refresh()
            self.pl.reset_camera()
            # geometry=True: initial install runs before plotter_ui exists, so the
            # first client connect must publish the scene (camera-only leaves a blank canvas).
            self._push(geometry=True, camera=True)

        refresh = lambda **_kw: self._refresh()  # noqa: E731
        for g in GROUPS:
            state.change(f"vis_{g}")(refresh)
        state.change("explode", "solo", "slice_phi", "slice_z", "preview_tally", "tally_mesh_index",
                     "demo_heating", "demo_opacity", "result_id", "result_quantity")(refresh)
        state.change("selected_opacity")(self._on_opacity)
        state.change("slice_phi_on", "slice_z_on")(self._on_slice_toggle)
        # File bytes are passed through a trigger attachment, not shared state.
        # Trame filters binary content out of synchronized file metadata.

    def _on_opacity(self, **_kwargs: Any) -> None:
        if self.state.selected_name:
            self.opacity_override[self.state.selected_name] = float(self.state.selected_opacity)
            self._refresh()

    def _on_slice_toggle(self, **_kwargs: Any) -> None:
        self._refresh()
        if self.state.slice_phi_on or self.state.slice_z_on:
            self.aim_at_cut()

    # -------------------------------------------------------------- materials
    def toggle_materials_panel(self) -> None:
        s = self.state
        opening = not bool(s.materials_panel_open)
        s.materials_panel_open = opening
        if opening:
            # List = materials used by the loaded model; composition only after a click.
            if self.model is not None:
                materials = sorted({c.material for c in self.model.components}, key=str.lower)
                self._publish_model_materials(materials)
            self._clear_material_detail()
            s.materials_status = (
                f"{len(s.materials_items)} in this model · click for OpenMC composition"
            )

    def select_material(self, key: str) -> None:
        """Show OpenMC composition for a model-legend material (resolved via catalog)."""
        s = self.state
        raw = str(key)
        s.materials_selected_key = raw
        s.materials_detail_label = material_label(raw)
        s.materials_status = "Loading composition…"
        try:
            self.materials_catalog.ensure_loaded()
            s.materials_source = self.materials_catalog.source_label
            record = self.materials_catalog.resolve(raw) or self.materials_catalog.get(raw)
        except Exception as exc:  # noqa: BLE001
            self._clear_material_detail(keep_key=raw)
            s.materials_detail_label = material_label(raw)
            s.materials_status = f"Could not load composition: {exc}"
            return
        if record is None:
            self._clear_material_detail(keep_key=raw)
            s.materials_detail_label = material_label(raw)
            s.materials_detail_density = "not in materials library"
            s.materials_status = f"No OpenMC definition matched “{raw}”"
            return
        detail = record.to_detail()
        dens = detail["density"]
        units = detail["density_units"]
        dens_label = f"{dens:g} {units}".strip() if dens is not None else "density unavailable"
        s.materials_detail_label = detail["label"]
        s.materials_detail_density = dens_label
        s.materials_detail_id = detail["openmc_id"] if detail["openmc_id"] is not None else "—"
        s.materials_composition = detail["composition"]
        s.materials_status = ""

    def _publish_model_materials(self, materials: list[str]) -> None:
        """Legend-order list of materials present on the current device model."""
        s = self.state
        s.materials_items = [
            {
                "key": m,
                "label": material_label(m),
                "color": _rgb_css(color_for_material(m)),
            }
            for m in materials
        ]
        s.materials_source = "materials used in this model"

    def _clear_material_detail(self, *, keep_key: str = "") -> None:
        s = self.state
        s.materials_selected_key = keep_key
        if not keep_key:
            s.materials_detail_label = ""
        s.materials_detail_density = ""
        s.materials_detail_id = "—"
        s.materials_composition = []

    # ---------------------------------------------------------------------- UI

def create_app(input_yaml: str | Path, *, n_theta: int = 96, server_name: str = "paratan-viewer",
               render_mode: str = "trame") -> Any:
    """Build the Trame server + plotter for a supported device YAML.

    ``render_mode``: ``trame`` (default, local WebGL with remote fallback), ``client`` (browser WebGL
    only) or ``server`` (JPEG stream from Python/VTK).
    """
    pv.global_theme.trame.still_ratio = 1.25
    pv.global_theme.trame.interactive_ratio = 0.4
    app = ViewerApp(input_yaml, n_theta=n_theta, server_name=server_name, render_mode=render_mode)
    # Controller callbacks can be weak references; retain the application for
    # the full server lifetime, including after on_server_ready completes.
    app.server._paratan_viewer = app
    return app.server

def run_app(input_yaml: str | Path, *, host: str = "127.0.0.1", port: int = 8080, n_theta: int = 96,
            open_browser: bool = True, render_mode: str = "trame", timeout: int = 2) -> None:
    """Serve the viewer until the last browser client disconnects.

    ``timeout`` is the idle reap delay in seconds after no clients remain
    (and also before the first client connects). Use ``0`` to keep the
    process alive until Ctrl-C.
    """
    server = create_app(input_yaml, n_theta=n_theta, render_mode=render_mode)
    server.start(
        host=host,
        port=port,
        open_browser=open_browser,
        exec_mode="main",
        timeout=timeout,
    )
