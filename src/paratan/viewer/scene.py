"""Persistent VTK scene controller. No Trame imports or widget callbacks."""
from __future__ import annotations
import numpy as np
import pyvista as pv
from src.paratan.viewer.model import MODES, MODE_SECTIONS, SLICE_WEDGE, GROUP_OFFSETS, OVERLAYS, BODY_STYLE, EDGE_COLOR, HEAT, SCALAR_BAR, _set_shown
from src.paratan.viewer.cyl_grid import CylGrid, demo_heating_field
from src.paratan.viewer.inspection import tally_configuration
from src.paratan.viewer.revolve import Section
from src.paratan.viewer.components import ComponentMesh
from src.paratan.viewer.results import synthetic_dataset
from src.paratan.viewer.manifest import component_id
from src.paratan.viewer.tally_overlay import describe_overlay, lift_boundary

class SceneController:
    def __init__(self, state, plotter, n_theta, push):
        self.state, self.pl, self.n_theta = state, plotter, n_theta
        self._push = push
        self.model = None
        self.body = {m: {} for m in MODES}
        self.edges = {m: {} for m in MODES}
        # Browser synchronization identifies VTK objects by native address.
        # Keep actors/properties/mappers stable across model replacements.
        self._body_pool = {m: {} for m in MODES}
        self._edge_pool = {m: {} for m in MODES}
        self.hide_body = False
        self.had_dynamic = False
        self._geometry_key = None
        self.datasets = {}
        self._demo_cache = {}
        self.active_dataset = None
        self.port_openings = None

    def install_model(self, model):
        pl = self.pl
        for registry in (*self.body.values(), *self.edges.values()):
            for actor in registry.values():
                pl.remove_actor(actor, reset_camera=False, render=False)
            registry.clear()
        self.model = model
        self.port_openings = None
        if model.device_label == 'Simple mirror' and (model.data.get('ports') or {}).get('ports_present'):
            from src.paratan.viewer.port_openings import PortOpenings
            self.port_openings = PortOpenings(model.data)
        self._geometry_key = None
        self.datasets.clear()
        self._demo_cache.clear()
        self.active_dataset = None
        self.state.hidden.clear()
        self.state.opacity.clear()
        self._ensure_mode('none')

    def _ensure_mode(self, mode):
        """Send a cutaway only on its first use, then retain its actors."""
        pl, model = self.pl, self.model
        for name, mesh in model.meshes[mode].items():
            if name in self.body[mode]:
                continue
            comp = model.by_name[name]
            if comp.group == 'ports':
                if not self.state.vis_ports:
                    continue
            key = f"{mode}|{name}"
            actor = self._body_pool[mode].get(name)
            if actor is None:
                actor = pl.add_mesh(mesh, name=key, color=comp.color, opacity=comp.opacity, pickable=True,
                                    render=False, reset_camera=False, **BODY_STYLE)
                self._body_pool[mode][name] = actor
            else:
                actor.GetMapper().SetInputData(mesh)
                actor.GetProperty().SetColor(*comp.color)
                actor.GetProperty().SetOpacity(comp.opacity)
                pl.add_actor(actor, name=key, reset_camera=False, render=False)
            actor.mapper.scalar_visibility = False  # solid material colour, never data-driven
            _set_shown(actor, False)
            self.body[mode][name] = actor
            lines = model.outlines[mode].get(name)
            if lines is not None:
                edge = self._edge_pool[mode].get(name)
                if edge is None:
                    edge = pl.add_mesh(lines, name=f"edge|{key}", color=EDGE_COLOR,
                                       opacity=min(1.0, 0.25 + 0.5 * comp.opacity), line_width=1.0,
                                       lighting=False, pickable=False, render=False, reset_camera=False)
                    self._edge_pool[mode][name] = edge
                else:
                    edge.GetMapper().SetInputData(lines)
                    edge.GetProperty().SetOpacity(min(1.0, 0.25 + 0.5 * comp.opacity))
                    pl.add_actor(edge, name=f"edge|{key}", reset_camera=False, render=False)
                _set_shown(edge, False)
                self.edges[mode][name] = edge

    @property
    def hidden(self):
        return self.state.hidden

    @property
    def opacity_override(self):
        return self.state.opacity

    def _comp(self, name: str | None = None) -> ComponentMesh | None:
        return self.model.by_name.get(self.state.selected_name if name is None else name)

    def _slice_comp(self) -> ComponentMesh | None:
        comp = self._comp()
        return comp if comp is not None and comp.profile is not None else None

    def _slicing(self) -> bool:
        return self._slice_comp() is not None and bool(self.state.slice_phi_on or self.state.slice_z_on)

    def _slice_section(self) -> Section:
        phi, z = (0.0, 2.0 * np.pi), (-np.inf, np.inf)
        if self.state.slice_phi_on:
            start = np.radians(float(self.state.slice_phi))
            phi = (start + SLICE_WEDGE, start + 2.0 * np.pi)
        if self.state.slice_z_on:
            z = (-np.inf, float(self.state.slice_z))
        return Section(phi, z)

    def _active_section(self) -> Section:
        return self._slice_section() if self._slicing() else MODE_SECTIONS[self.state.section_mode]

    def _effective_opacity(self, name: str) -> float:
        return float(self.opacity_override.get(name, self._comp(name).opacity))

    def _comp_visible(self, name: str) -> bool:
        s = self.state
        if not bool(getattr(s, f"vis_{self._comp(name).group}", True)) or name in self.hidden:
            return False
        if s.solo and s.selected_name:
            return name == s.selected_name
        return True

    def _offset(self, name: str) -> tuple[float, float, float]:
        scale = float(self.state.explode) * 80.0  # cm
        ox, oy, oz = GROUP_OFFSETS.get(self._comp(name).group, (0.0, 0.0, 0.0))
        return (ox * scale, oy * scale, oz * scale * (1.0 if self.model.centers[name][2] >= 0 else -1.0))

    def _overlay_actors(self):
        for key in OVERLAYS:
            actor = self.pl.renderer.actors.get(key)
            if actor is not None:
                yield actor

    def _sync_scene(self) -> None:
        """Ensure the requested cutaway, then sync visibility and appearance."""
        s = self.state
        mode, slicing, selected = s.section_mode, self._slicing(), s.selected_name
        if not slicing:
            self._ensure_mode(mode)
        for m in MODES:
            for name, actor in self.body[m].items():
                on = (m == mode and not slicing and self._comp_visible(name)
                      and self._effective_opacity(name) > 0
                      and not (name == selected and self.hide_body))
                _set_shown(actor, on)
                actor.SetPosition(*self._offset(name))
                prop = actor.GetProperty()
                opacity = self._effective_opacity(name)
                prop.SetOpacity(opacity)
                prop.SetAmbient(0.48 if name == selected else 0.3)
                edge = self.edges[m].get(name)
                if edge is not None:
                    _set_shown(edge, on)
                    edge.SetPosition(*self._offset(name))
        shown = bool(selected) and self._comp_visible(selected)
        offset = self._offset(selected) if selected else (0.0, 0.0, 0.0)
        extras = [self.pl.renderer.actors.get(k) for k in ("slice-body", "slice-edges")]
        for actor in (*extras, *self._overlay_actors()):
            if actor is not None:
                actor.SetVisibility(shown)
                actor.SetPosition(*offset)
                if self.hide_body or self.state.demo_heating or self.state.result_id:
                    actor.GetProperty().SetOpacity(float(self.state.demo_opacity))

    def _refresh(self) -> None:
        """Rebuild the small dynamic pieces (slice body, tally overlays) and sync the scene."""
        key = (self.model.fingerprint, self.state.selected_name, self.state.section_mode,
               self.state.slice_phi_on, self.state.slice_phi, self.state.slice_z_on, self.state.slice_z,
               self.state.preview_tally, self.state.demo_heating, self.state.tally_mesh_index,
               self.state.result_id, self.state.result_quantity)
        if key == self._geometry_key:
            self._sync_scene()
            self._push()
            return
        self._geometry_key = key
        pl = self.pl
        for key in ("slice-body", "slice-edges", *OVERLAYS):
            pl.remove_actor(key, reset_camera=False, render=False)
        for title in list(pl.scalar_bars.keys()):
            pl.remove_scalar_bar(title, render=False)
        self.hide_body = False
        self.state.slice_info = ""
        comp = self._slice_comp()
        if comp is not None and self._slicing():
            section = self._slice_section()
            body = comp.profile.revolve(section, n_theta=self.n_theta)
            lines = comp.profile.outline(section, n_theta=self.n_theta)
            if self.port_openings is not None and body is not None:
                trimmed = self.port_openings.trim(body, comp)
                if trimmed is not body and lines is not None:
                    lines = self.port_openings.trim_outline(lines)
                body = trimmed
            if body is not None and body.n_cells:
                actor = pl.add_mesh(body, name="slice-body", color=comp.color,
                                    opacity=self._effective_opacity(comp.name), pickable=False, **BODY_STYLE)
                actor.mapper.scalar_visibility = False
            if lines is not None and lines.n_points:
                pl.add_mesh(lines, name="slice-edges", color=EDGE_COLOR, opacity=0.6, line_width=1.0,
                            lighting=False, pickable=False)
        try:
            self._update_overlays()
        except (ValueError, IndexError) as exc:
            self.state.action_status = f'Results unavailable: {exc}' 
        self._sync_scene()
        dynamic = self._slicing() or self.hide_body or any(True for _ in self._overlay_actors())
        self._push(geometry=dynamic or self.had_dynamic)
        self.had_dynamic = dynamic

    def _grid_for_selection(self):
        """``(CylGrid, tally_entry)`` of the selected cylindrical component.

        The grid follows the chosen mesh tally; a component without one gets a
        single-bin grid so its slices still show the cross-section.
        """
        comp = self._slice_comp()
        if comp is None:
            return None, None
        _, entries = tally_configuration(self.model.data, comp)
        index = int(self.state.tally_mesh_index)
        entry = entries[index] if 0 <= index < len(entries) and entries[index]["kind"] == "mesh_tallies" else None
        self.active_dataset = self.datasets.get(self.state.result_id)
        if self.active_dataset is not None:
            expected = component_id(self.model.manifest.device, comp.name)
            if self.active_dataset.component_id != expected:
                raise ValueError('The selected result belongs to a different component')
            mesh = self.active_dataset.mesh
            if not np.allclose(mesh.origin[:2], [0, 0]):
                raise ValueError('Off-axis cylindrical result origins are not supported yet')
            if not np.allclose([mesh.phi[0], mesh.phi[-1]], [0, 2*np.pi]):
                raise ValueError('Partial-azimuth cylindrical results are not supported yet')
            grid = CylGrid(mesh.r, mesh.phi, mesh.z + mesh.origin[2], mesh.shape)
            return grid, {'kind': 'mesh_tallies', 'scores': [self.active_dataset.score]}
        grid = CylGrid.from_profile(comp.profile, entry.get("dimensions") if entry else [1, 1, 1])
        if self.state.demo_heating and entry and 'heating' in entry.get('scores', []):
            cache_key = (self.model.fingerprint, comp.name, grid.shape)
            if cache_key not in self._demo_cache:
                self._demo_cache[cache_key] = synthetic_dataset(component_id(self.model.manifest.device, comp.name), grid)
            self.active_dataset = self._demo_cache[cache_key]
        return grid, entry

    def _update_overlays(self) -> None:
        """Tally slice planes when a slice is on; otherwise the bin cage / heating volume."""
        s, pl = self.state, self.pl
        self.active_dataset = None
        s.heating_available = False
        s.tally_summary = {}
        try:
            grid, entry = self._grid_for_selection()
        except ValueError as exc:
            s.action_status = f"Preview unavailable: {exc}"
            return
        if grid is None:
            return
        if (s.preview_tally and entry is not None) or self.active_dataset is not None:
            s.tally_summary = describe_overlay(grid, self.active_dataset, s.result_quantity)
        s.heating_available = bool(entry and "heating" in entry.get("scores", []))
        section = self._active_section()
        heat = self.active_dataset is not None
        scalar = self.active_dataset.scalar_name(s.result_quantity) if heat else HEAT
        field = self.active_dataset.values(s.result_quantity) if heat else None
        finite = field[np.isfinite(field)] if heat else None
        if heat and not finite.size:
            s.action_status = 'No finite values for this quantity (zero-mean relative errors are undefined).'
            return
        clim = (float(finite.min()), float(finite.max())) if heat else None
        if heat and clim[0] == clim[1]:
            clim = (clim[0], clim[0] + max(abs(clim[0]) * .01, 1e-12))
        bar = dict(SCALAR_BAR, title=scalar)
        lift = 0.05  # cm: keeps a slice plane clear of the cap face it lies on
        planes, notes = [], []
        if s.slice_phi_on:
            angle = np.radians(float(s.slice_phi))
            planes.append(("tally-slice-phi", grid.plane_phi(angle, field, section, lift, name=scalar)))
            notes.append(f"φ = {float(s.slice_phi):.1f}° · {grid.bin_text_phi(angle)}")
        if s.slice_z_on:
            height = float(s.slice_z)
            planes.append(("tally-slice-z", grid.plane_z(height, field, section, lift, name=scalar)))
            notes.append(f"z = {height:.1f} cm · {grid.bin_text_z(height)}")
        planes = [(n, m) for n, m in planes if m is not None and m.n_cells]
        s.slice_info = "\n".join(notes)
        if planes:
            for i, (key, mesh) in enumerate(planes):
                if heat:
                    pl.add_mesh(mesh, name=key, scalars=scalar, preference="cell", cmap="inferno", clim=clim,
                                opacity=float(s.demo_opacity), lighting=False, show_edges=False, pickable=False,
                                show_scalar_bar=(i == 0), scalar_bar_args=bar)
                else:
                    pl.add_mesh(mesh, name=key, color="#facc15", opacity=0.55, lighting=False, pickable=False)
            if s.preview_tally:
                lines = pv.PolyData()
                if s.slice_phi_on:
                    lines = lines.merge(grid.lines_phi(np.radians(float(s.slice_phi)), section, 2 * lift))
                if s.slice_z_on:
                    lines = lines.merge(grid.lines_z(float(s.slice_z), section, 2 * lift))
                if lines.n_points:
                    pl.add_mesh(lines, name="tally-lines", color="#1e293b", line_width=1.2,
                                lighting=False, pickable=False)
            return
        if heat:
            volume = grid.volume(field, section, name=scalar)
            if volume is not None and volume.n_cells:
                volume = lift_boundary(volume, grid, section)
                pl.add_mesh(volume, name="tally-map", scalars=scalar, preference="cell", cmap="inferno",
                            clim=clim, opacity=float(s.demo_opacity), show_edges=False,
                            lighting=False, interpolate_before_map=False, pickable=False,
                            scalar_bar_args=bar)
        if s.preview_tally and entry is not None:
            cage = grid.cage(section)
            if cage.n_points:
                cage = lift_boundary(cage, grid, section, distance=.16)
                pl.add_mesh(cage, name="tally-grid", color="#334155" if heat else "#d97706",
                            opacity=.35 if heat else .85, line_width=1.1, lighting=False, pickable=False)
