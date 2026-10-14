"""Model preparation and render policy, independent of Trame."""
from __future__ import annotations
import hashlib
from collections import OrderedDict
from dataclasses import dataclass, field, fields
from copy import deepcopy
from math import isfinite
import numpy as np
import pyvista as pv
import yaml
from src.paratan.viewer.components import ComponentMesh, GROUP_LABELS
from src.paratan.viewer.revolve import FULL, Section, section_closed_mesh
from src.paratan.viewer.devices import adapter_for
from src.paratan.viewer.manifest import display_manifest, ModelManifest
GROUPS = tuple(GROUP_LABELS)
MODES = ("none", "half", "quarter")
MODE_LABELS = (("none", "Full"), ("half", "Half cut"), ("quarter", "Quarter cut"))
MODE_SECTIONS = {
    "none": FULL,
    "half": Section.without_wedge(0.0, np.pi),
    "quarter": Section.without_wedge(0.0, 0.5 * np.pi),
}
SLICE_WEDGE = 0.5 * np.pi  # material removed ahead of an azimuthal slice plane
MAX_UPLOAD_BYTES = 2_000_000
HEAT = "Demo heating [a.u.]"
OVERLAYS = ("tally-grid", "tally-map", "tally-slice-phi", "tally-slice-z", "tally-lines")
EDGE_COLOR = "#1f2937"
GROUP_OFFSETS = {
    "VV": (0.0, 0.0, 0.0), "FW": (0.0, 0.0, 0.0), "CC": (1.0, 0.0, 0.0), "LF": (0.0, 1.0, 0.0),
    "HF": (0.0, 0.0, 1.0), "ends": (0.0, 0.0, 1.5), "ports": (1.2, 0.0, 0.0),
    "IO": (1.0, 0.0, 0.5),
}
SCALAR_BAR = {"title": "Heating, demo [a.u.]", "color": "#334155", "fmt": "%.0f", "vertical": True,
              "position_x": 0.92, "position_y": 0.22, "height": 0.5, "width": 0.04, "n_labels": 5}
BODY_STYLE = dict(smooth_shading=False, ambient=0.3, diffuse=0.7, specular=0.08, specular_power=12,
                  show_scalar_bar=False)
_GEOMETRY_CACHE = OrderedDict()


def _geometry_key(comp, mode, n_theta):
    digest = hashlib.sha256()
    digest.update(comp.mesh.points.tobytes())
    digest.update(comp.mesh.faces.tobytes())
    if comp.profile is not None:
        for loop in comp.profile.loops:
            digest.update(loop.tobytes())
    return comp.name, mode, n_theta, digest.hexdigest()


@dataclass
class Model:
    """One loaded input file: components plus every precomputed cutaway."""

    name: str
    data: dict
    fingerprint: str
    components: list[ComponentMesh]
    by_name: dict[str, ComponentMesh]
    centers: dict[str, tuple[float, float, float]]
    meshes: dict[str, dict[str, pv.PolyData]]    # mode -> component -> mesh
    outlines: dict[str, dict[str, pv.PolyData]]  # mode -> component -> edge lines
    manifest: ModelManifest
    device_label: str


def _with_normals(mesh: pv.PolyData) -> pv.PolyData:
    """Revolved parts carry analytic normals; other closed meshes get crease-aware ones."""
    if "Normals" in mesh.point_data:
        return mesh
    # Keep the winding as built: VTK's "consistent" re-ordering flips capped clips inside out.
    return mesh.compute_normals(cell_normals=False, split_vertices=True, feature_angle=30.0,
                                consistent_normals=False, auto_orient_normals=False, inplace=False)



def build_model(text: str, name: str, n_theta: int) -> Model:
    """Parse YAML text and mesh every component for every cutaway mode."""
    data = yaml.safe_load(text)
    if not isinstance(data, dict):
        raise ValueError("the file is not a simple-mirror input (expected a YAML mapping)")
    adapter = adapter_for(data)
    components = adapter.build_components(data, n_theta)
    from src.paratan.viewer.port_openings import openings_for
    port_openings = openings_for(data)
    port_key = (
        hashlib.sha256(
            yaml.safe_dump(
                {
                    "ports": data.get("ports"),
                    "vacuum_vessel": data.get("vacuum_vessel"),
                    "blanket_inlet_outlet": data.get("blanket_inlet_outlet"),
                }
            ).encode()
        ).hexdigest()
        if port_openings
        else ""
    )
    meshes: dict[str, dict[str, pv.PolyData]] = {m: {} for m in MODES}
    outlines: dict[str, dict[str, pv.PolyData]] = {m: {} for m in MODES}
    for mode in MODES:
        section = MODE_SECTIONS[mode]
        for comp in components:
            cache_key = (*_geometry_key(comp, mode, n_theta),
                         port_key)
            cached = _GEOMETRY_CACHE.get(cache_key)
            if cached is not None:
                meshes[mode][comp.name], lines = cached
                if lines is not None:
                    outlines[mode][comp.name] = lines
                _GEOMETRY_CACHE.move_to_end(cache_key)
                continue
            if comp.profile is not None:
                mesh = comp.profile.revolve(section, n_theta=n_theta)
                lines = comp.profile.outline(section, n_theta=n_theta)
            else:
                try:  # clip the clean mesh: normal-split copies are not manifold
                    mesh = comp.mesh if section.is_full else section_closed_mesh(comp.mesh, section)
                except ValueError:
                    mesh = None
                lines = None
            if mesh is None or not mesh.n_cells:
                continue
            if port_openings is not None:
                trimmed = port_openings.trim(mesh, comp)
                if trimmed is not mesh:
                    mesh = trimmed
                    # Clip design outlines too, without drawing tessellation.
                    if lines is not None:
                        lines = port_openings.trim_outline(lines)
            if not mesh.n_cells:
                continue
            mesh = _with_normals(mesh.copy(deep=True))
            if mesh is None or not mesh.n_cells:
                continue
            mesh.GetPointData().SetActiveScalars(None)
            mesh.field_data["component_name"] = [comp.name]
            meshes[mode][comp.name] = mesh
            if lines is not None and lines.n_points:
                outlines[mode][comp.name] = lines
            _GEOMETRY_CACHE[cache_key] = (mesh, lines)
            while len(_GEOMETRY_CACHE) > 128:
                _GEOMETRY_CACHE.popitem(last=False)
    centers = {}
    for comp in components:
        x0, x1, y0, y1, z0, z1 = comp.mesh.bounds
        centers[comp.name] = (0.5 * (x0 + x1), 0.5 * (y0 + y1), 0.5 * (z0 + z1))
    return Model(name, data, hashlib.sha256(text.encode()).hexdigest(), components,
                 {c.name: c for c in components}, centers, meshes, outlines,
                 display_manifest(data, components),
                 {"tandem": "Tandem mirror", "ffh": "FFH"}.get(adapter.kind, "Simple mirror"))


HIDDEN_SCALE = 1e-6


def _set_shown(actor, shown: bool) -> None:
    """Show / hide a preloaded actor without dropping it from the browser's scene.

    trame-vtk serializes no geometry for invisible actors, so toggling visibility makes
    the browser re-download and rebuild meshes.  Collapsing the actor to a point keeps
    its geometry resident and turns every cut / visibility change into a transform update.
    """
    actor.SetVisibility(True)
    actor.SetScale(1.0 if shown else HIDDEN_SCALE)
    actor.SetPickable(shown)


def is_shown(actor) -> bool:
    return bool(actor.GetVisibility()) and actor.GetScale()[0] > HIDDEN_SCALE


# --- scene state (merged from state.py) ---


@dataclass
class SceneState:
    selected_name: str = ''
    has_selection: bool = False
    selected_opacity: float = 1.0
    solo: bool = False
    section_mode: str = 'none'
    explode: float = 0.0
    slice_phi_on: bool = False
    slice_phi: float = 0.0
    slice_z_on: bool = False
    slice_z: float = 0.0
    has_slice: bool = False
    slice_z_min: float = 0.0
    slice_z_max: float = 1.0
    slice_info: str = ''
    preview_tally: bool = False
    demo_heating: bool = False
    demo_opacity: float = 1.0
    heating_available: bool = False
    tally_summary: dict = field(default_factory=dict)
    tally_mesh_index: int = 0
    result_id: str = ''
    result_quantity: str = 'mean'
    groups: dict[str, bool] = field(default_factory=dict)
    hidden: set[str] = field(default_factory=set)
    opacity: dict[str, float] = field(default_factory=dict)
    camera: dict = field(default_factory=dict)
    action_status: str = ''

    def __getattr__(self, key):
        if key.startswith('vis_'):
            return self.groups.get(key[4:], True)
        raise AttributeError(key)

    def __setattr__(self, key, value):
        if key.startswith('vis_') and 'groups' in self.__dict__:
            self.groups[key[4:]] = bool(value)
        else:
            super().__setattr__(key, value)

    def snapshot(self):
        return deepcopy(self.__dict__)


SCENE_FIELDS = {f.name for f in fields(SceneState)} - {'groups', 'hidden', 'opacity', 'camera'}


class SceneActions:
    """The same actions can be used from Trame, scripts, tests, or another UI."""
    def __init__(self, state: SceneState):
        self.state = state
        self.components = set()

    def install(self, components, groups):
        self.components = set(components)
        self.state.groups = {g: True for g in groups}
        self.reset()

    def set(self, key, value):
        if key.startswith('vis_'):
            setattr(self.state, key, bool(value))
            return
        if key not in SCENE_FIELDS:
            raise ValueError(f'Unknown scene setting: {key}')
        if key in ('explode', 'selected_opacity', 'demo_opacity'):
            value = float(value)
            if not isfinite(value) or not 0 <= value <= 1:
                raise ValueError(f'{key} must be between 0 and 1')
        elif key in ('slice_phi', 'slice_z', 'slice_z_min', 'slice_z_max'):
            value = float(value)
            if not isfinite(value):
                raise ValueError(f'{key} must be finite')
        elif key == 'section_mode' and value not in ('none', 'half', 'quarter'):
            raise ValueError('Unknown section mode')
        elif key == 'result_quantity' and value not in ('mean', 'std_dev', 'relative_error'):
            raise ValueError('Unknown result quantity')
        elif key == 'tally_mesh_index':
            value = int(value)
            if value < 0:
                raise ValueError('Tally index must be nonnegative')
        setattr(self.state, key, value)

    def select(self, name):
        if name not in self.components:
            raise ValueError(f'Unknown component: {name}')
        self.state.selected_name = name
        self.state.has_selection = True
        self.state.result_id = ''

    def clear(self):
        s = self.state
        s.selected_name, s.has_selection, s.solo = '', False, False
        s.slice_phi_on = s.slice_z_on = s.has_slice = False
        s.result_id = ''

    def hide_selected(self):
        if self.state.selected_name:
            self.state.hidden.add(self.state.selected_name)

    def show_all(self):
        self.state.hidden.clear()
        self.state.groups = dict.fromkeys(self.state.groups, True)
        self.state.solo = False

    def reset(self):
        groups = dict.fromkeys(self.state.groups, True)
        fresh = SceneState(groups=groups)
        self.state.__dict__.clear()
        self.state.__dict__.update(fresh.__dict__)


class TrameStateBridge:
    """Trame is a projection of SceneState plus presentation-only fields.

    UI changes are validated at the boundary. Python actions publish back to
    Trame. Rendering never depends on Trame's event or state implementation.
    """
    def __init__(self, ui, scene, actions):
        object.__setattr__(self, 'ui', ui)
        object.__setattr__(self, 'scene', scene)
        object.__setattr__(self, 'actions', actions)

    def __getattr__(self, key):
        if key in SCENE_FIELDS or key.startswith('vis_'):
            return getattr(self.scene, key)
        return getattr(self.ui, key)

    def __setattr__(self, key, value):
        if key in SCENE_FIELDS or key.startswith('vis_'):
            self.actions.set(key, value)
        setattr(self.ui, key, value)

    def update(self, mapping):
        for key, value in mapping.items():
            setattr(self, key, value)

    def publish_scene(self):
        self.ui.update({key: getattr(self.scene, key) for key in SCENE_FIELDS})
        self.ui.update({f'vis_{g}': value for g, value in self.scene.groups.items()})

    def change(self, *names):
        def decorate(callback):
            def changed(**kwargs):
                before = self.scene.snapshot()
                try:
                    # One browser event can change multiple independently watched
                    # fields. Import the complete projection before any callback
                    # publishes it, otherwise an early callback can overwrite the
                    # other incoming changes with stale canonical values.
                    keys = SCENE_FIELDS | {f'vis_{g}' for g in self.scene.groups}
                    incoming = {key: getattr(self.ui, key) for key in keys}
                    for key, value in incoming.items():
                        self.actions.set(key, value)
                except (ValueError, TypeError) as exc:
                    self.scene.__dict__.clear()
                    self.scene.__dict__.update(before)
                    self.scene.action_status = str(exc)
                    self.publish_scene()
                    return
                callback(**kwargs)
                self.publish_scene()
            self.ui.change(*names)(changed)
            return changed
        return decorate
