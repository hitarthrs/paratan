"""Device-independent scene state and validated actions; no UI or renderer imports."""
from __future__ import annotations

from dataclasses import dataclass, field, fields
from copy import deepcopy
from math import isfinite


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
