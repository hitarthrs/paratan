"""Model preparation and render policy, independent of Trame."""
from __future__ import annotations
import hashlib
from collections import OrderedDict
from dataclasses import dataclass
import numpy as np
import pyvista as pv
import yaml
from src.paratan.viewer.components import ComponentMesh
from src.paratan.viewer.labels import GROUP_LABELS
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
    port_openings = None
    if adapter.kind == "simple_mirror" and (data.get("ports") or {}).get("ports_present"):
        from src.paratan.viewer.port_openings import PortOpenings
        port_openings = PortOpenings(data)
    port_key = hashlib.sha256(yaml.safe_dump({"ports":data.get("ports"), "vacuum_vessel":data.get("vacuum_vessel")}).encode()).hexdigest() if port_openings else ""
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
                 display_manifest(data, components), "Tandem mirror" if adapter.kind == "tandem" else "Simple mirror")


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
