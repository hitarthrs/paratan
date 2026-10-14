"""Shared viewer types for Paratan component meshes."""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any
import re

import pyvista as pv

from src.paratan.viewer.revolve import Profile


@dataclass
class ComponentMesh:
    """One named, pickable component in the viewer."""

    name: str
    group: str
    material: str
    mesh: pv.DataSet
    meta: dict[str, Any] = field(default_factory=dict)
    color: tuple[float, float, float] = (0.7, 0.7, 0.7)
    opacity: float = 1.0
    # (r, z) cross-section for axisymmetric parts; enables exact sections and tally slices.
    profile: Profile | None = None
    label: str = ""  # display name; ``name`` stays the stable identifier


# Material-role palette (RGB 0–1), aligned with FFH/SolidRayTrace vibes.
MATERIAL_COLORS: dict[str, tuple[float, float, float]] = {
    "vacuum": (0.85, 0.90, 0.95),
    "air": (0.92, 0.92, 0.94),
    "tungsten": (0.31, 0.31, 0.33),
    "FW_FNSF": (0.80, 0.31, 0.24),
    "BW_FNSF": (0.33, 0.52, 0.67),
    "eutectic_breeding_material": (0.84, 0.67, 0.07),
    "HT_Shield_filler": (0.28, 0.33, 0.41),
    "water_cooled_wc": (0.19, 0.23, 0.31),
    "stainless": (0.67, 0.67, 0.70),
    "ss316ln": (0.37, 0.66, 0.42),
    "Magnet_Winding_Pack_2": (0.26, 0.44, 0.88),
    "cooled_tungsten_boride": (0.31, 0.51, 0.36),
}


def color_for_material(name: str) -> tuple[float, float, float]:
    if name in MATERIAL_COLORS:
        return MATERIAL_COLORS[name]
    # Stable fallback from name hash
    import hashlib

    h = int.from_bytes(hashlib.sha256(name.encode()).digest()[:4], "big") % 360
    import colorsys

    r, g, b = colorsys.hsv_to_rgb(h / 360.0, 0.45, 0.75)
    return (r, g, b)


# --- display labels (merged from labels.py) ---

import re


GROUP_LABELS = {
    "VV": "Vacuum vessel",
    "FW": "First wall",
    "CC": "Central cell",
    "EP": "End plugs",
    "LF": "Low-field coils",
    "HF": "High-field coils",
    "ends": "End cells",
    "ports": "Ports",
    "IO": "Inlet / outlet",
}

_MATERIAL_LABELS = {
    "vacuum": "Vacuum",
    "air": "Air",
    "tungsten": "Tungsten",
    "stainless": "Stainless steel",
    "ss316ln": "SS 316LN",
    "FW_FNSF": "First wall (FNSF)",
    "BW_FNSF": "Back wall (FNSF)",
    "eutectic_breeding_material": "Eutectic breeder",
    "HT_Shield_filler": "HT shield filler",
    "water_cooled_wc": "Water-cooled WC",
    "Magnet_Winding_Pack_2": "Magnet winding pack",
    "cooled_tungsten_boride": "Cooled tungsten boride",
}

_META_LABELS = {
    "thickness_cm": "Thickness",
    "layer_index": "Layer number",
    "r_inner_cm": "Inner radius",
    "r_outer_cm": "Outer radius",
    "inner_radius_cm": "Inner radius",
    "outer_radius_cm": "Outer radius",
    "axial_length_cm": "Axial length",
    "z0_cm": "Axial centre",
    "z_min_cm": "Axial start",
    "z_max_cm": "Axial end",
    "bore_start_r_cm": "Bore start radius",
    "radial_thickness_cm": "Radial thickness",
    "central_radius_cm": "Central radius",
    "bottleneck_radius_cm": "Bottleneck radius",
    "shell_thickness_cm": "Shell thickness",
    "geometry_style": "Geometry style",
    "tilt_deg": "Tilt",
    "x_limits_cm": "X limits",
}


def humanize(text: str) -> str:
    """``first_wall`` -> ``First wall``; tokens that already carry capitals (FNSF) are kept."""
    words = text.replace("_", " ").split()
    words = [w if any(c.isupper() for c in w) else w.lower() for w in words]
    out = " ".join(words)
    return out[:1].upper() + out[1:]


def material_label(name: str) -> str:
    return _MATERIAL_LABELS.get(name) or humanize(name)


def meta_row(key: str, value) -> dict:
    """Inspector row: readable label, formatted value and unit."""
    label = _META_LABELS.get(key) or humanize(re.sub(r"_(cm|deg)$", "", key))
    unit = "cm" if key.endswith("_cm") else "°" if key.endswith("_deg") else ""
    if key == "layer_index":
        value = int(value) + 1  # shown 1-based, matching the component names
    if isinstance(value, (list, tuple)):
        text = " to ".join(f"{v:g}" for v in value)
    elif isinstance(value, (int, float)) and not isinstance(value, bool):
        text = f"{value:g}"
    else:
        text = humanize(str(value)) if key == "geometry_style" else str(value)
    return {"label": label, "value": text, "unit": unit}


def component_label(name: str, group: str, material: str, meta: dict) -> str:
    """Short display name; the group heading supplies the rest of the context."""
    if name == "vacuum_vessel_plasma":
        return "Plasma region"
    if name.startswith("vv_"):
        return humanize(name[3:])
    if group == "CC" and "layer_index" in meta:
        return f"Layer {int(meta['layer_index']) + 1} · {material_label(material)}"
    if m := re.fullmatch(r"lf_coil_(\d+)_(magnet|shield)", name):
        return f"Coil {int(m[1]) + 1} · {m[2].capitalize()}"
    if m := re.fullmatch(r"hf_coil_(left|right)_(magnet|shield|casing_(\d+))", name):
        part = f"Casing {int(m[3]) + 1}" if m[3] is not None else m[2].capitalize()
        return f"{m[1].capitalize()} · {part}"
    if m := re.fullmatch(r"end_cell_(left|right)_(shell|inner)", name):
        return f"{m[1].capitalize()} · {'Shell' if m[2] == 'shell' else 'Interior'}"
    if m := re.fullmatch(r"port_(left|right)_(core|wall_(\d+))", name):
        return f"{m[1].capitalize()} · {'Core' if m[2] == 'core' else 'Wall ' + m[3]}"
    if m := re.fullmatch(r"port_(left|right)_(pos|neg)_layer_(\d+)", name):
        layer = int(m[3])
        part = "Core" if layer == 0 else f"Wall {layer}"
        return f"{m[1].capitalize()} plug · {m[2]} · {part}"
    if m := re.fullmatch(r"(inlet|outlet)_(\d+)_(wall|flow)", name):
        part = "Wall" if m[3] == "wall" else "Flow"
        return f"{m[1].capitalize()} {m[2]} · {part}"
    return humanize(name)
