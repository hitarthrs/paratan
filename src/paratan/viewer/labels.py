"""Human-readable names for viewer components, groups, materials and metadata."""
from __future__ import annotations

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
    return humanize(name)
