"""Shared viewer types for Paratan component meshes."""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

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
