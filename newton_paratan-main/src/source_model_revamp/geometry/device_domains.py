"""ParaTAN component locations and source model axial domains

Convert geometry inputs in centimeters to machine coordinates in meters
Combine the layout with magnetic throats and selected material boundaries
"""
from __future__ import annotations
from collections.abc import Mapping, Sequence
from dataclasses import dataclass, replace
from typing import Any
import numpy as np
from source_model_revamp.constants import CM_TO_M

def _block(root: Mapping[str, Any], key: str) -> Mapping[str, Any]:
    """Read a required ParaTAN geometry mapping"""
    value = root.get(key)
    if not isinstance(value, Mapping):
        raise ValueError(f"ParaTAN geometry requires a {key!r} mapping")

    return value

def _positive(block: Mapping[str, Any], key: str, owner: str) -> float:
    """Read a required finite positive geometry dimension"""
    if key not in block:
        raise ValueError(f"{owner}.{key} is required to derive source model domains")
    value = float(block[key])
    if not np.isfinite(value) or value <= 0.0:
        raise ValueError(f"{owner}.{key} must be finite and positive")

    return value

def _sum_layer_thickness_cm(items: Any) -> float:
    """Sum thickness entries from a sequence of layer mappings in centimeters"""
    if not isinstance(items, Sequence) or isinstance(items, (str, bytes)):
        return 0.0

    return float(sum(float(item.get("thickness", 0.0)) for item in items if isinstance(item, Mapping)))

@dataclass(frozen=True)
class AxialDomain:
    """Closed axial interval in machine coordinates

    Bounds are in meters and source records how they were chosen
    """
    name: str
    z_min_m: float
    z_max_m: float
    source: str

    def __post_init__(self) -> None:
        """Require finite bounds and a positive interval length"""
        if not np.isfinite(self.z_min_m) or not np.isfinite(self.z_max_m):
            raise ValueError(f"{self.name} bounds must be finite")
        if self.z_max_m <= self.z_min_m:
            raise ValueError(f"{self.name} must have positive axial extent")

    @property
    def length_m(self) -> float:
        """Return the axial extent in meters"""
        return float(self.z_max_m - self.z_min_m)

    @property
    def midpoint_m(self) -> float:
        """Return the geometric midpoint of this interval"""
        return float(0.5 * (self.z_min_m + self.z_max_m))

    def contains(self, z_m: Any, *, atol_m: float = 0.0) -> np.ndarray:
        """Mark coordinates inside the closed interval with an optional boundary tolerance"""
        z = np.asarray(z_m, dtype=float)
        return (z >= self.z_min_m - atol_m) & (z <= self.z_max_m + atol_m)

    def as_dict(self) -> dict[str, Any]:
        """Return the interval bounds and their origin as metadata"""
        return {"name": self.name, "z_min_m": float(self.z_min_m), "z_max_m": float(self.z_max_m), "source": self.source}

@dataclass(frozen=True)
class ParaTANLayout:
    """Geometry derived component locations before magnetic throat fitting

    All coordinates and dimensions are in meters
    Conical domains are absent for perpendicular vessel transitions
    full_device_domain includes the outer end cell shell extents
    """
    midplane_z_m: float
    central_cell_plasma_domain: AxialDomain
    left_conical_vacuum_domain: AxialDomain | None
    right_conical_vacuum_domain: AxialDomain | None
    left_bottleneck_vacuum_domain: AxialDomain
    right_bottleneck_vacuum_domain: AxialDomain
    vacuum_vessel_domain: AxialDomain
    left_end_cell_vacuum_domain: AxialDomain
    right_end_cell_vacuum_domain: AxialDomain
    full_device_domain: AxialDomain
    central_cell_vacuum_radius_m: float
    bottleneck_vacuum_radius_m: float
    end_cell_vacuum_radius_m: float
    left_hf_coil_center_m: float
    right_hf_coil_center_m: float
    left_end_cell_center_m: float
    right_end_cell_center_m: float
    end_cell_shell_thickness_m: float
    hf_to_end_cell_gap_m: float

    def as_dict(self) -> dict[str, Any]: 
        """Return component locations and domains as nested metadata"""
        return {
            "midplane_z_m": float(self.midplane_z_m),
            "central_cell_plasma_domain": self.central_cell_plasma_domain.as_dict(),
            "left_conical_vacuum_domain": None if self.left_conical_vacuum_domain is None else self.left_conical_vacuum_domain.as_dict(),
            "right_conical_vacuum_domain": None if self.right_conical_vacuum_domain is None else self.right_conical_vacuum_domain.as_dict(),
            "left_bottleneck_vacuum_domain": self.left_bottleneck_vacuum_domain.as_dict(),
            "right_bottleneck_vacuum_domain": self.right_bottleneck_vacuum_domain.as_dict(),
            "vacuum_vessel_domain": self.vacuum_vessel_domain.as_dict(),
            "left_end_cell_vacuum_domain": self.left_end_cell_vacuum_domain.as_dict(),
            "right_end_cell_vacuum_domain": self.right_end_cell_vacuum_domain.as_dict(),
            "full_device_domain": self.full_device_domain.as_dict(),
            "central_cell_vacuum_radius_m": float(self.central_cell_vacuum_radius_m),
            "bottleneck_vacuum_radius_m": float(self.bottleneck_vacuum_radius_m),
            "end_cell_vacuum_radius_m": float(self.end_cell_vacuum_radius_m),
            "left_hf_coil_center_m": float(self.left_hf_coil_center_m),
            "right_hf_coil_center_m": float(self.right_hf_coil_center_m),
            "left_end_cell_center_m": float(self.left_end_cell_center_m),
            "right_end_cell_center_m": float(self.right_end_cell_center_m),
            "end_cell_shell_thickness_m": float(self.end_cell_shell_thickness_m),
            "hf_to_end_cell_gap_m": float(self.hf_to_end_cell_gap_m),
        }

@dataclass(frozen=True)
class DeviceDomains:
    """Axial supports used by the coupled source model

    The confined interval spans the magnetic throats
    Each expander spans a throat and the selected outer boundary
    Fusion, neutron source, and export domains span the selected device interval
    """
    midplane_z_m: float
    full_device_domain: AxialDomain
    central_cell_plasma_domain: AxialDomain
    confined_orbit_domain: AxialDomain
    left_expander_domain: AxialDomain
    right_expander_domain: AxialDomain
    fusion_domain: AxialDomain
    neutron_source_domain: AxialDomain
    openmc_export_domain: AxialDomain

    def as_dict(self) -> dict[str, Any]:
        """Return the named physics domains as nested metadata"""
        return {
            "midplane_z_m": float(self.midplane_z_m),
            "full_device_domain": self.full_device_domain.as_dict(),
            "central_cell_plasma_domain": self.central_cell_plasma_domain.as_dict(),
            "confined_orbit_domain": self.confined_orbit_domain.as_dict(),
            "left_expander_domain": self.left_expander_domain.as_dict(),
            "right_expander_domain": self.right_expander_domain.as_dict(),
            "fusion_domain": self.fusion_domain.as_dict(),
            "neutron_source_domain": self.neutron_source_domain.as_dict(),
            "openmc_export_domain": self.openmc_export_domain.as_dict(),
        }

def derive_paratan_layout(root: Mapping[str, Any]) -> ParaTANLayout:
    """Derive component intervals, radii, and coil centers from ParaTAN inputs

    Read dimensions in centimeters and return a layout in meters
    """
    vv = _block(root, "vacuum_vessel")
    cc = _block(root, "central_cell")
    hf = _block(root, "hf_coil")
    end = _block(root, "end_cell")
    shield = hf.get("shield")
    magnet = hf.get("magnet")
    if not isinstance(shield, Mapping) or not isinstance(magnet, Mapping):
        raise ValueError("hf_coil.shield and hf_coil.magnet mappings are required")
    geometry_style = str(vv.get("geometry_style") or "conical").strip().lower()
    if geometry_style not in {"conical", "perpendicular"}:
        raise ValueError("vacuum_vessel.geometry_style must be conical or perpendicular")
    midplane_cm = float(vv.get("axial_midplane", 0.0))
    if not np.isfinite(midplane_cm):
        raise ValueError("vacuum_vessel.axial_midplane must be finite")
    component_length_cm = _positive(cc, "axial_length", "central_cell")
    outer_axial_length_cm = _positive(vv, "outer_axial_length", "vacuum_vessel")
    vv_central_length_cm = _positive(vv, "central_axial_length", "vacuum_vessel")
    left_bottleneck_cm = _positive(vv, "left_bottleneck_length", "vacuum_vessel")
    right_bottleneck_cm = _positive(vv, "right_bottleneck_length", "vacuum_vessel")
    end_length_cm = _positive(end, "axial_length", "end_cell")
    end_shell_cm = _positive(end, "shell_thickness", "end_cell")
    end_diameter_cm = _positive(end, "diameter", "end_cell")
    end_gap_cm = float(end.get("hf_to_end_cell_gap_cm", 10.0))
    if not np.isfinite(end_gap_cm) or end_gap_cm < 0.0:
        raise ValueError("end_cell.hf_to_end_cell_gap_cm must be finite and nonnegative")
    casing_cm = _sum_layer_thickness_cm(hf.get("casing_layers", ()))
    shield_gap_cm = float(shield.get("shield_central_cell_gap", 0.0))
    shield_axial_values = shield.get("axial_thickness", 0.0)
    # A scalar gives equal end thicknesses while a pair gives inward then outward values
    if isinstance(shield_axial_values, Sequence) and not isinstance(shield_axial_values, (str, bytes)):
        if len(shield_axial_values) != 2:
            raise ValueError("hf_coil.shield.axial_thickness must contain central facing and outward facing values")
        shield_axial_inward_cm = float(shield_axial_values[0])
        shield_axial_outward_cm = float(shield_axial_values[1])
    else:
        shield_axial_inward_cm = float(shield_axial_values)
        shield_axial_outward_cm = float(shield_axial_values)
    if min(shield_axial_inward_cm, shield_axial_outward_cm) < 0.0 or not np.isfinite([shield_axial_inward_cm, shield_axial_outward_cm]).all():
        raise ValueError("hf_coil.shield axial thicknesses must be finite and nonnegative")
    magnet_axial_cm = _positive(magnet, "axial_thickness", "hf_coil.magnet")
    hf_offset_cm = float(hf.get("axial_offset_cm", hf.get("z0_offset_cm", hf.get("center_offset_cm", hf.get("outward_shift_cm", 0.0))),))
    if not np.isfinite(hf_offset_cm):
        raise ValueError("hf_coil.axial_offset_cm must be finite")
    cc_half_cm = 0.5 * component_length_cm
    vv_cc_offset_cm = max(0.5 * (vv_central_length_cm - component_length_cm), 0.0)
    # Build outward offsets from the midplane before converting positions to meters
    hf_center_rel_cm = (cc_half_cm + shield_gap_cm + shield_axial_inward_cm + casing_cm + 0.5 * magnet_axial_cm + vv_cc_offset_cm + hf_offset_cm)
    hf_outer_edge_rel_cm = hf_center_rel_cm + 0.5 * magnet_axial_cm + casing_cm + shield_axial_outward_cm
    end_center_rel_cm = hf_outer_edge_rel_cm + end_gap_cm + end_shell_cm + 0.5 * end_length_cm
    # ParaTAN uses different length inputs for the cylindrical center in each vessel style
    if geometry_style == "conical":
        central_plasma_length_cm = outer_axial_length_cm
        central_source = "registered vacuum_vessel.central_cylinder component (vacuum_vessel.outer_axial_length)"
        left_conical = AxialDomain("left_conical_vacuum_domain", (midplane_cm - 0.5 * vv_central_length_cm) * CM_TO_M, (midplane_cm - 0.5 * central_plasma_length_cm) * CM_TO_M, "vacuum_vessel conical transition component")
        right_conical = AxialDomain("right_conical_vacuum_domain", (midplane_cm + 0.5 * central_plasma_length_cm) * CM_TO_M, (midplane_cm + 0.5 * vv_central_length_cm) * CM_TO_M, "vacuum_vessel conical transition component")
    else:
        central_plasma_length_cm = vv_central_length_cm
        central_source = "registered vacuum_vessel.central_cylinder component (vacuum_vessel.central_axial_length)"
        left_conical = None
        right_conical = None
    central = AxialDomain("central_cell_plasma_domain", (midplane_cm - 0.5 * central_plasma_length_cm) * CM_TO_M, (midplane_cm + 0.5 * central_plasma_length_cm) * CM_TO_M, central_source)
    left_bottleneck = AxialDomain("left_bottleneck_vacuum_domain", (midplane_cm - 0.5 * vv_central_length_cm - left_bottleneck_cm) * CM_TO_M, (midplane_cm - 0.5 * vv_central_length_cm) * CM_TO_M, "vacuum_vessel left bottleneck component")
    right_bottleneck = AxialDomain("right_bottleneck_vacuum_domain", (midplane_cm + 0.5 * vv_central_length_cm) * CM_TO_M, (midplane_cm + 0.5 * vv_central_length_cm + right_bottleneck_cm) * CM_TO_M, "vacuum_vessel right bottleneck component")
    vv_domain = AxialDomain("vacuum_vessel_domain", (midplane_cm - 0.5 * vv_central_length_cm - left_bottleneck_cm) * CM_TO_M, (midplane_cm + 0.5 * vv_central_length_cm + right_bottleneck_cm) * CM_TO_M, "vacuum_vessel central length plus independent bottleneck lengths")
    left_center_cm = midplane_cm - end_center_rel_cm
    right_center_cm = midplane_cm + end_center_rel_cm
    left_vacuum = AxialDomain("left_end_cell_vacuum_domain", (left_center_cm - 0.5 * end_length_cm) * CM_TO_M, (left_center_cm + 0.5 * end_length_cm) * CM_TO_M, "end_cell placement and inner axial length")
    right_vacuum = AxialDomain("right_end_cell_vacuum_domain", (right_center_cm - 0.5 * end_length_cm) * CM_TO_M, (right_center_cm + 0.5 * end_length_cm) * CM_TO_M, "end_cell placement and inner axial length")
    full = AxialDomain("full_device_domain", min(vv_domain.z_min_m, left_vacuum.z_min_m - end_shell_cm * CM_TO_M), max(vv_domain.z_max_m, right_vacuum.z_max_m + end_shell_cm * CM_TO_M), "union of ParaTAN vacuum vessel and end cell material extents")

    return ParaTANLayout(
        midplane_z_m=midplane_cm * CM_TO_M,
        central_cell_plasma_domain=central,
        left_conical_vacuum_domain=left_conical,
        right_conical_vacuum_domain=right_conical,
        left_bottleneck_vacuum_domain=left_bottleneck,
        right_bottleneck_vacuum_domain=right_bottleneck,
        vacuum_vessel_domain=vv_domain,
        left_end_cell_vacuum_domain=left_vacuum,
        right_end_cell_vacuum_domain=right_vacuum,
        full_device_domain=full,
        central_cell_vacuum_radius_m=_positive(vv, "central_radius", "vacuum_vessel") * CM_TO_M,
        bottleneck_vacuum_radius_m=_positive(vv, "bottleneck_radius", "vacuum_vessel") * CM_TO_M,
        end_cell_vacuum_radius_m=0.5 * end_diameter_cm * CM_TO_M,
        left_hf_coil_center_m=(midplane_cm - hf_center_rel_cm) * CM_TO_M,
        right_hf_coil_center_m=(midplane_cm + hf_center_rel_cm) * CM_TO_M,
        left_end_cell_center_m=left_center_cm * CM_TO_M,
        right_end_cell_center_m=right_center_cm * CM_TO_M,
        end_cell_shell_thickness_m=end_shell_cm * CM_TO_M,
        hf_to_end_cell_gap_m=end_gap_cm * CM_TO_M,
    )

def attach_confined_domain(layout: ParaTANLayout, left_throat_z_m: float, right_throat_z_m: float, *, left_boundary_z_m: float | None = None, right_boundary_z_m: float | None = None) -> DeviceDomains:
    """Build physics domains from the layout, throats, and optional outer boundaries

    Omitted outer boundaries use the full layout extents
    """
    left = float(left_throat_z_m)
    right = float(right_throat_z_m)
    if not left < layout.midplane_z_m < right:
        raise ValueError("fitted throats must bracket the geometry derived midplane")
    if left <= layout.full_device_domain.z_min_m or right >= layout.full_device_domain.z_max_m:
        raise ValueError("fitted throats must lie strictly inside the full device domain")
    left_wall = layout.full_device_domain.z_min_m if left_boundary_z_m is None else float(left_boundary_z_m)
    right_wall = layout.full_device_domain.z_max_m if right_boundary_z_m is None else float(right_boundary_z_m)
    if not layout.full_device_domain.z_min_m <= left_wall < left:
        raise ValueError("selected left plasma boundary must lie between full device edge and left throat")
    if not right < right_wall <= layout.full_device_domain.z_max_m:
        raise ValueError("selected right plasma boundary must lie between right throat and full device edge")
    # Selected material interfaces bound the plasma domain within the full layout extents
    source_device = AxialDomain("full_device_domain", left_wall, right_wall, "selected plasma facing material interfaces bounding the accessible device")
    confined = AxialDomain("confined_orbit_domain", left, right, "independent fitted magnetic throats")
    left_expander = AxialDomain("left_expander_domain", left_wall, left, "selected left plasma facing boundary to left fitted throat")
    right_expander = AxialDomain("right_expander_domain", right, right_wall, "right fitted throat to selected right plasma facing boundary")
    fusion = AxialDomain("fusion_domain", left_wall, right_wall, "full local population coexistence support through both expanders")
    neutron = replace(fusion, name="neutron_source_domain")
    export = replace(source_device, name="openmc_export_domain")

    return DeviceDomains(
        midplane_z_m=layout.midplane_z_m,
        full_device_domain=source_device,
        central_cell_plasma_domain=layout.central_cell_plasma_domain,
        confined_orbit_domain=confined,
        left_expander_domain=left_expander,
        right_expander_domain=right_expander,
        fusion_domain=fusion,
        neutron_source_domain=neutron,
        openmc_export_domain=export,
    )
