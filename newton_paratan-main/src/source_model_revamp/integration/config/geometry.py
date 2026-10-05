"""Magnetic geometry, field fit, and plasma boundary configuration"""
from __future__ import annotations
from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field
from typing import Any
import numpy as np
from source_model_revamp.constants import pi, mu_0, CM_TO_M
from source_model_revamp.geometry.device_domains import ParaTANLayout, derive_paratan_layout
from source_model_revamp.geometry.plasma_boundaries import BoundaryOverride
from source_model_revamp.integration.config.common import canonical_model_name, _mapping, _check_unknown, _number, _integer, _sum_thickness_cm, _coerce_float_list
from source_model_revamp.integration.config.constants import *

@dataclass(frozen=True)
class FieldStrengthConfig:
    """
    Requested on axis magnetic field strengths
    
    `mirror_ratio = mirror_throat_B_T / midplane_B_T`
    """
    midplane_B_T: float
    mirror_throat_B_T: float

    @property
    def mirror_ratio(self) -> float:
        """Return `mirror_throat_B_T / midplane_B_T`"""
        return float(self.mirror_throat_B_T / self.midplane_B_T)

    @classmethod
    def from_mapping(cls, data: Mapping[str, Any], *, strict: bool = False) -> "FieldStrengthConfig":
        """Parse positive midplane and throat fields and require the throat field to exceed the midplane field"""
        allowed = {"midplane_B_T", "mirror_throat_B_T"}
        _check_unknown(data, allowed, "source_model.geometry.field_strength", strict)
        midplane = _number(data.get("midplane_B_T"), "geometry.field_strength.midplane_B_T", positive=True)
        throat = _number(data.get("mirror_throat_B_T"), "geometry.field_strength.mirror_throat_B_T", positive=True)
        if throat <= midplane:
            raise ValueError("geometry.field_strength.mirror_throat_B_T must be greater than midplane_B_T")

        return cls(midplane_B_T=midplane, mirror_throat_B_T=throat)

@dataclass(frozen=True)
class DerivedCircularCoilConfig:
    """
    Effective thin circular coil used by the on axis ParaTAN field fit
    
    `z_center_m` and `radius_m` are in m and `amp_turns` is the fitted effective current turn product
    """
    name: str
    z_center_m: float
    radius_m: float
    amp_turns: float
    group: str

def _paratan_central_outer_radius_cm(root: Mapping[str, Any]) -> float:
    """Return the ParaTAN central cell outer radius after vessel structure and central cell layers in cm"""
    vv = _mapping(root.get("vacuum_vessel"), "vacuum_vessel")
    cc = _mapping(root.get("central_cell"), "central_cell")
    radius_cm = float(vv.get("central_radius", 0.0)) + _sum_thickness_cm(vv.get("structure", {}))
    layers = cc.get("layers", ())
    if isinstance(layers, Sequence) and not isinstance(layers, (str, bytes)):
        radius_cm += sum(float(layer.get("thickness", 0.0)) for layer in layers if isinstance(layer, Mapping))

    return float(radius_cm)

def _paratan_bottleneck_outer_radius_cm(root: Mapping[str, Any]) -> float:
    """Return the ParaTAN bottleneck outer radius after vessel structure in cm"""
    vv = _mapping(root.get("vacuum_vessel"), "vacuum_vessel")

    return float(vv.get("bottleneck_radius", vv.get("central_radius", 0.0))) + _sum_thickness_cm(vv.get("structure", {}))

def _unit_paratan_coils_from_geometry(root: Mapping[str, Any], layout: ParaTANLayout | None = None) -> tuple[DerivedCircularCoilConfig, ...]:
    """
    Build unit current thin coil representatives from ParaTAN LF and HF geometry
    
    The returned coils retain the geometry derived centers and effective radii while setting `amp_turns = 1`
    """
    coils: list[DerivedCircularCoilConfig] = []
    lf = _mapping(root.get("lf_coil"), "lf_coil")
    hf = _mapping(root.get("hf_coil"), "hf_coil")
    inner_dimensions = _mapping(lf.get("inner_dimensions"), "lf_coil.inner_dimensions")
    shell = _mapping(lf.get("shell_thicknesses"), "lf_coil.shell_thicknesses")
    positions_cm = _coerce_float_list(lf.get("positions", ()))
    if not positions_cm:
        raise ValueError("lf_coil.positions must contain at least one axial position in cm")
    inner_radius_cm = _paratan_central_outer_radius_cm(root)
    midplane_cm = float(_mapping(root.get("vacuum_vessel"), "vacuum_vessel").get("axial_midplane", 0.0))
    radial_thickness_cm = float(inner_dimensions.get("radial_thickness", 0.0))
    shell_front_cm = float(shell.get("front", 0.0))
    lf_radius_cm = inner_radius_cm + shell_front_cm + 0.5 * radial_thickness_cm
    for i, z_cm in enumerate(positions_cm):
        coils.append(DerivedCircularCoilConfig(name=f"lf_coil_{i}", z_center_m=(midplane_cm + float(z_cm)) * CM_TO_M, radius_m=lf_radius_cm * CM_TO_M, amp_turns=1.0, group="lf_coil"))
    magnet = _mapping(hf.get("magnet"), "hf_coil.magnet")
    shield = _mapping(hf.get("shield"), "hf_coil.shield")
    casing_radial_cm = _sum_thickness_cm(hf.get("casing_layers", ()))
    shield_radial = shield.get("radial_thickness", 0.0)
    if isinstance(shield_radial, Sequence) and not isinstance(shield_radial, (str, bytes)):
        shield_radial_inner_cm = float(shield_radial[0]) if len(shield_radial) else 0.0
    else:
        shield_radial_inner_cm = float(shield_radial)
    magnet_inner_cm = _paratan_bottleneck_outer_radius_cm(root) + casing_radial_cm + shield_radial_inner_cm
    hf_radius_cm = magnet_inner_cm + 0.5 * float(magnet.get("radial_thickness", 0.0))
    resolved_layout = derive_paratan_layout(root) if layout is None else layout
    for name, z_m in (("hf_coil_left", resolved_layout.left_hf_coil_center_m), ("hf_coil_right", resolved_layout.right_hf_coil_center_m)):
        coils.append(DerivedCircularCoilConfig(name=name, z_center_m=float(z_m), radius_m=hf_radius_cm * CM_TO_M, amp_turns=1.0, group="hf_coil"))
        
    return tuple(coils)

def _coil_field_Bz_on_axis_T(z_m: np.ndarray, coils: Sequence[DerivedCircularCoilConfig]) -> np.ndarray:
    """Return the summed on axis `B_z` from the effective circular coils in T"""
    z = np.asarray(z_m, dtype=float)
    total = np.zeros_like(z, dtype=float)
    for coil in coils:
        R = float(coil.radius_m)
        dz = z - float(coil.z_center_m)
        total += mu_0 * float(coil.amp_turns) * R**2 / (2.0 * (R**2 + dz**2) ** 1.5)

    return total

def _scale_paratan_coils_to_field_strength( unit_coils: Sequence[DerivedCircularCoilConfig], field_strength: FieldStrengthConfig, search_z_min_m: float, search_z_max_m: float, reference_z_m: float, grid_points: int) -> tuple[tuple[DerivedCircularCoilConfig, ...], dict[str, Any]]:
    """
    Fit relative HF to LF current scaling to the requested mirror ratio and then scale all coils to the requested midplane field
    
    The mirror ratio fit is evaluated on the configured axial grid before the overall field amplitude is applied
    """
    if not unit_coils:
        raise ValueError("paratan_coil_field requires at least one derived coil")
    z_min = float(search_z_min_m)
    z_max = float(search_z_max_m)
    if not np.isfinite(z_min) or not np.isfinite(z_max) or z_max <= z_min:
        raise ValueError("ParaTAN coil field search bounds must be finite and increasing")
    z = np.linspace(z_min, z_max, int(grid_points))
    lf_unit = tuple(c for c in unit_coils if c.group == "lf_coil")
    hf_unit = tuple(c for c in unit_coils if c.group == "hf_coil")
    other_unit = tuple(c for c in unit_coils if c.group not in {"lf_coil", "hf_coil"})
    B_lf = _coil_field_Bz_on_axis_T(z, lf_unit) if lf_unit else np.zeros_like(z)
    B_hf = _coil_field_Bz_on_axis_T(z, hf_unit) if hf_unit else np.zeros_like(z)
    B_other = _coil_field_Bz_on_axis_T(z, other_unit) if other_unit else np.zeros_like(z)
    B_lf_ref = float(_coil_field_Bz_on_axis_T(np.array([reference_z_m]), lf_unit)[0]) if lf_unit else 0.0
    B_hf_ref = float(_coil_field_Bz_on_axis_T(np.array([reference_z_m]), hf_unit)[0]) if hf_unit else 0.0
    B_other_ref = float(_coil_field_Bz_on_axis_T(np.array([reference_z_m]), other_unit)[0]) if other_unit else 0.0
    target_ratio = field_strength.mirror_ratio
    def ratio_for(hf_scale: float) -> float:
        """Return the grid maximum to reference field ratio for one HF scale"""
        field = B_lf + float(hf_scale) * B_hf + B_other
        ref = B_lf_ref + float(hf_scale) * B_hf_ref + B_other_ref
        if ref <= 0.0:
            return 0.0
        return float(np.max(np.abs(field)) / abs(ref))
    # First determine the relative HF to LF scaling required by the mirror ratio
    low = 0.0
    high = 1.0
    while ratio_for(high) < target_ratio and high < 1.0e12:
        high *= 2.0
    if ratio_for(high) < target_ratio:
        raise ValueError("ParaTAN LF/HF coil geometry cannot reach requested mirror_throat_B_T on the fit grid")
    for _ in range(80):
        mid = 0.5 * (low + high)
        if ratio_for(mid) < target_ratio:
            low = mid
        else:
            high = mid
    hf_scale = high
    # Then apply one common amplitude scale so the reference field matches B0
    raw_ref = B_lf_ref + hf_scale * B_hf_ref + B_other_ref
    if raw_ref == 0.0:
        raise ValueError("Derived ParaTAN coil geometry produced zero midplane field before scaling")
    overall_scale = field_strength.midplane_B_T / raw_ref
    scaled: list[DerivedCircularCoilConfig] = []
    for coil in unit_coils:
        group_scale = hf_scale if coil.group == "hf_coil" else 1.0
        scaled.append(DerivedCircularCoilConfig(name=coil.name, z_center_m=coil.z_center_m, radius_m=coil.radius_m, amp_turns=float(overall_scale * group_scale), group=coil.group))
    B_scaled = _coil_field_Bz_on_axis_T(z, scaled)
    B_ref = float(_coil_field_Bz_on_axis_T(np.array([reference_z_m]), scaled)[0])
    imax = int(np.argmax(np.abs(B_scaled)))
    B_max = float(abs(B_scaled[imax]))
    metadata = {
        "coil_current_fit_model": FIELD_STRENGTH_FIT_MODEL,
        "coil_effective_amp_turns_are_internal_scaling": True,
        "requested_midplane_B_T": float(field_strength.midplane_B_T),
        "requested_mirror_throat_B_T": float(field_strength.mirror_throat_B_T),
        "requested_mirror_ratio": float(field_strength.mirror_ratio),
        "fitted_midplane_B_T": float(B_ref),
        "fitted_mirror_throat_B_T": float(B_max),
        "fitted_mirror_ratio": float(B_max / abs(B_ref)),
        "fitted_mirror_throat_z_m": float(z[imax]),
        "hf_to_lf_unit_amp_turn_scale": float(hf_scale),
        "overall_amp_turn_scale": float(overall_scale),
        "coil_field_reference_z_m": float(reference_z_m),
        "coil_field_fit_z_min_m": z_min,
        "coil_field_fit_z_max_m": z_max,
        "coil_field_grid_points": int(grid_points),
        "coil_field_shape_source": "paratan_lf_hf_geometry",
        "derived_lf_coil_z_m": [float(c.z_center_m) for c in scaled if c.group == "lf_coil"],
        "derived_hf_coil_z_m": [float(c.z_center_m) for c in scaled if c.group == "hf_coil"],
        "derived_lf_coil_radius_m": [float(c.radius_m) for c in scaled if c.group == "lf_coil"],
        "derived_hf_coil_radius_m": [float(c.radius_m) for c in scaled if c.group == "hf_coil"],
        "derived_effective_amp_turns": [float(c.amp_turns) for c in scaled],
    }

    return tuple(scaled), metadata

@dataclass(frozen=True)
class GeometryConfig:
    """
    Magnetic field and reduced flux tube geometry controls
    
    ParaTAN linked geometry derives the full device domain and effective coils from the root geometry mapping
    """
    model: str = "paratan_coil_field"
    half_length_m: float = 1.0
    plasma_radius_m: float = 0.10
    axial_bins: int = 12
    field_strength: FieldStrengthConfig = field(default_factory=FieldStrengthConfig)
    transition_width: float = 0.10
    shape_exponent: float = 1.30
    effective_coils: tuple[DerivedCircularCoilConfig, ...] = ()
    coil_field_grid_points: int = 1001
    coil_fit_metadata: dict[str, Any] = field(default_factory=dict)
    domain_mode: str = "geometry_linked"
    paratan_layout: ParaTANLayout | None = None

    @property
    def midplane_area_m2(self) -> float:
        """Return the circular midplane plasma area `π a²` in m²"""
        return float(pi * self.plasma_radius_m**2)

    @property
    def B0_T(self) -> float:
        """Return the configured midplane field `B0` in T"""
        return float(self.field_strength.midplane_B_T)

    @property
    def mirror_ratio(self) -> float:
        """Return the configured magnetic mirror ratio"""
        return float(self.field_strength.mirror_ratio)

    @classmethod
    def from_mapping(cls, data: Mapping[str, Any], *, root: Mapping[str, Any] | None = None, strict: bool = False) -> "GeometryConfig":
        """
        Parse the magnetic geometry and resolve the selected device domain
        
        For `paratan_coil_field`, the ParaTAN layout supplies the full axial domain and the LF and HF coil geometry is fitted to the requested field strengths
        """
        root = root or {}
        allowed = {"model", "domain_mode", "half_length_m", "manual_confined_half_length_m", "plasma_radius_m", "axial_bins", "field_strength", "transition_width", "shape_exponent", "coil_field_grid_points"}
        _check_unknown(data, allowed, "source_model.geometry", strict)
        model = canonical_model_name(data.get("model"), GEOMETRY_MODEL_ALIASES, "paratan_coil_field")
        field_strength = FieldStrengthConfig.from_mapping(_mapping(data.get("field_strength"), "source_model.geometry.field_strength"), strict=strict)
        default_domain_mode = "geometry_linked" if model == "paratan_coil_field" else "manual_benchmark"
        domain_mode = str(data.get("domain_mode", default_domain_mode))
        if domain_mode not in {"geometry_linked", "manual_benchmark"}:
            raise ValueError("geometry.domain_mode must be geometry_linked or manual_benchmark")
        if model == "paratan_coil_field" and domain_mode != "geometry_linked":
            raise ValueError("paratan_coil_field requires geometry.domain_mode: geometry_linked")
        manual_half = data.get("manual_confined_half_length_m", data.get("half_length_m"))
        half_length_m = _number(manual_half, "geometry.manual_confined_half_length_m", positive=True, default=1.0)
        plasma_radius_m = _number(data.get("plasma_radius_m"), "geometry.plasma_radius_m", positive=True, default=0.10)
        grid_points = _integer(data.get("coil_field_grid_points"), "geometry.coil_field_grid_points", minimum=3, default=1001)
        effective_coils: tuple[DerivedCircularCoilConfig, ...] = ()
        coil_fit_metadata: dict[str, Any] = {}
        paratan_layout: ParaTANLayout | None = None
        if model == "paratan_coil_field":
            paratan_layout = derive_paratan_layout(root)
            half_length_m = 0.5 * paratan_layout.full_device_domain.length_m
            unit_coils = _unit_paratan_coils_from_geometry(root, paratan_layout)
            effective_coils, coil_fit_metadata = _scale_paratan_coils_to_field_strength(unit_coils=unit_coils, field_strength=field_strength, search_z_min_m=paratan_layout.vacuum_vessel_domain.z_min_m, search_z_max_m=paratan_layout.vacuum_vessel_domain.z_max_m, reference_z_m=paratan_layout.midplane_z_m, grid_points=grid_points)

        return cls(
            model=model,
            half_length_m=half_length_m,
            plasma_radius_m=plasma_radius_m,
            axial_bins=_integer(data.get("axial_bins"), "geometry.axial_bins", minimum=1, default=12),
            field_strength=field_strength,
            transition_width=_number(data.get("transition_width"), "geometry.transition_width", positive=True, default=0.10),
            shape_exponent=_number(data.get("shape_exponent"), "geometry.shape_exponent", positive=True, default=1.30),
            effective_coils=effective_coils,
            coil_field_grid_points=grid_points,
            coil_fit_metadata=coil_fit_metadata,
            domain_mode=domain_mode,
            paratan_layout=paratan_layout,
        )

@dataclass(frozen=True)
class PlasmaBoundaryConfig:
    """Automatic plasma facing boundary detection controls plus explicit surface overrides"""
    detection_mode: str = "automatic_with_overrides"
    default_material_interface_role: str = "absorbing"
    default_electrical_model: str = "floating"
    overrides: tuple[BoundaryOverride, ...] = ()

    @classmethod
    def from_mapping(cls, data: Mapping[str, Any], *, strict: bool = False) -> "PlasmaBoundaryConfig":
        """Parse boundary defaults and convert each selector override into a BoundaryOverride"""
        allowed = {"detection_mode", "default_material_interface_role", "default_electrical_model", "overrides"}
        _check_unknown(data, allowed, "source_model.plasma_boundaries", strict)
        detection_mode = str(data.get("detection_mode", "automatic_with_overrides"))
        default_role = str(data.get("default_material_interface_role", "absorbing"))
        default_electrical = str(data.get("default_electrical_model", "floating"))
        if detection_mode != "automatic_with_overrides":
            raise ValueError("plasma_boundaries.detection_mode must be automatic_with_overrides")
        if default_role not in {"absorbing", "collector", "end_ring", "reflecting", "excluded"}:
            raise ValueError(f"unsupported default material-interface role {default_role!r}")
        if default_electrical not in {"floating", "grounded"}:
            raise ValueError("plasma_boundaries.default_electrical_model must be floating or grounded, prescribed bias requires a surface override with prescribed_bias_V")
        raw_overrides = data.get("overrides", ())
        if not isinstance(raw_overrides, Sequence) or isinstance(raw_overrides, (str, bytes)):
            raise ValueError("plasma_boundaries.overrides must be a sequence")
        overrides: list[BoundaryOverride] = []
        for index, raw in enumerate(raw_overrides):
            item = _mapping(raw, f"source_model.plasma_boundaries.overrides[{index}]")
            _check_unknown(item, {"selector", "particle_role", "electrical_model", "prescribed_bias_V"}, f"source_model.plasma_boundaries.overrides[{index}]", strict)
            selector = item.get("selector")
            if isinstance(selector, Mapping):
                _check_unknown(selector, {"surface_id", "component"}, f"source_model.plasma_boundaries.overrides[{index}].selector", strict)
                values = [str(value) for value in selector.values() if value is not None]
                if len(values) != 1:
                    raise ValueError("each boundary override selector must provide exactly one surface_id or component")
                selector = values[0]
            if selector is None:
                raise ValueError("each plasma boundary override requires selector")
            overrides.append(BoundaryOverride(selector=str(selector), particle_role=None if item.get("particle_role") is None else str(item["particle_role"]), electrical_model=None if item.get("electrical_model") is None else str(item["electrical_model"]), prescribed_bias_V=None if item.get("prescribed_bias_V") is None else _number(item["prescribed_bias_V"], "plasma_boundaries.prescribed_bias_V"),))
        
        return cls(
            detection_mode=detection_mode,
            default_material_interface_role=default_role,
            default_electrical_model=default_electrical,
            overrides=tuple(overrides),
        )
