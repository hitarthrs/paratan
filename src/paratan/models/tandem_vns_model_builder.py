"""Tandem VNS geometry with modular per-ID SCLL / DCLL test modules.

Build order::

    vacuum vessel → first wall → central cylinders (test envelope + end plugs)
    → test modules → HF coils → LF coils → end cells

Three input files::

    1. tandem / VNS machine YAML   — VV, FW, coils, end cells, test_region
    2. (same file) placement       — ``test_modules.layout`` + ``geometry_file``
    3. testing_module_geometry.yaml — ``defaults`` + per-ID ``modules`` solids

Example::

    path = Path("input_files/testing_module_input.yaml")
    with path.open() as f:
        data = yaml.safe_load(f)
    builder = build_tandem_vns_universe(data, input_path=path)
"""

from __future__ import annotations

from copy import deepcopy
from functools import reduce
from operator import or_
from pathlib import Path
from typing import Any

import openmc
import yaml

from src.paratan.geometry.core import annular_shell_region
from src.paratan.materials import material as default_materials
from src.paratan.models.tandem_model_builder import (
    TandemMachineBuilder,
    parse_machine_input,
)
from src.paratan.models.testing_module_builder import (
    BulkChannelRow,
    build_dcll_fnsf_module,
    build_scll_module,
)

_MODULE_ID_STRIDE = 10000
_MODULE_GAP_CELL_ID = 8900

_SCLL_DEFAULTS = {
    "n_fw_channels": 20,
    "channel_radial_thickness": 2.0,
    "front_thickness": 0.4,
    "back_thickness": 0.4,
    "axial_gap": 0.0,
    "axial_end_thickness": 0.0,
    "armor_thickness": 0.2,
    "divider_thickness": 0.6,
    "bulk_rows": (),
    "bulk_front_thickness": 0.5,
    "bulk_back_thickness": 1.0,
    "bulk_radial_gap": 1.0,
}

_DCLL_DEFAULTS = {
    "front_pbli_thickness": 37.15,
    "rear_pbli_thickness": 37.15,
    "n_axial_bins": 4,
    "fw_steel_front": 0.4,
    "fw_he_thickness": 2.0,
    "fw_steel_back": 0.4,
    "duct_steel_thickness": 0.5,
    "pbli_gap_thickness": 0.1,
    "fci_thickness": 0.5,
    "midplate_he_thickness": 0.8,
    "back_he_thickness": 0.8,
    "n_fw_he_channels": 12,
    "n_midplate_he_channels": 8,
    "n_back_he_channels": 8,
    "he_rib_angle_deg": 1.0,
    "he_side_wall_angle_deg": 2.0,
}


def load_testing_module_geometry(
    geometry_file: str | Path,
    *,
    base_dir: str | Path | None = None,
) -> dict[str, Any]:
    """Load the third-file module geometry YAML (``defaults`` + ``modules``)."""
    path = Path(geometry_file)
    if not path.is_absolute():
        search_roots = []
        if base_dir is not None:
            search_roots.append(Path(base_dir))
        # Repo input_files/ as a convenient fallback when only a basename is given.
        search_roots.append(Path(__file__).resolve().parents[3] / "input_files")
        search_roots.append(Path.cwd())
        for root in search_roots:
            candidate = root / path
            if candidate.is_file():
                path = candidate
                break
        else:
            raise FileNotFoundError(
                f"testing-module geometry file not found: {geometry_file!r} "
                f"(searched {[str(r) for r in search_roots]})"
            )
    if not path.is_file():
        raise FileNotFoundError(f"testing-module geometry file not found: {path}")
    with path.open() as stream:
        data = yaml.safe_load(stream) or {}
    if not isinstance(data, dict):
        raise ValueError(f"{path} must contain a YAML mapping")
    if "modules" not in data:
        raise ValueError(f"{path} must define a top-level 'modules' mapping")
    return data


def resolve_test_modules_block(
    test_modules: dict,
    *,
    base_dir: str | Path | None = None,
) -> dict:
    """Merge machine-file layout with optional third-file geometry.

    If ``geometry_file`` is set, load ``defaults`` / ``modules`` from that file
    and keep ``layout`` (and optional ``count``) from the machine YAML.
    Inline ``modules`` in the machine file remains supported.
    """
    if not isinstance(test_modules, dict):
        raise ValueError("central_cell.test_modules must be a mapping")

    merged = dict(test_modules)
    geometry_file = merged.pop("geometry_file", None)
    if geometry_file is None:
        return merged

    geom = load_testing_module_geometry(geometry_file, base_dir=base_dir)
    # Geometry file owns solids; machine file owns placement.
    if "defaults" in geom:
        merged["defaults"] = geom["defaults"]
    if "modules" in geom:
        if "modules" in merged:
            raise ValueError(
                "test_modules.modules cannot be set both inline and via geometry_file"
            )
        merged["modules"] = geom["modules"]
    # Allow geometry file to supply a default layout only if machine omits one.
    if "layout" not in merged and "layout" in geom:
        merged["layout"] = geom["layout"]
    return merged


def parse_vns_machine_input(
    input_data,
    material_ns=default_materials,
    *,
    base_dir: str | Path | None = None,
):
    """Parse tandem inputs; central blanket becomes a vacuum test envelope."""
    central = input_data["central_cell"]
    if "blanket" in central:
        raise ValueError(
            "VNS input uses central_cell.test_region, not central_cell.blanket"
        )

    test_region = central["test_region"]
    for key in ("axial_length", "radial_thickness"):
        value = test_region[key]
        if isinstance(value, bool) or not isinstance(value, (int, float)) or value <= 0:
            raise ValueError(f"central_cell.test_region.{key} must be positive")

    test_modules = resolve_test_modules_block(
        central["test_modules"], base_dir=base_dir
    )
    module_params = parse_test_module_params(test_modules)

    tandem_input = deepcopy(input_data)
    tandem_input["central_cell"]["blanket"] = {
        "axial_length": test_region["axial_length"],
        "layers": [
            {"thickness": test_region["radial_thickness"], "material": "vacuum"}
        ],
    }
    machine_params = parse_machine_input(tandem_input, material_ns)
    return machine_params, module_params


def _normalize_type(value: str) -> str:
    if value in ("SCLL",):
        return "SCLL"
    if value in ("DCLL", "DCLL_FNSF"):
        return "DCLL_FNSF"
    raise ValueError(f"module type must be SCLL or DCLL_FNSF; got {value!r}")


def _merge_dict(base: dict, override: dict | None) -> dict:
    merged = dict(base)
    if override:
        for key, value in override.items():
            if key == "type":
                continue
            merged[key] = value
    return merged


def _parse_bulk_rows(raw) -> tuple:
    if not raw:
        return ()
    return tuple(
        BulkChannelRow(
            n_channels=int(row["n_channels"]),
            radial_thickness=float(row["radial_thickness"]),
        )
        for row in raw
    )


def _finalize_scll(block: dict) -> dict:
    required = (
        "n_fw_channels",
        "channel_radial_thickness",
        "front_thickness",
        "back_thickness",
    )
    for key in required:
        if key not in block:
            raise ValueError(f"SCLL module missing required key {key!r}")
    return {
        "type": "SCLL",
        "n_fw_channels": int(block["n_fw_channels"]),
        "channel_radial_thickness": float(block["channel_radial_thickness"]),
        "front_thickness": float(block["front_thickness"]),
        "back_thickness": float(block["back_thickness"]),
        "axial_gap": float(block.get("axial_gap", 0.0)),
        "axial_end_thickness": float(block.get("axial_end_thickness", 0.0)),
        "armor_thickness": float(block.get("armor_thickness", 0.2)),
        "divider_thickness": float(block.get("divider_thickness", 0.6)),
        "bulk_rows": _parse_bulk_rows(block.get("bulk_rows", ())),
        "bulk_front_thickness": float(block.get("bulk_front_thickness", 0.5)),
        "bulk_back_thickness": float(block.get("bulk_back_thickness", 1.0)),
        "bulk_radial_gap": float(block.get("bulk_radial_gap", 1.0)),
    }


def _finalize_dcll(block: dict) -> dict:
    required = ("front_pbli_thickness", "rear_pbli_thickness")
    for key in required:
        if key not in block:
            raise ValueError(f"DCLL module missing required key {key!r}")
    return {
        "type": "DCLL_FNSF",
        "front_pbli_thickness": float(block["front_pbli_thickness"]),
        "rear_pbli_thickness": float(block["rear_pbli_thickness"]),
        "n_axial_bins": int(block.get("n_axial_bins", 4)),
        "fw_steel_front": float(block.get("fw_steel_front", 0.4)),
        "fw_he_thickness": float(block.get("fw_he_thickness", 2.0)),
        "fw_steel_back": float(block.get("fw_steel_back", 0.4)),
        "duct_steel_thickness": float(block.get("duct_steel_thickness", 0.5)),
        "pbli_gap_thickness": float(block.get("pbli_gap_thickness", 0.1)),
        "fci_thickness": float(block.get("fci_thickness", 0.5)),
        "midplate_he_thickness": float(block.get("midplate_he_thickness", 0.8)),
        "back_he_thickness": float(block.get("back_he_thickness", 0.8)),
        "n_fw_he_channels": int(block.get("n_fw_he_channels", 12)),
        "n_midplate_he_channels": int(block.get("n_midplate_he_channels", 8)),
        "n_back_he_channels": int(block.get("n_back_he_channels", 8)),
        "he_rib_angle_deg": float(block.get("he_rib_angle_deg", 1.0)),
        "he_side_wall_angle_deg": float(block.get("he_side_wall_angle_deg", 2.0)),
    }


def parse_test_module_params(test_modules: dict) -> dict:
    """Parse ID-keyed module specs + layout for the VNS builder.

    Supports optional ``defaults.SCLL`` / ``defaults.DCLL_FNSF`` (or ``defaults.dcll``)
    that are shallow-merged under each module before validation. Legacy single
    ``scll`` / ``dcll`` + type-string layout is still accepted and expanded into
    unique IDs.
    """
    if not isinstance(test_modules, dict):
        raise ValueError("central_cell.test_modules must be a mapping")

    # ---- modern path: modules: {ID: {type, ...}} + layout: [{left, right}, ...]
    if "modules" in test_modules:
        return _parse_modular_modules(test_modules)

    # ---- legacy path: shared scll/dcll blocks + type strings in layout
    return _parse_legacy_modules(test_modules)


def _parse_modular_modules(test_modules: dict) -> dict:
    defaults_raw = test_modules.get("defaults") or {}
    scll_defaults = _merge_dict(
        _SCLL_DEFAULTS,
        defaults_raw.get("SCLL") or defaults_raw.get("scll"),
    )
    dcll_defaults = _merge_dict(
        _DCLL_DEFAULTS,
        defaults_raw.get("DCLL_FNSF")
        or defaults_raw.get("DCLL")
        or defaults_raw.get("dcll"),
    )

    modules_raw = test_modules["modules"]
    if not isinstance(modules_raw, dict) or not modules_raw:
        raise ValueError("test_modules.modules must be a non-empty mapping of ID → spec")

    modules = {}
    for module_id, spec in modules_raw.items():
        if not isinstance(module_id, str) or not module_id.strip():
            raise ValueError(f"module ID must be a non-empty string; got {module_id!r}")
        if not isinstance(spec, dict) or "type" not in spec:
            raise ValueError(f"modules[{module_id!r}] must be a mapping with a type")
        module_type = _normalize_type(str(spec["type"]))
        if module_type == "SCLL":
            merged = _merge_dict(scll_defaults, spec)
            modules[module_id] = _finalize_scll(merged)
        else:
            merged = _merge_dict(dcll_defaults, spec)
            modules[module_id] = _finalize_dcll(merged)
        modules[module_id]["id"] = module_id

    layout_raw = test_modules.get("layout")
    if not layout_raw:
        raise ValueError("test_modules.layout is required when using modules:")
    if len(layout_raw) < 1:
        raise ValueError("test_modules.layout must have at least one station")

    layout = []
    used_ids = set()
    for i, station in enumerate(layout_raw):
        left_id = station.get("left")
        right_id = station.get("right")
        for side, mid in (("left", left_id), ("right", right_id)):
            if mid not in modules:
                raise ValueError(
                    f"test_modules.layout[{i}].{side}={mid!r} is not in modules"
                )
        if left_id == right_id:
            raise ValueError(
                f"test_modules.layout[{i}] reuses module ID {left_id!r} on both sides; "
                "each placement needs its own ID (duplicate the spec under a new key)"
            )
        if left_id in used_ids or right_id in used_ids:
            raise ValueError(
                f"test_modules.layout[{i}] reuses a module ID already placed; "
                "give each physical module a unique ID"
            )
        used_ids.add(left_id)
        used_ids.add(right_id)
        layout.append({"left": left_id, "right": right_id})

    count = 2 * len(layout)
    if "count" in test_modules and int(test_modules["count"]) != count:
        raise ValueError(
            f"test_modules.count={test_modules['count']} does not match "
            f"layout ({count} modules)"
        )

    unused = set(modules) - used_ids
    if unused:
        raise ValueError(
            f"modules defined but not placed in layout: {sorted(unused)}"
        )

    return {
        "count": count,
        "layout": layout,
        "modules": modules,
    }


def _parse_legacy_modules(test_modules: dict) -> dict:
    """Expand old shared-type layout into unique modular IDs."""
    count = test_modules.get("count")
    if isinstance(count, bool) or not isinstance(count, int) or count < 2:
        raise ValueError("central_cell.test_modules.count must be an integer >= 2")
    if count % 2 != 0:
        raise ValueError(
            "central_cell.test_modules.count must be even "
            "(left/right half-module pairs)"
        )

    n_stations = count // 2
    layout_raw = test_modules.get("layout")
    if layout_raw is None:
        layout_raw = [{"left": "SCLL", "right": "SCLL"} for _ in range(n_stations)]
    if len(layout_raw) != n_stations:
        raise ValueError(
            f"test_modules.layout must have {n_stations} stations; got {len(layout_raw)}"
        )

    scll_src = test_modules.get("scll") or test_modules.get("defaults", {}).get("SCLL")
    dcll_src = (
        test_modules.get("dcll")
        or test_modules.get("dcll_fnsf")
        or test_modules.get("defaults", {}).get("DCLL_FNSF")
    )
    scll_base = _finalize_scll(_merge_dict(_SCLL_DEFAULTS, scll_src))
    dcll_base = (
        _finalize_dcll(_merge_dict(_DCLL_DEFAULTS, dcll_src)) if dcll_src else None
    )

    modules = {}
    layout = []
    module_index = 0
    for station in layout_raw:
        station_ids = {}
        for side in ("left", "right"):
            module_type = _normalize_type(str(station[side]))
            module_id = f"M{module_index + 1}"
            if module_type == "SCLL":
                spec = dict(scll_base)
            else:
                if dcll_base is None:
                    raise ValueError(
                        "layout uses DCLL_FNSF but test_modules.dcll is missing"
                    )
                spec = dict(dcll_base)
            spec["id"] = module_id
            modules[module_id] = spec
            station_ids[side] = module_id
            module_index += 1
        layout.append(station_ids)

    return {"count": count, "layout": layout, "modules": modules}


def _estimate_scll_outer(r_fw: float, p: dict) -> float:
    fw = p["front_thickness"] + p["channel_radial_thickness"] + p["back_thickness"]
    bulk = 0.0
    if p["bulk_rows"]:
        bulk = (
            p["bulk_front_thickness"]
            + sum(row.radial_thickness for row in p["bulk_rows"])
            + (len(p["bulk_rows"]) - 1) * p["bulk_radial_gap"]
            + p["bulk_back_thickness"]
        )
    return r_fw + fw + p["divider_thickness"] + bulk


def _estimate_dcll_outer(r_inner: float, p: dict) -> float:
    fixed = (
        p["fw_steel_front"]
        + p["fw_he_thickness"]
        + p["fw_steel_back"]
        + 2
        * (
            p["duct_steel_thickness"]
            + p["pbli_gap_thickness"]
            + p["fci_thickness"]
        )
        + p["front_pbli_thickness"]
        + p["duct_steel_thickness"]
        + p["midplate_he_thickness"]
        + p["duct_steel_thickness"]
        + 2
        * (
            p["duct_steel_thickness"]
            + p["pbli_gap_thickness"]
            + p["fci_thickness"]
        )
        + p["rear_pbli_thickness"]
        + p["duct_steel_thickness"]
        + p["back_he_thickness"]
        + p["duct_steel_thickness"]
    )
    return r_inner + fixed


class TandemVNSMachineBuilder(TandemMachineBuilder):
    """Tandem machine with ID-keyed modular SCLL/DCLL test modules."""

    def __init__(self, *args, module_params, **kwargs):
        super().__init__(*args, **kwargs)
        self.module_params = module_params
        self.module_count = module_params["count"]
        self._universe.name = "TandemVNSMachine"
        self.test_region_cell = None
        self.module_builders = []
        self.module_builders_by_id = {}
        self.module_universe = None
        self._module_outer_radius = None

    def subtract_from_room(self, regions):
        super().subtract_from_room(regions)
        self._universe.cells[69].region = self._room_region

    def _facility_fw_outer_radius(self) -> float:
        base = float(self.fw_params["base_radii"]["central"])
        layers = self.fw_params["central_fw"]
        return base + sum(float(layer["thickness"]) for layer in layers)

    def _module_axial_layout(self):
        n_stations = self.module_count // 2
        cc_length = float(self.vv_params["central_cell_vv_parameters"].outer_length)
        module_length = cc_length / n_stations
        axial_centers = tuple(
            -cc_length / 2.0 + (i + 0.5) * module_length for i in range(n_stations)
        )
        return module_length, axial_centers

    def _estimate_max_module_outer(self) -> float:
        fw_outer = self._facility_fw_outer_radius()
        outers = []
        for spec in self.module_params["modules"].values():
            if spec["type"] == "SCLL":
                r_fw = fw_outer + spec["armor_thickness"]
                outers.append(_estimate_scll_outer(r_fw, spec))
            else:
                outers.append(_estimate_dcll_outer(fw_outer, spec))
        return max(outers)

    def sync_test_envelope_to_modules(self):
        """Resize the central vacuum envelope to the max module outer radius."""
        if self.fw_builder is None:
            raise RuntimeError("Call build_first_wall() before syncing the envelope")
        module_outer = self._estimate_max_module_outer()
        fw_outer = self._facility_fw_outer_radius()
        thickness = module_outer - fw_outer
        if thickness <= 0:
            raise ValueError(
                f"module outer {module_outer:.3f} cm must exceed FW outer {fw_outer:.3f} cm"
            )
        self.central_cyl_params["central"] = [
            {"thickness": thickness, "material": "vacuum"}
        ]
        self._module_outer_radius = module_outer
        return module_outer

    def build_central_cylinders(self):
        super().build_central_cylinders()
        self.test_region_cell = self.ccyl_builder.get_cells_by_region()["central"][0]
        self.test_region_cell.name = "central_test_region_placeholder"
        self._regions["test_region"] = self.ccyl_builder.get_regions_by_region()[
            "central"
        ][0]

    def _build_one_module(self, module_id, *, z0, orientation, module_length, module_index):
        spec = self.module_params["modules"][module_id]
        prefix = f"{module_id}_"
        offset = module_index * _MODULE_ID_STRIDE
        fw_outer = self._facility_fw_outer_radius()

        if spec["type"] == "SCLL":
            builder = build_scll_module(
                r_fw=fw_outer + spec["armor_thickness"],
                n_fw_channels=spec["n_fw_channels"],
                axial_length=module_length,
                channel_radial_thickness=spec["channel_radial_thickness"],
                front_thickness=spec["front_thickness"],
                back_thickness=spec["back_thickness"],
                z0=z0,
                orientation=orientation,
                axial_gap=spec["axial_gap"],
                axial_end_thickness=spec["axial_end_thickness"],
                armor_thickness=spec["armor_thickness"],
                divider_thickness=spec["divider_thickness"],
                bulk_rows=spec["bulk_rows"],
                bulk_front_thickness=spec["bulk_front_thickness"],
                bulk_back_thickness=spec["bulk_back_thickness"],
                bulk_radial_gap=spec["bulk_radial_gap"],
                material_ns=self.material_ns,
                cell_id_offset=offset,
                name_prefix=prefix,
            )
            r_inner = builder.r_armor()
        else:
            builder = build_dcll_fnsf_module(
                r_inner=fw_outer,
                axial_length=module_length,
                front_pbli_thickness=spec["front_pbli_thickness"],
                rear_pbli_thickness=spec["rear_pbli_thickness"],
                n_axial_bins=spec["n_axial_bins"],
                fw_steel_front=spec["fw_steel_front"],
                fw_he_thickness=spec["fw_he_thickness"],
                fw_steel_back=spec["fw_steel_back"],
                duct_steel_thickness=spec["duct_steel_thickness"],
                pbli_gap_thickness=spec["pbli_gap_thickness"],
                fci_thickness=spec["fci_thickness"],
                midplate_he_thickness=spec["midplate_he_thickness"],
                back_he_thickness=spec["back_he_thickness"],
                n_fw_he_channels=spec["n_fw_he_channels"],
                n_midplate_he_channels=spec["n_midplate_he_channels"],
                n_back_he_channels=spec["n_back_he_channels"],
                he_rib_angle_deg=spec["he_rib_angle_deg"],
                he_side_wall_angle_deg=spec["he_side_wall_angle_deg"],
                z0=z0,
                orientation=orientation,
                material_ns=self.material_ns,
                cell_id_offset=offset,
                name_prefix=prefix,
            )
            r_inner = builder.r_inner

        builder.module_id = module_id
        envelope = annular_shell_region(
            z0=z0,
            inner_radius=r_inner,
            radial_thickness=builder.module_outer_radius() - r_inner,
            axial_length=module_length,
            orientation=orientation,
        )
        return builder, envelope

    def build_test_modules(self):
        """Build each layout station from per-ID module specs."""
        if self.test_region_cell is None:
            raise RuntimeError(
                "Call build_central_cylinders() before build_test_modules()"
            )

        module_length, axial_centers = self._module_axial_layout()
        envelope_outer = float(
            self.ccyl_builder.get_radii_by_region()["central"][-1]
        )

        cells = []
        envelopes = []
        builders = []
        builders_by_id = {}
        module_index = 0
        for z0, station in zip(axial_centers, self.module_params["layout"]):
            for orientation in ("left", "right"):
                module_id = station[orientation]
                builder, envelope = self._build_one_module(
                    module_id,
                    z0=z0,
                    orientation=orientation,
                    module_length=module_length,
                    module_index=module_index,
                )
                if builder.module_outer_radius() > envelope_outer + 1e-6:
                    raise ValueError(
                        f"module {module_id!r} outer {builder.module_outer_radius():.3f} cm "
                        f"exceeds test envelope {envelope_outer:.3f} cm"
                    )
                cells.extend(builder.get_cells().values())
                envelopes.append(envelope)
                builders.append(builder)
                builders_by_id[module_id] = builder
                module_index += 1

        gap_cell = openmc.Cell(
            _MODULE_GAP_CELL_ID,
            name="space_between_test_modules",
            region=~reduce(or_, envelopes),
            fill=getattr(self.material_ns, "vacuum"),
        )
        cells.append(gap_cell)

        self.module_builders = builders
        self.module_builders_by_id = builders_by_id
        self.module_universe = openmc.Universe(name="VNS_TestModules", cells=cells)
        self.set_test_module_universe(self.module_universe)
        self._regions["test_modules"] = dict(builders_by_id)
        self._module_outer_radius = max(b.module_outer_radius() for b in builders)
        return self.module_universe

    def set_test_module_universe(self, universe):
        if self.test_region_cell is None:
            raise RuntimeError("Build central cylinders before installing modules")
        if not isinstance(universe, openmc.Universe):
            raise TypeError("Expected an OpenMC Universe")
        self.test_region_cell.fill = universe

    def build_lf_coils(self):
        """LF coils: central pack sits on the outermost testing-module radius."""
        if self._module_outer_radius is None:
            raise RuntimeError("Call build_test_modules() before build_lf_coils()")
        from src.paratan.models.tandem_model_builder import LFCoilBuilder

        self.lf_coil_builder = LFCoilBuilder(self.lf_coil_params)
        outer_radii = {
            key: self.ccyl_builder.get_radii_by_region()[key][-1]
            for key in ["left", "central", "right"]
        }
        outer_radii["central"] = float(self._module_outer_radius)
        self.lf_coil_builder.build_all_sections(
            outer_radii_by_region=outer_radii,
            midplanes=self.central_cyl_params["midplanes"],
            mat_ns=self.material_ns,
            add_cell_callback=self._add_cell,
        )
        self._regions["lf_coils"] = self.lf_coil_builder.get_regions_by_region()
        self.tally_builder.add_descriptors(self.lf_coil_builder.get_tally_descriptors())
        self.subtract_from_room(self.lf_coil_builder.get_full_region_list())

    def get_module_builders(self):
        return list(self.module_builders)

    def get_module_builder(self, module_id: str):
        return self.module_builders_by_id[module_id]

    def get_module_universe(self):
        return self.module_universe


def build_tandem_vns_universe(
    input_data,
    material_ns=default_materials,
    *,
    base_dir: str | Path | None = None,
    input_path: str | Path | None = None,
):
    """Build VNS: VV → FW → sync envelope → cylinders → modules → HF → LF → end.

    ``base_dir`` / ``input_path`` resolve ``test_modules.geometry_file`` relative
    to the machine YAML directory when a bare filename is given.
    """
    if base_dir is None and input_path is not None:
        base_dir = Path(input_path).resolve().parent
    machine_params, module_params = parse_vns_machine_input(
        input_data, material_ns, base_dir=base_dir
    )
    builder = TandemVNSMachineBuilder(
        *machine_params,
        material_ns=material_ns,
        module_params=module_params,
    )
    builder.build_vacuum_vessel()
    builder.build_first_wall()
    builder.sync_test_envelope_to_modules()
    builder.build_central_cylinders()
    builder.build_test_modules()
    builder.build_hf_coils()
    builder.build_lf_coils()
    builder.build_end_cells()
    return builder
