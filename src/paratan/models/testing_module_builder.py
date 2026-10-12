"""Testing-module builders for the tandem VNS central cell (classic style).

Same geometry as ``testing_module_model_builder``, restyled to match
``simple_model_builder`` / ``tandem_model_builder``: star imports, ``__init__``
classes, imperative ``build_*`` methods. Lengths are centimetres.

Run from the repository root::

    MPLCONFIGDIR=/tmp/mpl-vns .openmc-env/bin/python src/paratan/models/testing_module_builder.py
"""

import math

import numpy as np
import openmc

from src.paratan.geometry.core import *
from types import SimpleNamespace
import src.paratan.materials.material as m


# Geometry-only row spec; same name used in SCLL demos / THM notes.
BulkChannelRow = ChannelRowSpec

# Stable OpenMC cell IDs (easy to read in --geometry-debug / XML).
# SCLL 8xxx · DCLL 9xxx · debug "outside" uses 8999 / 9999.
SCLL_ID_ARMOR = 8001
SCLL_ID_FW_STRUCTURE = 8002
SCLL_ID_DIVIDER = 8003
SCLL_ID_BULK_STRUCTURE = 8004
SCLL_ID_FW_CHANNEL0 = 8100          # fw_channel_z{k} -> 8100 + k
SCLL_ID_BULK_ROW0 = 8200            # bulk_r{r}_z{k} -> 8200 + 100*r + k
SCLL_ID_OUTSIDE = 8999

DCLL_ID_LAYER0 = 9000               # solid radial shells, assigned in build order
DCLL_ID_FW_WEBS = 9081
DCLL_ID_MIDPLATE_WEBS = 9082
DCLL_ID_BACK_WEBS = 9083
DCLL_ID_FW_HE0 = 9100               # fw_helium_az{k} -> 9100 + k
DCLL_ID_MID_HE0 = 9200
DCLL_ID_BACK_HE0 = 9300
DCLL_ID_FRONT_PBLI0 = 9400          # front_pbli_z{k} -> 9400 + k
DCLL_ID_REAR_PBLI0 = 9500
DCLL_ID_OUTSIDE = 9999


def _cell(cell_id, name, region, fill, *, id_offset=0, name_prefix=""):
    """Create an OpenMC cell with an explicit, trackable ID."""
    return openmc.Cell(
        cell_id + id_offset,
        name=f"{name_prefix}{name}",
        region=region,
        fill=fill,
    )


class SCLLModuleBuilder:
    """
    Self-cooled PbLi half-annular test module.

    Radial build (plasma → outward):
      1. Armor ending at r_fw
      2. First-wall channel row
      3. Divider
      4. Optional bulk PbLi pack (per-row channel counts / thicknesses)

    Cell IDs (see module constants)::

        8001 armor · 8002 fw_structure · 8003 divider · 8004 bulk_structure
        81xx fw_channel_z{xx} · 82rrzz bulk_r{r}_z{zz}
    """

    def __init__(
        self,
        r_fw,
        n_fw_channels,
        axial_length,
        channel_radial_thickness,
        front_thickness,
        back_thickness,
        z0=0.0,
        orientation="right",
        axial_gap=0.0,
        axial_end_thickness=0.0,
        armor_thickness=0.2,
        divider_thickness=0.6,
        divider_material_name="tungsten",
        bulk_rows=(),
        bulk_front_thickness=0.5,
        bulk_back_thickness=1.0,
        bulk_radial_gap=1.0,
        bulk_axial_gap=None,
        bulk_axial_end_thickness=None,
        materials=None,
        material_ns=None,
        cell_id_offset=0,
        name_prefix="",
    ):
        self.r_fw = r_fw
        self.n_fw_channels = n_fw_channels
        self.axial_length = axial_length
        self.channel_radial_thickness = channel_radial_thickness
        self.front_thickness = front_thickness
        self.back_thickness = back_thickness
        self.z0 = z0
        self.orientation = orientation
        self.axial_gap = axial_gap
        self.axial_end_thickness = axial_end_thickness
        self.armor_thickness = armor_thickness
        self.divider_thickness = divider_thickness
        self.divider_material_name = divider_material_name
        self.bulk_rows = tuple(bulk_rows) if bulk_rows else ()
        self.bulk_front_thickness = bulk_front_thickness
        self.bulk_back_thickness = bulk_back_thickness
        self.bulk_radial_gap = bulk_radial_gap
        self.bulk_axial_gap = bulk_axial_gap
        self.bulk_axial_end_thickness = bulk_axial_end_thickness
        self.material_ns = material_ns if material_ns is not None else m
        self.materials = materials  # optional overrides dict / SimpleNamespace
        self.cell_id_offset = int(cell_id_offset)
        self.name_prefix = name_prefix or ""
        self.module_type = "SCLL"

        self._cells = {}
        self._channels = []
        self._universe = None
        self._extents = None
        self._fw_layout = {}
        self._bulk_layout = {}
        self._resolved_materials = None

        self._validate()

    def _validate(self):
        if self.r_fw <= 0:
            raise ValueError("r_fw must be positive")
        if self.armor_thickness <= 0:
            raise ValueError("armor_thickness must be positive")
        if self.divider_thickness <= 0:
            raise ValueError("divider_thickness must be positive")
        if not self.divider_material_name:
            raise ValueError("divider_material_name must be a non-empty name")
        if self.r_fw <= self.armor_thickness:
            raise ValueError("r_fw must exceed armor_thickness")
        if self.n_fw_channels < 1:
            raise ValueError("n_fw_channels must be >= 1")
        if self.axial_length <= 0:
            raise ValueError("axial_length must be positive")
        if min(self.channel_radial_thickness, self.front_thickness, self.back_thickness) <= 0:
            raise ValueError("channel and wall thicknesses must be positive")
        if self.axial_gap < 0 or self.axial_end_thickness < 0:
            raise ValueError("axial_gap and axial_end_thickness must be >= 0")
        if self.orientation not in ("left", "right"):
            raise ValueError("orientation must be 'left' or 'right'")
        if self.bulk_rows:
            if min(self.bulk_front_thickness, self.bulk_back_thickness) <= 0:
                raise ValueError("bulk front/back thicknesses must be positive")
            if self.bulk_radial_gap < 0:
                raise ValueError("bulk_radial_gap must be >= 0")
        if self.cell_id_offset < 0:
            raise ValueError("cell_id_offset must be >= 0")

    def _make_cell(self, cell_id, name, region, fill):
        return _cell(
            cell_id,
            name,
            region,
            fill,
            id_offset=self.cell_id_offset,
            name_prefix=self.name_prefix,
        )

    # --- derived radii / z (same stack as testing_module_model_builder) ---

    def _bulk_axial_gap_value(self):
        return self.axial_gap if self.bulk_axial_gap is None else self.bulk_axial_gap

    def _bulk_axial_end_value(self):
        if self.bulk_axial_end_thickness is None:
            return self.axial_end_thickness
        return self.bulk_axial_end_thickness

    def r_armor(self):
        return self.r_fw - self.armor_thickness

    def fw_radial_thickness(self):
        return self.front_thickness + self.channel_radial_thickness + self.back_thickness

    def fw_outer_radius(self):
        return self.r_fw + self.fw_radial_thickness()

    def r_divider(self):
        return self.fw_outer_radius()

    def divider_outer_radius(self):
        return self.r_divider() + self.divider_thickness

    def r_bulk(self):
        return self.divider_outer_radius()

    def bulk_radial_thickness(self):
        if not self.bulk_rows:
            return 0.0
        channel_sum = sum(row.radial_thickness for row in self.bulk_rows)
        gaps = (len(self.bulk_rows) - 1) * self.bulk_radial_gap
        return self.bulk_front_thickness + channel_sum + gaps + self.bulk_back_thickness

    def bulk_outer_radius(self):
        return self.r_bulk() + self.bulk_radial_thickness()

    def module_outer_radius(self):
        if self.bulk_rows:
            return self.bulk_outer_radius()
        return self.divider_outer_radius()

    def z_min(self):
        return self.z0 - self.axial_length / 2.0

    def z_max(self):
        return self.z0 + self.axial_length / 2.0

    def _resolve_materials(self):
        """Fill SCLL roles from overrides or material_ns defaults."""
        overrides = self.materials or {}
        if isinstance(overrides, dict):
            get = overrides.get
        else:
            get = lambda key, default=None: getattr(overrides, key, default)

        ns = self.material_ns
        return SimpleNamespace(
            armor=get("armor") or ns.tungsten,
            fw_structure=get("fw_structure") or ns.rafm_steel,
            fw_channel=get("fw_channel") or ns.LiPb_breeder,
            divider=get("divider") or ns.tungsten,
            bulk_structure=get("bulk_structure") or ns.rafm_steel,
            bulk_channel=get("bulk_channel") or ns.LiPb_breeder,
        )

    def build(self):
        """Construct regions, channel records, cells, universe, and extents."""
        mats = self._resolve_materials()
        self._resolved_materials = mats
        cells = {}
        channels = []

        # 1. Armor
        armor_region = annular_shell_region(
            z0=self.z0,
            inner_radius=self.r_armor(),
            radial_thickness=self.armor_thickness,
            axial_length=self.axial_length,
            orientation=self.orientation,
        )
        cells["armor"] = self._make_cell(
            SCLL_ID_ARMOR, "armor", armor_region, mats.armor
        )

        # 2. First-wall channels + structure
        fw_channels, fw_layout = annular_shell_channel_grid(
            z0=self.z0,
            module_inner_radius=self.r_fw,
            axial_length=self.axial_length,
            n_axial=self.n_fw_channels,
            n_radial=1,
            channel_radial_thickness=self.channel_radial_thickness,
            radial_gap=0.0,
            first_wall_thickness=self.front_thickness,
            outer_wall_thickness=self.back_thickness,
            axial_gap=self.axial_gap,
            axial_end_thickness=self.axial_end_thickness,
            orientation=self.orientation,
        )
        self._fw_layout = fw_layout

        fw_envelope = annular_shell_region(
            z0=self.z0,
            inner_radius=self.r_fw,
            radial_thickness=self.fw_radial_thickness(),
            axial_length=self.axial_length,
            orientation=self.orientation,
        )
        fw_structure = annular_shell_structure(fw_envelope, fw_channels)
        cells["fw_structure"] = self._make_cell(
            SCLL_ID_FW_STRUCTURE, "fw_structure", fw_structure, mats.fw_structure
        )

        r_ch_inner = float(fw_layout["row_inner_radii"][0])
        ch_axial = float(fw_layout["channel_axial_length"])
        for axial_index, (region, z_center) in enumerate(
            zip(fw_channels, fw_layout["z_centers"])
        ):
            zc = float(z_center)
            z_num = axial_index + 1
            name = f"fw_channel_z{z_num}"
            cell_id = SCLL_ID_FW_CHANNEL0 + z_num
            record = SimpleNamespace(
                name=f"{self.name_prefix}{name}",
                kind="fw",
                row=0,
                axial_index=axial_index,
                cell_id=cell_id + self.cell_id_offset,
                region=region,
                r_inner=r_ch_inner,
                r_outer=r_ch_inner + self.channel_radial_thickness,
                z_min=zc - ch_axial / 2.0,
                z_max=zc + ch_axial / 2.0,
                material=mats.fw_channel,
            )
            channels.append(record)
            cells[name] = self._make_cell(cell_id, name, region, mats.fw_channel)

        # 3. Divider
        divider_region = annular_shell_region(
            z0=self.z0,
            inner_radius=self.r_divider(),
            radial_thickness=self.divider_thickness,
            axial_length=self.axial_length,
            orientation=self.orientation,
        )
        cells["divider"] = self._make_cell(
            SCLL_ID_DIVIDER, "divider", divider_region, mats.divider
        )

        layer_names = ["armor_inner", "fw_inner", "fw_outer", "divider_outer"]
        r_interfaces = [
            self.r_armor(),
            self.r_fw,
            self.fw_outer_radius(),
            self.divider_outer_radius(),
        ]

        # 4. Bulk pack (optional)
        if self.bulk_rows:
            bulk_channels, bulk_layout = annular_shell_channel_rows(
                z0=self.z0,
                axial_length=self.axial_length,
                pack_inner_radius=self.r_bulk(),
                front_thickness=self.bulk_front_thickness,
                back_thickness=self.bulk_back_thickness,
                radial_gap=self.bulk_radial_gap,
                rows=self.bulk_rows,
                axial_gap=self._bulk_axial_gap_value(),
                axial_end_thickness=self._bulk_axial_end_value(),
                orientation=self.orientation,
            )
            self._bulk_layout = bulk_layout

            bulk_envelope = annular_shell_region(
                z0=self.z0,
                inner_radius=self.r_bulk(),
                radial_thickness=float(bulk_layout["pack_radial_thickness"]),
                axial_length=self.axial_length,
                orientation=self.orientation,
            )
            bulk_structure = annular_shell_structure(bulk_envelope, bulk_channels)
            cells["bulk_structure"] = self._make_cell(
                SCLL_ID_BULK_STRUCTURE,
                "bulk_structure",
                bulk_structure,
                mats.bulk_structure,
            )

            offset = 0
            for row_layout in bulk_layout["row_layouts"]:
                n_ch = int(row_layout["n_channels"])
                row_index = int(row_layout["row_index"])
                r_inner = float(row_layout["inner_radius"])
                r_outer = float(row_layout["outer_radius"])
                ch_len = float(row_layout["channel_axial_length"])
                z_centers = row_layout["z_centers"]
                for axial_index, (region, z_center) in enumerate(
                    zip(bulk_channels[offset : offset + n_ch], z_centers)
                ):
                    zc = float(z_center)
                    r_num = row_index + 1
                    z_num = axial_index + 1
                    name = f"bulk_r{r_num}_z{z_num}"
                    cell_id = SCLL_ID_BULK_ROW0 + 100 * r_num + z_num
                    record = SimpleNamespace(
                        name=f"{self.name_prefix}{name}",
                        kind="bulk",
                        row=row_index,
                        axial_index=axial_index,
                        cell_id=cell_id + self.cell_id_offset,
                        region=region,
                        r_inner=r_inner,
                        r_outer=r_outer,
                        z_min=zc - ch_len / 2.0,
                        z_max=zc + ch_len / 2.0,
                        material=mats.bulk_channel,
                    )
                    channels.append(record)
                    cells[name] = self._make_cell(cell_id, name, region, mats.bulk_channel)
                offset += n_ch

            layer_names = list(layer_names) + ["bulk_outer"]
            r_interfaces.append(self.bulk_outer_radius())

        self._cells = cells
        self._channels = channels
        self._universe = openmc.Universe(name="SCLLModule", cells=list(cells.values()))
        self._extents = SimpleNamespace(
            r_inner=self.r_armor(),
            r_outer=self.module_outer_radius(),
            z_min=self.z_min(),
            z_max=self.z_max(),
            orientation=self.orientation,
            r_interfaces=tuple(r_interfaces),
            layer_names=tuple(layer_names),
        )
        return self

    def _require_built(self):
        if self._universe is None:
            raise RuntimeError("Call SCLLModuleBuilder.build() before querying geometry")

    def get_channels(self):
        self._require_built()
        return list(self._channels)

    def get_fw_channels(self):
        return [c for c in self.get_channels() if c.kind == "fw"]

    def get_bulk_channels(self):
        return [c for c in self.get_channels() if c.kind == "bulk"]

    def get_cells(self):
        self._require_built()
        return dict(self._cells)

    def get_cell_id_map(self):
        """Return ``{cell_id: name}`` sorted by ID for debug / tallies."""
        self._require_built()
        return {cell.id: name for name, cell in sorted(
            self._cells.items(), key=lambda item: item[1].id
        )}

    def get_cells_by_role(self):
        """Group cells for tallies: structure roles vs channel kinds."""
        self._require_built()
        groups = {
            "armor": [],
            "fw_structure": [],
            "divider": [],
            "bulk_structure": [],
            "fw_channel": [],
            "bulk_channel": [],
        }
        for name, cell in self._cells.items():
            if name in groups:
                groups[name].append(cell)
            elif name.startswith("fw_channel"):
                groups["fw_channel"].append(cell)
            elif name.startswith("bulk_r"):
                groups["bulk_channel"].append(cell)
        return groups

    def get_universe(self):
        self._require_built()
        return self._universe

    def get_extents(self):
        self._require_built()
        return self._extents

    def get_materials(self):
        self._require_built()
        return self._resolved_materials

    def get_fw_layout(self):
        self._require_built()
        return dict(self._fw_layout)

    def get_bulk_layout(self):
        self._require_built()
        return dict(self._bulk_layout)

    def get_channel_extents(self, kind, axial_index, row=0):
        """
        Return cylindrical-mesh extents for one FW or bulk channel.

        kind: "fw" / "FW" or "bulk" / "BULK"
        axial_index: zero-based index along z within that row
        row: always 0 for FW; bulk radial row index otherwise

        Returns a SimpleNamespace with r/z/phi bounds plus helpers:
          cylindrical_grids(dimensions), to_cylindrical_mesh(dimensions)
        """
        self._require_built()
        kind_key = kind.strip().lower()
        if kind_key not in ("fw", "bulk"):
            raise ValueError("kind must be 'fw' or 'bulk'")
        if axial_index < 0:
            raise ValueError("axial_index must be >= 0")
        if row < 0:
            raise ValueError("row must be >= 0")
        if kind_key == "fw" and row != 0:
            raise ValueError("first-wall channels only have row=0")

        matches = [
            ch
            for ch in self._channels
            if ch.kind == kind_key and ch.row == row and ch.axial_index == axial_index
        ]
        if not matches:
            available = [
                (ch.kind, ch.row, ch.axial_index, ch.name)
                for ch in self._channels
                if ch.kind == kind_key
            ]
            raise KeyError(
                f"No {kind_key} channel with row={row}, axial_index={axial_index}. "
                f"Available: {available}"
            )
        ch = matches[0]

        # OpenMC requires phi in [0, 2pi]. Left half is one interval;
        # right half wraps across phi=0, so use full circle + CellFilter.
        if self.orientation == "right":
            phi_min, phi_max = 0.0, 2.0 * math.pi
        else:
            phi_min, phi_max = 0.5 * math.pi, 1.5 * math.pi

        def cylindrical_grids(dimensions=(10, 1, 10)):
            n_r, n_phi, n_z = dimensions
            if min(n_r, n_phi, n_z) < 1:
                raise ValueError("dimensions must be positive integers")
            origin = (0.0, 0.0, ch.z_min)
            r_grid = np.linspace(ch.r_inner, ch.r_outer, num=n_r + 1)
            phi_grid = np.linspace(phi_min, phi_max, num=n_phi + 1)
            z_grid = np.linspace(0.0, ch.z_max - ch.z_min, num=n_z + 1)
            return {
                "r_grid": r_grid,
                "phi_grid": phi_grid,
                "z_grid": z_grid,
                "origin": origin,
            }

        def to_cylindrical_mesh(dimensions=(10, 1, 10)):
            grids = cylindrical_grids(dimensions)
            return openmc.CylindricalMesh(
                r_grid=grids["r_grid"],
                phi_grid=grids["phi_grid"],
                z_grid=grids["z_grid"],
                origin=grids["origin"],
            )

        return SimpleNamespace(
            kind=ch.kind,
            row=ch.row,
            axial_index=ch.axial_index,
            name=ch.name,
            r_min=ch.r_inner,
            r_max=ch.r_outer,
            z_min=ch.z_min,
            z_max=ch.z_max,
            phi_min=phi_min,
            phi_max=phi_max,
            orientation=self.orientation,
            needs_cell_filter=(self.orientation == "right"),
            cylindrical_grids=cylindrical_grids,
            to_cylindrical_mesh=to_cylindrical_mesh,
        )

    def named_regions(self):
        self._require_built()
        return {
            name: cell.region
            for name, cell in self._cells.items()
            if cell.region is not None
        }


def build_scll_module(**kwargs):
    """Construct and build an SCLL testing module."""
    return SCLLModuleBuilder(**kwargs).build()


class DCLLFNSFModuleBuilder:
    """
    Half-annular, FNSF-inspired DCLL test-section CSG (cm).

    RAFM steel with axial He channels in FW, midplate, and back wall.
    Two PbLi ducts with RAFM walls, narrow PbLi gap, and SiC FCI.
    Front/rear PbLi cores are split into axial cells for tallies only
    (not physical baffles). Headers / turns / side walls still TBD.

    Cell IDs (see module constants)::

        90xx solid radial shells (build order) · 9081/82/83 He steel webs
        91xx fw He · 92xx midplate He · 93xx back He
        94xx front PbLi bins · 95xx rear PbLi bins
    """

    def __init__(
        self,
        r_inner,
        axial_length,
        front_pbli_thickness,
        rear_pbli_thickness,
        n_axial_bins=4,
        fw_steel_front=0.4,
        fw_he_thickness=2.0,
        fw_steel_back=0.4,
        duct_steel_thickness=0.5,
        pbli_gap_thickness=0.1,
        fci_thickness=0.5,
        midplate_he_thickness=0.8,
        back_he_thickness=0.8,
        n_fw_he_channels=12,
        n_midplate_he_channels=8,
        n_back_he_channels=8,
        he_rib_angle_deg=1.0,
        he_side_wall_angle_deg=2.0,
        z0=0.0,
        orientation="right",
        materials=None,
        material_ns=None,
        cell_id_offset=0,
        name_prefix="",
    ):
        self.r_inner = r_inner
        self.axial_length = axial_length
        self.front_pbli_thickness = front_pbli_thickness
        self.rear_pbli_thickness = rear_pbli_thickness
        self.n_axial_bins = n_axial_bins
        self.fw_steel_front = fw_steel_front
        self.fw_he_thickness = fw_he_thickness
        self.fw_steel_back = fw_steel_back
        self.duct_steel_thickness = duct_steel_thickness
        self.pbli_gap_thickness = pbli_gap_thickness
        self.fci_thickness = fci_thickness
        self.midplate_he_thickness = midplate_he_thickness
        self.back_he_thickness = back_he_thickness
        self.n_fw_he_channels = n_fw_he_channels
        self.n_midplate_he_channels = n_midplate_he_channels
        self.n_back_he_channels = n_back_he_channels
        self.he_rib_angle_deg = he_rib_angle_deg
        self.he_side_wall_angle_deg = he_side_wall_angle_deg
        self.z0 = z0
        self.orientation = orientation
        self.material_ns = material_ns if material_ns is not None else m
        self.materials = materials
        self.cell_id_offset = int(cell_id_offset)
        self.name_prefix = name_prefix or ""
        self.module_type = "DCLL_FNSF"

        self._cells = {}
        self._channels = []
        self._extents = None
        self._universe = None
        self._resolved_materials = None

        self._validate()
        if self.cell_id_offset < 0:
            raise ValueError("cell_id_offset must be >= 0")

    def _make_cell(self, cell_id, name, region, fill):
        return _cell(
            cell_id,
            name,
            region,
            fill,
            id_offset=self.cell_id_offset,
            name_prefix=self.name_prefix,
        )

    def module_outer_radius(self):
        self._require_built()
        return float(self._extents.r_outer)

    def _validate(self):
        lengths = (
            self.r_inner,
            self.axial_length,
            self.front_pbli_thickness,
            self.rear_pbli_thickness,
            self.fw_steel_front,
            self.fw_he_thickness,
            self.fw_steel_back,
            self.duct_steel_thickness,
            self.pbli_gap_thickness,
            self.fci_thickness,
            self.midplate_he_thickness,
            self.back_he_thickness,
        )
        if any(not math.isfinite(value) or value <= 0 for value in lengths):
            raise ValueError(
                "all radii, lengths, and layer thicknesses must be finite and positive"
            )
        if (
            isinstance(self.n_axial_bins, bool)
            or not isinstance(self.n_axial_bins, int)
            or self.n_axial_bins < 1
        ):
            raise ValueError("n_axial_bins must be a positive integer")
        for count in (
            self.n_fw_he_channels,
            self.n_midplate_he_channels,
            self.n_back_he_channels,
        ):
            if isinstance(count, bool) or not isinstance(count, int) or count < 1:
                raise ValueError("helium channel counts must be positive integers")
            if (
                2 * self.he_side_wall_angle_deg
                + (count - 1) * self.he_rib_angle_deg
                >= 180
            ):
                raise ValueError("helium ribs and side walls leave no open duct width")
        if self.orientation not in ("left", "right"):
            raise ValueError("orientation must be 'left' or 'right'")

    def _resolve_materials(self):
        overrides = self.materials or {}
        if isinstance(overrides, dict):
            get = overrides.get
        else:
            get = lambda key, default=None: getattr(overrides, key, default)

        ns = self.material_ns
        return SimpleNamespace(
            structure=get("structure") or ns.rafm_steel,
            helium=get("helium") or ns.helium_8mpa,
            pbli=get("pbli") or ns.LiPb_breeder,
            fci=get("fci") or ns.silicon_carbide_fci,
        )

    def build(self):
        """Build mutually exclusive material cells for the module envelope."""
        mats = self._resolve_materials()
        self._resolved_materials = mats
        cells = {}
        channels = []
        radius = self.r_inner
        interfaces = [radius]
        layer_names = ["inner"]
        z_min = self.z0 - self.axial_length / 2
        layer_id = DCLL_ID_LAYER0
        he_id0 = {
            "fw": DCLL_ID_FW_HE0,
            "midplate": DCLL_ID_MID_HE0,
            "back": DCLL_ID_BACK_HE0,
        }
        webs_id = {
            "fw": DCLL_ID_FW_WEBS,
            "midplate": DCLL_ID_MIDPLATE_WEBS,
            "back": DCLL_ID_BACK_WEBS,
        }
        pbli_id0 = {
            "front_pbli": DCLL_ID_FRONT_PBLI0,
            "rear_pbli": DCLL_ID_REAR_PBLI0,
        }

        def add_layer(name, thickness, fill):
            nonlocal radius, layer_id
            r0 = radius
            layer_id += 1
            region = annular_shell_region(
                z0=self.z0,
                inner_radius=r0,
                radial_thickness=thickness,
                axial_length=self.axial_length,
                orientation=self.orientation,
            )
            cells[name] = self._make_cell(layer_id, name, region, fill)
            radius += thickness
            interfaces.append(radius)
            layer_names.append(name)

        def add_helium_ducts(label, thickness, count):
            nonlocal radius
            r0 = radius
            envelope = annular_shell_region(
                z0=self.z0,
                inner_radius=r0,
                radial_thickness=thickness,
                axial_length=self.axial_length,
                orientation=self.orientation,
            )
            ducts, _ = half_annular_axial_channels(
                z0=self.z0,
                inner_radius=r0,
                radial_thickness=thickness,
                axial_length=self.axial_length,
                n_channels=count,
                rib_angle_deg=self.he_rib_angle_deg,
                side_wall_angle_deg=self.he_side_wall_angle_deg,
                orientation=self.orientation,
            )
            steel_name = f"{label}_steel_webs"
            cells[steel_name] = self._make_cell(
                webs_id[label],
                steel_name,
                annular_shell_structure(envelope, ducts),
                mats.structure,
            )
            for index, region in enumerate(ducts):
                az = index + 1
                name = f"{label}_helium_az{az}"
                cell_id = he_id0[label] + az
                cells[name] = self._make_cell(cell_id, name, region, mats.helium)
                channels.append(
                    SimpleNamespace(
                        name=f"{self.name_prefix}{name}",
                        kind=f"helium_{label}",
                        row=0,
                        axial_index=index,
                        cell_id=cell_id + self.cell_id_offset,
                        region=region,
                        r_inner=r0,
                        r_outer=r0 + thickness,
                        z_min=z_min,
                        z_max=z_min + self.axial_length,
                        material=mats.helium,
                    )
                )
            radius += thickness
            interfaces.append(radius)
            layer_names.append(label)

        def add_pbli_core(name, thickness, row):
            nonlocal radius
            r0 = radius
            dz = self.axial_length / self.n_axial_bins
            for index in range(self.n_axial_bins):
                center = z_min + (index + 0.5) * dz
                z_num = index + 1
                cell_name = f"{name}_z{z_num}"
                cell_id = pbli_id0[name] + z_num
                region = annular_shell_region(
                    z0=center,
                    inner_radius=r0,
                    radial_thickness=thickness,
                    axial_length=dz,
                    orientation=self.orientation,
                )
                cells[cell_name] = self._make_cell(cell_id, cell_name, region, mats.pbli)
                channels.append(
                    SimpleNamespace(
                        name=f"{self.name_prefix}{cell_name}",
                        kind=name,
                        row=row,
                        axial_index=index,
                        cell_id=cell_id + self.cell_id_offset,
                        region=region,
                        r_inner=r0,
                        r_outer=r0 + thickness,
                        z_min=center - dz / 2,
                        z_max=center + dz / 2,
                        material=mats.pbli,
                    )
                )
            radius += thickness
            interfaces.append(radius)
            layer_names.append(name)

        def add_duct(label, thickness, row):
            add_layer(f"{label}_steel_inner", self.duct_steel_thickness, mats.structure)
            add_layer(f"{label}_gap_inner", self.pbli_gap_thickness, mats.pbli)
            add_layer(f"{label}_fci_inner", self.fci_thickness, mats.fci)
            add_pbli_core(f"{label}_pbli", thickness, row)
            add_layer(f"{label}_fci_outer", self.fci_thickness, mats.fci)
            add_layer(f"{label}_gap_outer", self.pbli_gap_thickness, mats.pbli)
            add_layer(f"{label}_steel_outer", self.duct_steel_thickness, mats.structure)

        # Radial stack (plasma → out)
        add_layer("fw_steel_front", self.fw_steel_front, mats.structure)
        add_helium_ducts("fw", self.fw_he_thickness, self.n_fw_he_channels)
        add_layer("fw_steel_back", self.fw_steel_back, mats.structure)
        add_duct("front", self.front_pbli_thickness, 0)
        add_layer("midplate_steel_front", self.duct_steel_thickness, mats.structure)
        add_helium_ducts("midplate", self.midplate_he_thickness, self.n_midplate_he_channels)
        add_layer("midplate_steel_back", self.duct_steel_thickness, mats.structure)
        add_duct("rear", self.rear_pbli_thickness, 1)
        add_layer("back_steel_front", self.duct_steel_thickness, mats.structure)
        add_helium_ducts("back", self.back_he_thickness, self.n_back_he_channels)
        add_layer("back_steel_outer", self.duct_steel_thickness, mats.structure)

        self._cells = cells
        self._channels = channels
        self._extents = SimpleNamespace(
            r_inner=self.r_inner,
            r_outer=radius,
            z_min=z_min,
            z_max=z_min + self.axial_length,
            orientation=self.orientation,
            r_interfaces=tuple(interfaces),
            layer_names=tuple(layer_names),
        )
        self._universe = openmc.Universe(
            name="DCLL_FNSF_Module", cells=list(cells.values())
        )
        return self

    def _require_built(self):
        if self._universe is None:
            raise RuntimeError(
                "Call DCLLFNSFModuleBuilder.build() before querying geometry"
            )

    def get_channels(self):
        self._require_built()
        return list(self._channels)

    def get_cells(self):
        self._require_built()
        return dict(self._cells)

    def get_cell_id_map(self):
        """Return ``{cell_id: name}`` sorted by ID for debug / tallies."""
        self._require_built()
        return {cell.id: name for name, cell in sorted(
            self._cells.items(), key=lambda item: item[1].id
        )}

    def get_cells_by_role(self):
        """Group steel, helium, PbLi, and SiC cells for material tallies."""
        self._require_built()
        groups = {
            "structure": [],
            "helium": [],
            "pbli_core": [],
            "pbli_gap": [],
            "fci": [],
        }
        for name, cell in self._cells.items():
            if "steel" in name:
                role = "structure"
            elif "helium" in name:
                role = "helium"
            elif "_pbli_z" in name:
                role = "pbli_core"
            elif "_gap_" in name:
                role = "pbli_gap"
            else:
                role = "fci"
            groups[role].append(cell)
        return groups

    def get_universe(self):
        self._require_built()
        return self._universe

    def get_extents(self):
        self._require_built()
        return self._extents

    def get_materials(self):
        self._require_built()
        return self._resolved_materials

    def named_regions(self):
        self._require_built()
        return {
            name: cell.region
            for name, cell in self._cells.items()
            if cell.region is not None
        }


def build_dcll_fnsf_module(**kwargs):
    """Construct and build an FNSF-inspired dual-coolant test module."""
    return DCLLFNSFModuleBuilder(**kwargs).build()


if __name__ == "__main__":
    module = build_scll_module(
        r_fw=50.0,
        n_fw_channels=20,
        axial_length=100.0,  # 1 m
        channel_radial_thickness=2.0,
        front_thickness=0.4,
        back_thickness=0.4,
        orientation="right",
        axial_gap=0.4,
        axial_end_thickness=1.0,
        bulk_rows=(BulkChannelRow(n_channels=4, radial_thickness=80.0),),
        bulk_front_thickness=0.5,
        bulk_back_thickness=1.0,
        bulk_radial_gap=1.0,
    )
    extents = module.get_extents()
    fw0 = module.get_channel_extents("FW", axial_index=0)
    bulk0 = module.get_channel_extents("BULK", row=0, axial_index=0)
    print(f"universe cells: {len(module.get_universe().cells)}")
    print(
        f"module extents: r=[{extents.r_inner:.3f}, {extents.r_outer:.3f}] cm, "
        f"z=[{extents.z_min:.3f}, {extents.z_max:.3f}] cm"
    )
    print(
        f"FW[0]:  {fw0.name}  r=[{fw0.r_min:.3f}, {fw0.r_max:.3f}]  "
        f"z=[{fw0.z_min:.3f}, {fw0.z_max:.3f}]"
    )
    print(
        f"BULK[0,0]: {bulk0.name}  r=[{bulk0.r_min:.3f}, {bulk0.r_max:.3f}]  "
        f"z=[{bulk0.z_min:.3f}, {bulk0.z_max:.3f}]"
    )
    mesh = fw0.to_cylindrical_mesh((4, 1, 6))
    print(f"sample FW mesh dimension={mesh.dimension}, origin={mesh.origin}")
