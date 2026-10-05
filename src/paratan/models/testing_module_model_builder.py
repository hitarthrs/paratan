"""Testing-module model objects for the tandem VNS central cell.

Definitions live here. Geometry primitives stay in ``paratan.geometry.core``.
Lengths are centimetres.

Run from the paratan repository root::

    MPLCONFIGDIR=/tmp/mpl-vns .openmc-env/bin/python src/paratan/models/testing_module_model_builder.py
"""
from __future__ import annotations

from dataclasses import dataclass, field
import math
from pathlib import Path
import sys
from typing import Any

# Allow ``python path/to/testing_module_model_builder.py`` without an install.
_ROOT = Path(__file__).resolve().parents[3]  # .../paratan
if str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))

import numpy as np
import openmc

import src.paratan.geometry.core as core
import src.paratan.materials.material as default_materials

# Geometry-only row spec; re-exported under the SCLL name used in demos.
BulkChannelRow = core.ChannelRowSpec


@dataclass(frozen=True)
class ChannelRecord:
    """One flow channel owned by a built module."""

    name: str
    kind: str  # SCLL: "fw"/"bulk"; DCLL: helium or PbLi path name
    row: int
    axial_index: int
    region: openmc.Region
    r_inner: float
    r_outer: float
    z_min: float
    z_max: float
    material: openmc.Material | None = None


@dataclass(frozen=True)
class ModuleExtents:
    """Geometric extents for mesh / cell tally setup."""

    r_inner: float
    r_outer: float
    z_min: float
    z_max: float
    orientation: str
    r_interfaces: tuple[float, ...]
    layer_names: tuple[str, ...]


@dataclass(frozen=True)
class ChannelMeshExtents:
    """Bounds for one channel, ready for ``openmc.CylindricalMesh`` setup.

    ``r_min``/``r_max``/``z_min``/``z_max`` come from the channel packing.
    ``phi_min``/``phi_max`` are OpenMC-valid bounds in ``[0, 2π]``:

    - ``left`` (x < 0): ``[π/2, 3π/2]`` — matches the half-annulus exactly
    - ``right`` (x > 0): ``[0, 2π]`` — full azimuth (cannot span the +x half
      as one interval without wrapping). Pair the mesh filter with a
      ``CellFilter`` on the channel cell to restrict scoring to the half.

    Grids from :meth:`cylindrical_grids` use origin ``(0, 0, z_min)`` with
    mesh-local ``z_grid``, matching
    :func:`paratan.geometry.core.hollow_mesh_from_domain`.
    """

    kind: str
    row: int
    axial_index: int
    name: str
    r_min: float
    r_max: float
    z_min: float
    z_max: float
    phi_min: float
    phi_max: float
    orientation: str

    @property
    def origin(self) -> tuple[float, float, float]:
        return (0.0, 0.0, self.z_min)

    @property
    def r_grid_edges(self) -> tuple[float, float]:
        return (self.r_min, self.r_max)

    @property
    def z_grid_edges(self) -> tuple[float, float]:
        """Absolute lab-frame z edges (not mesh-relative)."""
        return (self.z_min, self.z_max)

    @property
    def phi_grid_bounds(self) -> tuple[float, float]:
        return (self.phi_min, self.phi_max)

    @property
    def needs_cell_filter(self) -> bool:
        """True when phi spans a full circle and the half-annulus must be clipped."""
        return self.orientation == "right"

    def cylindrical_grids(
        self,
        dimensions: tuple[int, int, int] = (10, 1, 10),
    ) -> dict[str, object]:
        """Build ``r_grid``, ``phi_grid``, ``z_grid``, and ``origin`` for a mesh.

        Parameters
        ----------
        dimensions
            ``(n_r, n_phi, n_z)`` bin counts. Defaults to one phi bin spanning
            ``phi_min``–``phi_max``.

        Returns
        -------
        dict
            Keys ``r_grid``, ``phi_grid``, ``z_grid``, ``origin`` suitable for
            ``openmc.CylindricalMesh(...)``.
        """
        n_r, n_phi, n_z = dimensions
        if min(n_r, n_phi, n_z) < 1:
            raise ValueError("dimensions must be positive integers")
        origin = self.origin
        r_grid = np.linspace(self.r_min, self.r_max, num=n_r + 1)
        phi_grid = np.linspace(self.phi_min, self.phi_max, num=n_phi + 1)
        # Mesh-local z: 0 at channel bottom → channel axial length at top.
        z_grid = np.linspace(0.0, self.z_max - self.z_min, num=n_z + 1)
        return {
            "r_grid": r_grid,
            "phi_grid": phi_grid,
            "z_grid": z_grid,
            "origin": origin,
        }

    def to_cylindrical_mesh(
        self,
        dimensions: tuple[int, int, int] = (10, 1, 10),
    ) -> openmc.CylindricalMesh:
        """Construct an ``openmc.CylindricalMesh`` covering this channel's r–z box."""
        grids = self.cylindrical_grids(dimensions)
        return openmc.CylindricalMesh(
            r_grid=grids["r_grid"],
            phi_grid=grids["phi_grid"],
            z_grid=grids["z_grid"],
            origin=grids["origin"],
        )


@dataclass
class SCLLModuleMaterials:
    """Material fills for SCLL roles. Defaults resolve from the materials module."""

    armor: openmc.Material | None = None
    fw_structure: openmc.Material | None = None
    fw_channel: openmc.Material | None = None
    divider: openmc.Material | None = None
    bulk_structure: openmc.Material | None = None
    bulk_channel: openmc.Material | None = None

    def resolve(self, material_ns: Any = default_materials) -> SCLLModuleMaterials:
        """Fill unset roles from ``material_ns`` defaults."""
        return SCLLModuleMaterials(
            armor=self.armor or material_ns.tungsten,
            fw_structure=self.fw_structure or material_ns.rafm_steel,
            fw_channel=self.fw_channel or material_ns.LiPb_breeder,
            divider=self.divider or material_ns.tungsten,
            bulk_structure=self.bulk_structure or material_ns.rafm_steel,
            bulk_channel=self.bulk_channel or material_ns.LiPb_breeder,
        )


@dataclass
class SCLLModule:
    """Self-cooled PbLi half-annular test module.

    Radial build order (plasma → outward):

    1. Solid armour ending at ``r_fw``.
    2. First-wall channel row.
    3. Solid divider.
    4. Bulk PbLi pack (per-row channel counts / thicknesses).

    Call :meth:`build` before querying channels, cells, universe, or extents.
    """

    r_fw: float
    n_fw_channels: int
    axial_length: float
    channel_radial_thickness: float
    front_thickness: float
    back_thickness: float
    z0: float = 0.0
    orientation: str = "right"
    axial_gap: float = 0.0
    axial_end_thickness: float = 0.0
    armor_thickness: float = 0.2
    divider_thickness: float = 0.6
    divider_material_name: str = "tungsten"
    bulk_rows: tuple[BulkChannelRow, ...] = ()
    bulk_front_thickness: float = 0.5
    bulk_back_thickness: float = 1.0
    bulk_radial_gap: float = 1.0
    bulk_axial_gap: float | None = None
    bulk_axial_end_thickness: float | None = None
    materials: SCLLModuleMaterials = field(default_factory=SCLLModuleMaterials)
    material_ns: Any = field(default=default_materials, repr=False)

    # Populated by build()
    _built: bool = field(default=False, init=False, repr=False)
    _resolved_materials: SCLLModuleMaterials | None = field(default=None, init=False, repr=False)
    _channels: list[ChannelRecord] = field(default_factory=list, init=False, repr=False)
    _cells: dict[str, openmc.Cell] = field(default_factory=dict, init=False, repr=False)
    _universe: openmc.Universe | None = field(default=None, init=False, repr=False)
    _extents: ModuleExtents | None = field(default=None, init=False, repr=False)
    _fw_layout: dict[str, object] = field(default_factory=dict, init=False, repr=False)
    _bulk_layout: dict[str, object] = field(default_factory=dict, init=False, repr=False)

    def __post_init__(self) -> None:
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

    @property
    def _bulk_axial_gap(self) -> float:
        return self.axial_gap if self.bulk_axial_gap is None else self.bulk_axial_gap

    @property
    def _bulk_axial_end(self) -> float:
        if self.bulk_axial_end_thickness is None:
            return self.axial_end_thickness
        return self.bulk_axial_end_thickness

    @property
    def r_armor(self) -> float:
        return self.r_fw - self.armor_thickness

    @property
    def fw_radial_thickness(self) -> float:
        return (
            self.front_thickness
            + self.channel_radial_thickness
            + self.back_thickness
        )

    @property
    def fw_outer_radius(self) -> float:
        return self.r_fw + self.fw_radial_thickness

    @property
    def r_divider(self) -> float:
        return self.fw_outer_radius

    @property
    def divider_outer_radius(self) -> float:
        return self.r_divider + self.divider_thickness

    @property
    def r_bulk(self) -> float:
        return self.divider_outer_radius

    @property
    def bulk_radial_thickness(self) -> float:
        if not self.bulk_rows:
            return 0.0
        channel_sum = sum(row.radial_thickness for row in self.bulk_rows)
        gaps = (len(self.bulk_rows) - 1) * self.bulk_radial_gap
        return self.bulk_front_thickness + channel_sum + gaps + self.bulk_back_thickness

    @property
    def bulk_outer_radius(self) -> float:
        return self.r_bulk + self.bulk_radial_thickness

    @property
    def module_outer_radius(self) -> float:
        if self.bulk_rows:
            return self.bulk_outer_radius
        return self.divider_outer_radius

    @property
    def z_min(self) -> float:
        return self.z0 - self.axial_length / 2.0

    @property
    def z_max(self) -> float:
        return self.z0 + self.axial_length / 2.0

    def _require_built(self) -> None:
        if not self._built:
            raise RuntimeError("Call SCLLModule.build() before querying built geometry")

    def build(self) -> SCLLModule:
        """Construct regions, channel records, cells, universe, and extents."""
        mats = self.materials.resolve(self.material_ns)
        self._resolved_materials = mats
        channels: list[ChannelRecord] = []
        cells: dict[str, openmc.Cell] = {}

        armor_region = core.annular_shell_region(
            z0=self.z0,
            inner_radius=self.r_armor,
            radial_thickness=self.armor_thickness,
            axial_length=self.axial_length,
            orientation=self.orientation,
        )
        cells["armor"] = openmc.Cell(name="armor", region=armor_region, fill=mats.armor)

        fw_channels, fw_layout = core.annular_shell_channel_grid(
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
        fw_envelope = core.annular_shell_region(
            z0=self.z0,
            inner_radius=self.r_fw,
            radial_thickness=self.fw_radial_thickness,
            axial_length=self.axial_length,
            orientation=self.orientation,
        )
        fw_structure = core.annular_shell_structure(fw_envelope, fw_channels)
        cells["fw_structure"] = openmc.Cell(
            name="fw_structure", region=fw_structure, fill=mats.fw_structure
        )

        r_ch_inner = float(fw_layout["row_inner_radii"][0])  # type: ignore[index]
        ch_axial = float(fw_layout["channel_axial_length"])  # type: ignore[arg-type]
        for axial_index, (region, z_center) in enumerate(
            zip(fw_channels, fw_layout["z_centers"])  # type: ignore[arg-type]
        ):
            zc = float(z_center)
            name = f"fw_channel_z{axial_index + 1}"
            record = ChannelRecord(
                name=name,
                kind="fw",
                row=0,
                axial_index=axial_index,
                region=region,
                r_inner=r_ch_inner,
                r_outer=r_ch_inner + self.channel_radial_thickness,
                z_min=zc - ch_axial / 2.0,
                z_max=zc + ch_axial / 2.0,
                material=mats.fw_channel,
            )
            channels.append(record)
            cells[name] = openmc.Cell(name=name, region=region, fill=mats.fw_channel)

        divider_region = core.annular_shell_region(
            z0=self.z0,
            inner_radius=self.r_divider,
            radial_thickness=self.divider_thickness,
            axial_length=self.axial_length,
            orientation=self.orientation,
        )
        cells["divider"] = openmc.Cell(name="divider", region=divider_region, fill=mats.divider)

        layer_names = ["armor_inner", "fw_inner", "fw_outer", "divider_outer"]
        r_interfaces = [
            self.r_armor,
            self.r_fw,
            self.fw_outer_radius,
            self.divider_outer_radius,
        ]

        if self.bulk_rows:
            bulk_channels, bulk_layout = core.annular_shell_channel_rows(
                z0=self.z0,
                axial_length=self.axial_length,
                pack_inner_radius=self.r_bulk,
                front_thickness=self.bulk_front_thickness,
                back_thickness=self.bulk_back_thickness,
                radial_gap=self.bulk_radial_gap,
                rows=self.bulk_rows,
                axial_gap=self._bulk_axial_gap,
                axial_end_thickness=self._bulk_axial_end,
                orientation=self.orientation,
            )
            self._bulk_layout = bulk_layout
            bulk_envelope = core.annular_shell_region(
                z0=self.z0,
                inner_radius=self.r_bulk,
                radial_thickness=float(bulk_layout["pack_radial_thickness"]),
                axial_length=self.axial_length,
                orientation=self.orientation,
            )
            bulk_structure = core.annular_shell_structure(bulk_envelope, bulk_channels)
            cells["bulk_structure"] = openmc.Cell(
                name="bulk_structure", region=bulk_structure, fill=mats.bulk_structure
            )

            offset = 0
            for row_layout in bulk_layout["row_layouts"]:  # type: ignore[assignment]
                assert isinstance(row_layout, dict)
                n_ch = int(row_layout["n_channels"])
                row_index = int(row_layout["row_index"])
                r_inner = float(row_layout["inner_radius"])
                r_outer = float(row_layout["outer_radius"])
                ch_len = float(row_layout["channel_axial_length"])
                z_centers = row_layout["z_centers"]
                assert isinstance(z_centers, tuple)
                for axial_index, (region, z_center) in enumerate(
                    zip(bulk_channels[offset : offset + n_ch], z_centers)
                ):
                    zc = float(z_center)
                    name = f"bulk_r{row_index + 1}_z{axial_index + 1}"
                    record = ChannelRecord(
                        name=name,
                        kind="bulk",
                        row=row_index,
                        axial_index=axial_index,
                        region=region,
                        r_inner=r_inner,
                        r_outer=r_outer,
                        z_min=zc - ch_len / 2.0,
                        z_max=zc + ch_len / 2.0,
                        material=mats.bulk_channel,
                    )
                    channels.append(record)
                    cells[name] = openmc.Cell(
                        name=name, region=region, fill=mats.bulk_channel
                    )
                offset += n_ch

            layer_names = (*layer_names, "bulk_outer")
            r_interfaces.append(self.bulk_outer_radius)

        self._channels = channels
        self._cells = cells
        self._universe = openmc.Universe(name="SCLLModule", cells=list(cells.values()))
        self._extents = ModuleExtents(
            r_inner=self.r_armor,
            r_outer=self.module_outer_radius,
            z_min=self.z_min,
            z_max=self.z_max,
            orientation=self.orientation,
            r_interfaces=tuple(r_interfaces),
            layer_names=tuple(layer_names),
        )
        self._built = True
        return self

    def get_channels(self) -> list[ChannelRecord]:
        self._require_built()
        return list(self._channels)

    def get_fw_channels(self) -> list[ChannelRecord]:
        return [c for c in self.get_channels() if c.kind == "fw"]

    def get_bulk_channels(self) -> list[ChannelRecord]:
        return [c for c in self.get_channels() if c.kind == "bulk"]

    def get_cells(self) -> dict[str, openmc.Cell]:
        self._require_built()
        return dict(self._cells)

    def get_cells_by_role(self) -> dict[str, list[openmc.Cell]]:
        """Group cells for tallies: structure roles vs channel kinds."""
        self._require_built()
        groups: dict[str, list[openmc.Cell]] = {
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

    def get_universe(self) -> openmc.Universe:
        self._require_built()
        assert self._universe is not None
        return self._universe

    def get_extents(self) -> ModuleExtents:
        self._require_built()
        assert self._extents is not None
        return self._extents

    def get_materials(self) -> SCLLModuleMaterials:
        self._require_built()
        assert self._resolved_materials is not None
        return self._resolved_materials

    def get_fw_layout(self) -> dict[str, object]:
        self._require_built()
        return dict(self._fw_layout)

    def get_bulk_layout(self) -> dict[str, object]:
        self._require_built()
        return dict(self._bulk_layout)

    def get_channel_extents(
        self,
        kind: str,
        *,
        axial_index: int,
        row: int = 0,
    ) -> ChannelMeshExtents:
        """Return cylindrical-mesh extents for one FW or bulk channel.

        Parameters
        ----------
        kind
            ``\"fw\"`` / ``\"FW\"`` or ``\"bulk\"`` / ``\"BULK\"``.
        axial_index
            Zero-based index along z within that row (increasing z).
        row
            Radial row index. Always ``0`` for the first wall; for bulk this
            selects which radial channel row.

        Returns
        -------
        ChannelMeshExtents
            ``r_min``/``r_max``/``z_min``/``z_max`` from the channel packing,
            plus half-annulus ``phi`` bounds from ``orientation``. Use
            :meth:`ChannelMeshExtents.cylindrical_grids` or
            :meth:`ChannelMeshExtents.to_cylindrical_mesh` for tally setup.
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

        # OpenMC requires phi ∈ [0, 2π]. Left half is one contiguous interval;
        # right half wraps across φ=0, so use a full circle and CellFilter.
        if self.orientation == "right":
            phi_min, phi_max = 0.0, 2.0 * math.pi
        else:
            phi_min, phi_max = 0.5 * math.pi, 1.5 * math.pi

        return ChannelMeshExtents(
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
        )

    def named_regions(self) -> dict[str, openmc.Region]:
        """Region map for plotting (void fills). Requires :meth:`build`."""
        self._require_built()
        return {name: cell.region for name, cell in self._cells.items() if cell.region is not None}


def build_scll_module(**kwargs: Any) -> SCLLModule:
    """Construct and build an :class:`SCLLModule` from keyword arguments."""
    return SCLLModule(**kwargs).build()


@dataclass
class DCLLFNSFMaterials:
    """Fills for the two-coolant FNSF-inspired module."""

    structure: openmc.Material | None = None
    helium: openmc.Material | None = None
    pbli: openmc.Material | None = None
    fci: openmc.Material | None = None

    def resolve(self, material_ns: Any = default_materials) -> DCLLFNSFMaterials:
        return DCLLFNSFMaterials(
            structure=self.structure or material_ns.rafm_steel,
            helium=self.helium or material_ns.helium_8mpa,
            pbli=self.pbli or material_ns.LiPb_breeder,
            fci=self.fci or material_ns.silicon_carbide_fci,
        )


@dataclass
class DCLLFNSFModule:
    """Half-annular, FNSF-inspired DCLL test-section CSG, in centimetres.

    RAFM steel contains axial He channels in the first wall, midplate, and
    back wall. Two PbLi ducts have RAFM walls, a narrow PbLi gap,
    and a non-structural SiC flow-channel insert on each radial side. The
    front/rear PbLi cores are split into touching axial cells for tallies:
    these cuts are not physical baffles. Headers, the front-to-rear turn,
    azimuthal side walls for the PbLi ducts remain to be designed. The He
    channels are annular-sector approximations of discrete ducts.
    """

    r_inner: float
    axial_length: float
    front_pbli_thickness: float
    rear_pbli_thickness: float
    n_axial_bins: int = 4
    fw_steel_front: float = 0.4
    fw_he_thickness: float = 2.0
    fw_steel_back: float = 0.4
    duct_steel_thickness: float = 0.5
    pbli_gap_thickness: float = 0.1
    fci_thickness: float = 0.5
    midplate_he_thickness: float = 0.8
    back_he_thickness: float = 0.8
    n_fw_he_channels: int = 12
    n_midplate_he_channels: int = 8
    n_back_he_channels: int = 8
    he_rib_angle_deg: float = 1.0
    he_side_wall_angle_deg: float = 2.0
    z0: float = 0.0
    orientation: str = "right"
    materials: DCLLFNSFMaterials = field(default_factory=DCLLFNSFMaterials)
    material_ns: Any = field(default=default_materials, repr=False)

    _built: bool = field(default=False, init=False, repr=False)
    _cells: dict[str, openmc.Cell] = field(default_factory=dict, init=False, repr=False)
    _channels: list[ChannelRecord] = field(default_factory=list, init=False, repr=False)
    _extents: ModuleExtents | None = field(default=None, init=False, repr=False)
    _universe: openmc.Universe | None = field(default=None, init=False, repr=False)
    _resolved_materials: DCLLFNSFMaterials | None = field(default=None, init=False, repr=False)

    def __post_init__(self) -> None:
        lengths = (
            self.r_inner, self.axial_length, self.front_pbli_thickness,
            self.rear_pbli_thickness, self.fw_steel_front,
            self.fw_he_thickness, self.fw_steel_back,
            self.duct_steel_thickness, self.pbli_gap_thickness,
            self.fci_thickness, self.midplate_he_thickness,
            self.back_he_thickness,
        )
        if any(not math.isfinite(value) or value <= 0 for value in lengths):
            raise ValueError("all radii, lengths, and layer thicknesses must be finite and positive")
        if isinstance(self.n_axial_bins, bool) or not isinstance(self.n_axial_bins, int) or self.n_axial_bins < 1:
            raise ValueError("n_axial_bins must be a positive integer")
        for count in (self.n_fw_he_channels, self.n_midplate_he_channels, self.n_back_he_channels):
            if isinstance(count, bool) or not isinstance(count, int) or count < 1:
                raise ValueError("helium channel counts must be positive integers")
            if 2 * self.he_side_wall_angle_deg + (count - 1) * self.he_rib_angle_deg >= 180:
                raise ValueError("helium ribs and side walls leave no open duct width")
        if self.orientation not in ("left", "right"):
            raise ValueError("orientation must be 'left' or 'right'")

    def build(self) -> DCLLFNSFModule:
        """Build mutually exclusive material cells for the module envelope."""
        mats = self.materials.resolve(self.material_ns)
        cells: dict[str, openmc.Cell] = {}
        channels: list[ChannelRecord] = []
        radius = self.r_inner
        interfaces = [radius]
        layer_names = ["inner"]
        z_min = self.z0 - self.axial_length / 2

        def add_layer(name: str, thickness: float, fill: openmc.Material | None) -> None:
            nonlocal radius
            r0 = radius
            region = core.annular_shell_region(
                z0=self.z0, inner_radius=r0, radial_thickness=thickness,
                axial_length=self.axial_length, orientation=self.orientation,
            )
            cells[name] = openmc.Cell(name=name, region=region, fill=fill)
            radius += thickness
            interfaces.append(radius)
            layer_names.append(name)

        def add_helium_ducts(label: str, thickness: float, count: int) -> None:
            nonlocal radius
            r0 = radius
            envelope = core.annular_shell_region(
                z0=self.z0, inner_radius=r0, radial_thickness=thickness,
                axial_length=self.axial_length, orientation=self.orientation,
            )
            ducts, _ = core.half_annular_axial_channels(
                z0=self.z0, inner_radius=r0, radial_thickness=thickness,
                axial_length=self.axial_length, n_channels=count,
                rib_angle_deg=self.he_rib_angle_deg,
                side_wall_angle_deg=self.he_side_wall_angle_deg,
                orientation=self.orientation,
            )
            steel_name = f"{label}_steel_webs"
            cells[steel_name] = openmc.Cell(
                name=steel_name,
                region=core.annular_shell_structure(envelope, ducts),
                fill=mats.structure,
            )
            for index, region in enumerate(ducts):
                name = f"{label}_helium_az{index + 1}"
                cells[name] = openmc.Cell(name=name, region=region, fill=mats.helium)
                channels.append(ChannelRecord(
                    name=name, kind=f"helium_{label}", row=0, axial_index=index,
                    region=region, r_inner=r0, r_outer=r0 + thickness,
                    z_min=z_min, z_max=z_min + self.axial_length,
                    material=mats.helium,
                ))
            radius += thickness
            interfaces.append(radius)
            layer_names.append(label)

        def add_pbli_core(name: str, thickness: float, row: int) -> None:
            nonlocal radius
            r0 = radius
            dz = self.axial_length / self.n_axial_bins
            for index in range(self.n_axial_bins):
                center = z_min + (index + 0.5) * dz
                cell_name = f"{name}_z{index + 1}"
                region = core.annular_shell_region(
                    z0=center, inner_radius=r0, radial_thickness=thickness,
                    axial_length=dz, orientation=self.orientation,
                )
                cells[cell_name] = openmc.Cell(name=cell_name, region=region, fill=mats.pbli)
                channels.append(ChannelRecord(
                    name=cell_name, kind=name, row=row, axial_index=index,
                    region=region, r_inner=r0, r_outer=r0 + thickness,
                    z_min=center - dz / 2, z_max=center + dz / 2,
                    material=mats.pbli,
                ))
            radius += thickness
            interfaces.append(radius)
            layer_names.append(name)

        def add_duct(label: str, thickness: float, row: int) -> None:
            add_layer(f"{label}_steel_inner", self.duct_steel_thickness, mats.structure)
            add_layer(f"{label}_gap_inner", self.pbli_gap_thickness, mats.pbli)
            add_layer(f"{label}_fci_inner", self.fci_thickness, mats.fci)
            add_pbli_core(f"{label}_pbli", thickness, row)
            add_layer(f"{label}_fci_outer", self.fci_thickness, mats.fci)
            add_layer(f"{label}_gap_outer", self.pbli_gap_thickness, mats.pbli)
            add_layer(f"{label}_steel_outer", self.duct_steel_thickness, mats.structure)

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
        self._resolved_materials = mats
        self._extents = ModuleExtents(
            r_inner=self.r_inner, r_outer=radius, z_min=z_min,
            z_max=z_min + self.axial_length, orientation=self.orientation,
            r_interfaces=tuple(interfaces), layer_names=tuple(layer_names),
        )
        self._universe = openmc.Universe(name="DCLL_FNSF_Module", cells=list(cells.values()))
        self._built = True
        return self

    def _require_built(self) -> None:
        if not self._built:
            raise RuntimeError("Call DCLLFNSFModule.build() before querying geometry")

    def get_channels(self) -> list[ChannelRecord]:
        self._require_built()
        return list(self._channels)

    def get_cells(self) -> dict[str, openmc.Cell]:
        self._require_built()
        return dict(self._cells)

    def get_cells_by_role(self) -> dict[str, list[openmc.Cell]]:
        """Group steel, helium, PbLi, and SiC cells for material tallies."""
        self._require_built()
        groups: dict[str, list[openmc.Cell]] = {
            "structure": [], "helium": [], "pbli_core": [],
            "pbli_gap": [], "fci": [],
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

    def get_universe(self) -> openmc.Universe:
        self._require_built()
        assert self._universe is not None
        return self._universe

    def get_extents(self) -> ModuleExtents:
        self._require_built()
        assert self._extents is not None
        return self._extents

    def get_materials(self) -> DCLLFNSFMaterials:
        self._require_built()
        assert self._resolved_materials is not None
        return self._resolved_materials

    def named_regions(self) -> dict[str, openmc.Region]:
        self._require_built()
        return {name: cell.region for name, cell in self._cells.items() if cell.region is not None}


def build_dcll_fnsf_module(**kwargs: Any) -> DCLLFNSFModule:
    """Construct and build an FNSF-inspired dual-coolant test module."""
    return DCLLFNSFModule(**kwargs).build()


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
    print(f"FW[0]:  {fw0.name}  r=[{fw0.r_min:.3f}, {fw0.r_max:.3f}]  "
          f"z=[{fw0.z_min:.3f}, {fw0.z_max:.3f}]")
    print(f"BULK[0,0]: {bulk0.name}  r=[{bulk0.r_min:.3f}, {bulk0.r_max:.3f}]  "
          f"z=[{bulk0.z_min:.3f}, {bulk0.z_max:.3f}]")
    mesh = fw0.to_cylindrical_mesh((4, 1, 6))
    print(f"sample FW mesh dimension={mesh.dimension}, origin={mesh.origin}")
