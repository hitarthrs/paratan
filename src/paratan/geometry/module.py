"""Test-module CSG builders for the tandem VNS.

Lengths are centimetres. The SCLL module builds armour → first wall → divider →
bulk PbLi rows. Bulk rows are specified individually (channel count and radial
thickness per row), not as a uniform n×m grid.
"""
from __future__ import annotations

from dataclasses import dataclass
from functools import reduce
from operator import or_
from pathlib import Path
import sys

# Allow `python module.py` from this directory without an editable install.
_SRC = Path(__file__).resolve().parents[2]  # .../paratan/src
if str(_SRC) not in sys.path:
    sys.path.insert(0, str(_SRC))

import matplotlib.pyplot as plt
from matplotlib.colors import to_rgb
import openmc
import paratan.geometry.core as core


@dataclass(frozen=True)
class BulkChannelRow:
    """One radial row in the bulk PbLi pack.

    ``n_channels`` axial channels of radial height ``radial_thickness`` tile the
    module axial length (with the pack's axial gaps / end walls). Channel axial
    length is derived: more channels → shorter segments.
    """

    n_channels: int
    radial_thickness: float

    def __post_init__(self) -> None:
        if self.n_channels < 1:
            raise ValueError("BulkChannelRow.n_channels must be >= 1")
        if self.radial_thickness <= 0:
            raise ValueError("BulkChannelRow.radial_thickness must be positive")


@dataclass(frozen=True)
class SCLLModule:
    """Self-cooled PbLi half-annular test module.

    Radial build order (plasma → outward):

    1. Solid armour of ``armor_thickness`` ending at ``r_fw``.
    2. First-wall layer from ``r_fw`` with one radial row of ``n_fw_channels``
       axial channels, bounded by steel front and back thicknesses.
    3. Solid divider of ``divider_thickness`` immediately behind the first wall
       (``divider_material`` label; tungsten for now).
    4. Bulk PbLi pack: each entry in ``bulk_rows`` is one radial row with its
       own channel count and radial thickness.
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
    armor_thickness: float = 0.2  # 2 mm, cm
    divider_thickness: float = 0.6  # 6 mm, cm
    divider_material: str = "tungsten"
    bulk_rows: tuple[BulkChannelRow, ...] = ()
    bulk_front_thickness: float = 0.5
    bulk_back_thickness: float = 1.0
    bulk_radial_gap: float = 1.0
    bulk_axial_gap: float | None = None
    bulk_axial_end_thickness: float | None = None

    def __post_init__(self) -> None:
        if self.r_fw <= 0:
            raise ValueError("r_fw must be positive")
        if self.armor_thickness <= 0:
            raise ValueError("armor_thickness must be positive")
        if self.divider_thickness <= 0:
            raise ValueError("divider_thickness must be positive")
        if not self.divider_material:
            raise ValueError("divider_material must be a non-empty name")
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
            axial_gap = self._bulk_axial_gap
            axial_end = self._bulk_axial_end_thickness
            if axial_gap < 0 or axial_end < 0:
                raise ValueError("bulk axial gap/end thickness must be >= 0")

    @property
    def _bulk_axial_gap(self) -> float:
        return self.axial_gap if self.bulk_axial_gap is None else self.bulk_axial_gap

    @property
    def _bulk_axial_end_thickness(self) -> float:
        if self.bulk_axial_end_thickness is None:
            return self.axial_end_thickness
        return self.bulk_axial_end_thickness

    @property
    def r_armor(self) -> float:
        """Inner radius of the armour (plasma-facing surface)."""
        return self.r_fw - self.armor_thickness

    @property
    def fw_radial_thickness(self) -> float:
        """Outer radius of the first-wall layer minus ``r_fw``."""
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
        """Inner radius of the divider (outer radius of the first wall)."""
        return self.fw_outer_radius

    @property
    def divider_outer_radius(self) -> float:
        return self.r_divider + self.divider_thickness

    @property
    def r_bulk(self) -> float:
        """Inner radius of the bulk pack (outer radius of the divider)."""
        return self.divider_outer_radius

    @property
    def bulk_radial_thickness(self) -> float:
        """Full radial thickness of the bulk pack, or 0 if no rows."""
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
        """Outermost radius of the layers built so far."""
        if self.bulk_rows:
            return self.bulk_outer_radius
        return self.divider_outer_radius

    def armor_region(self) -> openmc.Region:
        """Solid half-annular armour immediately inside the first wall."""
        return core.annular_shell_region(
            z0=self.z0,
            inner_radius=self.r_armor,
            radial_thickness=self.armor_thickness,
            axial_length=self.axial_length,
            orientation=self.orientation,
        )

    def first_wall_envelope(self) -> openmc.Region:
        """Half-annular region occupied by the entire first-wall layer."""
        return core.annular_shell_region(
            z0=self.z0,
            inner_radius=self.r_fw,
            radial_thickness=self.fw_radial_thickness,
            axial_length=self.axial_length,
            orientation=self.orientation,
        )

    def first_wall_channels(self) -> tuple[tuple[openmc.Region, ...], dict[str, object]]:
        """One radial row of ``n_fw_channels`` axial channels inside the first wall.

        Returns
        -------
        channels
            Length-``n_fw_channels`` tuple of channel regions, ordered by increasing z.
        layout
            Packing metadata from :func:`core.annular_shell_channel_grid`.
        """
        return core.annular_shell_channel_grid(
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

    def first_wall_structure(self) -> openmc.Region:
        """Steel volume of the first wall: envelope minus the channel voids."""
        channels, _ = self.first_wall_channels()
        return self.first_wall_envelope() & ~reduce(or_, channels)

    def divider_region(self) -> openmc.Region:
        """Solid half-annular divider immediately behind the first wall."""
        return core.annular_shell_region(
            z0=self.z0,
            inner_radius=self.r_divider,
            radial_thickness=self.divider_thickness,
            axial_length=self.axial_length,
            orientation=self.orientation,
        )

    def bulk_envelope(self) -> openmc.Region | None:
        """Half-annular region occupied by the entire bulk pack, or None if empty."""
        if not self.bulk_rows:
            return None
        return core.annular_shell_region(
            z0=self.z0,
            inner_radius=self.r_bulk,
            radial_thickness=self.bulk_radial_thickness,
            axial_length=self.axial_length,
            orientation=self.orientation,
        )

    def bulk_channels(self) -> tuple[tuple[openmc.Region, ...], dict[str, object]]:
        """Build bulk PbLi channels from per-row specs.

        Returns
        -------
        channels
            Flat tuple in row-major order: radial row slowest, axial index fastest.
        layout
            Packing metadata including per-row z centers and channel axial lengths.
        """
        if not self.bulk_rows:
            return (), {
                "row_inner_radii": (),
                "row_layouts": (),
                "bulk_radial_thickness": 0.0,
            }

        axial_gap = self._bulk_axial_gap
        axial_end = self._bulk_axial_end_thickness
        channels: list[openmc.Region] = []
        row_inner_radii: list[float] = []
        row_layouts: list[dict[str, object]] = []

        r_inner = self.r_bulk + self.bulk_front_thickness
        for row_index, row in enumerate(self.bulk_rows):
            z_centers, channel_axial_length = core.axial_channel_centers_and_length(
                z0=self.z0,
                axial_length=self.axial_length,
                n_axial=row.n_channels,
                axial_gap=axial_gap,
                axial_end_thickness=axial_end,
            )
            row_inner_radii.append(r_inner)
            for z_center in z_centers:
                channels.append(
                    core.annular_shell_region(
                        z0=z_center,
                        inner_radius=r_inner,
                        radial_thickness=row.radial_thickness,
                        axial_length=channel_axial_length,
                        orientation=self.orientation,
                    )
                )
            row_layouts.append(
                {
                    "row_index": row_index,
                    "n_channels": row.n_channels,
                    "radial_thickness": row.radial_thickness,
                    "inner_radius": r_inner,
                    "z_centers": z_centers,
                    "channel_axial_length": channel_axial_length,
                }
            )
            r_inner += row.radial_thickness
            if row_index < len(self.bulk_rows) - 1:
                r_inner += self.bulk_radial_gap

        layout: dict[str, object] = {
            "row_inner_radii": tuple(row_inner_radii),
            "row_layouts": tuple(row_layouts),
            "bulk_radial_thickness": self.bulk_radial_thickness,
            "bulk_front_thickness": self.bulk_front_thickness,
            "bulk_back_thickness": self.bulk_back_thickness,
            "bulk_radial_gap": self.bulk_radial_gap,
            "axial_gap": axial_gap,
            "axial_end_thickness": axial_end,
        }
        return tuple(channels), layout

    def bulk_structure(self) -> openmc.Region | None:
        """Steel volume of the bulk pack: envelope minus channel voids."""
        envelope = self.bulk_envelope()
        if envelope is None:
            return None
        channels, _ = self.bulk_channels()
        if not channels:
            return envelope
        return envelope & ~reduce(or_, channels)

    def named_regions(self) -> dict[str, openmc.Region]:
        """Named regions for plotting / later cell assembly."""
        fw_channels, _ = self.first_wall_channels()
        named: dict[str, openmc.Region] = {
            "armor": self.armor_region(),
            "fw_structure": self.first_wall_structure(),
            "divider": self.divider_region(),
        }
        for index, region in enumerate(fw_channels):
            named[f"fw_channel_z{index + 1}"] = region

        if self.bulk_rows:
            bulk_channels, layout = self.bulk_channels()
            structure = self.bulk_structure()
            if structure is not None:
                named["bulk_structure"] = structure
            offset = 0
            for row_layout in layout["row_layouts"]:
                assert isinstance(row_layout, dict)
                n_ch = int(row_layout["n_channels"])
                row_index = int(row_layout["row_index"])
                for axial_index, region in enumerate(bulk_channels[offset : offset + n_ch]):
                    named[f"bulk_r{row_index + 1}_z{axial_index + 1}"] = region
                offset += n_ch
        return named

    # Back-compat alias used by earlier demos.
    def named_first_wall_regions(self) -> dict[str, openmc.Region]:
        return self.named_regions()


def plot_regions(
    named_regions: dict[str, openmc.Region],
    *,
    outer_radius: float,
    axial_length: float,
    colors: dict[str, str | tuple[int, int, int]],
    output_path: Path | str,
    title: str = "regions",
    xy_cut_z: float = 0.0,
    radial_span: float | None = None,
    xz_origin: tuple[float, float, float] | None = None,
    xy_origin: tuple[float, float, float] | None = None,
) -> Path:
    """Plot multiple named regions as side-by-side xz and xy cuts.

    ``xy_cut_z`` sets the z of the xy slice so it can pass through a channel
    instead of an axial gap / end wall.
    """
    cells = [
        openmc.Cell(name=name, region=region, fill=None)
        for name, region in named_regions.items()
    ]
    union = reduce(or_, named_regions.values())
    void = openmc.Cell(name="outside", region=~union, fill=None)
    universe = openmc.Universe(cells=[*cells, void])

    cell_colors = {cell: colors[cell.name] for cell in cells}
    cell_colors[void] = "white"

    if radial_span is None:
        radial_span = 2.2 * outer_radius
    axial_span = 1.4 * axial_length
    if xz_origin is None:
        xz_origin = (0.0, 0.0, 0.0)
    if xy_origin is None:
        xy_origin = (0.0, 0.0, xy_cut_z)

    fig, axes = plt.subplots(1, 2, figsize=(12, 5.5))
    plot_specs = (
        ("xz", xz_origin, (radial_span, axial_span), "x [cm]", "z [cm]", None),
        ("xy", xy_origin, (radial_span, radial_span), "x [cm]", "y [cm]", xy_cut_z),
    )
    openmc_exec = str(Path(sys.executable).parent / "openmc")
    for ax, (basis, origin, width, xlabel, ylabel, cut_z) in zip(axes, plot_specs):
        universe.plot(
            basis=basis,
            origin=origin,
            width=width,
            pixels=(800, 800),
            color_by="cell",
            colors=cell_colors,
            axes=ax,
            openmc_exec=openmc_exec,
        )
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        if cut_z is None:
            ax.set_title(f"{title} ({basis})")
        else:
            ax.set_title(f"{title} ({basis}, z={cut_z:.2f} cm)")

    output = Path(output_path)
    fig.tight_layout()
    fig.savefig(output, bbox_inches="tight", dpi=150)
    plt.close(fig)
    return output


if __name__ == "__main__":
    module = SCLLModule(
        r_fw=50.0,
        n_fw_channels=20,
        axial_length=60.0,
        channel_radial_thickness=2.0,
        front_thickness=0.4,
        back_thickness=0.4,
        orientation="right",
        axial_gap=0.4,
        axial_end_thickness=1.0,
        armor_thickness=0.2,  # 2 mm
        divider_thickness=0.6,  # 6 mm
        divider_material="tungsten",
        # Per-row bulk: give channel count + radial thickness for each row.
        bulk_rows=(
            BulkChannelRow(n_channels=4, radial_thickness=80.0),
        ),
        bulk_front_thickness=0.5,
        bulk_back_thickness=1.0,
        bulk_radial_gap=1.0,
    )
    fw_channels, fw_layout = module.first_wall_channels()
    bulk_channels, bulk_layout = module.bulk_channels()
    named_regions = module.named_regions()

    colors: dict[str, str | tuple[int, int, int]] = {
        "armor": tuple(int(255 * c) for c in to_rgb("dimgray")),
        "fw_structure": "steelblue",
        "divider": tuple(int(255 * c) for c in to_rgb("gray")),
        "bulk_structure": tuple(int(255 * c) for c in to_rgb("slategray")),
    }
    channel_cmap = plt.cm.YlOrBr
    for index in range(module.n_fw_channels):
        name = f"fw_channel_z{index + 1}"
        frac = (index + 1) / (module.n_fw_channels + 1)
        rgba = channel_cmap(0.30 + 0.55 * frac)
        colors[name] = tuple(int(255 * c) for c in to_rgb(rgba))

    bulk_cmap = plt.cm.YlGn
    for name in named_regions:
        if name.startswith("bulk_r"):
            # Color by radial row index embedded in the name.
            row_token = name.split("_")[1]  # r1, r2, ...
            row_index = int(row_token[1:]) - 1
            frac = (row_index + 1) / (len(module.bulk_rows) + 1)
            rgba = bulk_cmap(0.25 + 0.55 * frac)
            colors[name] = tuple(int(255 * c) for c in to_rgb(rgba))

    z_centers = fw_layout["z_centers"]
    xy_cut_z = float(min(z_centers, key=abs))

    wall_mid = 0.5 * (module.r_armor + module.module_outer_radius)
    wall_span = max(12.0, 2.5 * (module.module_outer_radius - module.r_armor))
    path = plot_regions(
        named_regions,
        outer_radius=module.module_outer_radius,
        axial_length=module.axial_length,
        colors=colors,
        title=(
            f"SCLL · armor + FW + divider + bulk "
            f"({', '.join(str(r.n_channels) for r in module.bulk_rows)} ch/row)"
        ),
        output_path=Path(__file__).with_name("scll_first_wall_xz_xy.png"),
        xy_cut_z=xy_cut_z,
        radial_span=wall_span,
        xz_origin=(wall_mid, 0.0, 0.0),
        xy_origin=(wall_mid, 0.0, xy_cut_z),
    )
    print(f"r_armor: {module.r_armor:.3f} cm")
    print(f"r_fw: {module.r_fw:.3f} cm")
    print(f"fw outer / divider outer: {module.fw_outer_radius:.3f} / {module.divider_outer_radius:.3f} cm")
    print(f"bulk outer / module outer: {module.bulk_outer_radius:.3f} / {module.module_outer_radius:.3f} cm")
    print(f"n_fw_channels: {module.n_fw_channels}")
    for row_layout in bulk_layout["row_layouts"]:
        print(
            f"bulk row {row_layout['row_index'] + 1}: "
            f"n={row_layout['n_channels']}, "
            f"dr={row_layout['radial_thickness']:.3f} cm, "
            f"r_inner={row_layout['inner_radius']:.3f} cm, "
            f"ch_axial={row_layout['channel_axial_length']:.3f} cm"
        )
    print(f"Wrote {path}")
