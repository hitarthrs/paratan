"""Plot demos for SCLL and DCLL_FNSF testing modules.

Module definitions live in
``src.paratan.models.testing_module_model_builder``. Geometry primitives live in
``paratan.geometry.core``.

Run from the paratan repository root::

    MPLCONFIGDIR=/tmp/mpl-vns .openmc-env/bin/python src/paratan/geometry/module.py --type DCLL_FNSF
"""
from __future__ import annotations

import argparse
from functools import reduce
from operator import or_
from pathlib import Path
import sys

_SRC = Path(__file__).resolve().parents[2]  # .../paratan/src
_ROOT = Path(__file__).resolve().parents[3]  # .../paratan
if str(_SRC) not in sys.path:
    sys.path.insert(0, str(_SRC))
if str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))

import matplotlib.pyplot as plt
from matplotlib.colors import to_rgb
from matplotlib.patches import Patch
import openmc

from src.paratan.models.testing_module_model_builder import (
    BulkChannelRow,
    build_dcll_fnsf_module,
    build_scll_module,
)


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
    legend_colors: dict[str, str] | None = None,
) -> Path:
    """Plot multiple named regions as side-by-side xz and xy cuts."""
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
    if legend_colors:
        fig.legend(
            handles=[Patch(facecolor=color, label=label) for label, color in legend_colors.items()],
            loc="lower center", ncol=len(legend_colors), frameon=False,
        )
        fig.tight_layout(rect=(0, 0.07, 1, 1))
    else:
        fig.tight_layout()
    fig.savefig(output, bbox_inches="tight", dpi=150)
    plt.close(fig)
    return output


# Material-role palette for SCLL plots (aligned with the DCLL figure style).
SCLL_LEGEND = {
    "tungsten armor/divider": "#2f2f37",
    "RAFM structure": "#8a9aab",
    "FW PbLi": "#d4a017",
    "bulk PbLi": "#2f6f4e",
}


def scll_material_colors(
    named_regions: dict[str, openmc.Region],
) -> dict[str, tuple[int, int, int]]:
    """Color SCLL cells by material role for legend-style plots."""
    hex_by_name: dict[str, str] = {}
    for name in named_regions:
        if name in ("armor", "divider"):
            hex_by_name[name] = SCLL_LEGEND["tungsten armor/divider"]
        elif name in ("fw_structure", "bulk_structure"):
            hex_by_name[name] = SCLL_LEGEND["RAFM structure"]
        elif name.startswith("fw_channel"):
            hex_by_name[name] = SCLL_LEGEND["FW PbLi"]
        elif name.startswith("bulk_r"):
            hex_by_name[name] = SCLL_LEGEND["bulk PbLi"]
        else:
            hex_by_name[name] = "#bbbbbb"
    return {
        name: tuple(round(255 * value) for value in to_rgb(color))
        for name, color in hex_by_name.items()
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Plot a testing-module cross section")
    parser.add_argument("--type", choices=("SCLL", "DCLL_FNSF"), default="SCLL")
    module_type = parser.parse_args().type
    if module_type == "DCLL_FNSF":
        module = build_dcll_fnsf_module(
            r_inner=50.0, axial_length=60.0,
            front_pbli_thickness=30.0, rear_pbli_thickness=30.0,
            n_axial_bins=4, orientation="right",
        )
        named_regions = module.named_regions()
        extents = module.get_extents()
        colors = {}
        for name in named_regions:
            if "steel" in name:
                colors[name] = "#6d8194"
            elif "helium" in name:
                colors[name] = "#4cbfdf"
            elif "fci" in name:
                colors[name] = "#e2a354"
            elif "pbli" in name or "gap" in name:
                colors[name] = "#83b66a" if name.startswith("front") else "#4d965b"
        colors = {
            name: tuple(round(255 * value) for value in to_rgb(color))
            for name, color in colors.items()
        }
        wall_mid = 0.5 * (extents.r_inner + extents.r_outer)
        path = plot_regions(
            named_regions,
            outer_radius=extents.r_outer,
            axial_length=module.axial_length,
            colors=colors,
            output_path=Path(__file__).with_name("dcll_fnsf_xz_xy.png"),
            title="DCLL_FNSF: He-cooled RAFM / SiC FCI / front-return PbLi",
            xy_cut_z=0.0,
            radial_span=max(12.0, 2.5 * (extents.r_outer - extents.r_inner)),
            xz_origin=(wall_mid, 5.0, 0.0),
            xy_origin=(wall_mid, 0.0, 0.0),
            legend_colors={
                "RAFM steel": "#6d8194", "helium": "#4cbfdf",
                "SiC insert": "#e2a354", "front PbLi": "#83b66a",
                "return PbLi": "#4d965b",
            },
        )
        print(f"DCLL_FNSF extents: r=[{extents.r_inner:.3f}, {extents.r_outer:.3f}] cm")
        print(f"PbLi axial scoring cells: {2 * module.n_axial_bins}; "
              f"axial He channels: {sum(name.startswith(('fw_helium_', 'midplate_helium_', 'back_helium_')) for name in named_regions)}")
        print(f"Wrote {path}")
        sys.exit(0)

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
        armor_thickness=0.2,
        divider_thickness=0.6,
        bulk_rows=(BulkChannelRow(n_channels=4, radial_thickness=80.0),),
        bulk_front_thickness=0.5,
        bulk_back_thickness=1.0,
        bulk_radial_gap=1.0,
    )
    named_regions = module.named_regions()
    fw_layout = module.get_fw_layout()
    bulk_layout = module.get_bulk_layout()
    extents = module.get_extents()
    colors = scll_material_colors(named_regions)

    z_centers = fw_layout["z_centers"]
    xy_cut_z = float(min(z_centers, key=abs))  # type: ignore[arg-type]

    wall_mid = 0.5 * (extents.r_inner + extents.r_outer)
    wall_span = max(12.0, 2.5 * (extents.r_outer - extents.r_inner))
    path = plot_regions(
        named_regions,
        outer_radius=extents.r_outer,
        axial_length=module.axial_length,
        colors=colors,
        title="SCLL: W armor / RAFM structure / PbLi FW + bulk",
        output_path=Path(__file__).with_name("scll_first_wall_xz_xy.png"),
        xy_cut_z=xy_cut_z,
        radial_span=wall_span,
        xz_origin=(wall_mid, 5.0, 0.0),
        xy_origin=(wall_mid, 0.0, xy_cut_z),
        legend_colors=SCLL_LEGEND,
    )
    print(f"extents: r=[{extents.r_inner:.3f}, {extents.r_outer:.3f}] cm, "
          f"z=[{extents.z_min:.3f}, {extents.z_max:.3f}] cm")
    print(f"r_interfaces ({extents.layer_names}): {extents.r_interfaces}")
    print(f"n_fw_channels: {len(module.get_fw_channels())}")
    print(f"n_bulk_channels: {len(module.get_bulk_channels())}")
    for row_layout in bulk_layout.get("row_layouts", ()):
        assert isinstance(row_layout, dict)
        print(
            f"bulk row {row_layout['row_index'] + 1}: "
            f"n={row_layout['n_channels']}, "
            f"dr={row_layout['radial_thickness']:.3f} cm, "
            f"r_inner={row_layout['inner_radius']:.3f} cm, "
            f"ch_axial={row_layout['channel_axial_length']:.3f} cm"
        )
    print(f"Wrote {path}")
