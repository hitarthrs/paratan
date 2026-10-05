"""Plot four copies of geometry/module.py's CSG in the tandem VNS.

Run from the paratan repository root:
    MPLCONFIGDIR=/tmp/mpl-vns .openmc-env/bin/python examples/plot_tandem_vns_modules.py

Geometry plotting only: no particle transport, source, or model.run(). Module
cells are void-filled and colored by their intended role for this layout study.
"""
from pathlib import Path
from functools import reduce
from operator import or_
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import to_rgb
from matplotlib.patches import Patch
import openmc
import yaml

from src.paratan.geometry.module import core
from src.paratan.models.tandem_vns_model_builder import build_tandem_vns_universe


def main():
    with (ROOT / "input_files/tandem_vns_parametric_input.yaml").open() as stream:
        data = yaml.safe_load(stream)
    data["central_cell"]["test_modules"]["count"] = 4
    builder = build_tandem_vns_universe(data)
    machine = builder.get_universe()
    colors = {cell: "#b9c2ca" for cell in machine.cells.values()}
    colors[machine.cells[69]] = "white"
    for cell in builder.vv_builder.get_cells().values():
        colors[cell] = "#f1f5f8"
    for cells in builder.fw_builder.get_cells_by_region().values():
        for cell in cells:
            colors[cell] = "#455a64"
    for cell in machine.cells.values():
        if 4100 <= cell.id < 4200:
            colors[cell] = "#7553a1"

    inner_radius = (
        data["vacuum_vessel"]["central_cell"]["central_radius"]
        + sum(layer["thickness"] for layer in data["first_wall"]["central_cell"]["layers"])
        + 0.5
    )
    # LF shielding starts at the central test envelope's outer radius.
    outer_radius = float(builder.ccyl_builder.get_radii_by_region()["central"][-1])
    channel_radial_thickness = (outer_radius - inner_radius - 2.0 - 3.0 - 2 * 1.5) / 3
    if channel_radial_thickness <= 0:
        raise ValueError("Insufficient space between the first wall and LF shielding")
    # In core.redefined_vacuum_vessel_region the constant-radius cylinder
    # extends to +/- outer_axial_length / 2; the cones begin there.
    cc_length = data["vacuum_vessel"]["central_cell"]["outer_axial_length"]
    module_length = cc_length / 2
    axial_centers = (-cc_length / 4, cc_length / 4)
    module_colors = ("#238b8e", "#e28b26", "#507dcc", "#ce567a")
    cells, envelopes, layouts = [], [], []
    for z0 in axial_centers:
        for side in ("left", "right"):
            module_index = len(envelopes)
            channels, layout = core.annular_shell_channel_grid(
                z0=z0, module_inner_radius=inner_radius, axial_length=module_length,
                n_axial=4, n_radial=3, channel_radial_thickness=channel_radial_thickness,
                radial_gap=1.5, first_wall_thickness=2.0,
                outer_wall_thickness=3.0, axial_gap=0.4,
                axial_end_thickness=1.0, orientation=side,
            )
            envelope = core.annular_shell_region(
                z0=z0, inner_radius=inner_radius,
                radial_thickness=layout["module_radial_thickness"],
                axial_length=module_length, orientation=side,
            )
            structure = openmc.Cell(
                name=f"{side}_{z0:g}_module_structure",
                region=envelope & ~reduce(or_, channels),
            )
            cells.append(structure)
            colors[structure] = "#34424f"
            for index, region in enumerate(channels):
                cell = openmc.Cell(name=f"{side}_{z0:g}_pbli_{index}", region=region)
                cells.append(cell)
                colors[cell] = module_colors[module_index]
            envelopes.append(envelope)
            layouts.append(layout)
    background = openmc.Cell(name="space_between_test_modules", region=~reduce(or_, envelopes))
    cells.append(background)
    colors[background] = "#fcfcfc"
    builder.set_test_module_universe(openmc.Universe(cells=cells))
    colors = {cell: tuple(round(255 * v) for v in to_rgb(color)) for cell, color in colors.items()}

    output = ROOT / "plots"
    output.mkdir(exist_ok=True)
    # Actual OpenMC CSG slices, with each leaf cell colored explicitly.
    fig, axes = plt.subplots(1, 3, figsize=(17, 10), gridspec_kw={"width_ratios": [0.85, 1.15, 1.3]})
    cut_z = float(layouts[2]["z_centers"][1])
    specs = [
        ("xz", (0, 0, 0), (600, 3500), (550, 2200), "Whole tandem VNS", "x [cm]", "z [cm]"),
        ("xz", (0, 0, 0), (460, 660), (1500, 2200), "Central cylinder: four modules", "x [cm]", "z [cm]"),
        ("xy", (0, 0, cut_z), (460, 460), (1600, 1600), f"Transverse section at z = {cut_z:.1f} cm", "x [cm]", "y [cm]"),
    ]
    for ax, (basis, origin, width, pixels, title, xlabel, ylabel) in zip(axes, specs):
        machine.plot(basis=basis, origin=origin, width=width, pixels=pixels,
                     color_by="cell", colors=colors, axes=ax,
                     openmc_exec=str(Path(sys.executable).parent / "openmc"))
        ax.set(title=title, xlabel=xlabel, ylabel=ylabel)
    for station, z0 in enumerate(axial_centers):
        for sign, side in ((-1, "L"), (1, "R")):
            module_number = 2 * station + (1 if side == "L" else 2)
            axes[1].text(sign * (inner_radius + outer_radius) / 2, z0,
                         f"M{module_number}", ha="center", va="center", fontsize=12,
                         bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.85})
    axes[2].axvline(0, color="#999999", linestyle=":", linewidth=0.7)
    fig.suptitle("Tandem VNS — four half-annular PbLi test modules\nOpenMC CSG geometry slices", fontsize=17)
    fig.legend(handles=[*[Patch(color=color, label=f"M{i+1}: {side}, {span} cm")
                          for i, (color, side, span) in enumerate(zip(
                              module_colors, ("left", "right", "left", "right"),
                              (f"z = {-cc_length/2:g} to 0", f"z = {-cc_length/2:g} to 0",
                               f"z = 0 to {cc_length/2:g}", f"z = 0 to {cc_length/2:g}")))],
                        Patch(color="#34424f", label="Module walls and ribs"),
                        Patch(color="#455a64", label="Facility first wall"),
                        Patch(color="#7553a1", label="LF coils")],
               loc="lower center", ncol=3, frameon=False)
    fig.tight_layout(rect=(0, 0.09, 1, 0.94))
    path = output / "tandem_vns_four_modules.png"
    fig.savefig(path, dpi=180, bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote {path}")
    print(f"4 modules; 48 PbLi cells; inner/outer radius {inner_radius}/{outer_radius} cm")
    print(f"Each PbLi radial row: {channel_radial_thickness:.3f} cm")
    print(f"Each module length: {module_length:g} cm; combined coverage: +/-{cc_length/2:g} cm")


if __name__ == "__main__":
    main()
