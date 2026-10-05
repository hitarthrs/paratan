"""Build an SCLL module via the model builder and plot xz/xy slices.

Uses ``testing_module_model_builder`` for geometry. Coloring matches the
DCLL material-legend style (tungsten / RAFM / PbLi).

Run from the paratan repository root::

    MPLCONFIGDIR=/tmp/mpl-vns .openmc-env/bin/python tests/plot_scll_module_builder.py
"""
from __future__ import annotations

from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.paratan.geometry.module import (
    SCLL_LEGEND,
    plot_regions,
    scll_material_colors,
)
from src.paratan.models.testing_module_model_builder import (
    BulkChannelRow,
    build_scll_module,
)


def main() -> Path:
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

    assert len(module.get_fw_channels()) == 20
    assert len(module.get_bulk_channels()) == 4
    assert len(module.get_universe().cells) == 28
    extents = module.get_extents()
    fw0 = module.get_channel_extents("FW", axial_index=0)
    bulk0 = module.get_channel_extents("BULK", row=0, axial_index=0)
    assert fw0.r_min < fw0.r_max and fw0.z_min < fw0.z_max
    assert bulk0.r_min < bulk0.r_max and bulk0.z_min < bulk0.z_max

    named_regions = module.named_regions()
    fw_layout = module.get_fw_layout()
    colors = scll_material_colors(named_regions)

    z_centers = fw_layout["z_centers"]
    xy_cut_z = float(min(z_centers, key=abs))  # type: ignore[arg-type]
    wall_mid = 0.5 * (extents.r_inner + extents.r_outer)
    wall_span = max(12.0, 2.5 * (extents.r_outer - extents.r_inner))

    output_dir = ROOT / "plots"
    output_dir.mkdir(exist_ok=True)
    path = plot_regions(
        named_regions,
        outer_radius=extents.r_outer,
        axial_length=module.axial_length,
        colors=colors,
        title="SCLL: W armor / RAFM structure / PbLi FW + bulk",
        output_path=output_dir / "scll_module_builder_xz_xy.png",
        xy_cut_z=xy_cut_z,
        radial_span=wall_span,
        xz_origin=(wall_mid, 5.0, 0.0),
        xy_origin=(wall_mid, 0.0, xy_cut_z),
        legend_colors=SCLL_LEGEND,
    )

    print(f"extents: r=[{extents.r_inner:.3f}, {extents.r_outer:.3f}] cm, "
          f"z=[{extents.z_min:.3f}, {extents.z_max:.3f}] cm")
    print(f"FW[0] mesh box: r=[{fw0.r_min:.3f}, {fw0.r_max:.3f}] "
          f"z=[{fw0.z_min:.3f}, {fw0.z_max:.3f}]")
    print(f"BULK[0,0] mesh box: r=[{bulk0.r_min:.3f}, {bulk0.r_max:.3f}] "
          f"z=[{bulk0.z_min:.3f}, {bulk0.z_max:.3f}]")
    print(f"Wrote {path}")
    return path


if __name__ == "__main__":
    main()
