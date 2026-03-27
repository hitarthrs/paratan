"""
Quick trial build — constructs the simple mirror geometry and shows
the cross-section plot so you can visually confirm it looks right.

Run from the paratan2 directory:
    python trial_build.py
"""

import matplotlib.pyplot as plt
import openmc

from paratan2.config.models import SimpleMachineConfig
from paratan2.models.simple_machine_builder import SimpleMachineBuilder, _load_materials

YAML = "/home/hrshah3/paratan/paratan/compare_hf_study_runs/hpc_run/simple_parametric_input_new.yaml"
OUTPUT_DIR = "/tmp/paratan2_trial"

import os
os.makedirs(OUTPUT_DIR, exist_ok=True)

# --- Parse config ---
print("Loading config...")
config = SimpleMachineConfig.from_yaml(YAML)
print(f"  VV central radius:    {config.vacuum_vessel.central_radius} cm")
print(f"  VV bottleneck radius: {config.vacuum_vessel.bottleneck_radius} cm")
print(f"  CC axial length:      {config.central_cell.axial_length} cm")
print(f"  LF coil positions:    {config.lf_coil.positions}")
print(f"  HF magnet thickness:  {config.hf_coil.magnet.axial_thickness} cm")

# --- Build geometry ---
print("\nBuilding geometry...")
m = _load_materials()
builder = SimpleMachineBuilder(config, m)
builder.build()

geometry = builder.get_geometry()
print(f"  Total cells: {len(builder.all_cells)}")

# --- Plot ---
print("\nGenerating cross-section plot...")
fig, axes = plt.subplots(1, 2, figsize=(18, 8))

# XZ plane (side view — the interesting one for a mirror machine)
geometry.root_universe.plot(
    basis="xz",
    origin=(0, 0, 0),
    width=(900, 2800),
    pixels=(1800, 2800),
    color_by="material",
    axes=axes[0],
)
axes[0].set_title("Side view (XZ)", fontsize=13)
axes[0].set_xlabel("Z (cm)")
axes[0].set_ylabel("X (cm)")

# XY plane (end-on cross section at midplane)
geometry.root_universe.plot(
    basis="xy",
    origin=(0, 0, 0),
    width=(700, 700),
    pixels=(800, 800),
    color_by="material",
    axes=axes[1],
)
axes[1].set_title("End-on view (XY) at midplane", fontsize=13)
axes[1].set_xlabel("X (cm)")
axes[1].set_ylabel("Y (cm)")

plt.suptitle("Paratan 2 — Simple Mirror Geometry Trial", fontsize=15, fontweight="bold")
plt.tight_layout()

plot_path = os.path.join(OUTPUT_DIR, "trial_geometry.png")
plt.savefig(plot_path, dpi=150, bbox_inches="tight")
print(f"\nPlot saved to: {plot_path}")
plt.show()

print("\nDone.")
