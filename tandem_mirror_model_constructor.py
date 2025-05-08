import openmc
import numpy as np
import matplotlib.pyplot as plt
import material as m
from geometry_lib import *
from source_lib import *
import yaml
from tandem_geometry import *

# ------------------ Load Input ------------------ #
with open("input_files/tandem_parametric_input.yaml", "r") as f:
    input_data = yaml.safe_load(f)

vv_params, fw_params, cyl_params, lf_coil_params, hf_coil_params, end_cell_params = parse_machine_input(input_data, m)

# ------------------ Build Machine ------------------ #
builder = TandemMachineBuilder(vv_params, fw_params, cyl_params, lf_coil_params, hf_coil_params, end_cell_params, m)
builder.build()

universe_machine = builder.get_universe()
geometry = openmc.Geometry([openmc.Cell(fill=universe_machine)])

# ------------------ Tally Setup ------------------ #
# 3. Get tallies and export to XML
tallies = builder.get_all_tallies()
print(tallies)

# ------------------ Plot Geometry ------------------ #
geometry.root_universe.plot(
    basis='xz', 
    width=(1000, 5100), 
    pixels=(700, 700), 
    color_by='material'
)
plt.savefig('tandem_mirror_cross_section.png', bbox_inches="tight")

# ------------------ Write XMLs ------------------ #
geometry.export_to_xml()
