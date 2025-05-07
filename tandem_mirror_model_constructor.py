import openmc
import numpy as np
import matplotlib.pyplot as plt
import material as m
from geometry_lib import *
from source_lib import *
import yaml
from tandem_geometry import *


with open('input_files/tandem_parametric_input.yaml', 'r') as f:
    input_data = yaml.safe_load(f)


################### Universe that contains the machine #####################

universe_machine = openmc.Universe(786)
with open("input_files/tandem_parametric_input.yaml", "r") as f:
    input_data = yaml.safe_load(f)

vv_params, fw_params, cyl_params, lf_coil_params, hf_coil_params, end_cell_params = parse_machine_input(input_data, m)

builder = TandemMachineBuilder(vv_params, fw_params, cyl_params, lf_coil_params, hf_coil_params, end_cell_params, m)
builder.build()

universe_machine = builder.get_universe()

# ________________________________________ FINALIZE MODEL CROSS-SECTION ________________________________________ #

# Create the OpenMC geometry with the machine universe as the root
geometry = openmc.Geometry([openmc.Cell(fill=universe_machine)])

# Generate and save the cross-sectional plot of the model
geometry.root_universe.plot(
    basis='xz', 
    width=(1000, 5100), 
    pixels=(700, 700), 
    color_by='material', 
    #openmc_exec='/opt/openmc/bin/openmc'  # Specify the OpenMC executable path
)

# Save the generated plot to the results directory
plt.savefig('tandem_mirror_cross_section.png', bbox_inches="tight")

