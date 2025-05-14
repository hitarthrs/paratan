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

# ------------------ Plot Geometry ------------------ #
geometry.root_universe.plot(
    basis='xz', 
    width=(1000, 5100), 
    pixels=(1000, 1200), 
    color_by='material'
)
plt.savefig('tandem_mirror_cross_section.png', bbox_inches="tight")


materials = m.materials
# ------------------ Write XMLs ------------------ #
geometry.merge_surfaces = True
geometry.export_to_xml()



with open('input_files/source_information.yaml', 'r') as f:
    source_data = yaml.safe_load(f)

# openmc_source = load_source_from_yaml('input_files/source_information.yaml')

source = openmc.Source()

source.space = openmc.stats.CylindricalIndependent(
    r=openmc.stats.Discrete(
        [0], [1.0]
    ),
    phi=openmc.stats.Uniform(a=0.0, b=2 * np.pi),
    z=openmc.stats.Uniform(a=-600.0, b=600.0),
          # axis-aligned at center
)

source.angle = openmc.stats.Isotropic()
source.energy = openmc.stats.Discrete([14e6], [1.0]) 


settings = openmc.Settings()
settings.run_mode = "fixed source"
settings.particles = int(source_data['settings']['particles_per_batch'])
settings.batches = source_data['settings']['batches']
settings.output = {'tallies': False}
settings.statepoint = {
    'batches': [1] + list(range(source_data['settings']['statepoint_frequency'], settings.batches, source_data['settings']['statepoint_frequency'])) + [settings.batches]
}
settings.weight_windows_on = source_data['settings']['weight_windows']
settings.weight_window_checkpoints = {'collision': True, 'surface': True}
wwg = openmc.WeightWindowGenerator(openmc.RegularMesh(), 
                                   [0, 14e6], 
                                   'neutron', 
                                   'magic', 
                                   max_realizations=25, 
                                   update_interval=1, 
                                   on_the_fly=True)
settings.weight_windows_generator = [wwg]
settings.photon_transport = source_data['settings']['photon_transport']
settings.source = source

# Export the finalized geometry to an OpenMC XML file
settings.export_to_xml("settings.xml")

tallies = openmc.Tallies(tallies)

model = openmc.Model(geometry, materials, settings, tallies)

model.geometry.merge_surfaces = True
model.run(geometry_debug = True)


