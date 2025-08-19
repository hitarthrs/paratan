import openmc
import numpy as np
import matplotlib.pyplot as plt
import material as m
from src.paratan.geometry.core import *
from src.paratan.source.core import *
import yaml
from simple_model_builder import *

# ------------------ Load Input ------------------ #
with open("input_files/simple_parametric_input.yaml", "r") as f:
    input_data = yaml.safe_load(f)

vv_params, cc_params, lf_params, hf_params, end_params = parse_simple_machine_input(input_data, m)

# ------------------ Build Machine ------------------ #
builder = SimpleMachineBuilder(
    vv_params, cc_params, lf_params, hf_params, end_params, m
)

builder.build_vacuum_vessel()
builder.build_central_cell()
builder.build_lf_coils()
builder.build_hf_coils()
builder.build_end_cells()

# ------------------ Geometry and Universe ------------------ #
universe_machine = builder.get_universe()
geometry = openmc.Geometry([openmc.Cell(fill=universe_machine)])
geometry.merge_surfaces = True

# ------------------ Plot Geometry ------------------ #
geometry.root_universe.plot(
    basis='xz',
    width=(1000, 2800),
    pixels=(700, 700),
    color_by='material'
)
plt.savefig('simple_mirror_modular_cross_section.png', bbox_inches="tight")

# ------------------ Tallies ------------------ #
tallies = openmc.Tallies(builder.get_all_tallies())

# ------------------ Materials ------------------ #
materials = m.materials

# ------------------ Source ------------------ #
with open('input_files/source_information.yaml', 'r') as f:
    source_data = yaml.safe_load(f)

source = openmc.Source()
source.space = openmc.stats.CylindricalIndependent(
    r=openmc.stats.Discrete([0], [1.0]),
    phi=openmc.stats.Uniform(a=0.0, b=2 * np.pi),
    z=openmc.stats.Uniform(a=-600.0, b=600.0),
)
source.angle = openmc.stats.Isotropic()
source.energy = openmc.stats.Discrete([14e6], [1.0])

# ------------------ Settings ------------------ #
settings = openmc.Settings()
settings.run_mode = "fixed source"
settings.particles = int(source_data['settings']['particles_per_batch'])
settings.batches = source_data['settings']['batches']
settings.output = {'tallies': False}

freq = source_data['settings']['statepoint_frequency']
settings.statepoint = {
    'batches': [1] + list(range(freq, settings.batches, freq)) + [settings.batches]
}

settings.weight_windows_on = source_data['settings']['weight_windows']
settings.weight_window_checkpoints = {'collision': True, 'surface': True}

# Placeholder weight window generator setup (adjust mesh as needed)
wwg = openmc.WeightWindowGenerator(
    openmc.RegularMesh(),
    energy_bounds=[0, 14e6],
    particle_type='neutron',
    method='magic',
    max_realizations=25,
    update_interval=1,
    on_the_fly=True
)
settings.weight_windows_generator = [wwg]

settings.photon_transport = source_data['settings']['photon_transport']
settings.source = source

# ------------------ Final Model ------------------ #
model = openmc.Model(geometry, materials, settings, tallies)
model.export_to_xml('model_xml_files')

# ------------------ Run Simulation ------------------ #
model.run(geometry_debug=False) 