import numpy as np
import openmc
from mirror_paratan.geometry.core import hollow_mesh_from_domain

def strings_to_openmc_filters(filter_strings: list[str], filter_cache: dict[str, object] | None = None):
    filters = []
    cache = {} if filter_cache is None else filter_cache

    for filter_string in filter_strings:
        if filter_string not in cache:
            if filter_string == "photon_filter":
                cache[filter_string] = openmc.ParticleFilter("photon")
            elif filter_string == "neutron_filter":
                cache[filter_string] = openmc.ParticleFilter("neutron")
            elif filter_string == "fast_energies_filter":
                edges = np.array([1e5, 20e6])
                cache[filter_string] = openmc.EnergyFilter(edges)
            elif filter_string == "thermal_energies_filter":
                edges = np.logspace(np.log10(1e-3), np.log10(1e3), 501)
                cache[filter_string] = openmc.EnergyFilter(edges)
            elif filter_string == "vitaminj_filter":
                group_bounds = openmc.mgxs.GROUP_STRUCTURES["VITAMIN-J-175"]
                cache[filter_string] = openmc.EnergyFilter(group_bounds)
            else:
                raise ValueError(f"Unknown filter string: {filter_string}")
        filters.append(cache[filter_string])

    return filters


def generate_cell_tallies_from_region(input_data: dict, region_name: str, cell: openmc.Cell):
    """
    Given a YAML-loaded dictionary, a region name like 'lf_coil_tallies',
    and an OpenMC cell, generate the corresponding list of tallies.
    """
    cell_tally_list = input_data.get(region_name, {}).get("cell_tallies", [])
    mesh_tally_list = input_data.get(region_name, {}).get("mesh_tallies", [])
    tallies = []
    filter_cache: dict[str, object] = {}

    for i, entry in enumerate(cell_tally_list):
        scores = entry.get("scores", [])
        filter_strings = entry.get("filters", [])
        nuclides = entry.get("nuclides", [])

        filters = [openmc.CellFilter(cell)] + strings_to_openmc_filters(filter_strings, filter_cache)

        tally = openmc.Tally(name=f"{region_name}_cell_tally_{i+1}")
        tally.filters = filters
        tally.scores = scores

        if nuclides:
            tally.nuclides = nuclides

        tallies.append(tally)

    for i, entry in enumerate(mesh_tally_list):
        scores = entry.get("scores", [])
        filter_strings = entry.get("filters", [])
        nuclides = entry.get("nuclides", [])
        dimensions = entry.get("dimensions", [])
        print(dimensions)

        cylindrical_mesh = hollow_mesh_from_domain(domain = cell,dimensions=dimensions)


        filters = [openmc.MeshFilter(cylindrical_mesh)] + strings_to_openmc_filters(filter_strings, filter_cache)

        tally = openmc.Tally(name=f"{region_name}_mesh_tally_{i+1}")
        tally.filters = filters
        tally.scores = scores

        if nuclides:
            tally.nuclides = nuclides

        tallies.append(tally)

    return tallies
class TallyBuilder:
    """
    Collects tally descriptors and generates OpenMC Tally objects.
    Supports both cell and mesh tallies.
    """

    def __init__(self):
        self._tallies = []
        self._named_filter_cache: dict[str, object] = {}

    def add_descriptors(self, descriptors):
        for desc in descriptors:
            self._add_cell_tallies(desc)
            self._add_mesh_tallies(desc)

    def _add_cell_tallies(self, desc):
        for i, entry in enumerate(desc.get("cell_tallies", [])):
            filters = [openmc.CellFilter(desc["cell"])] + strings_to_openmc_filters(entry.get("filters", []), self._named_filter_cache)
            tally = openmc.Tally(name=f"{desc['type']}_{desc['location']}_{desc['description']}_cell_tally_{i+1}")
            tally.filters = filters
            tally.scores = entry.get("scores", [])
            if "nuclides" in entry:
                tally.nuclides = entry["nuclides"]
            self._tallies.append(tally)

    def _add_mesh_tallies(self, desc):
        for i, entry in enumerate(desc.get("mesh_tallies", [])):
            mesh = hollow_mesh_from_domain(desc["cell"], entry["dimensions"])
            filters = [openmc.MeshFilter(mesh)] + strings_to_openmc_filters(entry.get("filters", []), self._named_filter_cache)
            tally = openmc.Tally(name=f"{desc['type']}_{desc['location']}_{desc['description']}_mesh_tally_{i+1}")
            tally.filters = filters
            tally.scores = entry.get("scores", [])
            if "nuclides" in entry:
                tally.nuclides = entry["nuclides"]
            self._tallies.append(tally)

    def get_tallies(self):
        return self._tallies
# with open('parametric_input.yaml', 'r') as f:
#     input_data = yaml.safe_load(f)

# lf_coil_shell, lf_coil_inner = hollow_cylinder_with_shell(
#         2,      # Coil center position
#         10,  # Outer reference radius
#         15,    # Radial thickness of the coil
#         15,  # Inner axial length
#         3,  # Shell front thickness
#         3,   # Shell back thickness
#         3   # Shell axial thickness
#     )

# lf_coil_chell_cell = openmc.Cell(region = lf_coil_shell)

# tallies = generate_cell_tallies_from_region(input_data, 'lf_coil_tallies', lf_coil_chell_cell)
