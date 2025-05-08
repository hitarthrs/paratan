import openmc
from tallies_lib import hollow_mesh_from_domain, strings_to_openmc_filters

class TallyBuilder:
    """
    Collects tally descriptors and generates OpenMC Tally objects.
    Supports both cell and mesh tallies.
    """

    def __init__(self):
        self._tallies = []

    def add_descriptors(self, descriptors):
        for desc in descriptors:
            self._add_cell_tallies(desc)
            self._add_mesh_tallies(desc)

    def _add_cell_tallies(self, desc):
        for i, entry in enumerate(desc.get("cell_tallies", [])):
            filters = [openmc.CellFilter(desc["cell"])] + strings_to_openmc_filters(entry.get("filters", []))
            tally = openmc.Tally(name=f"{desc['type']}_{desc['location']}_cell_{i+1}")
            tally.filters = filters
            tally.scores = entry.get("scores", [])
            if "nuclides" in entry:
                tally.nuclides = entry["nuclides"]
            self._tallies.append(tally)

    def _add_mesh_tallies(self, desc):
        for i, entry in enumerate(desc.get("mesh_tallies", [])):
            mesh = hollow_mesh_from_domain(desc["cell"], entry["dimensions"])
            filters = [openmc.MeshFilter(mesh)] + strings_to_openmc_filters(entry.get("filters", []))
            if desc['type'] == "hf_coil":
                tally = openmc.Tally(name=f"{desc['type']}_{desc['location']}_{desc['direction']}_mesh_{i+1}")
            else:
                tally = openmc.Tally(name=f"{desc['type']}_{desc['location']}_mesh_{i+1}")
            tally.filters = filters
            tally.scores = entry.get("scores", [])
            if "nuclides" in entry:
                tally.nuclides = entry["nuclides"]
            self._tallies.append(tally)

    def get_tallies(self):
        return self._tallies