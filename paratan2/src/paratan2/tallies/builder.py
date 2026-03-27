import openmc
from paratan2.tallies.base import hollow_mesh_from_domain, strings_to_openmc_filters


class TallyBuilder:
    """
    Consumes tally descriptor dicts from component builders and produces
    openmc.Tally objects. Fully score-agnostic — whatever is in the descriptor
    gets passed straight through to OpenMC.
    """

    def __init__(self):
        self._tallies: list[openmc.Tally] = []

    def add_descriptors(self, descriptors: list[dict]) -> None:
        for desc in descriptors:
            self._add_cell_tallies(desc)
            self._add_mesh_tallies(desc)

    def _add_cell_tallies(self, desc: dict) -> None:
        for i, entry in enumerate(desc.get("cell_tallies", [])):
            filters = [openmc.CellFilter(desc["cell"])] + strings_to_openmc_filters(
                entry.get("filters", [])
            )
            tally = openmc.Tally(
                name=f"{desc['type']}_{desc['location']}_{desc['description']}_cell_tally_{i+1}"
            )
            tally.filters = filters
            tally.scores = entry.get("scores", [])
            if entry.get("nuclides"):
                tally.nuclides = entry["nuclides"]
            self._tallies.append(tally)

    def _add_mesh_tallies(self, desc: dict) -> None:
        for i, entry in enumerate(desc.get("mesh_tallies", [])):
            mesh = hollow_mesh_from_domain(desc["cell"], entry.get("dimensions", [10, 10, 10]))
            filters = [openmc.MeshFilter(mesh)] + strings_to_openmc_filters(
                entry.get("filters", [])
            )
            tally = openmc.Tally(
                name=f"{desc['type']}_{desc['location']}_{desc['description']}_mesh_tally_{i+1}"
            )
            tally.filters = filters
            tally.scores = entry.get("scores", [])
            if entry.get("nuclides"):
                tally.nuclides = entry["nuclides"]
            self._tallies.append(tally)

    def get_tallies(self) -> list[openmc.Tally]:
        return self._tallies
