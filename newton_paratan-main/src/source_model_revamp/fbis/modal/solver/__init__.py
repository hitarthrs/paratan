"""Modal FBIS solver entry points"""
from source_model_revamp.fbis.modal.solver.species_solver import solve_fast_ion_species
from source_model_revamp.fbis.modal.solver.system_solver import solve_fast_ion_system

__all__ = ["solve_fast_ion_species", "solve_fast_ion_system"]
