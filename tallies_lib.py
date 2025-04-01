import numpy as np
import pandas as pd
import openmc
import parametric_input as param
import yaml
from geometry_lib import *

def hollow_mesh_from_domain(domain, dimensions= [10, 10, 10], phi_grid_bounds=(0.0, 2 * np.pi)):
    """
    Generate a cylindrical mesh overs a hollow region defined by an OpenMC region.
    
    Parameters:
        domain (openmc.Region/ openmc.Cell): The domain to bound and mesh (not necessarily hollow).
        dimensions (tuple): Number of divisions in (r, phi, z), i.e., (nr, nphi, nz).
        phi_grid_bounds (tuple): Angular bounds in radians for phi. Default is (0, 2π).
    
    Returns:
        openmc.CylindricalMesh: A cylindrical mesh over the hollow region.
    """
    # Get the bounding box of the region
    bounding_box = domain.bounding_box
    
    # Determine max radial extent from bounding box corners
    max_radius = max(
        bounding_box[0][0],  # x-min
        bounding_box[0][1],  # y-min
        bounding_box[1][0],  # x-max
        bounding_box[1][1]   # y-max
    )
    
    # Create outer bounding cylindrical surfaces
    outer_cylinder = openmc.ZCylinder(r=max_radius)
    lower_z = openmc.ZPlane(bounding_box[0][2])
    upper_z = openmc.ZPlane(bounding_box[1][2])
    
    outer_region = -outer_cylinder & +lower_z & -upper_z

    # Check if the type of domain is an openmc.Cell
    if type(domain) == openmc.Cell:
        region = domain.region
    
    # Subtract the original region to define hollow space
    hollow_region = outer_region & ~region
    
    # Extract all surfaces in the resulting region
    surfaces = hollow_region.get_surfaces()
    
    # Find all z-cylindrical surfaces and collect their radii
    radii = [
        surface.coefficients['r']
        for surface in surfaces.values()
        if surface.type == 'z-cylinder'
    ]
    
    # Set inner radius based on smallest detected cylindrical surface
    if radii:
        min_radius = min(radii)
    else:
        min_radius = 0.0  # fallback if no cylinders are found
    
    # Build the r, phi, z grids
    r_grid = np.linspace(min_radius, max_radius, num=dimensions[0] + 1)
    phi_grid = np.linspace(phi_grid_bounds[0], phi_grid_bounds[1], num=dimensions[1] + 1)
    z_grid = np.linspace(bounding_box[0][2], bounding_box[1][2], num=dimensions[2] + 1)


    origin = (bounding_box.center[0], bounding_box.center[1], z_grid[0])

    z_grid -= origin[2]

    # Construct and return the cylindrical mesh

    cyl_mesh = openmc.CylindricalMesh(r_grid=r_grid, phi_grid=phi_grid, z_grid=z_grid, origin=origin)
    
    return cyl_mesh

def strings_to_openmc_filters(filter_strings: list[str]):
    filters = []

    for filter_string in filter_strings:
        if filter_string == "photon_filter":
            filters.append(openmc.ParticleFilter("photon"))

        elif filter_string == "neutron_filter":
            filters.append(openmc.ParticleFilter("neutron"))

        elif filter_string == "fast_energies_filter":
            edges = np.logspace(np.log10(1e5), np.log10(20e6), 501)  # eV
            filters.append(openmc.EnergyFilter(edges))

        elif filter_string == "thermal_energies_filter":
            edges = np.logspace(np.log10(1e-3), np.log10(1e3), 501)  # eV
            filters.append(openmc.EnergyFilter(edges))

        elif filter_string == "vitaminj_filter":
            group_bounds = openmc.mgxs.GROUP_STRUCTURES["VITAMIN-J"]
            filters.append(openmc.EnergyFilter(group_bounds[::-1]))  # descending

        else:
            raise ValueError(f"Unknown filter string: {filter_string}")

    return filters


cell = openmc.Cell()
def generate_cell_tallies_from_region(input_data: dict, region_name: str, cell: openmc.Cell):
    """
    Given a YAML-loaded dictionary, a region name like 'lf_coil_tallies',
    and an OpenMC cell, generate the corresponding list of tallies.
    """
    cell_tally_list = input_data.get(region_name, {}).get("cell_tallies", [])
    mesh_tally_list = input_data.get(region_name, {}).get("mesh_tallies", [])
    tallies = []

    for i, entry in enumerate(cell_tally_list):
        scores = entry.get("scores", [])
        filter_strings = entry.get("filters", [])
        nuclides = entry.get("nuclides", [])

        filters = [openmc.CellFilter(cell)] + strings_to_openmc_filters(filter_strings)

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


        filters = [openmc.MeshFilter(cylindrical_mesh)] + strings_to_openmc_filters(filter_strings)

        tally = openmc.Tally(name=f"{region_name}_mesh_tally_{i+1}")
        tally.filters = filters
        tally.scores = scores

        if nuclides:
            tally.nuclides = nuclides

        tallies.append(tally)

    return tallies

with open('parametric_input.yaml', 'r') as f:
    input_data = yaml.safe_load(f)

lf_coil_shell, lf_coil_inner = hollow_cylinder_with_shell(
        2,      # Coil center position
        10,  # Outer reference radius
        15,    # Radial thickness of the coil
        15,  # Inner axial length
        3,  # Shell front thickness
        3,   # Shell back thickness
        3   # Shell axial thickness
    )

lf_coil_chell_cell = openmc.Cell(region = lf_coil_shell)

tallies = generate_cell_tallies_from_region(input_data, 'lf_coil_tallies', lf_coil_chell_cell)

