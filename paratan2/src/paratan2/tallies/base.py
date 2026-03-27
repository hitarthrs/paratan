import numpy as np
import openmc


def hollow_mesh_from_domain(domain, dimensions=(10, 10, 10), phi_grid_bounds=(0.0, 2 * np.pi)):
    """
    Generate a cylindrical mesh over a region defined by an OpenMC cell or region.

    Parameters
    ----------
    domain : openmc.Cell or openmc.Region
    dimensions : tuple of int — (nr, nphi, nz)
    phi_grid_bounds : tuple of float

    Returns
    -------
    openmc.CylindricalMesh
    """
    bounding_box = domain.bounding_box

    region = domain.region if isinstance(domain, openmc.Cell) else domain

    surfaces = region.get_surfaces()
    radii = [
        s.coefficients["r"]
        for s in surfaces.values()
        if s.type == "z-cylinder"
    ]

    if radii:
        max_radius = max(radii)
    else:
        max_radius = max(
            abs(bounding_box[0][0]),
            abs(bounding_box[0][1]),
            abs(bounding_box[1][0]),
            abs(bounding_box[1][1]),
        )

    outer_cylinder = openmc.ZCylinder(r=max_radius)
    lower_z = openmc.ZPlane(bounding_box[0][2])
    upper_z = openmc.ZPlane(bounding_box[1][2])
    outer_region = -outer_cylinder & +lower_z & -upper_z
    hollow_region = outer_region & ~region

    hollow_surfaces = hollow_region.get_surfaces()
    hollow_radii = [
        s.coefficients["r"]
        for s in hollow_surfaces.values()
        if s.type == "z-cylinder"
    ]
    min_radius = min(hollow_radii) if hollow_radii else 0.0

    r_grid = np.linspace(min_radius, max_radius, dimensions[0] + 1)
    phi_grid = np.linspace(phi_grid_bounds[0], phi_grid_bounds[1], dimensions[1] + 1)
    z_grid = np.linspace(bounding_box[0][2], bounding_box[1][2], dimensions[2] + 1)

    origin = (bounding_box.center[0], bounding_box.center[1], z_grid[0])
    z_grid -= origin[2]

    return openmc.CylindricalMesh(
        r_grid=r_grid, phi_grid=phi_grid, z_grid=z_grid, origin=origin
    )


def strings_to_openmc_filters(filter_strings: list[str]) -> list:
    filters = []
    for fs in filter_strings:
        if fs == "photon_filter":
            filters.append(openmc.ParticleFilter("photon"))
        elif fs == "neutron_filter":
            filters.append(openmc.ParticleFilter("neutron"))
        elif fs == "fast_energies_filter":
            filters.append(openmc.EnergyFilter([1e5, 20e6]))
        elif fs == "thermal_energies_filter":
            edges = np.logspace(np.log10(1e-3), np.log10(1e3), 501)
            filters.append(openmc.EnergyFilter(edges))
        elif fs == "vitaminj_filter":
            group_bounds = openmc.mgxs.GROUP_STRUCTURES["VITAMIN-J-175"]
            filters.append(openmc.EnergyFilter(group_bounds))
        else:
            raise ValueError(f"Unknown filter string: {fs!r}")
    return filters
