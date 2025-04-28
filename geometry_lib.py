import openmc
import numpy as np
import openmc.lib
import matplotlib.pyplot as plt

def sample_source(model, n_samples=1000, seed=1):
    """ Initializes the OpenMC problem and samples source sites

    Parameters
    ----------
    model : openmc.Model
        The model object used for the source. Assumed to be a fixed source problem.

    n_samples: Integral
        The number of source sites to sample. Default is 1000.

    seed: Integral
        The random number seed used in OpenMC. Default is 1.
    """
    with TemporaryDirectory() as d:
        model.export_to_model_xml()
        try:
            openmc.lib.init()
            openmc.lib.settings.seed = seed
            openmc.lib.simulation_init()
            sites = openmc.lib.sample_external_source(n_samples)
        finally:
            openmc.lib.finalize()

    return sites

def plot_sites(ax, sites, basis='xy', arrows=False):

    """Add sites to the matplotlib axes

    Parameters
    ----------
    ax : matplotlib.axes object
        Generally, the axes returned from `openmc.Universe.plot`

    sites : Iterable of openmc.SourceParticle
        The source sites to plot

    basis : One of ('xy, 'yz', 'xz')
        The basis of the plot. Default is 'xy'.

    arrows : boolmodel.export_to_xml('os.abs')
        Whether or not to plot the sites as arrows or dots.
    """

    basis_indices = {'xy': [0, 1],
                     'yz': [1, 2],
                     'xz': [0, 2]}

    indices = basis_indices[basis]

    for site in sites:
        x, y = np.asarray(site.r)[indices]
        u, v = np.asarray(site.u)[indices]
        if arrows:
            ax.arrow(x, y, u, v, head_width=0.1)
        else:
            ax.plot(x, y,  marker='o', markerfacecolor='blue')


def plot_geom_and_source(model, source_sites, cell_colors=None):

    """Plots the geometry and source sites.

    Two plots are produced:

        1. A plot in the Y-Z plane with source sites plotted as arrows

        2. A plot in the X-Y plane with source sites plotted as dots to show the spatial distributions.
    """
    for basis in ('yz', 'xy'):
        ax = model.geometry.root_universe.plot(basis=basis,
                                        pixels=(600, 600),
                                        colors=cell_colors,
                                        legend=cell_colors is not None)

        plot_sites(ax, source_sites, basis=basis, arrows=basis == 'xz')

        plt.show()

def vacuum_vessel_region(first_plane_distance, outermost_plane_distance, central_radius, outer_radius, angle_degrees):

    """Generates an OpenMC region in the shape of a typical vacuum vessel shape. Generates a symmetrical shape around the midplane with the following structure:
    1. A central cylindrical section.
    2. The cylinder tapers outward into a cone.
    3. Beyond the cone, a thin cylindrical segment extends outward.

    This sequence of cylinder-cone-thin cylinder is mirrored on both sides of the midplane, creating an inverted-hourglass-like geometry.

    Parameters
    ------------

    first_plane_distance (cm): Distance of the end of the first cylindrical section from the midplane
    outermost_plane distance (cm) : Distance of the end of the last cylindrical section from the midplane
    central_radius (cm) : Radius of the first (central) cylindrical section
    outer_radius (cm) : Radius of the outer cylindrical section
    angle_degrees (degrees) : Angle between the slant surface of the cone and the plane perpendicular to the central axis of the cone

    """ 

    # if (first_plane_distance >= second_plane_distance | second_plane_distance >= outermost_plane_distance):
    #     raise ValueError("Check distances of planes from the midplane. Ideal order is first_plane_distance < second_plane_distance < outermost_plane_distance.")
    
    if (central_radius < outer_radius):
        raise ValueError("The radius of central cylinder should not be lesser than the radius of the outer cylindrical section.")
    
    angle = (np.pi/180)*angle_degrees

    second_plane_distance = first_plane_distance + (central_radius - outer_radius) * np.tan(angle)

    central_cylinder_left_plane = openmc.ZPlane(-first_plane_distance)
    central_cylinder_right_plane = openmc.ZPlane(first_plane_distance)
    central_cell_cylinder = openmc.ZCylinder(r = central_radius)
    central_cylinder = -central_cell_cylinder & +central_cylinder_left_plane & -central_cylinder_right_plane
    
    left_outer_cylinder_1 = openmc.ZPlane(-outermost_plane_distance)
    right_outer_cylinder_1 = openmc.ZPlane(-second_plane_distance)
    
    left_outer_cylinder_2 = openmc.ZPlane(second_plane_distance)
    right_outer_cylinder_2 = openmc.ZPlane(outermost_plane_distance)
    
    outer_cylinder = openmc.ZCylinder(r = outer_radius)
    
    outer_cylinders_region = -outer_cylinder & ((+left_outer_cylinder_1 & -right_outer_cylinder_1)| (+left_outer_cylinder_2 & -right_outer_cylinder_2))
    
    left_cone = openmc.model.ZConeOneSided(x0=0.0, y0=0.0, z0= -(central_radius/np.tan(angle)+first_plane_distance), r2=(np.tan(angle))**2)
    right_cone = openmc.model.ZConeOneSided(x0=0.0, y0=0.0, z0=central_radius/np.tan(angle)+first_plane_distance, r2=(np.tan(angle))**2, up=False)
    
    left_cone_region = -left_cone & -central_cylinder_left_plane
    right_cone_region = -right_cone & +central_cylinder_right_plane
    
    vessel_region = left_cone_region | right_cone_region | outer_cylinders_region | central_cylinder
    
    return vessel_region

def hollow_cylindrical_region(z0, innermost_radius, inner_radial_thickness, axial_length):

    """ Generates an OpenMC region in the shape of a hollow cylindrical region. 
    This geometry is symmetrical around the midplane (z = z0) and consists of the following structure:

    1. A central hollow cylindrical section with a defined axial and radial thickness
    2. Axial boundaries for the region allowing for controlled extensions along the z-axis

    Parameters
    ------------
    z0 (cm): Axial position of the central midplane (z = z0) for the geometry.
    innermost_radius (cm): Radius of the innermost boundary of the shell.
    inner_radial_thickness (cm): Radial thickness of the inner cylindrical section.
    axial_length (cm): Axial thickness of the inner cylindrical section

    Returns
    ------------
    hollowed_cylinder: A hollow cylindrical region
    """    
    left_boundary = openmc.ZPlane(z0-axial_length/2)
    right_boundary = openmc.ZPlane(z0+axial_length/2)
    
    outer_cylinder = openmc.ZCylinder(r=innermost_radius+inner_radial_thickness)
    inner_cylinder = openmc.ZCylinder(r=innermost_radius)
    
    hollowed_cylinder = +left_boundary & -right_boundary & (-outer_cylinder & +inner_cylinder)
    
    return hollowed_cylinder


def hollow_cylinder_with_shell(z0, innermost_radius, inner_radial_thickness,inner_axial_thickness, shell_thickness_front, shell_thickness_back, shell_axial_thickness):

    """ Generates an OpenMC region in the shape of a cylindrical shell with distinct inner and shell regions. 
    This geometry is symmetrical around the midplane (z = z0) and consists of the following structure:
    
    1. A central inner cylindrical section with a defined axial and radial thickness.
    2. An outer cylindrical shell surrounding the inner cylinder, with adjustable front and back thicknesses.
    3. Axial boundaries for both inner and outer regions, allowing for a controlled extension along the z-axis.
    
    This configuration creates a hollow cylindrical shell structure, with the inner and shell regions 
    mirroring each other symmetrically around the midplane.

    Parameters
    ------------
    z0 (cm): Axial position of the central midplane (z = z0) for the geometry.
    innermost_radius (cm): Radius of the innermost boundary of the shell.
    inner_radial_thickness (cm): Radial thickness of the inner cylindrical section.
    inner_axial_thickness (cm): Axial thickness of the inner cylindrical section.
    shell_thickness_front (cm): Radial thickness of the shell on the front (inner) side.
    shell_thickness_back (cm): Radial thickness of the shell on the back (outer) side.
    shell_axial_thickness (cm): Additional axial thickness of the outer boundary along the z-axis.

    Returns
    ------------
    tuple:
        shell_region : The annular region representing the outer shell, located between the inner and outer cylindrical sections.
        inner_region : The inner cylindrical region defined within the innermost boundaries.
    """
        
    inner_boundary_left = openmc.ZPlane(z0-inner_axial_thickness/2)
    outer_boundary_left = openmc.ZPlane(z0-inner_axial_thickness/2-shell_axial_thickness)
    
    inner_boundary_right = openmc.ZPlane(z0+inner_axial_thickness/2)
    outer_boundary_right = openmc.ZPlane(z0+inner_axial_thickness/2+shell_axial_thickness)
    
    inner_cylinder_front = openmc.ZCylinder(r=innermost_radius)
    inner_cylinder_back = openmc.ZCylinder(r=innermost_radius+shell_thickness_front)
    
    outer_cylinder_front = openmc.ZCylinder(r=innermost_radius+shell_thickness_front+inner_radial_thickness)
    outer_cylinder_back = openmc.ZCylinder(r=innermost_radius+shell_thickness_front+inner_radial_thickness+shell_thickness_back)
    
    inner_region = (+inner_cylinder_back & -outer_cylinder_front) & -inner_boundary_right & +inner_boundary_left
    
    outer_region = (+inner_cylinder_front & -outer_cylinder_back) & -outer_boundary_right & +outer_boundary_left
    
    shell_region = outer_region & ~inner_region
    
    return shell_region , inner_region

def cylinder_with_shell(z0, innermost_radius, inner_axial_thickness, shell_thickness):
    """ Generates an OpenMC region representing a cylindrical structure with an outer shell of uniform radial and axial thickness. 
    This geometry is symmetrical around the midplane (z = z0) and consists of the following structure:
    
    1. An inner cylindrical region that is hollow (no material inside), bounded axially and radially.
    2. An outer cylindrical shell surrounding the inner hollow cylinder, with adjustable shell thickness along both radial and axial directions.
    
    This configuration creates a hollow cylinder with a defined outer shell around it

    Parameters
    ------------
    z0 (cm) : Axial position of the central midplane (z = z0) for the geometry.
    innermost_radius (cm) : Radius of the innermost boundary of the hollowed cylinder.
    inner_axial_thickness (cm) : Axial thickness of the inner cylindrical section.
    shell_thickness (cm) : Radial and axial thickness of the surrounding shell.

    Returns
    ------------
    tuple:
        outer_shell_region : The annular region representing the shell around the cylinder.
        inner_cylinder_region : The region of the inner cylinder within the defined boundaries.
    """
    inner_boundary_left = openmc.ZPlane(z0-inner_axial_thickness/2)
    outer_boundary_left = openmc.ZPlane(z0-inner_axial_thickness/2-shell_thickness)
    
    inner_boundary_right = openmc.ZPlane(z0+inner_axial_thickness/2)
    outer_boundary_right = openmc.ZPlane(z0+inner_axial_thickness/2+shell_thickness)
    
    inner_cylinder = openmc.ZCylinder(r=innermost_radius)
    outer_cylinder = openmc.ZCylinder(r=innermost_radius+shell_thickness)
    
    inner_cylinder_region = +inner_boundary_left & -inner_boundary_right & -inner_cylinder
    outer_cylinder_region = +outer_boundary_left & -outer_boundary_right & -outer_cylinder
    
    outer_shell_region = outer_cylinder_region & ~inner_cylinder_region
    
    return outer_shell_region, inner_cylinder_region

def cylindrical_region_no_outer_surface(z0, innermost_radius, inner_radial_thickness, inner_axial_thickness, shell_thickness):

    """Generates an OpenMC region representing a cylindrical shell without a defined outer surface.
    This geometry is symmetrical around the midplane (z = z0) and consists of the following structure:
    
    1. An inner cylindrical section with specified axial and radial thickness.
    2. A shell region surrounding the inner section, extending both radially and axially.
    3. The outer surface of the shell is open (not enclosed by an additional boundary surface).
    
    This configuration is useful for simulations where only a partial shell is needed around a central cylindrical core,
    without the requirement for a fully enclosed outer boundary.

    Parameters
    ------------
    z0 (cm) : Axial position of the central midplane (z = z0) for the geometry.
    innermost_radius (cm) : Radius of the innermost boundary of the cylindrical region.
    inner_radial_thickness (cm) : Radial thickness of the inner cylindrical section.
    inner_axial_thickness (cm) : Axial thickness of the inner cylindrical section.
    shell_thickness (cm): Radial and axial thickness of the surrounding shell.

    Returns
    ------------
    shell_region : The region representing the shell around the inner cylindrical section, without a defined outer boundary.
    """

    inner_boundary_left = openmc.ZPlane(z0-inner_axial_thickness/2)
    outer_boundary_left = openmc.ZPlane(z0-inner_axial_thickness/2-shell_thickness)
    
    inner_boundary_right = openmc.ZPlane(z0+inner_axial_thickness/2)
    outer_boundary_right = openmc.ZPlane(z0+inner_axial_thickness/2+shell_thickness)
    
    inner_cylinder_front = openmc.ZCylinder(r=innermost_radius)
    inner_cylinder_back = openmc.ZCylinder(r=innermost_radius+shell_thickness)
    
    outer_cylinder_front = openmc.ZCylinder(r=innermost_radius+shell_thickness+inner_radial_thickness)
    
    inner_region = (+inner_cylinder_back & -outer_cylinder_front) & -inner_boundary_right & +inner_boundary_left
    
    outer_region = (+inner_cylinder_front & -outer_cylinder_front) & -outer_boundary_right & +outer_boundary_left
    
    shell_region = outer_region & ~inner_region
    
    return shell_region

def nested_cylindrical_shells(z0, innermost_radius, inner_radial_thickness,inner_axial_thickness, layer_front_thickness, layer_back_thickness, layer_axial_thickness):
    """Generates a set of nested OpenMC regions representing cylindrical shells with user-defined inner and outer boundaries.
    This geometry is symmetrical around the midplane (z = z0) and consists of multiple layers structured as follows:
    
    1. An innermost cylindrical core with specified axial and radial thickness.
    2. Multiple shell regions surrounding the core, each with specified radial and axial thicknesses.
    3. Each shell region is defined by distinct inner and outer boundaries, extending outward from the core.
       The inner radius specifies the radius of the inner side of the outermost shell.

    This configuration is useful for simulations requiring layered, concentric cylindrical geometries, such as 
    neutron transport or radiation shielding studies, where each layer represents a material with different properties.

    Parameters
    ------------
    z0 (cm) : Axial position of the central midplane (z = z0) for the geometry (in cm).
        
    inner_radius (cm) : Radius of the inner side of the outermost shell (in cm).
        
    inner_radial_thickness (cm) : Radial thickness of the innermost cylindrical core (in cm).
        
    inner_axial_thickness (cm) : Axial thickness of the innermost cylindrical core (in cm).
        
    layer_front_thickness  (list of float [cm]) : Radial thicknesses for the front side of each successive shell layer (in cm).
        
    layer_back_thickness (list of float [cm]) : Radial thicknesses for the back side of each successive shell layer (in cm).
        
    layer_axial_thickness (list of float [cm]) : Axial thicknesses for each successive shell layer (in cm).

    Returns
    ------------
    regions : list of openmc.Region
        A list of OpenMC region objects, where each region represents a nested cylindrical shell with distinct boundaries.
    
    Notes
    -----
    This function constructs a set of nested cylindrical shells, each bounded by radial and axial surfaces.
    The outermost boundary is defined only up to the last specified layer's thickness, without an enclosing outer boundary.

    """
    
    inner_left_layer = z0-inner_axial_thickness/2
    left_boundaries = [openmc.ZPlane(inner_left_layer)]

    inner_right_layer = z0+inner_axial_thickness/2
    right_boundaries = [openmc.ZPlane(inner_right_layer)]

    for i in range(len(layer_axial_thickness)):
        inner_left_layer -= layer_axial_thickness[i]
        inner_right_layer += layer_axial_thickness[i]

        left_plane = openmc.ZPlane(inner_left_layer)
        right_plane = openmc.ZPlane(inner_right_layer)


        left_boundaries.append(left_plane)
        right_boundaries.append(right_plane)


    regions = []

    inner_cylinders_radii_array = np.cumsum(layer_front_thickness[::-1])[::-1]+innermost_radius
    inner_cylinders_radii_array = np.append(inner_cylinders_radii_array, inner_cylinders_radii_array[-1]-layer_front_thickness[-1])


    outer_cylinders_radii_array = np.cumsum(layer_back_thickness)+inner_cylinders_radii_array[0]+inner_radial_thickness
    outer_cylinders_radii_array = np.append(outer_cylinders_radii_array, outer_cylinders_radii_array[-1]+layer_back_thickness[-1])


    innermost_cylinder_front = openmc.ZCylinder(r = inner_cylinders_radii_array[0])
    innermost_cylinder_back = openmc.ZCylinder(r = outer_cylinders_radii_array[0])

    innermost_region = (-innermost_cylinder_back & +innermost_cylinder_front) & +left_boundaries[0] & -right_boundaries[0]
    regions.append(innermost_region)


    for i in range(1,len(inner_cylinders_radii_array)):

        old_cylinder_front = openmc.ZCylinder(r = inner_cylinders_radii_array[i-1])
        old_cylinder_back = openmc.ZCylinder(r = outer_cylinders_radii_array[i-1])
        
        new_cylinder_front = openmc.ZCylinder(r = inner_cylinders_radii_array[i])
        new_cylinder_back = openmc.ZCylinder(r = outer_cylinders_radii_array[i])

        old_cylinder = (-old_cylinder_back & +old_cylinder_front) & +left_boundaries[i-1] & -right_boundaries[i-1]
        new_cylinder = (-new_cylinder_back & +new_cylinder_front) & +left_boundaries[i] & -right_boundaries[i]
        
        new_region = new_cylinder & ~old_cylinder

        regions.append(new_region) 

        
    return regions


def hollow_mesh_from_domain(region, dimensions= [10, 10, 10], phi_grid_bounds=(0.0, 2 * np.pi)):
    """
    Generate a cylindrical mesh overs a hollow region defined by an OpenMC region.
    
    Parameters:
        region (openmc.Region): The region to bound and mesh (not necessarily hollow).
        dimensions (tuple): Number of divisions in (r, phi, z), i.e., (nr, nphi, nz).
        phi_grid_bounds (tuple): Angular bounds in radians for phi. Default is (0, 2π).
    
    Returns:
        openmc.CylindricalMesh: A cylindrical mesh over the hollow region.
    """
    # Get the bounding box of the region
    bounding_box = region.bounding_box
    
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

def redefined_vacuum_vessel_region(outer_axial_length, central_axial_length, central_radius, bottleneck_radius, left_bottleneck_length, right_bottleneck_length, axial_midplane = 0.0):

    """
    Generates an OpenMC region representing a typical vacuum vessel shape for tandem mirror devices.

    The geometry is symmetric about an axial midplane and consists of:
      1. A central cylindrical section.
      2. A conical taper connecting to an outer bottleneck cylindrical section.
      3. Thin outer cylindrical sections extending outward.

    Parameters
    ----------
    outer_axial_length : float
        Total axial length (cm) of the outer cylindrical segments beyond the cone sections.
    central_axial_length : float
        Total axial length (cm) of the central cylindrical section.
    central_radius : float
        Radius (cm) of the central cylindrical section.
    bottleneck_radius : float
        Radius (cm) of the outer (bottleneck) cylindrical sections.
    left_bottleneck_length : float
        Axial length (cm) of the bottleneck section on the left side of the midplane.
    right_bottleneck_length : float
        Axial length (cm) of the bottleneck section on the right side of the midplane.
    axial_midplane : float, optional
        Z-coordinate (cm) of the midplane. Default is 0.0.

    Returns
    -------
    openmc.Region
        An OpenMC region object representing the full vacuum vessel shape.

    Raises
    ------
    ValueError
        If any input length or radius is non-positive.
        If central_radius is less than bottleneck_radius.
    """

    # --- Input validation ---
    if outer_axial_length <= 0 or central_axial_length <= 0:
        raise ValueError("Axial lengths must be positive.")
    if central_radius <= 0 or bottleneck_radius <= 0:
        raise ValueError("Radii must be positive.")
    if left_bottleneck_length <= 0 or right_bottleneck_length <= 0:
        raise ValueError("Bottleneck lengths must be positive.")
    if central_radius < bottleneck_radius:
        raise ValueError("Central radius must be greater than bottleneck radius.")

    # --- Set important z-positions ---
    first_plane_distance = outer_axial_length / 2.0
    second_plane_distance = central_axial_length / 2.0
    
    angle = np.arctan((central_axial_length - outer_axial_length)/(central_radius - bottleneck_radius))

    first_plane_distance = outer_axial_length/2
    second_plane_distance = central_axial_length/2

    right_outermost_plane_distance = second_plane_distance + right_bottleneck_length
    left_outermost_plane_distance = -second_plane_distance - left_bottleneck_length

    central_cylinder_left_plane = openmc.ZPlane(z0 = axial_midplane -first_plane_distance)
    central_cylinder_right_plane = openmc.ZPlane(z0 = axial_midplane + first_plane_distance)
    central_cell_cylinder = openmc.ZCylinder(r = central_radius)
    central_cylinder = -central_cell_cylinder & +central_cylinder_left_plane & -central_cylinder_right_plane
    
    left_outer_cylinder_1 = openmc.ZPlane(z0 = axial_midplane + left_outermost_plane_distance)
    right_outer_cylinder_1 = openmc.ZPlane(z0 = axial_midplane - second_plane_distance)
    
    left_outer_cylinder_2 = openmc.ZPlane(z0 = axial_midplane + second_plane_distance)
    right_outer_cylinder_2 = openmc.ZPlane(z0 = axial_midplane + right_outermost_plane_distance)
    
    outer_cylinder = openmc.ZCylinder(r = bottleneck_radius)
    
    outer_cylinders_region = -outer_cylinder & ((+left_outer_cylinder_1 & -right_outer_cylinder_1)| (+left_outer_cylinder_2 & -right_outer_cylinder_2))
    
    left_cone = openmc.model.ZConeOneSided(x0=0.0, y0=0.0, z0= axial_midplane - (central_radius/np.tan(angle)+first_plane_distance), r2=(np.tan(angle))**2)
    right_cone = openmc.model.ZConeOneSided(x0=0.0, y0=0.0, z0=axial_midplane + (central_radius/np.tan(angle)+first_plane_distance), r2=(np.tan(angle))**2, up=False)
    
    left_cone_region = -left_cone & -central_cylinder_left_plane
    right_cone_region = -right_cone & +central_cylinder_right_plane
    
    # --- Full vessel region ---
    vessel_region = left_cone_region | right_cone_region | outer_cylinders_region | central_cylinder
    
    return vessel_region

def make_tandem_vacuum_vessel_regions(end_plug_vv_parameters, central_cell_vv_parameters, central_cell_end_plug_separation_distance, bottleneck_radius, end_axial_distance):

    end_plug_vv_central_radius = end_plug_vv_parameters.radius
    end_plug_vv_central_axis_length = end_plug_vv_parameters.central_axis_length
    end_plug_vv_outer_length = end_plug_vv_parameters.outer_length
    
    central_cell_central_radius = central_cell_vv_parameters.radius
    central_cell_vv_central_axis_length = central_cell_vv_parameters.central_axis_length
    central_cell_vv_outer_length = central_cell_vv_parameters.outer_length

    right_end_plug_axial_midplane = central_cell_vv_central_axis_length/2 + central_cell_end_plug_separation_distance + end_plug_vv_central_axis_length/2
    left_end_plug_axial_midplane = -central_cell_vv_central_axis_length/2 - central_cell_end_plug_separation_distance - end_plug_vv_central_axis_length/2

    central_cell_vv_region = redefined_vacuum_vessel_region(central_cell_vv_outer_length, central_cell_vv_central_axis_length, central_cell_central_radius, bottleneck_radius, central_cell_end_plug_separation_distance/2, central_cell_end_plug_separation_distance/2, axial_midplane = 0.0)

    right_end_plug_vv_region = redefined_vacuum_vessel_region(end_plug_vv_outer_length, end_plug_vv_central_axis_length, end_plug_vv_central_radius, bottleneck_radius, central_cell_end_plug_separation_distance/2, end_axial_distance, axial_midplane = right_end_plug_axial_midplane)
    left_end_plug_vv_region = redefined_vacuum_vessel_region(end_plug_vv_outer_length, end_plug_vv_central_axis_length, end_plug_vv_central_radius, bottleneck_radius, end_axial_distance, central_cell_end_plug_separation_distance/2, axial_midplane = left_end_plug_axial_midplane)


    return central_cell_vv_region, [left_end_plug_vv_region, right_end_plug_vv_region]