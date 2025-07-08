
import openmc
import pytest
from types import SimpleNamespace
from function_libraries.geometry_lib import *
from function_libraries.tallies_lib import *

import openmc

import re

def normalize_structure(region_str):
    """Normalize a region string: (+X -X ...) while preserving parens and ops"""
    region_str = re.sub(r'(?<![\d\-\+])(\d+)', r'+\1', region_str)  # add '+' to bare numbers
    return re.sub(r'[\+\-]\d+', lambda m: m.group()[0] + 'X', region_str)

def extract_surface_ids(region_str):
    """Extract signed IDs like '+3', '-4', etc., even if '+' was missing originally."""
    region_str = re.sub(r'(?<![\d\-\+])(\d+)', r'+\1', region_str)
    return re.findall(r'[\+\-]\d+', region_str)

def compare_region_structures_and_surfaces(r1, r2):
    s1 = str(r1)
    s2 = str(r2)

    struct1 = normalize_structure(s1)
    struct2 = normalize_structure(s2)

    if struct1 != struct2:
        print("Region structure mismatch:")
        print(f"  expected: {struct2}")
        print(f"  actual  : {struct1}")
        return False

    ids1 = extract_surface_ids(s1)
    ids2 = extract_surface_ids(s2)

    if len(ids1) != len(ids2):
        print("Surface ID count mismatch")
        return False

    # Surface lookup by unsigned ID only
    surfaces1 = r1.get_surfaces()
    surfaces2 = r2.get_surfaces()

    for sid1, sid2 in zip(ids1, ids2):
        id1 = int(sid1[1:])  # remove sign
        id2 = int(sid2[1:])

        s1 = surfaces1[id1]
        s2 = surfaces2[id2]

        if not s1.is_equal(s2):
            print(f"Surfaces {sid1} and {sid2} not geometrically equal")
            return False

    return True



def test_hollow_cylindrical_region_structure():
    z0 = 0.0
    r = 10.0
    dr = 2.0
    dz = 20.0

    region = hollow_cylindrical_region(z0, r, dr, dz)

    assert isinstance(region, openmc.Region), "Should return an OpenMC Region object"

    # Check that each component surface is used
    left = openmc.ZPlane(z0 - dz/2)
    right = openmc.ZPlane(z0 + dz/2)
    inner = openmc.ZCylinder(r=r)
    outer = openmc.ZCylinder(r=r + dr)

    expected_region = +left & -right & (-outer & +inner)

    # Compare structure
    assert compare_region_structures_and_surfaces(region, expected_region)

def test_hollow_cylinder_with_shell():
    z0 = 0.0
    r0 = 10.0
    dr_inner = 2.0
    dz_inner = 20.0
    shell_front = 1.0
    shell_back = 1.5
    dz_shell = 5.0

    # Get actual output
    shell_region, inner_region = hollow_cylinder_with_shell(
        z0, r0, dr_inner, dz_inner, shell_front, shell_back, dz_shell
    )

    # Define surfaces manually
    zl_inner = openmc.ZPlane(z0 - dz_inner / 2)
    zl_outer = openmc.ZPlane(z0 - dz_inner / 2 - dz_shell)
    zr_inner = openmc.ZPlane(z0 + dz_inner / 2)
    zr_outer = openmc.ZPlane(z0 + dz_inner / 2 + dz_shell)

    r_front = openmc.ZCylinder(r=r0)
    r_inner = openmc.ZCylinder(r=r0 + shell_front)
    r_outer = openmc.ZCylinder(r=r0 + shell_front + dr_inner)
    r_back = openmc.ZCylinder(r=r0 + shell_front + dr_inner + shell_back)

    expected_inner = (+r_inner & -r_outer) & -zr_inner & +zl_inner
    expected_outer = (+r_front & -r_back) & -zr_outer & +zl_outer
    expected_shell = expected_outer & ~expected_inner

    assert compare_region_structures_and_surfaces(shell_region, expected_shell)
    assert compare_region_structures_and_surfaces(inner_region, expected_inner)

def test_cylinder_with_shell():
    z0 = 0.0
    r0 = 10.0
    dz_inner = 20.0
    shell_thickness = 5.0

    shell_region, inner_region = cylinder_with_shell(z0, r0, dz_inner, shell_thickness)

    zl_inner = openmc.ZPlane(z0 - dz_inner / 2)
    zl_outer = openmc.ZPlane(z0 - dz_inner / 2 - shell_thickness)
    zr_inner = openmc.ZPlane(z0 + dz_inner / 2)
    zr_outer = openmc.ZPlane(z0 + dz_inner / 2 + shell_thickness)

    r_inner = openmc.ZCylinder(r=r0)
    r_outer = openmc.ZCylinder(r=r0 + shell_thickness)

    expected_inner = +zl_inner & -zr_inner & -r_inner
    expected_outer = +zl_outer & -zr_outer & -r_outer
    expected_shell = expected_outer & ~expected_inner

    assert compare_region_structures_and_surfaces(shell_region, expected_shell)
    assert compare_region_structures_and_surfaces(inner_region, expected_inner)

def test_cylindrical_region_no_outer_surface():
    z0 = 0.0
    r0 = 10.0
    dr_inner = 2.0
    dz_inner = 20.0
    shell_thickness = 5.0

    shell_region_actual = cylindrical_region_no_outer_surface(
        z0, r0, dr_inner, dz_inner, shell_thickness
    )

    # Manual surfaces
    zl_inner = openmc.ZPlane(z0 - dz_inner / 2)
    zl_outer = openmc.ZPlane(z0 - dz_inner / 2 - shell_thickness)
    zr_inner = openmc.ZPlane(z0 + dz_inner / 2)
    zr_outer = openmc.ZPlane(z0 + dz_inner / 2 + shell_thickness)

    r_inner = openmc.ZCylinder(r=r0)
    r_inner_plus_shell = openmc.ZCylinder(r=r0 + shell_thickness)
    r_outer = openmc.ZCylinder(r=r0 + shell_thickness + dr_inner)

    # Rebuild expected region logic
    inner_region = (+r_inner_plus_shell & -r_outer) & -zr_inner & +zl_inner
    outer_region = (+r_inner & -r_outer) & -zr_outer & +zl_outer
    expected_shell = outer_region & ~inner_region

    # Compare structure and surface geometry
    assert compare_region_structures_and_surfaces(shell_region_actual, expected_shell)

def test_redefined_vacuum_vessel_region():
    # ---- PARAMETERS ----
    outer_axial = 10.0
    central_axial = 20.0
    central_r = 15.0
    bottleneck_r = 10.0
    left_bottleneck_len = 5.0
    right_bottleneck_len = 5.0
    z0 = 0.0

    # ---- ACTUAL OUTPUT ----
    region_actual, region_actual_components = redefined_vacuum_vessel_region(
        outer_axial, central_axial, central_r, bottleneck_r,
        left_bottleneck_len, right_bottleneck_len, axial_midplane=z0
    )

    # ---- EXPECTED REGION ----
    # Build all planes
    first_plane = outer_axial / 2
    second_plane = central_axial / 2
    rightmost = second_plane + right_bottleneck_len
    leftmost = -second_plane - left_bottleneck_len

    z_left = openmc.ZPlane(z0 - first_plane)
    z_right = openmc.ZPlane(z0 + first_plane)
    cyl = openmc.ZCylinder(r=central_r)
    central = -cyl & +z_left & -z_right

    # Outer sections
    z_l1 = openmc.ZPlane(z0 + leftmost)
    z_l2 = openmc.ZPlane(z0 - second_plane)
    z_r1 = openmc.ZPlane(z0 + second_plane)
    z_r2 = openmc.ZPlane(z0 + rightmost)
    outer_cyl = openmc.ZCylinder(r=bottleneck_r)

    left_outer = -outer_cyl & (+z_l1 & -z_l2)
    right_outer = -outer_cyl & (+z_r1 & -z_r2)

    outer = left_outer | right_outer

    # Cones
    angle = np.arctan(2 * (central_r - bottleneck_r) / (central_axial - outer_axial))
    z_cone_left = z0 - (central_r / np.tan(angle) + first_plane)
    z_cone_right = z0 + (central_r / np.tan(angle) + first_plane)

    cone_l = openmc.model.ZConeOneSided(0.0, 0.0, z_cone_left, r2=np.tan(angle)**2)
    cone_r = openmc.model.ZConeOneSided(0.0, 0.0, z_cone_right, r2=np.tan(angle)**2, up=False)

    cone_l_region = -cone_l & -z_left & +z_l1
    cone_r_region = -cone_r & +z_right & -z_r2

    # Union full
    region_expected = cone_l_region | cone_r_region | outer | central

    assert compare_region_structures_and_surfaces(region_actual, region_expected)


def test_hollow_mesh_from_domain():
    z0 = 0.0
    r0 = 10.0
    dr = 2.0
    dz = 20.0

    # Create a basic hollow cylindrical region
    region = openmc.Cell(region=hollow_cylindrical_region(z0, r0, dr, dz))

    mesh = hollow_mesh_from_domain(region, dimensions=(4, 8, 6))

    assert isinstance(mesh, openmc.CylindricalMesh)

    # Bounding box center for origin
    bbox = region.bounding_box
    center_x = (bbox[0][0] + bbox[1][0]) / 2
    center_y = (bbox[0][1] + bbox[1][1]) / 2
    center_z = (bbox[0][2] + bbox[1][2]) / 2
    expected_origin = (center_x, center_y, bbox[0][2])

    # Check origin
    assert np.allclose(mesh.origin, expected_origin)

    # Check grid lengths
    assert len(mesh.r_grid) == 5  # 4 + 1
    assert len(mesh.phi_grid) == 9
    assert len(mesh.z_grid) == 7

    # Check that r_grid starts at expected inner radius
    # Should match r0 from original inner cylinder
    assert np.isclose(mesh.r_grid[0], r0)

    # Check that r_grid ends at max bounding cylinder radius
    max_radius = max(bbox[0][0], bbox[0][1], bbox[1][0], bbox[1][1])
    assert np.isclose(mesh.r_grid[-1], max_radius)

def test_strings_to_openmc_filters():
    filters = strings_to_openmc_filters([
        "neutron_filter",
        "photon_filter",
        "fast_energies_filter",
        "thermal_energies_filter",
        "vitaminj_filter"
    ])

    assert len(filters) == 5

    # Check types
    assert isinstance(filters[0], openmc.ParticleFilter)
    assert filters[0] == openmc.ParticleFilter("neutron")

    assert isinstance(filters[1], openmc.ParticleFilter)
    assert filters[1] == openmc.ParticleFilter("photon")

    assert isinstance(filters[2], openmc.EnergyFilter)
    assert np.isclose(filters[2].values[0], 1e5, rtol=1e-2)
    assert np.isclose(filters[2].values[-1], 20e6, rtol=1e-2)
    assert len(filters[2].values) == 501

    assert isinstance(filters[3], openmc.EnergyFilter)
    assert np.isclose(filters[3].values[0], 1e-3, rtol=1e-2)
    assert np.isclose(filters[3].values[-1], 1e3, rtol=1e-2)
    assert len(filters[3].values) == 501

