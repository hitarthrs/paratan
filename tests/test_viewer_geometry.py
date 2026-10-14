"""Regression checks for viewer dimensions, materials, and layer volumes."""
from pathlib import Path

import numpy as np
import math
import unittest

from src.paratan.viewer.simple_mirror_meshes import (
    build_simple_mirror_components,
    load_simple_mirror_yaml,
)


def model_input():
    path = Path(__file__).resolve().parents[1] / "input_files/simple_parametric_input_new.yaml"
    return load_simple_mirror_yaml(path)


def check_ports_have_core_correct_materials_and_x_plane_end_faces(model_input):
    model_input["ports"] = {
        "ports_present": True, "inner_radius": 25, "port_x_limits": 250,
        "port_layers_thicknesses": [2, 5],
        "port_layers_materials": ["water_cooled_wc", "Magnet_Winding_Pack_2"],
    }
    parts = build_simple_mirror_components(model_input, n_theta=72)
    for side in ("left", "right"):
        ports = [c for c in parts if c.name.startswith(f"port_{side}_")]
        assert [c.material for c in ports] == ["vacuum", "water_cooled_wc", "Magnet_Winding_Pack_2"]
        for c, ri, ro in zip(ports, [0, 25, 27], [25, 27, 32]):
            assert np.allclose(c.mesh.bounds[:2], [-250, 250])
            assert c.mesh.n_open_edges == 0
            expected = np.pi * (ro**2 - ri**2) * 500 / np.sin(np.pi / 4)
            # Independent cross-section integration of the part excluded by
            # the central plasma cylinder. Ports now contain two outer arms.
            nodes, weights = np.polynomial.legendre.leggauss(32)
            radii = (nodes+1)*(ro-ri)/2+ri
            angles = np.arange(256)*2*np.pi/256
            lengths = 2*np.sqrt(model_input['vacuum_vessel']['central_radius']**2-
                                (radii[:,None]*np.sin(angles))**2)/np.sin(np.pi/4)
            removed = np.sum(weights*radii*lengths.mean(axis=1))*(ro-ri)/2*2*np.pi
            expected -= removed
            assert c.mesh.volume == Approx(expected, rel=0.01)


def check_hf_nested_layers_match_openmc_dimensions_and_volumes(model_input):
    parts = {c.name: c for c in build_simple_mirror_components(model_input, n_theta=72)}
    hf = model_input["hf_coil"]
    casing = [layer["thickness"] for layer in hf["casing_layers"]]
    front = casing + [hf["shield"]["radial_thickness"][0]]
    back = casing + [hf["shield"]["radial_thickness"][1]]
    axial = casing + [hf["shield"]["axial_thickness"][0]]
    vv = model_input["vacuum_vessel"]
    base = vv["bottleneck_radius"] + sum(v["thickness"] for v in vv["structure"].values())
    ri = base + hf["shield"]["radial_gap_before_casing"] + sum(front)
    ro = ri + hf["magnet"]["radial_thickness"]
    height = hf["magnet"]["axial_thickness"]
    for side in ("left", "right"):
        magnet = parts[f"hf_coil_{side}_magnet"]
        assert magnet.meta["bore_start_r_cm"] == Approx(ri)
        inner, outer, length = ri, ro, height
        labels = [f"casing_{j}" for j in range(len(casing))] + ["shield"]
        for j, label in enumerate(labels):
            old_volume = np.pi * (outer**2 - inner**2) * length
            inner -= front[j]
            outer += back[j]
            length += 2 * axial[j]
            expected = np.pi * (outer**2 - inner**2) * length - old_volume
            mesh = parts[f"hf_coil_{side}_{label}"].mesh
            assert mesh.bounds[1] == Approx(outer)
            assert mesh.bounds[5] - mesh.bounds[4] == Approx(length)
            assert mesh.volume == Approx(expected, rel=0.01)


def _perpendicular_vessel_volume(vv, thickness=0.0):
    """Closed-form volume matching simple_mirror_volumes._profile('perpendicular')."""
    mid = vv.get("axial_midplane", 0.0)
    half = vv["central_axial_length"] / 2.0
    left = mid - half - vv["left_bottleneck_length"]
    right = mid + half + vv["right_bottleneck_length"]
    h = half + thickness
    rc = vv["central_radius"] + thickness
    rb = vv["bottleneck_radius"] + thickness
    return math.pi * (
        rb**2 * ((mid - h) - left) + rc**2 * (2 * h) + rb**2 * (right - (mid + h))
    )


def check_perpendicular_vessel_and_breeder_volumes(data):
    data = dict(data)
    data["vacuum_vessel"] = dict(data["vacuum_vessel"])
    data["vacuum_vessel"]["geometry_style"] = "perpendicular"
    parts = {c.name: c for c in build_simple_mirror_components(data, n_theta=72)}
    vv = data["vacuum_vessel"]
    plasma = parts["vacuum_vessel_plasma"]
    assert plasma.meta["geometry_style"] == "perpendicular"
    assert plasma.mesh.volume == Approx(_perpendicular_vessel_volume(vv), rel=0.01)

    previous = _perpendicular_vessel_volume(vv)
    thickness = 0.0
    for name, layer in vv["structure"].items():
        thickness += layer["thickness"]
        total = _perpendicular_vessel_volume(vv, thickness)
        assert parts[f"vv_{name}"].mesh.volume == Approx(total - previous, rel=0.01)
        previous = total

    # Breeder is a CC cylinder minus the outer stepped vessel, within CC length.
    cc = data["central_cell"]
    wall = sum(layer["thickness"] for layer in vv["structure"].values())
    rc = vv["central_radius"] + wall
    end = cc["axial_length"] / 2.0
    # With default inputs the CC lies inside the thickened central cylinder.
    assert end <= vv["central_axial_length"] / 2.0 + wall
    ro = rc + cc["layers"][0]["thickness"]
    expected = math.pi * (ro**2 - rc**2) * 2 * end
    breeder = parts[f"cc_layer_0_{cc['layers'][0]['material']}"]
    assert breeder.mesh.volume == Approx(expected, rel=0.01)


class Approx:
    def __init__(self, expected, rel=1e-6):
        self.expected, self.rel = expected, rel

    def __eq__(self, actual):
        return math.isclose(actual, self.expected, rel_tol=self.rel, abs_tol=1e-12)


class ViewerGeometryTests(unittest.TestCase):
    def test_ports(self):
        check_ports_have_core_correct_materials_and_x_plane_end_faces(model_input())

    def test_hf_layers(self):
        check_hf_nested_layers_match_openmc_dimensions_and_volumes(model_input())

    def test_perpendicular_vacuum_vessel(self):
        check_perpendicular_vessel_and_breeder_volumes(model_input())

    def test_unknown_geometry_style_is_rejected(self):
        data = model_input()
        data["vacuum_vessel"]["geometry_style"] = "axisymmetric"
        with self.assertRaises(ValueError):
            build_simple_mirror_components(data, n_theta=24)


if __name__ == "__main__":
    unittest.main()
