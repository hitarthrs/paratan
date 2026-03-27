"""
Neutron source models for simple mirror machines.

VolumetricSource is the primary source for a simple mirror — it approximates
the plasma fusion neutron distribution as a central cylinder plus stepped
conical sections matching the VV geometry.

All other source types (Uniform, 1D, 2D) are also available here.
"""

import numpy as np
import openmc
import yaml


# ---------------------------------------------------------------------------
# Power → neutron rate conversion
# ---------------------------------------------------------------------------

def _power_to_strength(power_mw: float, energy_mev: float = 14.1) -> float:
    """Convert fusion power (MW) to neutron emission rate (n/s)."""
    return power_mw * 0.8 * 1e6 / (energy_mev * 1e6 * 1.602e-19)


# ---------------------------------------------------------------------------
# Core volumetric approximation function
# ---------------------------------------------------------------------------

def volumetric_source_approximation(
    vacuum_vessel_axial_length: float,
    vacuum_vessel_outer_axial_length: float,
    vacuum_vessel_central_radius: float,
    throat_radius: float,
    z_origin: float = 0.0,
    conical_sources: int = 5,
) -> list[openmc.IndependentSource]:
    """
    Approximate a volumetric plasma source as:
      - one central cylindrical source (uniform in z, power-law in r)
      - N stepped conical sections on each side of the central cylinder

    Source strengths are volume-weighted so they sum to 1.0 (relative).
    Scale the returned source list using `source.strength` after calling this.

    Parameters
    ----------
    vacuum_vessel_axial_length : float
        Total axial length of the VV central section (cm).
    vacuum_vessel_outer_axial_length : float
        Axial length of the outer (bottleneck) cylinder section (cm).
    vacuum_vessel_central_radius : float
        Radius of the central VV cylinder (cm).
    throat_radius : float
        Radius at the VV bottleneck / mirror throat (cm).
    z_origin : float
        Z-coordinate of the machine midplane (cm).
    conical_sources : int
        Number of stepped conical slices per side.

    Returns
    -------
    list[openmc.IndependentSource]
    """
    phi_dist = openmc.stats.Uniform(0, 2 * np.pi)
    energy_dist = openmc.stats.Discrete([14.1e6], [1])

    # --- Central cylindrical source ---
    z_dist_central = openmc.stats.Uniform(
        -(vacuum_vessel_outer_axial_length / 2) + z_origin,
        (vacuum_vessel_outer_axial_length / 2) + z_origin,
    )
    r_dist_central = openmc.stats.PowerLaw(0.001, vacuum_vessel_central_radius, -2)

    # --- Conical sections ---
    angle_cone = np.arctan(
        2 * (vacuum_vessel_central_radius - throat_radius)
        / (vacuum_vessel_axial_length - vacuum_vessel_outer_axial_length)
    )
    step_z = (vacuum_vessel_axial_length - vacuum_vessel_outer_axial_length) / (2 * conical_sources)
    step_r = np.tan(angle_cone) * step_z

    z_dists_conical = []
    r_dists_conical = []
    volumes_conical = []

    r_current = vacuum_vessel_central_radius
    for i in range(conical_sources):
        z_left_start = z_origin - vacuum_vessel_outer_axial_length / 2 - (i + 1) * step_z
        z_left_end = z_left_start + step_z
        z_right_start = z_origin + vacuum_vessel_outer_axial_length / 2 + i * step_z
        z_right_end = z_right_start + step_z

        z_mix = openmc.stats.Mixture(
            [0.5, 0.5],
            [openmc.stats.Uniform(z_left_start, z_left_end),
             openmc.stats.Uniform(z_right_start, z_right_end)],
        )
        z_dists_conical.append(z_mix)

        r_slice = r_current - (i + 1) * step_r
        r_dists_conical.append(openmc.stats.PowerLaw(0.001, r_slice, -2))
        volumes_conical.append(2 * np.pi * (r_slice ** 2) * step_z)

    vol_central = np.pi * (vacuum_vessel_central_radius ** 2) * vacuum_vessel_outer_axial_length
    total_volume = vol_central + np.sum(volumes_conical)

    # --- Build source objects ---
    sources = []

    central = openmc.IndependentSource()
    central.particle = "neutron"
    central.angle = openmc.stats.Isotropic()
    central.energy = energy_dist
    central.space = openmc.stats.CylindricalIndependent(r_dist_central, phi_dist, z_dist_central)
    central.strength = vol_central / total_volume
    sources.append(central)

    for i in range(conical_sources):
        src = openmc.IndependentSource()
        src.particle = "neutron"
        src.angle = openmc.stats.Isotropic()
        src.energy = energy_dist
        src.space = openmc.stats.CylindricalIndependent(
            r_dists_conical[i], phi_dist, z_dists_conical[i]
        )
        src.strength = volumes_conical[i] / total_volume
        sources.append(src)

    return sources


# ---------------------------------------------------------------------------
# Source classes
# ---------------------------------------------------------------------------

class Source:
    """Base class — converts MW power to absolute neutron emission rate."""

    def __init__(self, power_output: float, energy_mev: float = 14.1):
        self.strength = _power_to_strength(power_output, energy_mev)

    def describe(self):
        raise NotImplementedError


class UniformSource(Source):
    def __init__(self, power_output: float, length: float, radius: float,
                 z_origin: float = 0.0, energy: float = 14.1e6):
        super().__init__(power_output)
        self.length = length
        self.radius = radius
        self.z_origin = z_origin
        self.energy = energy

    def create_openmc_source(self) -> openmc.IndependentSource:
        source = openmc.IndependentSource()
        source.particle = "neutron"
        source.angle = openmc.stats.Isotropic()
        source.energy = openmc.stats.Discrete([self.energy], [1])
        source.space = openmc.stats.CylindricalIndependent(
            openmc.stats.PowerLaw(0.001, self.radius, -2),
            openmc.stats.Uniform(0, 2 * np.pi),
            openmc.stats.Uniform(-self.length / 2 + self.z_origin,
                                  self.length / 2 + self.z_origin),
        )
        return source

    def describe(self):
        print(f"UniformSource: strength={self.strength:.3e} n/s, "
              f"length={self.length} cm, radius={self.radius} cm")


class VolumetricSource(Source):
    def __init__(self, power_output: float, vacuum_vessel_axial_length: float,
                 vacuum_vessel_outer_axial_length: float, vacuum_vessel_central_radius: float,
                 throat_radius: float, z_origin: float = 0.0, conical_sources: int = 5):
        super().__init__(power_output)
        self.vacuum_vessel_axial_length = vacuum_vessel_axial_length
        self.vacuum_vessel_outer_axial_length = vacuum_vessel_outer_axial_length
        self.vacuum_vessel_central_radius = vacuum_vessel_central_radius
        self.throat_radius = throat_radius
        self.z_origin = z_origin
        self.conical_sources = conical_sources

    def create_openmc_source(self) -> list[openmc.IndependentSource]:
        """Returns a list of normalized sources (strengths sum to 1.0)."""
        return volumetric_source_approximation(
            self.vacuum_vessel_axial_length,
            self.vacuum_vessel_outer_axial_length,
            self.vacuum_vessel_central_radius,
            self.throat_radius,
            self.z_origin,
            self.conical_sources,
        )

    def describe(self):
        print(f"VolumetricSource: strength={self.strength:.3e} n/s, "
              f"central_radius={self.vacuum_vessel_central_radius} cm, "
              f"throat_radius={self.throat_radius} cm")


# ---------------------------------------------------------------------------
# YAML loader
# ---------------------------------------------------------------------------

def load_source_from_yaml(yaml_file: str):
    """
    Build an OpenMC source (or list of sources) from a source_information.yaml file.

    Returns a single source or list depending on source type.
    """
    with open(yaml_file) as f:
        data = yaml.safe_load(f)

    src_data = data["source"]
    source_type = src_data["type"]
    power_output = src_data["power_output"]

    if source_type == "Uniform":
        length = src_data["uniform"]["length"]
        radius = src_data["uniform"]["radius"]
        z_origin = src_data["uniform"].get("z_origin", 0.0)
        return UniformSource(power_output, length, radius, z_origin).create_openmc_source()

    elif source_type == "Volumetric":
        v = src_data["volumetric"]
        return VolumetricSource(
            power_output,
            v["central_cell_axial_length"],
            v["central_cell_outer_axial_length"],
            v["central_cell_radius"],
            v["throat_radius"],
            v.get("z_origin", 0.0),
            v.get("conical_sources", 5),
        ).create_openmc_source()

    else:
        raise ValueError(f"Unsupported source type: {source_type!r}")
