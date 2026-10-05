"""D and T seed or reference profiles on an axial grid

Normalize a shared shape at the magnetic midplane and integrate its inventory
The integration layer supplies the magnetic startup shape for stationary NBI runs
"""
from __future__ import annotations
from dataclasses import dataclass
from typing import Any
import numpy as np
from source_model_revamp.geometry.device_domains import AxialDomain

# The integration layer applies this role label after constructing the startup table
MAGNETICALLY_MAPPED_STARTUP_PROFILE_MODEL = "magnetically_mapped_loss_cone_depleted_maxwellian"

@dataclass(frozen=True)
class BackgroundProfileSpecification:
    """Shared axial shape controls parsed from configuration

    Polynomial coefficients run from constant term to highest power
    Tabulated coordinates use machine positions in meters
    """
    model: str = "uniform_central_cell"
    analytic_polynomial_coefficients: tuple[float, ...] = ()
    tabulated_z_m: tuple[float, ...] = ()
    tabulated_scale: tuple[float, ...] = ()

@dataclass(frozen=True)
class IonSeedReferenceProfiles:
    """Geometry owned seed or reference ion profiles

    Profile arrays share the z_m grid and inventories contain particle counts
    The shape is normalized at the exact magnetic midplane
    model records the profile role while specification defines its evaluation
    """
    z_m: np.ndarray
    shape: np.ndarray
    deuterium_density_m3: np.ndarray
    tritium_density_m3: np.ndarray
    ion_temperature_keV: np.ndarray
    support_domain: AxialDomain
    specification: BackgroundProfileSpecification
    model: str
    magnetic_midplane_z_m: float
    deuterium_midplane_density_m3: float
    tritium_midplane_density_m3: float
    shape_midplane_value_before_normalization: float
    symmetry_relative_error: float
    symmetry_tolerance: float
    inventories: dict[str, float]

    def normalized_shape_at(self, z_m: Any) -> float | np.ndarray:
        """Evaluate the exact normalized shared shape in machine coordinates"""
        coordinates = np.asarray(z_m, dtype=float)
        if np.any(~np.isfinite(coordinates)):
            raise ValueError("background profile evaluation coordinates must be finite")
        values = (_raw_shape_values(coordinates, self.support_domain, self.specification) / self.shape_midplane_value_before_normalization)
        if np.ndim(z_m) == 0:
            return float(values)

        return np.asarray(values, dtype=float)

    def deuterium_density_m3_at(self, z_m: Any) -> float | np.ndarray:
        """Evaluate the seed or reference deuterium density in machine coordinates"""
        return self.deuterium_midplane_density_m3 * self.normalized_shape_at(z_m)

    def tritium_density_m3_at(self, z_m: Any) -> float | np.ndarray:
        """Evaluate the seed or reference tritium density in machine coordinates"""
        return self.tritium_midplane_density_m3 * self.normalized_shape_at(z_m)

    def positive_charge_density_m3_at(self, z_m: Any) -> float | np.ndarray:
        """Return the sum of singly charged D and T number densities at given positions"""
        return self.positive_charge_midplane_density_m3 * self.normalized_shape_at(z_m)

    @property
    def positive_charge_density_m3(self) -> np.ndarray:
        """Return D plus T number density on the stored grid"""
        return self.deuterium_density_m3 + self.tritium_density_m3

    @property
    def positive_charge_midplane_density_m3(self) -> float:
        """Return D plus T number density at the exact magnetic midplane"""
        return self.deuterium_midplane_density_m3 + self.tritium_midplane_density_m3

    def as_metadata(self) -> dict[str, Any]:
        """Return profile arrays, normalization checks, and inventories for run metadata"""
        return {
            "central_cell_plasma_region_source": self.support_domain.source,
            "central_cell_plasma_region_bounds": [self.support_domain.z_min_m, self.support_domain.z_max_m],
            "background_profile_model": self.model,
            "background_profile_magnetic_midplane_z_m": self.magnetic_midplane_z_m,
            "background_profile_shape": self.shape,
            "background_profile_shape_midplane_value": 1.0,
            "background_profile_shape_midplane_value_before_normalization": self.shape_midplane_value_before_normalization,
            "background_profile_midplane_normalized": True,
            "background_profile_finite": True,
            "background_profile_nonnegative": True,
            "background_profile_symmetry_relative_error": self.symmetry_relative_error,
            "background_profile_symmetry_tolerance": self.symmetry_tolerance,
            "background_profile_symmetric_within_tolerance": True,
            "deuterium_background_density_profile_m3": self.deuterium_density_m3,
            "tritium_background_density_profile_m3": self.tritium_density_m3,
            "background_positive_charge_density_profile_m3": self.positive_charge_density_m3,
            "background_deuterium_midplane_density_m3": self.deuterium_midplane_density_m3,
            "background_tritium_midplane_density_m3": self.tritium_midplane_density_m3,
            "background_positive_charge_midplane_density_m3": self.positive_charge_midplane_density_m3,
            "background_temperature_profiles": {"ion_temperature_keV": self.ion_temperature_keV},
            "background_inventory_by_species": dict(self.inventories),
            "background_zero_outside_support": True,
            "background_population_scope": "startup_seed_initialization_only" if self.model == MAGNETICALLY_MAPPED_STARTUP_PROFILE_MODEL else "prescribed_reference_profile_only",
            "background_electron_density_profile_available": False,
        }

def _raw_shape_values(z_m: np.ndarray, support: AxialDomain, specification: BackgroundProfileSpecification) -> np.ndarray:
    """Evaluate the selected shape and set values outside its support to zero"""
    z = np.asarray(z_m, dtype=float)
    inside = support.contains(z)
    result = np.zeros_like(z, dtype=float)
   
    if specification.model == "uniform_central_cell":
        result[inside] = 1.0
  
    elif specification.model == "analytic_axial_profile":
        coefficients = np.asarray(specification.analytic_polynomial_coefficients, dtype=float)
        if coefficients.size == 0 or np.any(~np.isfinite(coefficients)):
            raise ValueError("analytic_axial_profile requires finite polynomial_coefficients")
       
        # Map the support endpoints to negative one and positive one for the polynomial
        coordinate = 2.0 * (z[inside] - support.midpoint_m) / support.length_m
        result[inside] = np.polynomial.polynomial.polyval(coordinate, coefficients)
   
    elif specification.model == "tabulated_axial_profile":
        table_z = np.asarray(specification.tabulated_z_m, dtype=float)
        table_scale = np.asarray(specification.tabulated_scale, dtype=float)
        if (table_z.size < 2 or table_z.size != table_scale.size or np.any(~np.isfinite(table_z)) or np.any(~np.isfinite(table_scale)) or np.any(table_scale < 0.0) or np.any(np.diff(table_z) <= 0.0)):
            raise ValueError("tabulated_axial_profile requires equal length, finite, increasing z_m and nonnegative scale arrays")
        result[inside] = np.interp(z[inside], table_z, table_scale, left=0.0, right=0.0)
    else:
        raise ValueError(f"unsupported background profile model {specification.model!r}")
    result[~inside] = 0.0

    return result

def _analytic_extrema_coordinates(support: AxialDomain, specification: BackgroundProfileSpecification) -> np.ndarray:
    """Return support endpoints and interior polynomial extrema in machine coordinates"""
    if specification.model != "analytic_axial_profile":
        return np.empty(0, dtype=float)
   
    coefficients = np.asarray(specification.analytic_polynomial_coefficients, dtype=float)
    derivative = np.polynomial.polynomial.polyder(coefficients)
   
    if derivative.size <= 1 and np.all(derivative == 0.0):
        roots = np.empty(0, dtype=float)
    else:
        complex_roots = np.polynomial.polynomial.polyroots(derivative)
        roots = np.real(complex_roots[np.abs(np.imag(complex_roots)) <= 1.0e-12])
        roots = roots[(roots > -1.0) & (roots < 1.0)]
    coordinate = np.concatenate((np.asarray([-1.0, 1.0]), roots))

    return support.midpoint_m + 0.5 * support.length_m * coordinate

def _normalized_shape(z_m: np.ndarray, support: AxialDomain, specification: BackgroundProfileSpecification, *, magnetic_midplane_z_m: float, symmetry_tolerance: float) -> tuple[np.ndarray, float, float]:
    """Normalize at the exact magnetic midplane and check nonnegativity and symmetry

    Return the sampled shape, its original midplane value, and the symmetry error
    """
    midpoint = float(magnetic_midplane_z_m)
    tolerance = float(symmetry_tolerance)
    if not np.isfinite(midpoint):
        raise ValueError("magnetic_midplane_z_m must be finite")
    if not bool(support.contains(np.asarray([midpoint]))[0]):
        raise ValueError("background profile support must contain the magnetic midplane")
    if not np.isfinite(tolerance) or tolerance < 0.0:
        raise ValueError("background profile symmetry tolerance must be finite and nonnegative")
    if specification.model == "tabulated_axial_profile":
        table_z = np.asarray(specification.tabulated_z_m, dtype=float)
        if table_z.size < 2 or midpoint < table_z[0] or midpoint > table_z[-1]:
            raise ValueError("tabulated_axial_profile must cover the magnetic midplane")
    midpoint_value = float(_raw_shape_values(np.asarray([midpoint]), support, specification)[0])
    if not np.isfinite(midpoint_value) or midpoint_value <= 0.0:
        raise ValueError("background profile shape must be finite and positive at the magnetic midplane")
    # Include polynomial extrema to catch negative values between sampled coordinates
    validation_z = np.concatenate((np.linspace(support.z_min_m, support.z_max_m, 4097), np.asarray(z_m, dtype=float), np.asarray([midpoint]), _analytic_extrema_coordinates(support, specification)))
    raw_validation = _raw_shape_values(validation_z, support, specification)
    if np.any(~np.isfinite(raw_validation)) or np.any(raw_validation < 0.0):
        raise ValueError("background profile shape must be finite and nonnegative")

    radius = max(midpoint - support.z_min_m, support.z_max_m - midpoint)
    # Compare the zero extended shape on both sides of the magnetic midplane
    symmetry_z = midpoint + np.linspace(-radius, radius, 4097)
    if specification.model == "tabulated_axial_profile":
        table_z = np.asarray(specification.tabulated_z_m, dtype=float)
        symmetry_z = np.concatenate((symmetry_z, table_z, 2.0 * midpoint - table_z))
    direct = _raw_shape_values(symmetry_z, support, specification) / midpoint_value
    mirrored = _raw_shape_values(2.0 * midpoint - symmetry_z, support, specification) / midpoint_value
    # This average sets the error scale without changing the returned profile
    symmetric = 0.5 * (direct + mirrored)
    scale = max(float(np.max(np.abs(symmetric))), np.finfo(float).tiny)
    symmetry_error = float(np.max(np.abs(direct - mirrored)) / scale)
    if symmetry_error > tolerance:
        raise ValueError("background profile is not symmetric about the magnetic midplane within " f"the configured tolerance: relative_error={symmetry_error:.6e}, tolerance={tolerance:.6e}")
    normalized = _raw_shape_values(np.asarray(z_m, dtype=float), support, specification) / midpoint_value
    if np.any(~np.isfinite(normalized)) or np.any(normalized < 0.0):
        raise ValueError("normalized background profile shape must be finite and nonnegative")
    normalized[~support.contains(z_m)] = 0.0

    return normalized, midpoint_value, symmetry_error

def build_ion_seed_reference_profiles(z_m: Any, cell_volumes_m3: Any, support: AxialDomain, specification: BackgroundProfileSpecification, *, magnetic_midplane_z_m: float, deuterium_midplane_density_m3: float, tritium_midplane_density_m3: float, ion_temperature_keV: float, symmetry_tolerance: float = 1.0e-2) -> IonSeedReferenceProfiles:
    """Build D and T profiles with the requested magnetic midplane densities

    Use matching axial coordinates and cell volumes to calculate particle inventories
    """
    z = np.asarray(z_m, dtype=float)
    volumes = np.asarray(cell_volumes_m3, dtype=float)
    if z.ndim != 1 or volumes.shape != z.shape:
        raise ValueError("background profile coordinates and cell volumes must be matching 1D arrays")
    if z.size < 2 or np.any(~np.isfinite(z)) or np.any(np.diff(z) <= 0.0):
        raise ValueError("background profile coordinates must be finite and strictly increasing")
    if np.any(~np.isfinite(volumes)) or np.any(volumes <= 0.0):
        raise ValueError("background profile cell volumes must be finite and positive")

    deuterium_midpoint = float(deuterium_midplane_density_m3)
    tritium_midpoint = float(tritium_midplane_density_m3)
    ion_temperature = float(ion_temperature_keV)
    if (not np.isfinite(deuterium_midpoint) or not np.isfinite(tritium_midpoint) or deuterium_midpoint < 0.0 or tritium_midpoint < 0.0):
        raise ValueError("background ion midplane densities must be finite and nonnegative")
    if not np.isfinite(ion_temperature) or ion_temperature <= 0.0:
        raise ValueError("background ion temperature must be finite and positive")

    shape, raw_midpoint, symmetry_error = _normalized_shape(z, support, specification, magnetic_midplane_z_m=magnetic_midplane_z_m, symmetry_tolerance=symmetry_tolerance)
    deuterium = deuterium_midpoint * shape
    tritium = tritium_midpoint * shape
    # Store the temperature only where the shared shape is positive
    temperature = np.where(shape > 0.0, ion_temperature, 0.0)
    inventories = {"deuterium": float(np.sum(deuterium * volumes)), "tritium": float(np.sum(tritium * volumes))}

    return IonSeedReferenceProfiles(
        z_m=z,
        shape=shape,
        deuterium_density_m3=deuterium,
        tritium_density_m3=tritium,
        ion_temperature_keV=temperature,
        support_domain=support,
        specification=specification,
        model=specification.model,
        magnetic_midplane_z_m=float(magnetic_midplane_z_m),
        deuterium_midplane_density_m3=deuterium_midpoint,
        tritium_midplane_density_m3=tritium_midpoint,
        shape_midplane_value_before_normalization=raw_midpoint,
        symmetry_relative_error=symmetry_error,
        symmetry_tolerance=float(symmetry_tolerance),
        inventories=inventories,
    )
