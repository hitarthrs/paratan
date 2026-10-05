"""Startup seed or maintained reference ion and electron closure configuration"""
from __future__ import annotations
from collections.abc import Mapping
from dataclasses import dataclass, field
from typing import Any
from source_model_revamp.geometry.background_profiles import BackgroundProfileSpecification
from source_model_revamp.integration.plasma_closure import NBI_SUPPORTED_STATIONARY_PLASMA_CLOSURE_MODEL
from source_model_revamp.integration.config.common import canonical_model_name, _mapping, _check_unknown, _number, _array
from source_model_revamp.integration.config.constants import *

@dataclass(frozen=True)
class PlasmaClosureConfig:
    """
    Plasma closure inputs for background D and T, electron initial guesses, and axial profile shape
    
    Densities are in m⁻³ and temperatures are in keV
    """
    model: str = NBI_SUPPORTED_STATIONARY_PLASMA_CLOSURE_MODEL
    background_deuterium_midplane_density_m3: float = 1.0e20
    background_tritium_midplane_density_m3: float = 1.0e20
    background_ion_temperature_keV: float = 10.0
    electron_density_mode: str = QUASINEUTRAL_ELECTRON_DENSITY_MODE
    electron_density_initial_guess_m3: float = 2.0e20
    electron_density_initial_guess_supplied: bool = False
    electron_density_initial_guess_valid: bool = True
    electron_density_initial_guess_source: str = "positive_charge_seed_default"
    benchmark_prescribed_electron_density_m3: float | None = None
    electron_temperature_initial_guess_keV: float = 10.0
    background_profile: BackgroundProfileSpecification = field(default_factory=BackgroundProfileSpecification)

    @classmethod
    def from_mapping(cls, data: Mapping[str, Any], *, strict: bool = False) -> "PlasmaClosureConfig":
        """
        Parse the plasma closure and reject obsolete input blocks
        
        The electron density input is retained as an initial guess while the background profile can be uniform, analytic, or tabulated
        """
        if "initial_target" in data:
            raise ValueError("plasma_closure.initial_target is obsolete; use plasma_closure.background_ions and plasma_closure.electrons")
        if "fuel_mix" in data:
            raise ValueError("plasma_closure.fuel_mix is obsolete; use plasma_closure.background_ions")
        allowed = {
            "model",
            "background_ions",
            "electrons",
            "background_profile",
        }
        _check_unknown(data, allowed, "source_model.plasma_closure", strict)
        model = canonical_model_name(data.get("model"), PLASMA_CLOSURE_MODEL_ALIASES, NBI_SUPPORTED_STATIONARY_PLASMA_CLOSURE_MODEL)
        background_ions = _mapping(data.get("background_ions"), "source_model.plasma_closure.background_ions")
        electrons = _mapping(data.get("electrons"), "source_model.plasma_closure.electrons")
        profile = _mapping(data.get("background_profile"), "source_model.plasma_closure.background_profile")
        _check_unknown(background_ions, {"deuterium_midplane_density_m3", "tritium_midplane_density_m3", "ion_temperature_keV"}, "source_model.plasma_closure.background_ions", strict)
        _check_unknown(electrons, {"density_mode", "density_initial_guess_m3", "temperature_initial_guess_keV", "benchmark_prescribed_density_m3"}, "source_model.plasma_closure.electrons", strict)
        _check_unknown(profile, {"model", "polynomial_coefficients", "z_m", "scale"}, "source_model.plasma_closure.background_profile", strict)
        profile_model = str(profile.get("model", "uniform_central_cell"))
        supported_profile_models = {
            "uniform_central_cell",
            "analytic_axial_profile",
            "tabulated_axial_profile",
        }
        if profile_model not in supported_profile_models:
            raise ValueError("plasma_closure.background_profile.model must be one of " f"{sorted(supported_profile_models)}")
        coefficients = _array(profile.get("polynomial_coefficients"), "plasma_closure.background_profile.polynomial_coefficients", default=(1.0,)) if profile_model == "analytic_axial_profile" else ()
        table_z = _array(profile.get("z_m"), "plasma_closure.background_profile.z_m") if profile_model == "tabulated_axial_profile" else ()
        table_scale = _array(profile.get("scale"), "plasma_closure.background_profile.scale", nonnegative=True) if profile_model == "tabulated_axial_profile" else ()
        if profile_model == "tabulated_axial_profile" and len(table_z) != len(table_scale):
            raise ValueError("plasma_closure.background_profile z_m and scale must have equal length")
        if profile_model == "tabulated_axial_profile" and (len(table_z) < 2 or any(right <= left for left, right in zip(table_z[:-1], table_z[1:], strict=True))):
            raise ValueError("plasma_closure.background_profile z_m must be strictly increasing")

        deuterium_density = _number(background_ions.get("deuterium_midplane_density_m3"), "plasma_closure.background_ions.deuterium_midplane_density_m3", nonnegative=True, default=1.0e20)
        tritium_density = _number(background_ions.get("tritium_midplane_density_m3"), "plasma_closure.background_ions.tritium_midplane_density_m3", nonnegative=True, default=1.0e20)
        ion_temperature = _number(
            background_ions.get("ion_temperature_keV"),
            "plasma_closure.background_ions.ion_temperature_keV",
            positive=True,
            default=10.0)
        electron_density_mode = canonical_model_name(electrons.get("density_mode"), ELECTRON_DENSITY_MODE_ALIASES, QUASINEUTRAL_ELECTRON_DENSITY_MODE)
        # Start the electron guess from the singly charged positive ion seed when available
        default_density_guess = deuterium_density + tritium_density
        if default_density_guess <= 0.0:
            default_density_guess = 2.0e20
        density_guess_input = electrons.get("density_initial_guess_m3")
        density_guess_supplied = density_guess_input is not None
        density_guess = _number(density_guess_input, "plasma_closure.electrons.density_initial_guess_m3", positive=True, default=default_density_guess)
        density_guess_valid = True
        if density_guess_supplied:
            density_guess_source = 'plasma_closure_electrons_density_initial_guess_m3'
        elif model == NBI_SUPPORTED_STATIONARY_PLASMA_CLOSURE_MODEL:
            density_guess_source = 'startup_positive_charge_seed_default'
        else:
            density_guess_source = 'prescribed_reference_charge_default'

        return cls(
            model=model,
            background_deuterium_midplane_density_m3=deuterium_density,
            background_tritium_midplane_density_m3=tritium_density,
            background_ion_temperature_keV=ion_temperature,
            electron_density_mode=electron_density_mode,
            electron_density_initial_guess_m3=density_guess,
            electron_density_initial_guess_supplied=density_guess_supplied,
            electron_density_initial_guess_valid=density_guess_valid,
            electron_density_initial_guess_source=density_guess_source,
            benchmark_prescribed_electron_density_m3=(
                _number(electrons.get("benchmark_prescribed_density_m3"), "plasma_closure.electrons.benchmark_prescribed_density_m3", positive=True)
                if electrons.get("benchmark_prescribed_density_m3") is not None
                else (default_density_guess if model == "egedal_beam_plasma_quasineutral" else None)
            ),
            electron_temperature_initial_guess_keV=_number(electrons.get("temperature_initial_guess_keV"), "plasma_closure.electrons.temperature_initial_guess_keV", positive=True),
            background_profile=BackgroundProfileSpecification(model=profile_model, analytic_polynomial_coefficients=tuple(coefficients), tabulated_z_m=tuple(table_z), tabulated_scale=tuple(table_scale)),
        )

    @property
    def uses_startup_seed_only(self) -> bool:
        """Return whether the configured D and T background is used only to initialize the solved state"""
        return self.model == NBI_SUPPORTED_STATIONARY_PLASMA_CLOSURE_MODEL

    @property
    def startup_deuterium_midplane_density_m3(self) -> float:
        """Return the configured deuterium startup seed density in m⁻³"""
        return self.background_deuterium_midplane_density_m3

    @property
    def startup_tritium_midplane_density_m3(self) -> float:
        """Return the configured tritium startup seed density in m⁻³"""
        return self.background_tritium_midplane_density_m3

    @property
    def startup_ion_temperature_keV(self) -> float:
        """Return the configured startup ion temperature in keV"""
        return self.background_ion_temperature_keV

    @property
    def background_positive_charge_midplane_density_m3(self) -> float:
        """Return the singly charged D plus T positive charge number density at the midplane in m⁻³"""
        return self.background_deuterium_midplane_density_m3 + self.background_tritium_midplane_density_m3
