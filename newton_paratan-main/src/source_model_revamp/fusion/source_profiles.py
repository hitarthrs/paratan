"""Fusion reaction, neutron source, and power profile assembly"""
from __future__ import annotations
from dataclasses import dataclass
from collections.abc import Sequence
import numpy as np
from numpy.typing import ArrayLike
from source_model_revamp.geometry.axial_grid_geometry import volume_integral
from source_model_revamp.fusion.reactions import FusionReaction

@dataclass(frozen=True)
class FusionSourceProfile:
    """Local and volume integrated source quantities for one reaction branch

    Profile arrays have shape (n_z,) with rates in m^−3 s^−1 and powers in W/m^3
    Integrated totals are available only when cell_volumes_m3 is supplied
    """
    coordinate: np.ndarray
    reaction: FusionReaction
    reaction_rate_density_m3_s: np.ndarray
    neutron_source_density_m3_s: np.ndarray
    fusion_power_density_W_m3: np.ndarray
    neutron_power_density_W_m3: np.ndarray
    charged_product_power_density_W_m3: np.ndarray
    cell_volumes_m3: np.ndarray | None = None
    total_reaction_rate_s: float | None = None
    total_neutron_rate_s: float | None = None
    total_fusion_power_W: float | None = None
    total_neutron_power_W: float | None = None
    total_charged_product_power_W: float | None = None

@dataclass(frozen=True)
class MultiReactionFusionSourceProfile:
    """Collection of reaction profiles and their summed rate and power profiles"""
    reaction_profiles: tuple[FusionSourceProfile, ...]
    total_reaction_rate_density_m3_s: np.ndarray
    total_neutron_source_density_m3_s: np.ndarray
    total_fusion_power_density_W_m3: np.ndarray
    total_neutron_power_density_W_m3: np.ndarray
    total_charged_product_power_density_W_m3: np.ndarray
    total_reaction_rate_s: float | None
    total_neutron_rate_s: float | None
    total_fusion_power_W: float | None
    total_neutron_power_W: float | None
    total_charged_product_power_W: float | None

def _as_1d_profile(name: str, values: ArrayLike) -> np.ndarray:
    """Return a finite nonempty one dimensional profile"""
    array = np.asarray(values, dtype=float)
    if array.ndim != 1:
        raise ValueError(f"{name} must be a 1D profile")
    if array.size == 0:
        raise ValueError(f"{name} must not be empty")
    if not np.all(np.isfinite(array)):
        raise ValueError(f"{name} must contain only finite values")
   
    return array

def fusion_source_profile_from_rate_density(coordinate: ArrayLike, reaction: FusionReaction, reaction_rate_density_m3_s: ArrayLike, cell_volumes_m3: ArrayLike | None = None) -> FusionSourceProfile:
    """Build neutron and power profiles from one reaction rate density

    FusionReaction supplies the branch energy and neutron yield used for each conversion
    """
    z = _as_1d_profile("coordinate", coordinate)
    rate = _as_1d_profile("reaction_rate_density_m3_s", reaction_rate_density_m3_s)
    if rate.size != z.size:
        raise ValueError("reaction_rate_density_m3_s must have one value per coordinate")
    if np.any(rate < 0.0):
        raise ValueError("reaction_rate_density_m3_s must be nonnegative")
    neutron_density = reaction.neutron_yield_per_reaction * rate
    fusion_power = reaction.total_energy_J * rate
    neutron_power = reaction.neutron_yield_per_reaction * reaction.neutron_energy_J * rate
    charged_power = reaction.charged_product_energy_J * rate
    volumes = None
    total_reaction_rate = None
    total_neutron_rate = None
    total_fusion_power = None
    total_neutron_power = None
    total_charged_power = None
    if cell_volumes_m3 is not None:
        volumes = _as_1d_profile("cell_volumes_m3", cell_volumes_m3)
        if volumes.size != z.size:
            raise ValueError("cell_volumes_m3 must have one value per coordinate")
        if np.any(volumes <= 0.0):
            raise ValueError("cell_volumes_m3 must be positive")
        total_reaction_rate = volume_integral(rate, volumes)
        total_neutron_rate = volume_integral(neutron_density, volumes)
        total_fusion_power = volume_integral(fusion_power, volumes)
        total_neutron_power = volume_integral(neutron_power, volumes)
        total_charged_power = volume_integral(charged_power, volumes)

    return FusionSourceProfile(
        coordinate=z,
        reaction=reaction,
        reaction_rate_density_m3_s=rate,
        neutron_source_density_m3_s=neutron_density,
        fusion_power_density_W_m3=fusion_power,
        neutron_power_density_W_m3=neutron_power,
        charged_product_power_density_W_m3=charged_power,
        cell_volumes_m3=volumes,
        total_reaction_rate_s=total_reaction_rate,
        total_neutron_rate_s=total_neutron_rate,
        total_fusion_power_W=total_fusion_power,
        total_neutron_power_W=total_neutron_power,
        total_charged_product_power_W=total_charged_power,
    )

def combine_fusion_source_profiles(reaction_profiles: Sequence[FusionSourceProfile]) -> MultiReactionFusionSourceProfile:
    """Sum compatible FusionSourceProfile objects into total rate and power profiles"""
    profiles = tuple(reaction_profiles)
    if not profiles:
        raise ValueError("at least one FusionSourceProfile is required")
    shapes = {p.reaction_rate_density_m3_s.shape for p in profiles}
    if len(shapes) != 1:
        raise ValueError("all reaction profiles must have the same profile shape")
    total_rate_density = np.sum([p.reaction_rate_density_m3_s for p in profiles], axis=0)
    total_neutron_density = np.sum([p.neutron_source_density_m3_s for p in profiles], axis=0)
    total_power = np.sum([p.fusion_power_density_W_m3 for p in profiles], axis=0)
    total_neutron_power = np.sum([p.neutron_power_density_W_m3 for p in profiles], axis=0)
    total_charged_power = np.sum([p.charged_product_power_density_W_m3 for p in profiles], axis=0)
    def maybe_sum(attribute: str) -> float | None:
        """Sum one integrated quantity only when every profile provides it"""
        values = [getattr(p, attribute) for p in profiles]
        if any(value is None for value in values):
            return None
        return float(np.sum(values))

    return MultiReactionFusionSourceProfile(
        reaction_profiles=profiles,
        total_reaction_rate_density_m3_s=total_rate_density,
        total_neutron_source_density_m3_s=total_neutron_density,
        total_fusion_power_density_W_m3=total_power,
        total_neutron_power_density_W_m3=total_neutron_power,
        total_charged_product_power_density_W_m3=total_charged_power,
        total_reaction_rate_s=maybe_sum("total_reaction_rate_s"),
        total_neutron_rate_s=maybe_sum("total_neutron_rate_s"),
        total_fusion_power_W=maybe_sum("total_fusion_power_W"),
        total_neutron_power_W=maybe_sum("total_neutron_power_W"),
        total_charged_product_power_W=maybe_sum("total_charged_product_power_W"),
    )