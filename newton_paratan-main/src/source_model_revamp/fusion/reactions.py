"""Fusion reaction branch metadata and unit conversions"""
from __future__ import annotations
from dataclasses import dataclass
from scipy.constants import elementary_charge

EV_TO_J = elementary_charge
KEV_TO_J = 1.0e3 * EV_TO_J
MEV_TO_J = 1.0e6 * EV_TO_J
J_TO_KEV = 1.0 / KEV_TO_J
J_TO_MEV = 1.0 / MEV_TO_J
BARN_TO_M2 = 1.0e-28
M2_TO_BARN = 1.0e28
CM3_TO_M3 = 1.0e-6
M3_TO_CM3 = 1.0e6

@dataclass(frozen=True)
class FusionReaction:
    """Bookkeeping metadata for one fusion reaction branch

    Energies are stored in J and neutron_yield_per_reaction is dimensionless
    """
    key: str
    name: str
    reactant_a: str
    reactant_b: str
    total_energy_J: float
    neutron_energy_J: float = 0.0
    neutron_yield_per_reaction: float = 0.0

    @property
    def identical_species(self) -> bool:
        """True if both reactants are the same species label"""
        return self.reactant_a == self.reactant_b
    @property
    def same_species_pair_counting_factor(self) -> float:
        """Species only pair counting factor, 1/2 for identical species"""
        return 0.5 if self.identical_species else 1.0
    @property
    def charged_product_energy_J(self) -> float:
        """Fusion energy not carried by neutrons"""
        return self.total_energy_J - self.neutron_yield_per_reaction * self.neutron_energy_J

DT_NEUTRON = FusionReaction(key="dt_n", name="D(T,n)alpha", reactant_a="D", reactant_b="T", total_energy_J=17.589 * MEV_TO_J, neutron_energy_J=14.1 * MEV_TO_J, neutron_yield_per_reaction=1.0)
DD_NEUTRON = FusionReaction(key="dd_n", name="D(D,n)He3", reactant_a="D", reactant_b="D", total_energy_J=3.269 * MEV_TO_J, neutron_energy_J=2.45 * MEV_TO_J, neutron_yield_per_reaction=1.0)
DD_PROTON = FusionReaction(key="dd_p", name="D(D,p)T", reactant_a="D", reactant_b="D", total_energy_J=4.033 * MEV_TO_J, neutron_energy_J=0.0, neutron_yield_per_reaction=0.0)

ACTIVE_REACTIONS = (DT_NEUTRON, DD_NEUTRON, DD_PROTON)
REACTIONS_BY_KEY = {reaction.key: reaction for reaction in ACTIVE_REACTIONS}
REACTIONS_BY_NAME = {reaction.name: reaction for reaction in ACTIVE_REACTIONS}

def energy_J_from_keV(energy_keV: float):
    """Convert keV to joules"""
    return float(energy_keV) * KEV_TO_J

def energy_keV_from_J(energy_J: float):
    """Convert joules to keV"""
    return float(energy_J) * J_TO_KEV

def energy_J_from_MeV(energy_MeV: float):
    """Convert MeV to joules"""
    return float(energy_MeV) * MEV_TO_J

def cross_section_m2_from_barns(cross_section_barns):
    """Convert barns to m^2"""
    return cross_section_barns * BARN_TO_M2

def cross_section_barns_from_m2(cross_section_m2):
    """Convert m^2 to barns"""
    return cross_section_m2 * M2_TO_BARN

def reactivity_m3_s_from_cm3_s(reactivity_cm3_s):
    """Convert cm^3/s to m^3/s"""
    return reactivity_cm3_s * CM3_TO_M3

def reactivity_cm3_s_from_m3_s(reactivity_m3_s):
    """Convert m^3/s to cm^3/s"""
    return reactivity_m3_s * M3_TO_CM3

def reaction_from_key(key: str) -> FusionReaction:
    """Return supported reaction metadata from a compact reaction key"""
    try:
        return REACTIONS_BY_KEY[str(key).strip().lower()]
    except KeyError as exc:
        raise ValueError(f"Unsupported active fusion reaction key {key!r}") from exc