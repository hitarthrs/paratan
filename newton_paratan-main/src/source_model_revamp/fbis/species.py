"""Canonical charged ion species data"""
from __future__ import annotations
from dataclasses import dataclass
from types import MappingProxyType
import numpy as np
from source_model_revamp.constants import PROTON_MASS_KG, DEUTERON_MASS_KG, TRITON_MASS_KG, ALPHA_PARTICLE_MASS_KG

@dataclass(frozen=True)
class IonSpecies:
    """Immutable ion species definition"""
    species_id: str
    name: str
    symbol: str
    mass_kg: float
    charge_number: float
    nuclear_identity: str

    def __post_init__(self) -> None:
        if not self.species_id or not self.name or not self.symbol or not self.nuclear_identity:
            raise ValueError("ion species text fields must be nonempty")
        if not np.isfinite(self.mass_kg) or self.mass_kg <= 0.0:
            raise ValueError("ion species mass_kg must be positive and finite")
        if not np.isfinite(self.charge_number) or self.charge_number <= 0.0:
            raise ValueError("ion species charge_number must be positive and finite")

PROTON = IonSpecies(species_id="hydrogen", name="proton", symbol="H", mass_kg=PROTON_MASS_KG, charge_number=1.0, nuclear_identity="H-1")
DEUTERON = IonSpecies(species_id="deuterium", name="deuteron", symbol="D", mass_kg=DEUTERON_MASS_KG, charge_number=1.0, nuclear_identity="H-2")
TRITON = IonSpecies(species_id="tritium", name="triton", symbol="T", mass_kg=TRITON_MASS_KG, charge_number=1.0, nuclear_identity="H-3")
ALPHA_PARTICLE = IonSpecies(species_id="alpha", name="alpha_particle", symbol="alpha", mass_kg=ALPHA_PARTICLE_MASS_KG, charge_number=2.0, nuclear_identity="He-4")

ION_SPECIES_BY_ID = MappingProxyType({species.species_id: species for species in (PROTON, DEUTERON, TRITON, ALPHA_PARTICLE)})

def ion_species(species_id: str) -> IonSpecies:
    """Return a canonical ion species by identifier"""
    key = str(species_id).strip().lower()
    try:
        return ION_SPECIES_BY_ID[key]
    except KeyError as exc:
        raise ValueError(f"unsupported ion species {species_id!r}") from exc
