"""Represented deuterium and tritium fusion populations and allowed reaction pairs"""
from __future__ import annotations
from dataclasses import dataclass
from types import MappingProxyType
import numpy as np
from source_model_revamp.fbis.species import DEUTERON, TRITON, IonSpecies
from source_model_revamp.fbis.velocity_space_grid import PitchGrid, SpeedGrid
from source_model_revamp.fusion.reactions import DD_NEUTRON, DD_PROTON, DT_NEUTRON, FusionReaction

FAST_D = "fast_D"
FAST_T = "fast_T"

@dataclass(frozen=True)
class FusionReactantPopulation:
    """One represented full device fusion reactant population

    local_distribution_z_v_pitch stores f(z, v, ξ) with shape (n_z, n_speed, n_pitch)
    full_device_population and includes_directed_lost_population record spatial coverage
    """
    population_id: str
    species: IonSpecies
    population_kind: str
    speed_grid: SpeedGrid
    pitch_grid: PitchGrid
    local_distribution_z_v_pitch: np.ndarray
    full_device_population: bool
    includes_directed_lost_population: bool

    def __post_init__(self) -> None:
        """Validate the population identity and local distribution shape"""
        distribution = np.asarray(self.local_distribution_z_v_pitch, dtype=float)
        if distribution.ndim != 3:
            raise ValueError("fusion population distribution must be a 3D axial speed pitch array")
        expected = (distribution.shape[0], self.speed_grid.centers_m_s.size, self.pitch_grid.centers.size)
        if distribution.shape != expected:
            raise ValueError("fusion population distribution does not match its axial and velocity grids")
        if np.any(~np.isfinite(distribution)) or np.any(distribution < 0.0):
            raise ValueError("fusion population distribution must be finite and nonnegative")
        expected_identity = {FAST_D: (DEUTERON.species_id, 'fast'), FAST_T: (TRITON.species_id, 'fast')}
        try:
            species_id, kind = expected_identity[self.population_id]
        except KeyError as exc:
            raise ValueError(f"unsupported fusion population {self.population_id!r}") from exc
        if self.species.species_id != species_id or self.population_kind != kind:
            raise ValueError("fusion population identity does not match its species and kind")
        object.__setattr__(self, "local_distribution_z_v_pitch", distribution)

@dataclass(frozen=True)
class FusionPopulationPairDefinition:
    """One allowed represented population pair and its active reaction branches"""
    population_a_id: str
    population_b_id: str
    component_kind: str
    reactions: tuple[FusionReaction, ...]
    family: str

    @property
    def identical_population(self) -> bool:
        """True when both reactants are drawn from the same represented population"""
        return self.population_a_id == self.population_b_id

FUSION_POPULATION_SPECIES = MappingProxyType({FAST_D: DEUTERON, FAST_T: TRITON})
FUSION_POPULATION_DISPLAY_LABELS = MappingProxyType({FAST_D: 'beam-born kinetic deuterium', FAST_T: 'beam-born kinetic tritium'})

def fusion_population_display_label(population_id: str) -> str:
    """Return the display label for a represented population id"""
    return FUSION_POPULATION_DISPLAY_LABELS.get(str(population_id), str(population_id))

FUSION_POPULATION_PAIR_DEFINITIONS = (FusionPopulationPairDefinition(FAST_D, FAST_D, 'fast_fast', (DD_NEUTRON, DD_PROTON), 'fast_fast'), FusionPopulationPairDefinition(FAST_D, FAST_T, 'fast_fast', (DT_NEUTRON,), 'fast_fast'))

def active_population_pair_definitions(*, include_fast_fast: bool) -> tuple[FusionPopulationPairDefinition, ...]:
    """Return allowed pair definitions selected by the fast fast family control"""
    enabled = {'fast_fast': bool(include_fast_fast)}
    return tuple(definition for definition in FUSION_POPULATION_PAIR_DEFINITIONS if enabled[definition.family])