"""
Correlated neutron event bank assembly from authoritative fusion populations

The stage constructs event specifications for each neutron producing fusion component and applies evaluated center of momentum angular laws through the neutron events package
"""
from __future__ import annotations
from dataclasses import dataclass
from collections.abc import Mapping
import numpy as np
from source_model_revamp.fbis.species import DEUTERON, TRITON
from source_model_revamp.fusion.populations import FusionReactantPopulation
from source_model_revamp.integration.full_device_populations import FullDevicePopulationGrid
from source_model_revamp.integration.pipeline_types import FusionComponent
from source_model_revamp.integration.config.root import SourceModelRunConfig
from source_model_revamp.neutrons.events import CorrelatedNeutronEventBank, CorrelatedNeutronEventComponentSpec, build_correlated_neutron_event_bank, event_bank_source_rate_matrix_s
from source_model_revamp.neutrons.kinematics import neutron_kinematics_for_reaction

@dataclass(frozen=True)
class CorrelatedNeutronEventStageResult:
    """
    Optional correlated event bank plus its diagnostic axial energy histogram
    
    `source_rate_matrix_s` stores physical neutron rate in axial and energy bins when an event bank is built
    """
    event_bank: CorrelatedNeutronEventBank | None
    source_rate_matrix_s: np.ndarray | None
    out_of_range_rate_s: float | None

def _component_pair_gyroangles(component: FusionComponent) -> int:
    """Return the relative gyro angle quadrature count used by the fusion pair kernel for one component"""
    if component.pair_kernel is None:
        raise ValueError(f"fusion component {component.label!r} has no pair kernel")
 
    return int(component.pair_kernel.num_gyroangle_points)

def _component_event_spec(*, component: FusionComponent, populations_by_id: Mapping[str, FusionReactantPopulation]) -> CorrelatedNeutronEventComponentSpec:
    """
    Build one correlated neutron event specification from the fusion component and its reactant populations
    
    D T uses canonical deuteron first ordering
    Identical D D pairs use random projectile exchange symmetry before equivalent incident energy and angular law evaluation
    """
    if component.pair_kernel is None:
        raise ValueError(f"fusion component {component.label!r} has no pair kernel")
    try:
        population_a = populations_by_id[component.reactant_a_population_id]
        population_b = populations_by_id[component.reactant_b_population_id]
    except KeyError as exc:
        raise ValueError(f"fusion component {component.label!r} references an unavailable reactant population") from exc
  
    identical_population = component.reactant_a_population_id == component.reactant_b_population_id
  
    if component.reaction.key == "dt_n":
        if population_a.species.species_id != DEUTERON.species_id or population_b.species.species_id != TRITON.species_id:
            raise ValueError("DT correlated events require canonical deuteron first population ordering")
        projectile_assignment = "reactant_a_deuteron_projectile"
    elif identical_population:
        projectile_assignment = "random_exchange_symmetrized"
    else:
        if population_a.species.species_id != DEUTERON.species_id or population_b.species.species_id != DEUTERON.species_id:
            raise ValueError("DD correlated events require deuteron reactants")
        projectile_assignment = "reactant_a_deuteron_projectile"
   
    return CorrelatedNeutronEventComponentSpec(
        label=component.label,
        kind=component.kind,
        reaction=component.reaction,
        reaction_kinematics=neutron_kinematics_for_reaction(component.reaction),
        speed_grid_a=population_a.speed_grid,
        pitch_grid_a=population_a.pitch_grid,
        distributions_a_z_v_xi=population_a.local_distribution_z_v_pitch,
        mass_a_kg=population_a.species.mass_kg,
        speed_grid_b=population_b.speed_grid,
        pitch_grid_b=population_b.pitch_grid,
        distributions_b_z_v_xi=population_b.local_distribution_z_v_pitch,
        mass_b_kg=population_b.species.mass_kg,
        physical_rate_density_m3_s=np.asarray(component.pair_kernel.rate_density_m3_s, dtype=float),
        identical_population=identical_population,
        projectile_assignment_model=projectile_assignment,
        num_gyroangle_points=_component_pair_gyroangles(component),
    )

def build_correlated_neutron_event_stage(*, config: SourceModelRunConfig, neutron_components: tuple[FusionComponent, ...], populations_by_id: Mapping[str, FusionReactantPopulation], population_grid: FullDevicePopulationGrid, energy_edges_J: np.ndarray) -> CorrelatedNeutronEventStageResult:
    """
    Build the correlated event bank when the evaluated ENDF B VIII.1 center of momentum angular model is selected
    
    Events retain sampled source position, reactant velocities, evaluated angular state, and relativistic lab neutron energy and direction
    """
    ncfg = config.neutrons
    if ncfg.angular_model != "endf_b_viii1_evaluated_cm" or not neutron_components:
        return CorrelatedNeutronEventStageResult(None, None, None)
   
    event_specs = tuple(_component_event_spec(component=component, populations_by_id=populations_by_id) for component in neutron_components)
    event_bank = build_correlated_neutron_event_bank(
        event_specs,
        z_edges_m=np.asarray(population_grid.z_edges_m, dtype=float),
        cell_volumes_m3=np.asarray(population_grid.cell_volumes_m3, dtype=float),
        radial_inner_radius_m_by_z=np.asarray(population_grid.radial_inner_radius_m_by_z, dtype=float),
        radial_outer_radius_m_by_z=np.asarray(population_grid.radial_outer_radius_m_by_z, dtype=float),
        radial_outer_radius_m_by_edge=np.asarray(population_grid.radial_outer_radius_m_by_edge, dtype=float),
        event_count=ncfg.correlated_event_count,
        random_seed=ncfg.correlated_event_seed,
        radial_source_profile_model=ncfg.radial_source_profile_model,
        source_connected_radial_rho_limit_edges=np.asarray(population_grid.source_connected_radial_rho_max_edges, dtype=float),
    )
    source_rate_matrix_s, out_of_range_rate_s = event_bank_source_rate_matrix_s(event_bank, energy_edges_J, np.asarray(population_grid.z_centers_m).size, allow_out_of_range=ncfg.allow_out_of_range)
   
    return CorrelatedNeutronEventStageResult(event_bank=event_bank, source_rate_matrix_s=source_rate_matrix_s, out_of_range_rate_s=out_of_range_rate_s)
