"""
Geometry linked neutral beam attenuation for the mirror source model

attenuation equation is the cell integrated Beer-Lambert update for a primary neutral rate N0 along path coordinate s:
    dN0/ds = -κ_tot(s) N0
    κ_tot = κ_ion + κ_cx + κ_other

For constant coefficients in a path segment Δs:
    N_out = N_in * exp(-κ_tot * Δs)
    ΔN_c = N_in * (κ_c/κ_tot) * [1 - exp(-κ_tot * Δs)]

where c is ionization, charge exchange, or other neutral loss  
Ionization and primary neutral charge exchange both create a fast beam ion
"""
from __future__ import annotations
from dataclasses import dataclass
import numpy as np
from numpy.typing import ArrayLike
from source_model_revamp.beam.beam_path_geometry import BeamPathOnAxisymmetricGrid
from source_model_revamp.fbis.beam_source_definition import component_powers_from_fractions_W, particle_birth_rate_from_power_s

@dataclass(frozen=True)
class NeutralBeamAttenuation:
    """
    Cell integrated attenuation result for one neutral beam energy component
    
    Axial rate coefficient arrays use m⁻¹ and cell event arrays use s⁻¹
    Rate density arrays use m⁻³ s⁻¹ after division by the full flux tube cell volume
    neutral_rate_faces_s follows traversed beam cells and therefore has one more entry than cell_indices_along_beam rather than one entry per axial grid face
    """
    coordinate_centers_m: np.ndarray
    coordinate_faces_m: np.ndarray
    path_lengths_m: np.ndarray
    cell_volumes_m3: np.ndarray
    B_tilde_centers: np.ndarray
    reference_area_m2: float
    input_neutral_rate_s: float
    neutral_rate_faces_s: np.ndarray
    ionization_rate_per_m: np.ndarray
    charge_exchange_rate_per_m: np.ndarray
    other_loss_rate_per_m: np.ndarray
    total_attenuation_rate_per_m: np.ndarray
    cell_ionization_rates_s: np.ndarray
    cell_charge_exchange_rates_s: np.ndarray
    cell_other_loss_rates_s: np.ndarray
    ionization_rate_density_m3_s: np.ndarray
    charge_exchange_rate_density_m3_s: np.ndarray
    other_loss_rate_density_m3_s: np.ndarray
    fast_ion_birth_rate_density_m3_s: np.ndarray
    net_fueling_rate_density_m3_s: np.ndarray
    electron_ionization_rate_per_m: np.ndarray
    deuterium_ionization_rate_per_m: np.ndarray
    tritium_ionization_rate_per_m: np.ndarray
    deuterium_charge_exchange_rate_per_m: np.ndarray
    tritium_charge_exchange_rate_per_m: np.ndarray
    cell_electron_ionization_rates_s: np.ndarray
    cell_deuterium_ionization_rates_s: np.ndarray
    cell_tritium_ionization_rates_s: np.ndarray
    cell_deuterium_charge_exchange_rates_s: np.ndarray
    cell_tritium_charge_exchange_rates_s: np.ndarray
    electron_ionization_rate_density_m3_s: np.ndarray
    deuterium_ionization_rate_density_m3_s: np.ndarray
    tritium_ionization_rate_density_m3_s: np.ndarray
    deuterium_charge_exchange_rate_density_m3_s: np.ndarray
    tritium_charge_exchange_rate_density_m3_s: np.ndarray

@dataclass(frozen=True)
class MultiEnergyNeutralBeamAttenuation:
    """
    Attenuation results for all energy components of one configured beam
    
    component_results preserves one NeutralBeamAttenuation per energy component
    Total rate density arrays are sums over components on the shared axial grid
    """
    component_results: tuple[NeutralBeamAttenuation, ...]
    component_energies_J: np.ndarray
    component_powers_W: np.ndarray
    component_input_rates_s: np.ndarray
    total_input_power_W: float
    total_input_neutral_rate_s: float
    total_fast_ion_birth_rate_density_m3_s: np.ndarray
    total_net_fueling_rate_density_m3_s: np.ndarray
    total_ionization_rate_density_m3_s: np.ndarray
    total_charge_exchange_rate_density_m3_s: np.ndarray
    total_other_loss_rate_density_m3_s: np.ndarray
    total_electron_ionization_rate_density_m3_s: np.ndarray
    total_deuterium_ionization_rate_density_m3_s: np.ndarray
    total_tritium_ionization_rate_density_m3_s: np.ndarray
    total_deuterium_charge_exchange_rate_density_m3_s: np.ndarray
    total_tritium_charge_exchange_rate_density_m3_s: np.ndarray
SUPPORTED_ATTENUATION_SOURCE_CHANNELS = frozenset({"fast_ion_birth", "ionization", "charge_exchange", "net_fueling", "other_loss",})

def normalize_attenuation_source_channel(channel: str = "fast_ion_birth") -> str:
    """Normalize and validate a supported attenuation source channel name"""
    normalized = str(channel).strip().lower().replace("-", "_")
    if normalized not in SUPPORTED_ATTENUATION_SOURCE_CHANNELS:
        allowed = ", ".join(sorted(SUPPORTED_ATTENUATION_SOURCE_CHANNELS))
        raise ValueError(f"Unsupported attenuation source channel {channel!r}, use one of: {allowed}")
    
    return normalized

def attenuation_channel_rate_density(attenuation: NeutralBeamAttenuation, channel: str = "fast_ion_birth") -> np.ndarray:
    """Return the selected attenuation channel rate density [m^-3 s^-1]"""
    channel_norm = normalize_attenuation_source_channel(channel)
    if channel_norm == "fast_ion_birth":
        return np.asarray(attenuation.fast_ion_birth_rate_density_m3_s, dtype=float)
    if channel_norm == "ionization":
        return np.asarray(attenuation.ionization_rate_density_m3_s, dtype=float)
    if channel_norm == "charge_exchange":
        return np.asarray(attenuation.charge_exchange_rate_density_m3_s, dtype=float)
    if channel_norm == "net_fueling":
        return np.asarray(attenuation.net_fueling_rate_density_m3_s, dtype=float)
    if channel_norm == "other_loss":
        return np.asarray(attenuation.other_loss_rate_density_m3_s, dtype=float)
    
    raise AssertionError("validated attenuation channel was not handled")

def _finite_array(name: str, values: ArrayLike) -> np.ndarray:
    array = np.asarray(values, dtype=float)
    if not np.all(np.isfinite(array)):
        raise ValueError(f"{name} must contain only finite values")
    
    return array

def _nonnegative_array(name: str, values: ArrayLike) -> np.ndarray:
    array = _finite_array(name, values)
    if np.any(array < 0.0):
        raise ValueError(f"{name} must be non negative")
    
    return array

def attenuation_rate_per_m(target_density_m3: ArrayLike, cross_section_m2: ArrayLike) -> np.ndarray:
    """Effective line attenuation rate: κ = n * σ [m^-1]"""
    density = _nonnegative_array("target_density_m3", target_density_m3)
    sigma = _nonnegative_array("cross_section_m2", cross_section_m2)
    return density * sigma

def sum_attenuation_rates_per_m(rate_terms_per_m: ArrayLike) -> np.ndarray:
    """Sum one or more attenuation rate contributions over the first axis"""
    terms = _nonnegative_array("rate_terms_per_m", rate_terms_per_m)

    return np.sum(terms, axis=0)

def total_attenuation_rate_per_m(ionization_rate_per_m: ArrayLike, charge_exchange_rate_per_m: ArrayLike, other_loss_rate_per_m: ArrayLike | float = 0.0) -> np.ndarray:
    """Return κ_tot = κ_ion + κ_cx + κ_other"""
    k_ion = _nonnegative_array("ionization_rate_per_m", ionization_rate_per_m)
    k_cx = _nonnegative_array("charge_exchange_rate_per_m", charge_exchange_rate_per_m)
    k_other = _nonnegative_array("other_loss_rate_per_m", other_loss_rate_per_m) + np.zeros_like(k_ion)
    if k_cx.shape != k_ion.shape or k_other.shape != k_ion.shape:
        raise ValueError("attenuation rate arrays must have matching shapes")
    
    return k_ion + k_cx + k_other

def neutral_transmission_fraction(total_attenuation_rate_per_m: ArrayLike, path_lengths_m: ArrayLike) -> np.ndarray:
    """Return exp(-κ_tot Δs) for each path cell"""
    k_total = _nonnegative_array("total_attenuation_rate_per_m", total_attenuation_rate_per_m)
    ds = _nonnegative_array("path_lengths_m", path_lengths_m)
    if k_total.shape != ds.shape:
        raise ValueError("total_attenuation_rate_per_m and path_lengths_m must have matching shapes")
    
    return np.exp(-k_total * ds)

def channel_deposition_fraction(channel_rate_per_m: ArrayLike, total_rate_per_m: ArrayLike, path_lengths_m: ArrayLike) -> np.ndarray:
    """Fraction of incoming primary neutrals removed by one channel in each cell"""
    k_channel = _nonnegative_array("channel_rate_per_m", channel_rate_per_m)
    k_total = _nonnegative_array("total_rate_per_m", total_rate_per_m)
    ds = _nonnegative_array("path_lengths_m", path_lengths_m)
    if k_channel.shape != k_total.shape or ds.shape != k_total.shape:
        raise ValueError("channel_rate_per_m, total_rate_per_m, and path_lengths_m must have matching shapes")
    removed = 1.0 - np.exp(-k_total * ds)
    channel_share = np.divide(k_channel, k_total, out=np.zeros_like(k_total, dtype=float), where=k_total > 0.0)

    return channel_share * removed

def attenuate_neutral_beam_on_path_grid(input_neutral_rate_s: float, path_geometry: BeamPathOnAxisymmetricGrid, ionization_rate_per_m: ArrayLike, charge_exchange_rate_per_m: ArrayLike, other_loss_rate_per_m: ArrayLike | float = 0.0, *, electron_ionization_rate_per_m: ArrayLike, deuterium_ionization_rate_per_m: ArrayLike, tritium_ionization_rate_per_m: ArrayLike, deuterium_charge_exchange_rate_per_m: ArrayLike, tritium_charge_exchange_rate_per_m: ArrayLike) -> NeutralBeamAttenuation:
    """Cell integrated attenuation of a primary neutral beam along a straight beam path"""
    input_rate = float(input_neutral_rate_s)
    if not np.isfinite(input_rate) or input_rate < 0.0:
        raise ValueError("input_neutral_rate_s must be non-negative and finite")
    geom = path_geometry
    centers = np.asarray(geom.coordinate_centers_m, dtype=float)
    B_tilde = np.asarray(geom.B_tilde_centers, dtype=float)
    k_ion = _nonnegative_array("ionization_rate_per_m", ionization_rate_per_m) + np.zeros_like(centers)
    k_cx = _nonnegative_array("charge_exchange_rate_per_m", charge_exchange_rate_per_m) + np.zeros_like(centers)
    k_other = _nonnegative_array("other_loss_rate_per_m", other_loss_rate_per_m) + np.zeros_like(centers)
    k_e_ion = _nonnegative_array("electron_ionization_rate_per_m", electron_ionization_rate_per_m) + np.zeros_like(centers)
    k_D_ion = _nonnegative_array("deuterium_ionization_rate_per_m", deuterium_ionization_rate_per_m) + np.zeros_like(centers)
    k_T_ion = _nonnegative_array("tritium_ionization_rate_per_m", tritium_ionization_rate_per_m) + np.zeros_like(centers)
    k_D_cx = _nonnegative_array("deuterium_charge_exchange_rate_per_m", deuterium_charge_exchange_rate_per_m) + np.zeros_like(centers)
    k_T_cx = _nonnegative_array("tritium_charge_exchange_rate_per_m", tritium_charge_exchange_rate_per_m) + np.zeros_like(centers)
    if not np.allclose(k_ion, k_e_ion + k_D_ion + k_T_ion, rtol=1.0e-12, atol=0.0):
        raise ValueError("resolved ionization rates must sum to the aggregate ionization rate")
    if not np.allclose(k_cx, k_D_cx + k_T_cx, rtol=1.0e-12, atol=0.0):
        raise ValueError("resolved charge exchange rates must sum to the aggregate charge exchange rate")
    k_total = total_attenuation_rate_per_m(k_ion, k_cx, k_other)
    path_lengths = _nonnegative_array("effective_path_lengths_m", geom.effective_path_lengths_m)
    volumes = _nonnegative_array("cell_volumes_m3", geom.cell_volumes_m3)
    if centers.shape != path_lengths.shape or centers.shape != volumes.shape:
        raise ValueError("path geometry arrays have inconsistent axial dimensions")
    neutral_path_faces = np.empty(geom.cell_indices_along_beam.size + 1, dtype=float)
    ionization_cells = np.zeros(centers.size, dtype=float)
    charge_exchange_cells = np.zeros(centers.size, dtype=float)
    other_loss_cells = np.zeros(centers.size, dtype=float)
    electron_ionization_cells = np.zeros(centers.size, dtype=float)
    deuterium_ionization_cells = np.zeros(centers.size, dtype=float)
    tritium_ionization_cells = np.zeros(centers.size, dtype=float)
    deuterium_charge_exchange_cells = np.zeros(centers.size, dtype=float)
    tritium_charge_exchange_cells = np.zeros(centers.size, dtype=float)
    neutral_path_faces[0] = input_rate
    incoming = input_rate
    for j, i in enumerate(np.asarray(geom.cell_indices_along_beam, dtype=int)):
        ds = float(path_lengths[i])
        if ds <= 0.0:
            neutral_path_faces[j + 1] = incoming
            continue
        kt = float(k_total[i])
        if kt <= 0.0:
            neutral_path_faces[j + 1] = incoming
            continue
        removed_fraction = 1.0 - np.exp(-kt * ds)
        ionization_cells[i] = incoming * (float(k_ion[i]) / kt) * removed_fraction
        charge_exchange_cells[i] = incoming * (float(k_cx[i]) / kt) * removed_fraction
        other_loss_cells[i] = incoming * (float(k_other[i]) / kt) * removed_fraction
        electron_ionization_cells[i] = incoming * (float(k_e_ion[i]) / kt) * removed_fraction
        deuterium_ionization_cells[i] = incoming * (float(k_D_ion[i]) / kt) * removed_fraction
        tritium_ionization_cells[i] = incoming * (float(k_T_ion[i]) / kt) * removed_fraction
        deuterium_charge_exchange_cells[i] = incoming * (float(k_D_cx[i]) / kt) * removed_fraction
        tritium_charge_exchange_cells[i] = incoming * (float(k_T_cx[i]) / kt) * removed_fraction
        incoming *= np.exp(-kt * ds)
        neutral_path_faces[j + 1] = incoming
    ionization_density = ionization_cells / volumes
    charge_exchange_density = charge_exchange_cells / volumes
    other_loss_density = other_loss_cells / volumes
    fast_ion_birth_density = ionization_density + charge_exchange_density
    net_fueling_density = ionization_density
    electron_ionization_density = electron_ionization_cells / volumes
    deuterium_ionization_density = deuterium_ionization_cells / volumes
    tritium_ionization_density = tritium_ionization_cells / volumes
    deuterium_charge_exchange_density = deuterium_charge_exchange_cells / volumes
    tritium_charge_exchange_density = tritium_charge_exchange_cells / volumes

    return NeutralBeamAttenuation(
        coordinate_centers_m=centers,
        coordinate_faces_m=np.asarray(geom.coordinate_faces_m, dtype=float),
        path_lengths_m=path_lengths,
        cell_volumes_m3=volumes,
        B_tilde_centers=B_tilde,
        reference_area_m2=float(geom.reference_area_m2),
        input_neutral_rate_s=input_rate,
        neutral_rate_faces_s=neutral_path_faces,
        ionization_rate_per_m=k_ion,
        charge_exchange_rate_per_m=k_cx,
        other_loss_rate_per_m=k_other,
        total_attenuation_rate_per_m=k_total,
        cell_ionization_rates_s=ionization_cells,
        cell_charge_exchange_rates_s=charge_exchange_cells,
        cell_other_loss_rates_s=other_loss_cells,
        ionization_rate_density_m3_s=ionization_density,
        charge_exchange_rate_density_m3_s=charge_exchange_density,
        other_loss_rate_density_m3_s=other_loss_density,
        fast_ion_birth_rate_density_m3_s=fast_ion_birth_density,
        net_fueling_rate_density_m3_s=net_fueling_density,
        electron_ionization_rate_per_m=k_e_ion,
        deuterium_ionization_rate_per_m=k_D_ion,
        tritium_ionization_rate_per_m=k_T_ion,
        deuterium_charge_exchange_rate_per_m=k_D_cx,
        tritium_charge_exchange_rate_per_m=k_T_cx,
        cell_electron_ionization_rates_s=electron_ionization_cells,
        cell_deuterium_ionization_rates_s=deuterium_ionization_cells,
        cell_tritium_ionization_rates_s=tritium_ionization_cells,
        cell_deuterium_charge_exchange_rates_s=deuterium_charge_exchange_cells,
        cell_tritium_charge_exchange_rates_s=tritium_charge_exchange_cells,
        electron_ionization_rate_density_m3_s=electron_ionization_density,
        deuterium_ionization_rate_density_m3_s=deuterium_ionization_density,
        tritium_ionization_rate_density_m3_s=tritium_ionization_density,
        deuterium_charge_exchange_rate_density_m3_s=deuterium_charge_exchange_density,
        tritium_charge_exchange_rate_density_m3_s=tritium_charge_exchange_density,
    )

def total_ionization_rate_s(attenuation: NeutralBeamAttenuation) -> float:
    """Total impact ionization event rate"""
    return float(np.sum(attenuation.cell_ionization_rates_s))

def total_charge_exchange_rate_s(attenuation: NeutralBeamAttenuation) -> float:
    """Total primary charge exchange event rate"""
    return float(np.sum(attenuation.cell_charge_exchange_rates_s))

def total_electron_ionization_rate_s(attenuation: NeutralBeamAttenuation) -> float:
    """Total electron impact ionization event rate"""
    return float(np.sum(attenuation.cell_electron_ionization_rates_s))

def total_deuterium_ionization_rate_s(attenuation: NeutralBeamAttenuation) -> float:
    """Total thermal deuterium impact ionization event rate"""
    return float(np.sum(attenuation.cell_deuterium_ionization_rates_s))

def total_tritium_ionization_rate_s(attenuation: NeutralBeamAttenuation) -> float:
    """Total thermal tritium impact ionization event rate"""
    return float(np.sum(attenuation.cell_tritium_ionization_rates_s))

def total_deuterium_charge_exchange_rate_s(attenuation: NeutralBeamAttenuation) -> float:
    """Total charge exchange event rate on thermal deuterium"""
    return float(np.sum(attenuation.cell_deuterium_charge_exchange_rates_s))

def total_tritium_charge_exchange_rate_s(attenuation: NeutralBeamAttenuation) -> float:
    """Total charge exchange event rate on thermal tritium"""
    return float(np.sum(attenuation.cell_tritium_charge_exchange_rates_s))

def total_other_loss_rate_s(attenuation: NeutralBeamAttenuation) -> float:
    """Total other primary neutral loss rate"""
    return float(np.sum(attenuation.cell_other_loss_rates_s))

def total_fast_ion_birth_rate_s(attenuation: NeutralBeamAttenuation) -> float:
    """Fast ion birth from primary ionization plus primary charge exchange"""
    return total_ionization_rate_s(attenuation) + total_charge_exchange_rate_s(attenuation)

def total_net_fueling_rate_s(attenuation: NeutralBeamAttenuation) -> float:
    """Net plasma fueling from impact ionization only"""
    return total_ionization_rate_s(attenuation)

def shine_through_rate_s(attenuation: NeutralBeamAttenuation) -> float:
    """Primary neutral rate exiting the final path cell"""
    return float(attenuation.neutral_rate_faces_s[-1])

def shine_through_fraction(attenuation: NeutralBeamAttenuation) -> float:
    """Fraction of primary neutral rate exiting the retained beam path"""
    if attenuation.input_neutral_rate_s <= 0.0:
        return 0.0
    return shine_through_rate_s(attenuation) / float(attenuation.input_neutral_rate_s)

def total_attenuated_fraction(attenuation: NeutralBeamAttenuation) -> float:
    """Fraction of primary neutral rate removed in the path domain"""
    return 1.0 - shine_through_fraction(attenuation)

def channel_fraction_of_input(attenuation: NeutralBeamAttenuation, cell_rates_s: ArrayLike) -> float:
    """Integrated channel rate divided by input neutral rate"""
    if attenuation.input_neutral_rate_s <= 0.0:
        return 0.0
    rates = np.asarray(cell_rates_s, dtype=float)
    return float(np.sum(rates) / attenuation.input_neutral_rate_s)

def ionization_fraction_of_input(attenuation: NeutralBeamAttenuation) -> float:
    """Impact ionization fraction of the input primary neutral rate"""
    return channel_fraction_of_input(attenuation, attenuation.cell_ionization_rates_s)

def charge_exchange_fraction_of_input(attenuation: NeutralBeamAttenuation) -> float:
    """Charge exchange fraction of the input primary neutral rate"""
    return channel_fraction_of_input(attenuation, attenuation.cell_charge_exchange_rates_s)

def other_loss_fraction_of_input(attenuation: NeutralBeamAttenuation) -> float:
    """Other loss fraction of the input primary neutral rate"""
    return channel_fraction_of_input(attenuation, attenuation.cell_other_loss_rates_s)

def charge_exchange_fraction_of_primary_interactions(attenuation: NeutralBeamAttenuation) -> float:
    """CX fraction among ionization + CX primary interactions"""
    cx = total_charge_exchange_rate_s(attenuation)
    ion = total_ionization_rate_s(attenuation)
    denom = cx + ion
    return 0.0 if denom <= 0.0 else float(cx / denom)

def integrate_rate_density_over_volume_s(rate_density_m3_s: ArrayLike, attenuation: NeutralBeamAttenuation) -> float:
    """Integral of a channel rate density over attenuation cell volumes"""
    density = np.asarray(rate_density_m3_s, dtype=float)
    return float(np.sum(density * attenuation.cell_volumes_m3))

def attenuation_conservation_residual_s(attenuation: NeutralBeamAttenuation) -> float:
    """Input - shine through - sum(all channel losses)"""
    removed = total_ionization_rate_s(attenuation) + total_charge_exchange_rate_s(attenuation) + total_other_loss_rate_s(attenuation)
    return float(attenuation.input_neutral_rate_s - shine_through_rate_s(attenuation) - removed)

def _component_rate_array(rate_per_m: ArrayLike, component_index: int) -> np.ndarray:
    """Return cell rates for one component from scalar, 1D, or 2D input"""
    rates = np.asarray(rate_per_m, dtype=float)
    if rates.ndim == 0:
        return rates
    if rates.ndim == 1:
        return rates
    if rates.ndim == 2:
        return rates[component_index]
    
    raise ValueError("component attenuation rates must be scalar, 1D, or component-by-cell 2D arrays")

def multi_energy_attenuation_from_power_fractions_on_path(total_power_W: float, component_energies_J: ArrayLike, power_fractions: ArrayLike, path_geometry: BeamPathOnAxisymmetricGrid, ionization_rate_per_m: ArrayLike, charge_exchange_rate_per_m: ArrayLike, other_loss_rate_per_m: ArrayLike | float = 0.0, *, electron_ionization_rate_per_m: ArrayLike, deuterium_ionization_rate_per_m: ArrayLike, tritium_ionization_rate_per_m: ArrayLike, deuterium_charge_exchange_rate_per_m: ArrayLike, tritium_charge_exchange_rate_per_m: ArrayLike) -> MultiEnergyNeutralBeamAttenuation:
    """Attenuate multiple energy components of one beam along a straight beam path"""
    energies = _nonnegative_array("component_energies_J", component_energies_J)
    if np.any(energies <= 0.0):
        raise ValueError("component_energies_J must be positive")
    powers = component_powers_from_fractions_W(float(total_power_W), power_fractions)
    if powers.shape != energies.shape:
        raise ValueError("power_fractions must match component_energies_J")
    input_rates = powers / energies
    components = []
    for i in range(energies.size):
        components.append(attenuate_neutral_beam_on_path_grid(input_neutral_rate_s=float(input_rates[i]), path_geometry=path_geometry, ionization_rate_per_m=_component_rate_array(ionization_rate_per_m, i), charge_exchange_rate_per_m=_component_rate_array(charge_exchange_rate_per_m, i), other_loss_rate_per_m=_component_rate_array(other_loss_rate_per_m, i), electron_ionization_rate_per_m=_component_rate_array(electron_ionization_rate_per_m, i), deuterium_ionization_rate_per_m=_component_rate_array(deuterium_ionization_rate_per_m, i), tritium_ionization_rate_per_m=_component_rate_array(tritium_ionization_rate_per_m, i), deuterium_charge_exchange_rate_per_m=_component_rate_array(deuterium_charge_exchange_rate_per_m, i), tritium_charge_exchange_rate_per_m=_component_rate_array(tritium_charge_exchange_rate_per_m, i)))
    total_fast = np.sum([component.fast_ion_birth_rate_density_m3_s for component in components], axis=0)
    total_fuel = np.sum([component.net_fueling_rate_density_m3_s for component in components], axis=0)
    total_ion = np.sum([component.ionization_rate_density_m3_s for component in components], axis=0)
    total_cx = np.sum([component.charge_exchange_rate_density_m3_s for component in components], axis=0)
    total_other = np.sum([component.other_loss_rate_density_m3_s for component in components], axis=0)
    total_electron_ionization = np.sum([component.electron_ionization_rate_density_m3_s for component in components], axis=0)
    total_deuterium_ionization = np.sum([component.deuterium_ionization_rate_density_m3_s for component in components], axis=0)
    total_tritium_ionization = np.sum([component.tritium_ionization_rate_density_m3_s for component in components], axis=0)
    total_deuterium_charge_exchange = np.sum([component.deuterium_charge_exchange_rate_density_m3_s for component in components], axis=0)
    total_tritium_charge_exchange = np.sum([component.tritium_charge_exchange_rate_density_m3_s for component in components], axis=0)

    return MultiEnergyNeutralBeamAttenuation(
        component_results=tuple(components),
        component_energies_J=energies,
        component_powers_W=powers,
        component_input_rates_s=input_rates,
        total_input_power_W=float(total_power_W),
        total_input_neutral_rate_s=float(np.sum(input_rates)),
        total_fast_ion_birth_rate_density_m3_s=total_fast,
        total_net_fueling_rate_density_m3_s=total_fuel,
        total_ionization_rate_density_m3_s=total_ion,
        total_charge_exchange_rate_density_m3_s=total_cx,
        total_other_loss_rate_density_m3_s=total_other,
        total_electron_ionization_rate_density_m3_s=total_electron_ionization,
        total_deuterium_ionization_rate_density_m3_s=total_deuterium_ionization,
        total_tritium_ionization_rate_density_m3_s=total_tritium_ionization,
        total_deuterium_charge_exchange_rate_density_m3_s=total_deuterium_charge_exchange,
        total_tritium_charge_exchange_rate_density_m3_s=total_tritium_charge_exchange,
    )

def total_fast_ion_birth_rate_multi_s(multi_result: MultiEnergyNeutralBeamAttenuation) -> float:
    """Total fast ion birth rate over all components"""
    return float(np.sum([total_fast_ion_birth_rate_s(component) for component in multi_result.component_results]))

def total_net_fueling_rate_multi_s(multi_result: MultiEnergyNeutralBeamAttenuation) -> float:
    """Total net fueling rate over all components"""
    return float(np.sum([total_net_fueling_rate_s(component) for component in multi_result.component_results]))

def shine_through_power_W(multi_result: MultiEnergyNeutralBeamAttenuation) -> float:
    """Primary neutral beam power exiting the retained beam path"""
    return float(np.sum([shine_through_rate_s(component) * multi_result.component_energies_J[i] for i, component in enumerate(multi_result.component_results)]))

def shine_through_power_fraction(multi_result: MultiEnergyNeutralBeamAttenuation) -> float:
    """Primary neutral beam power shine through fraction"""
    if multi_result.total_input_power_W <= 0.0:
        return 0.0
    return shine_through_power_W(multi_result) / multi_result.total_input_power_W

def input_rate_from_power_W(power_W: float, energy_J: float) -> float:
    """Return Ndot = P/E for beam attenuation inputs"""
    return float(particle_birth_rate_from_power_s(power_W=power_W, energy_J=energy_J))