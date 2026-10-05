"""
Coulomb collision parameters for FBIS model

Main variables
    T_e: electron temperature as an energy in joules
    n_e: electron density in m^-3
    m_f: fast ion mass in kg
    z_f: fast ion charge number, q_f = z_f * e
    z_j: background ion charge number
    n_j: background ion density in m^-3
    m_j: background ion mass in kg
    ln_Lambda_e: Coulomb logarithm for fast-ion/electron slowing down
    ln_Lambda_i: Coulomb logarithm for fast-ion/ion drag and scattering

Reference formulas
    τ_s = 3(2π)^(3/2) * eps0^2 * T_e^(3/2) * m_f / (e^4 * z_f^2 * n_e * ln_Λ_e * sqrt(m_e))

    v_c = [ (3 sqrt(π) / 4) * (m_e * ln_Λ_i / (n_e * ln_Λ_e)) * sum_j(z_j^2 * n_j / m_j) ]^(1/3) * sqrt(2 T_e / m_e)

    E_c = 0.5 * m_f * v_c^2

    z_eff = sum_j(z_j^2 * n_j) / n_e

    β_m = z_eff * m_i / (2m_f)
"""
from __future__ import annotations
from dataclasses import dataclass
import numpy as np
from numpy.typing import ArrayLike
from source_model_revamp.constants import ELECTRON_CHARGE_C, VACUUM_PERMITTIVITY_F_PER_M, ELECTRON_MASS_KG, EV_TO_J, HBAR_J_S, KEV_TO_J, J_TO_KEV
from source_model_revamp.fbis.species import DEUTERON, TRITON, IonSpecies
from source_model_revamp.fbis.pairwise_collisions import PairCoulombLogState, PairwiseCollisionRecord, PairwiseCollisionScreeningState, PairwiseCollisionState, ReducedEq59CollisionProjection, pairwise_collision_metadata

def energy_J_from_eV(energy_eV: ArrayLike):
    """E[J] = E[eV] * e"""
    return np.asarray(energy_eV, dtype=float) * EV_TO_J

def energy_J_from_keV(energy_keV: ArrayLike):
    """E[J] = E[keV] * 1000 * e"""
    return np.asarray(energy_keV, dtype=float) * KEV_TO_J

def energy_eV_from_J(energy_J: ArrayLike):
    """E[eV] = E[J] / e"""
    return np.asarray(energy_J, dtype=float) / EV_TO_J

def energy_keV_from_J(energy_J: ArrayLike):
    """E[keV] = E[J] / (1000 * e)"""
    return np.asarray(energy_J, dtype=float) * J_TO_KEV

def thermal_speed_m_s(temperature_energy_J: ArrayLike, mass_kg: float):
    """v_th = sqrt(2T / m)"""
    T = np.asarray(temperature_energy_J, dtype=float)
    return np.sqrt(2.0 * T / mass_kg)

def electron_thermal_speed_m_s(electron_temperature_J: ArrayLike):
    """v_th,e = sqrt(2T_e / m_e)"""
    return thermal_speed_m_s(temperature_energy_J=electron_temperature_J, mass_kg=ELECTRON_MASS_KG)

def ion_charge_density_m3(ion_densities_m3: ArrayLike, ion_charge_numbers: ArrayLike):
    """sum_j z_j n_j"""
    densities = np.asarray(ion_densities_m3, dtype=float)
    charges = np.asarray(ion_charge_numbers, dtype=float)
    return np.sum(charges * densities, axis=0)

def quasineutral_electron_density_m3(ion_densities_m3: ArrayLike, ion_charge_numbers: ArrayLike):
    """n_e = sum_j z_j n_j for quasineutral ions and electrons"""
    return ion_charge_density_m3(ion_densities_m3=ion_densities_m3, ion_charge_numbers=ion_charge_numbers)

def effective_charge(ion_densities_m3: ArrayLike, ion_charge_numbers: ArrayLike, electron_density_m3: ArrayLike):
    """z_eff = sum_j z_j^2 n_j / n_e"""
    densities = np.asarray(ion_densities_m3, dtype=float)
    charges = np.asarray(ion_charge_numbers, dtype=float)
    n_e = np.asarray(electron_density_m3, dtype=float)
  
    return np.sum(charges**2 * densities, axis=0) / n_e

def z2_density_over_mass_sum(ion_densities_m3: ArrayLike, ion_charge_numbers: ArrayLike, ion_masses_kg: ArrayLike):
    """sum_j z_j^2 n_j / m_j"""
    densities = np.asarray(ion_densities_m3, dtype=float)
    charges = np.asarray(ion_charge_numbers, dtype=float)
    masses = np.asarray(ion_masses_kg, dtype=float)
  
    return np.sum(charges**2 * densities / masses, axis=0)

def _require_positive_finite_scalar(value: float | None, name: str) -> float:
    """Return a positive finite scalar and reject hidden defaults"""
    if value is None:
        raise ValueError(f"{name} must be supplied by the backend collision parameter state")
    scalar = float(value)
    if not np.isfinite(scalar) or scalar <= 0.0:
        raise ValueError(f"{name} must be positive and finite")
 
    return scalar

def _require_positive_finite_array(value: ArrayLike, name: str) -> np.ndarray:
    arr = np.asarray(value, dtype=float)
    if not np.all(np.isfinite(arr)) or np.any(arr <= 0.0):
        raise ValueError(f"{name} must contain only positive finite values")
 
    return arr

def _require_nonnegative_finite_array(value: ArrayLike, name: str) -> np.ndarray:
    arr = np.asarray(value, dtype=float)
    if not np.all(np.isfinite(arr)) or np.any(arr < 0.0):
        raise ValueError(f"{name} must contain only nonnegative finite values")
  
    return arr

def spitzer_slowing_down_time_s(electron_temperature_J: ArrayLike, electron_density_m3: ArrayLike, fast_ion_mass_kg: float, fast_ion_charge_number: float = 1.0, coulomb_log_electron: float | None = None):
    """τ_s = 3(2*π)^(3/2) * eps0^2 * T_e^(3/2) * m_f / (e^4 * z_f^2 * n_e * ln_Λ_e sqrt(m_e))"""
    T_e = _require_positive_finite_array(electron_temperature_J, "electron_temperature_J")
    n_e = _require_positive_finite_array(electron_density_m3, "electron_density_m3")
    ln_e = _require_positive_finite_scalar(coulomb_log_electron, "coulomb_log_electron")
    numerator = (3.0 * (2.0 * np.pi) ** 1.5 * VACUUM_PERMITTIVITY_F_PER_M**2 * T_e**1.5 * fast_ion_mass_kg)
    denominator = (ELECTRON_CHARGE_C**4 * fast_ion_charge_number**2 * n_e * ln_e * np.sqrt(ELECTRON_MASS_KG))

    return numerator / denominator

def critical_velocity_m_s(electron_temperature_J: ArrayLike, electron_density_m3: ArrayLike, ion_densities_m3: ArrayLike, ion_charge_numbers: ArrayLike, ion_masses_kg: ArrayLike, coulomb_log_electron: float | None = None, coulomb_log_ion: float | None = None):
    """v_c = [ (3 sqrt(π) / 4) * (m_e * ln_Λ_i / (n_e * ln_Λ_e)) * sum_j(z_j^2 * n_j / m_j) ]^(1/3) * sqrt(2 T_e / m_e)"""
    T_e = _require_positive_finite_array(electron_temperature_J, "electron_temperature_J")
    n_e = _require_positive_finite_array(electron_density_m3, "electron_density_m3")
    ln_e = _require_positive_finite_scalar(coulomb_log_electron, "coulomb_log_electron")
    ln_i = _require_positive_finite_scalar(coulomb_log_ion, "coulomb_log_ion")
    ion_sum = z2_density_over_mass_sum(ion_densities_m3=ion_densities_m3, ion_charge_numbers=ion_charge_numbers, ion_masses_kg=ion_masses_kg)
    mass_density_factor = ((3.0 * np.sqrt(np.pi) / 4.0) * (ELECTRON_MASS_KG * ln_i) / (n_e * ln_e) * ion_sum)

    return mass_density_factor ** (1.0 / 3.0) * electron_thermal_speed_m_s(T_e)

def critical_energy_J(critical_velocity_m_s: ArrayLike, fast_ion_mass_kg: float):
    """E_c = 0.5 * m_f * v_c^2"""
    v_c = np.asarray(critical_velocity_m_s, dtype=float)
  
    return 0.5 * fast_ion_mass_kg * v_c**2

def critical_energy_keV(critical_velocity_m_s: ArrayLike, fast_ion_mass_kg: float):
    """E_c[keV] = 0.5 * m_f * v_c^2 / (1000 e)"""
    return energy_keV_from_J(critical_energy_J(critical_velocity_m_s=critical_velocity_m_s, fast_ion_mass_kg=fast_ion_mass_kg,))

def beta_m(effective_charge_value: ArrayLike, background_ion_mass_kg: float, fast_ion_mass_kg: float):
    """beta_m = z_eff * m_i / (2 m_f)"""
    z_eff = np.asarray(effective_charge_value, dtype=float)
  
    return z_eff * background_ion_mass_kg / (2.0 * fast_ion_mass_kg)

def beta_m_from_ion_composition(ion_densities_m3: ArrayLike, ion_charge_numbers: ArrayLike, electron_density_m3: ArrayLike, background_ion_mass_kg: float, fast_ion_mass_kg: float):
    """beta_m = [sum_j z_j^2 n_j / n_e] m_i / (2 m_f)"""
    z_eff = effective_charge(ion_densities_m3=ion_densities_m3, ion_charge_numbers=ion_charge_numbers, electron_density_m3=electron_density_m3)
  
    return beta_m(effective_charge_value=z_eff, background_ion_mass_kg=background_ion_mass_kg, fast_ion_mass_kg=fast_ion_mass_kg)

@dataclass(frozen=True)
class CoulombLogState:
    """backend Coulomb logs for the plasma closure state"""
    model: str
    electron_fast_ion: float
    ion_fast_ion: float
    electron_electron: float
    ion_self: float
    debye_length_m: float
    electron_fast_ion_b_min_m: float
    ion_fast_ion_b_min_m: float
    electron_electron_b_min_m: float
    ion_self_b_min_m: float
    electron_density_m3: float
    ion_density_m3: float
    electron_temperature_J: float
    ion_temperature_J: float
    temperature_closure_model: str

@dataclass(frozen=True)
class FBISCollisionParameterState:
    """Derived collision parameters used by the modal FBIS backend"""
    fast_ion_species: IonSpecies
    coulomb_logs: CoulombLogState
    spitzer_slowing_down_time_s: float
    critical_velocity_m_s: float
    critical_energy_J: float
    beta_m: float
    electron_electron_collision_frequency_s: float
    pairwise_collision_state: PairwiseCollisionState
    thermal_ion_densities_m3: np.ndarray | None = None
    thermal_ion_masses_kg: np.ndarray | None = None
    thermal_ion_charge_numbers: np.ndarray | None = None
    fast_ion_density_m3: float = 0.0
    multispecies_beta_m_applicable: bool = True
    fast_ion_collision_decomposition_applicable: bool | None = None
    thermal_screening_composition_applicable: bool = True
    ion_fast_coulomb_log_applicable: bool = True
    collision_physics_applicable: bool = True
    collision_scope: str = "single_thermal_background_species"

@dataclass(frozen=True)
class CollisionModelApplicability:
    """Temperature independent collision model applicability"""
    multispecies_beta_m_applicable: bool
    fast_ion_collision_decomposition_applicable: bool | None
    thermal_screening_composition_applicable: bool
    ion_fast_coulomb_log_applicable: bool
    collision_physics_applicable: bool | None

def _nbi_supported_collision_scope(collision_scope: str) -> bool:
    """Return whether the collision state belongs to the NBI supported closure"""
    return str(collision_scope).startswith("nbi_supported_fast_ion_self_background")

def collision_model_applicability(*, ion_densities_m3: ArrayLike, ion_masses_kg: ArrayLike, fast_ion_density_present: bool | None, collision_scope: str = "single_thermal_background_species") -> CollisionModelApplicability:
    """Evaluate collision applicability without a temperature"""
    thermal_densities = _require_nonnegative_finite_array(ion_densities_m3, "ion_densities_m3")
    thermal_masses = _require_positive_finite_array(ion_masses_kg, "ion_masses_kg")
    if thermal_densities.shape != thermal_masses.shape:
        raise ValueError("thermal ion density and mass arrays must have matching shapes")
    active_masses = thermal_masses[thermal_densities > 0.0]
    if active_masses.size == 0:
        raise ValueError("at least one field ion density must be positive")
    single_field_species = np.unique(active_masses).size == 1
    pure_fast_benchmark = str(collision_scope).startswith("pure_fast_") and str(collision_scope).endswith("_benchmark")
    beta_applicable = bool(single_field_species or pure_fast_benchmark)
    ion_log_applicable = bool(single_field_species or pure_fast_benchmark)
    if fast_ion_density_present is None:
        hybrid_applicable = None
    elif not fast_ion_density_present:
        hybrid_applicable = None
    else:
        hybrid_applicable = bool(pure_fast_benchmark)
    if hybrid_applicable is None and fast_ion_density_present is None:
        production_applicable = None
    else:
        production_applicable = bool(beta_applicable and ion_log_applicable and hybrid_applicable is not False and not pure_fast_benchmark)

    return CollisionModelApplicability(
        multispecies_beta_m_applicable=beta_applicable,
        fast_ion_collision_decomposition_applicable=hybrid_applicable,
        thermal_screening_composition_applicable=True,
        ion_fast_coulomb_log_applicable=ion_log_applicable,
        collision_physics_applicable=production_applicable,
    )

def debye_length_m(electron_density_m3: ArrayLike, electron_temperature_J: ArrayLike, ion_densities_m3: ArrayLike, ion_charge_numbers: ArrayLike, ion_temperature_J: ArrayLike) -> np.ndarray:
    """Multispecies Debye length using temperatures as energies in joules"""
    n_e = _require_positive_finite_array(electron_density_m3, "electron_density_m3")
    T_e = _require_positive_finite_array(electron_temperature_J, "electron_temperature_J")
    n_i = _require_nonnegative_finite_array(ion_densities_m3, "ion_densities_m3")
    z_i = np.asarray(ion_charge_numbers, dtype=float)
    T_i = _require_positive_finite_array(ion_temperature_J, "ion_temperature_J")
    if n_i.shape != z_i.shape:
        raise ValueError("ion_densities_m3 and ion_charge_numbers must have matching shapes")
    ion_term = np.sum((z_i**2) * n_i, axis=0) / T_i
    inv_lambda2 = ELECTRON_CHARGE_C**2 / VACUUM_PERMITTIVITY_F_PER_M * (n_e / T_e + ion_term)

    return 1.0 / np.sqrt(inv_lambda2)

def _pair_b_min_m(charge_a: float, charge_b: float, reduced_mass_kg: float, relative_energy_J: float) -> tuple[float, float, float]:
    energy = _require_positive_finite_scalar(relative_energy_J, "relative_energy_J")
    mu = _require_positive_finite_scalar(reduced_mass_kg, "reduced_mass_kg")
    zprod = abs(float(charge_a) * float(charge_b))
    b90 = zprod * ELECTRON_CHARGE_C**2 / (4.0 * np.pi * VACUUM_PERMITTIVITY_F_PER_M * energy)
    bq = HBAR_J_S / np.sqrt(2.0 * mu * energy)

    return float(max(b90, bq)), float(b90), float(bq)

def coulomb_log_debye_classical_quantum(debye_length_m_value: float, b_min_m: float) -> float:
    """collision Coulomb log with a log(1 + lambda_D^2 / b_min^2) form"""
    lam = _require_positive_finite_scalar(debye_length_m_value, "debye_length_m")
    bmin = _require_positive_finite_scalar(b_min_m, "b_min_m")
    value = 0.5 * np.log1p((lam / bmin) ** 2)
    if not np.isfinite(value) or value <= 0.0:
        raise ValueError("computed Coulomb logarithm must be positive and finite")
    
    return float(value)

def build_representative_coulomb_log_state(*, electron_density_m3: float, electron_temperature_J: float, ion_densities_m3: ArrayLike, ion_charge_numbers: ArrayLike, ion_masses_kg: ArrayLike, ion_temperature_J: float, fast_ion_species: IonSpecies, representative_field_ion_mass_kg: float | None = None, representative_field_ion_charge_number: float | None = None, model: str = "closure_debye_classical_quantum") -> CoulombLogState:
    """Compute scalar Coulomb logs from the current plasma closure state"""
    n_e = _require_positive_finite_scalar(electron_density_m3, "electron_density_m3")
    T_e = _require_positive_finite_scalar(electron_temperature_J, "electron_temperature_J")
    n_i_arr = _require_nonnegative_finite_array(ion_densities_m3, "ion_densities_m3")
    z_i_arr = np.asarray(ion_charge_numbers, dtype=float)
    m_i_arr = _require_positive_finite_array(ion_masses_kg, "ion_masses_kg")
    T_i = _require_positive_finite_scalar(ion_temperature_J, "ion_temperature_J")
    if n_i_arr.shape != z_i_arr.shape or n_i_arr.shape != m_i_arr.shape:
        raise ValueError("ion density, charge, and mass arrays must have the same shape")
    ion_density_total = float(np.sum(n_i_arr))
    lam = float(debye_length_m(n_e, T_e, n_i_arr, z_i_arr, T_i))
    if ion_density_total <= 0.0:
        raise ValueError("thermal ion density must be positive for Coulomb logs")
    active = n_i_arr > 0.0
    active_masses = np.unique(m_i_arr[active])
    active_charges = np.unique(np.abs(z_i_arr[active]))
    if representative_field_ion_mass_kg is None:
        if active_masses.size != 1:
            raise ValueError("mixed field ions require an explicit supported representative species")
        representative_field_ion_mass_kg = float(active_masses[0])
    if representative_field_ion_charge_number is None:
        if active_charges.size != 1:
            raise ValueError("mixed field ion charges require an explicit representative species")
        representative_field_ion_charge_number = float(active_charges[0])
    field_ion_mass = _require_positive_finite_scalar(representative_field_ion_mass_kg, "representative_field_ion_mass_kg")
    field_ion_charge = _require_positive_finite_scalar(representative_field_ion_charge_number, "representative_field_ion_charge_number")
    reduced_e_fast = ELECTRON_MASS_KG * fast_ion_species.mass_kg / (ELECTRON_MASS_KG + fast_ion_species.mass_kg)
    reduced_i_fast = field_ion_mass * fast_ion_species.mass_kg / (field_ion_mass + fast_ion_species.mass_kg)
    reduced_ee = 0.5 * ELECTRON_MASS_KG
    reduced_ii = 0.5 * field_ion_mass
    bmin_e_fast, _, _ = _pair_b_min_m(1.0, abs(fast_ion_species.charge_number), reduced_e_fast, T_e)
    bmin_i_fast, _, _ = _pair_b_min_m(field_ion_charge, abs(fast_ion_species.charge_number), reduced_i_fast, T_i)
    bmin_ee, _, _ = _pair_b_min_m(1.0, 1.0, reduced_ee, T_e)
    bmin_ii, _, _ = _pair_b_min_m(field_ion_charge, field_ion_charge, reduced_ii, T_i)

    return CoulombLogState(
        model=model,
        electron_fast_ion=coulomb_log_debye_classical_quantum(lam, bmin_e_fast),
        ion_fast_ion=coulomb_log_debye_classical_quantum(lam, bmin_i_fast),
        electron_electron=coulomb_log_debye_classical_quantum(lam, bmin_ee),
        ion_self=coulomb_log_debye_classical_quantum(lam, bmin_ii),
        debye_length_m=lam,
        electron_fast_ion_b_min_m=bmin_e_fast,
        ion_fast_ion_b_min_m=bmin_i_fast,
        electron_electron_b_min_m=bmin_ee,
        ion_self_b_min_m=bmin_ii,
        electron_density_m3=n_e,
        ion_density_m3=ion_density_total,
        electron_temperature_J=T_e,
        ion_temperature_J=T_i,
        temperature_closure_model="fixed_plasma_closure_state_option_C",
    )

def electron_electron_collision_frequency_s(electron_density_m3: ArrayLike, electron_temperature_J: ArrayLike, coulomb_log_electron: float | None) -> np.ndarray:
    """Electron electron collision scale used by Egedal Eq 66"""
    n_e = _require_positive_finite_array(electron_density_m3, "electron_density_m3")
    T_e = _require_positive_finite_array(electron_temperature_J, "electron_temperature_J")
    ln_e = _require_positive_finite_scalar(coulomb_log_electron, "coulomb_log_electron")
    v_t = np.sqrt(2.0 * T_e / ELECTRON_MASS_KG)

    return n_e * ELECTRON_CHARGE_C**4 * ln_e / (4.0 * np.pi * VACUUM_PERMITTIVITY_F_PER_M**2 * ELECTRON_MASS_KG**2 * v_t**3)

def _hydrogen_isotope_for_mass(mass_kg: float) -> IonSpecies:
    """Return the canonical singly charged D or T species for a mass"""
    mass = _require_positive_finite_scalar(mass_kg, "field_ion_mass_kg")
    tolerance = 128.0 * np.finfo(float).eps * max(mass, DEUTERON.mass_kg, TRITON.mass_kg)
    if abs(mass - DEUTERON.mass_kg) <= tolerance:
        return DEUTERON
    if abs(mass - TRITON.mass_kg) <= tolerance:
        return TRITON
    raise ValueError("pairwise collisions support only singly charged deuterium and tritium field ions")

def _pair_coulomb_log_state(*, test_mass_kg: float, test_charge_number: float, field_mass_kg: float, field_charge_number: float, relative_energy_J: float, relative_energy_model: str, screening_length_m: float, model: str, value_override: float | None = None, active_in_current_eq59_operator: bool = False, applicability: bool = True, limitation: str | None = None) -> PairCoulombLogState:
    """Build one detailed pair Coulomb logarithm state"""
    reduced_mass = test_mass_kg * field_mass_kg / (test_mass_kg + field_mass_kg)
    b_min, b_90, b_quantum = _pair_b_min_m(test_charge_number, field_charge_number, reduced_mass, relative_energy_J)
    computed_value = coulomb_log_debye_classical_quantum(screening_length_m, b_min)
    value = computed_value if value_override is None else _require_positive_finite_scalar(value_override, "pair_coulomb_log_value")
  
    return PairCoulombLogState(model=model, value=value, screening_length_m=float(screening_length_m), reduced_mass_kg=float(reduced_mass), relative_energy_J=float(relative_energy_J), relative_energy_model=str(relative_energy_model), b_90_m=float(b_90), b_quantum_m=float(b_quantum), b_min_m=float(b_min), reference_identity="repository_47_closure_debye_classical_quantum", active_in_current_eq59_operator=bool(active_in_current_eq59_operator), applicability=bool(applicability), limitation=limitation)

def maxwellian_pair_coulomb_log_state(*, test_mass_kg: float, test_charge_number: float, test_temperature_J: float, field_mass_kg: float, field_charge_number: float, field_temperature_J: float, screening_length_m: float, model: str = "thermal_pair_debye_classical_quantum") -> PairCoulombLogState:
    """Build a thermal pair Coulomb log from the mean relative center of mass energy"""
    test_temperature = _require_positive_finite_scalar(test_temperature_J, "test_temperature_J")
    field_temperature = _require_positive_finite_scalar(field_temperature_J, "field_temperature_J")
    test_mass = _require_positive_finite_scalar(test_mass_kg, "test_mass_kg")
    field_mass = _require_positive_finite_scalar(field_mass_kg, "field_mass_kg")
    mean_relative_energy = 1.5 * (field_mass * test_temperature + test_mass * field_temperature) / (test_mass + field_mass)
   
    return _pair_coulomb_log_state(test_mass_kg=test_mass, test_charge_number=test_charge_number, field_mass_kg=field_mass, field_charge_number=field_charge_number, relative_energy_J=mean_relative_energy, relative_energy_model="mean_Maxwellian_relative_center_of_mass_energy", screening_length_m=screening_length_m, model=model)

def _legacy_fast_self_coulomb_log_state(*, fast_ion_species: IonSpecies, screening_length_m: float, legacy_shared_ion_coulomb_log: float, active_in_current_eq59_operator: bool) -> PairCoulombLogState:
    """Represent the preserved shared ion log without inventing a fast self relative energy"""
    return PairCoulombLogState(model="legacy_shared_ion_log_without_fast_self_relative_energy", value=float(legacy_shared_ion_coulomb_log), screening_length_m=float(screening_length_m), reduced_mass_kg=0.5 * fast_ion_species.mass_kg, relative_energy_J=None, relative_energy_model="unavailable_non_Maxwellian_fast_self_relative_energy", b_90_m=None, b_quantum_m=None, b_min_m=None, reference_identity="repository_47_representative_ion_log_preserved_for_Eq59", active_in_current_eq59_operator=bool(active_in_current_eq59_operator), applicability=False, limitation="Pass 5 preserves the legacy shared ion log without inventing a non Maxwellian fast self relative energy convention")

def _rosenbluth_pair_weights(*, field_density_m3: float, test_charge_number: float, field_charge_number: float, test_mass_kg: float, field_mass_kg: float, coulomb_log: float) -> tuple[float, float]:
    """Return the exact Killeen g and h potential weights"""
    density = float(field_density_m3)
    charge_ratio_squared = (float(field_charge_number) / float(test_charge_number)) ** 2
    g_weight = density * charge_ratio_squared * float(coulomb_log)
    h_weight = g_weight * (float(test_mass_kg) + float(field_mass_kg)) / float(field_mass_kg)
   
    return float(g_weight), float(h_weight)

def _legacy_critical_velocity_cubed_contribution(*, electron_density_m3: float, electron_temperature_J: float, electron_coulomb_log: float, legacy_shared_ion_coulomb_log: float, field_density_m3: float, field_charge_number: float, field_mass_kg: float) -> float:
    """Return one exact algebraic contribution to the repository 47 critical velocity cubed"""
    prefactor = (3.0 * np.sqrt(np.pi) * ELECTRON_MASS_KG / (4.0 * float(electron_density_m3) * float(electron_coulomb_log)))
    electron_speed_cubed = (2.0 * float(electron_temperature_J) / ELECTRON_MASS_KG) ** 1.5
    field_term = float(field_charge_number) ** 2 * float(field_density_m3) * float(legacy_shared_ion_coulomb_log) / float(field_mass_kg)
  
    return float(prefactor * electron_speed_cubed * field_term)

def _build_pairwise_collision_state(*, fast_ion_species: IonSpecies, electron_density_m3: float, electron_temperature_J: float, ion_densities_m3: np.ndarray, ion_charge_numbers: np.ndarray, ion_masses_kg: np.ndarray, ion_temperature_J: float, fast_ion_density_m3: float, collision_scope: str, logs: CoulombLogState, spitzer_slowing_down_time_s_value: float, critical_velocity_m_s_value: float, beta_m_value: float, representative_field_ion_mass_kg: float, applicability: CollisionModelApplicability) -> PairwiseCollisionState:
    """Build the authoritative ordered collision identity state"""
    if fast_ion_species not in (DEUTERON, TRITON) or abs(float(fast_ion_species.charge_number) - 1.0) > 128.0 * np.finfo(float).eps:
        raise ValueError("pairwise collisions support only singly charged fast deuterium and tritium")
    thermal_density_by_species = {DEUTERON.species_id: 0.0, TRITON.species_id: 0.0}
    thermal_charge_by_species = {DEUTERON.species_id: 1.0, TRITON.species_id: 1.0}
    input_species: list[IonSpecies] = []
    for density, charge, mass in zip(ion_densities_m3, ion_charge_numbers, ion_masses_kg, strict=True):
        species = _hydrogen_isotope_for_mass(float(mass))
        if abs(float(charge) - 1.0) > 128.0 * np.finfo(float).eps:
            raise ValueError("pairwise collisions require singly charged positive deuterium and tritium field ions")
        if species.species_id in (item.species_id for item in input_species):
            raise ValueError("thermal ion composition must not duplicate canonical D or T entries")
        input_species.append(species)
        thermal_density_by_species[species.species_id] = float(density)
        thermal_charge_by_species[species.species_id] = float(charge)
    pure_fast_benchmark = str(collision_scope).startswith("pure_fast_") and str(collision_scope).endswith("_benchmark")
    nbi_supported = _nbi_supported_collision_scope(collision_scope)
    if pure_fast_benchmark or nbi_supported:
        thermal_density_by_species[DEUTERON.species_id] = 0.0
        thermal_density_by_species[TRITON.species_id] = 0.0
    test_label = f"fast_{fast_ion_species.symbol}"
    pair_ids = {"electrons": f"{test_label}<-electrons", DEUTERON.species_id: f"{test_label}<-thermal_D", TRITON.species_id: f"{test_label}<-thermal_T", "fast_self": f"{test_label}<-fast_{fast_ion_species.symbol}"}
    electron_log = _pair_coulomb_log_state(test_mass_kg=fast_ion_species.mass_kg, test_charge_number=fast_ion_species.charge_number, field_mass_kg=ELECTRON_MASS_KG, field_charge_number=-1.0, relative_energy_J=electron_temperature_J, relative_energy_model="electron_temperature_energy", screening_length_m=logs.debye_length_m, model=logs.model, value_override=logs.electron_fast_ion, active_in_current_eq59_operator=True, applicability=True, limitation="The full electron g and h potentials remain outside the reduced Egedal Eq 59 ordering")
    electron_g, electron_h = _rosenbluth_pair_weights(field_density_m3=electron_density_m3, test_charge_number=fast_ion_species.charge_number, field_charge_number=-1.0, test_mass_kg=fast_ion_species.mass_kg, field_mass_kg=ELECTRON_MASS_KG, coulomb_log=electron_log.value)
    records: dict[str, PairwiseCollisionRecord] = {}
    records[pair_ids["electrons"]] = PairwiseCollisionRecord(pair_id=pair_ids["electrons"], test_species_id=fast_ion_species.species_id, field_population_id="electrons", field_species_id="electron", field_distribution_model="isotropic_Maxwellian_analytic_reduced_role", field_mass_kg=ELECTRON_MASS_KG, field_charge_number=-1.0, field_density_m3=float(electron_density_m3), field_temperature_J_or_none=float(electron_temperature_J), field_distribution_reference_or_none="analytic_Maxwellian_electron_background", coulomb_log_state=electron_log, rosenbluth_g_weight_m3=electron_g, rosenbluth_h_weight_m3=electron_h, legacy_critical_velocity_cubed_contribution_m3_s3_or_none=None, operator_role="analytic_electron_slowing_time_tau_s", active_in_current_eq59_operator=True, active_operator_terms=("analytic_electron_drag_normalization",), applicability=True, limitation="The full electron pair potential is stored but the active reduced operator uses the analytic Spitzer slowing normalization")
    legacy_contributions_by_species = {DEUTERON.species_id: 0.0, TRITON.species_id: 0.0}
    for density, charge, mass, species in zip(ion_densities_m3, ion_charge_numbers, ion_masses_kg, input_species, strict=True):
        contribution = _legacy_critical_velocity_cubed_contribution(electron_density_m3=electron_density_m3, electron_temperature_J=electron_temperature_J, electron_coulomb_log=logs.electron_fast_ion, legacy_shared_ion_coulomb_log=logs.ion_fast_ion, field_density_m3=float(density), field_charge_number=float(charge), field_mass_kg=float(mass))
        legacy_contributions_by_species[species.species_id] += contribution
    for species in (DEUTERON, TRITON):
        density = thermal_density_by_species[species.species_id]
        charge = thermal_charge_by_species[species.species_id]
        active = bool(density > 0.0 and not pure_fast_benchmark)
        pair_log = _pair_coulomb_log_state(test_mass_kg=fast_ion_species.mass_kg, test_charge_number=fast_ion_species.charge_number, field_mass_kg=species.mass_kg, field_charge_number=charge, relative_energy_J=ion_temperature_J, relative_energy_model="maxwellian_reference_ion_temperature_energy", screening_length_m=logs.debye_length_m, model=logs.model, active_in_current_eq59_operator=active, applicability=True, limitation=None)
        g_weight, h_weight = _rosenbluth_pair_weights(field_density_m3=density, test_charge_number=fast_ion_species.charge_number, field_charge_number=charge, test_mass_kg=fast_ion_species.mass_kg, field_mass_kg=species.mass_kg, coulomb_log=pair_log.value)
        contribution = None if pure_fast_benchmark else legacy_contributions_by_species[species.species_id]
        records[pair_ids[species.species_id]] = PairwiseCollisionRecord(pair_id=pair_ids[species.species_id], test_species_id=fast_ion_species.species_id, field_population_id=f"thermal_{species.species_id}", field_species_id=species.species_id, field_distribution_model="isotropic_Maxwellian_reference_species", field_mass_kg=species.mass_kg, field_charge_number=charge, field_density_m3=density, field_temperature_J_or_none=float(ion_temperature_J), field_distribution_reference_or_none=f"maxwellian_reference_{species.symbol}_Maxwellian", coulomb_log_state=pair_log, rosenbluth_g_weight_m3=g_weight, rosenbluth_h_weight_m3=h_weight, legacy_critical_velocity_cubed_contribution_m3_s3_or_none=contribution, operator_role="pairwise_Maxwellian_ion_collision_terms", active_in_current_eq59_operator=active, active_operator_terms=("ion_drag", "energy_diffusion", "pitch_scattering"), applicability=True, limitation=None)
    fast_density = float(fast_ion_density_m3)
    fast_self_active = fast_density > 0.0
    fast_self_log = _legacy_fast_self_coulomb_log_state(fast_ion_species=fast_ion_species, screening_length_m=logs.debye_length_m, legacy_shared_ion_coulomb_log=logs.ion_fast_ion, active_in_current_eq59_operator=fast_self_active)
    fast_self_g, fast_self_h = _rosenbluth_pair_weights(field_density_m3=fast_density, test_charge_number=fast_ion_species.charge_number, field_charge_number=fast_ion_species.charge_number, test_mass_kg=fast_ion_species.mass_kg, field_mass_kg=fast_ion_species.mass_kg, coulomb_log=fast_self_log.value)
    fast_self_contribution = float(sum(legacy_contributions_by_species.values())) if pure_fast_benchmark else None
    records[pair_ids["fast_self"]] = PairwiseCollisionRecord(pair_id=pair_ids["fast_self"], test_species_id=fast_ion_species.species_id, field_population_id="fast_self", field_species_id=fast_ion_species.species_id, field_distribution_model="non_Maxwellian_species_first_mode_Rosenbluth_state", field_mass_kg=fast_ion_species.mass_kg, field_charge_number=fast_ion_species.charge_number, field_density_m3=fast_density, field_temperature_J_or_none=None, field_distribution_reference_or_none=f"{fast_ion_species.species_id}:modal_eq59_first_mode_rosenbluth_state", coulomb_log_state=fast_self_log, rosenbluth_g_weight_m3=fast_self_g, rosenbluth_h_weight_m3=fast_self_h, legacy_critical_velocity_cubed_contribution_m3_s3_or_none=fast_self_contribution, operator_role="pairwise_non_Maxwellian_fast_self_collision_terms", active_in_current_eq59_operator=fast_self_active, active_operator_terms=("ion_drag", "energy_diffusion", "pitch_scattering"), applicability=False, limitation="The shared legacy ion log is active without a non Maxwellian fast self relative energy prescription")
    representative_species = _hydrogen_isotope_for_mass(representative_field_ion_mass_kg)
    active_thermal_species = tuple(species.species_id for species in (DEUTERON, TRITON) if thermal_density_by_species[species.species_id] > 0.0)
    mixed_background = len(active_thermal_species) > 1
    projection_limitations: list[str] = []
    if mixed_background:
        projection_limitations.append("the representative beta_m and shared ion Coulomb log are unqualified for mixed D and T field masses")
    if applicability.fast_ion_collision_decomposition_applicable is False:
        projection_limitations.append("the reduced fast ion collision model is not a complete multispecies Landau operator" if nbi_supported else "the thermal background plus non Maxwellian fast self reduction is not a complete multispecies Landau operator")
    if pure_fast_benchmark:
        projection_limitations.append("the pure fast projection is a regression benchmark rather than a qualified stationary operating point")
    if nbi_supported:
        projection_limitations.append("the NBI supported scalar Coulomb log uses the existing representative ion reference while the active Eq 59 ion terms use kinetic self and external fast field Rosenbluth states")
    projection_limitation = "; ".join(projection_limitations) if projection_limitations else None
    projection_source_pair_ids = tuple(record.pair_id for record in records.values() if record.active_in_current_eq59_operator)
    projection = ReducedEq59CollisionProjection(spitzer_slowing_down_time_s=float(spitzer_slowing_down_time_s_value), critical_velocity_m_s=float(critical_velocity_m_s_value), beta_m=float(beta_m_value), electron_coulomb_log=float(logs.electron_fast_ion), legacy_shared_ion_coulomb_log=float(logs.ion_fast_ion), projection_model="legacy_scalar_Egedal_Eq14_reference_and_hot_Eq59_diagnostic", source_pair_ids=projection_source_pair_ids, representative_field_ion_species_id=representative_species.species_id, mixed_background_projection_used=mixed_background, qualified=bool(applicability.collision_physics_applicable), limitation=projection_limitation)
    if pure_fast_benchmark:
        screening_ids = ("fast_self",)
    elif nbi_supported:
        screening_ids = tuple(f"fast_{species.species_id}" for species in input_species)
    else:
        screening_ids = tuple(f"thermal_{species.species_id}" for species in input_species)
    screening = PairwiseCollisionScreeningState(model=logs.model, screening_length_m=float(logs.debye_length_m), electron_density_m3=float(electron_density_m3), electron_temperature_J=float(electron_temperature_J), ion_population_ids=screening_ids, ion_densities_m3=tuple(float(value) for value in ion_densities_m3), ion_charge_numbers=tuple(float(value) for value in ion_charge_numbers), ion_temperature_J=float(ion_temperature_J), fast_self_included_in_screening=bool(pure_fast_benchmark or nbi_supported), reference_identity="nbi_supported_kinetic_iteration_representative_Debye_screening_state" if nbi_supported else "repository_47_multispecies_Debye_screening_state")
    missing_cross_pair = f"fast_{fast_ion_species.symbol}<-fast_{TRITON.symbol if fast_ion_species.species_id == DEUTERON.species_id else DEUTERON.symbol}:not_in_base_pairwise_state_requires_external_system_coupling"
    temperature_convention = "kinetic_D_T_iteration_densities_with_existing_representative_ion_energy_for_Coulomb_log_and_fast_self_distribution_reference" if nbi_supported else "maxwellian_reference_D_and_T_use_the_shared_input_ion_temperature_and_fast_self_uses_a_distribution_reference"
    applicability_scope = "ordered_pair_state_active_for_NBI_supported_fast_self_and_external_fast_cross_species_Eq59_with_documented_limitations" if nbi_supported else "ordered_pair_state_active_for_Egedal_Eq59_thermal_D_thermal_T_and_fast_self_with_documented_limitations"
  
    return PairwiseCollisionState(test_species=fast_ion_species, records_by_pair_id=records, screening_state=screening, density_convention="confined_physical_volume_average_scalar_density_for_each_field_population", temperature_convention=temperature_convention, reduced_eq59_projection=projection, missing_pair_physics=(missing_cross_pair,), applicability=applicability_scope)

def rebuild_fbis_collision_parameter_state_with_fast_density(state: FBISCollisionParameterState, fast_ion_density_m3: float) -> FBISCollisionParameterState:
    """Rebuild the pairwise state after the solved fast species density changes"""
    projection = state.pairwise_collision_state.reduced_eq59_projection
    representative_species = DEUTERON if projection.representative_field_ion_species_id == DEUTERON.species_id else TRITON
    if state.thermal_ion_densities_m3 is None or state.thermal_ion_masses_kg is None or state.thermal_ion_charge_numbers is None:
        raise ValueError("collision state thermal composition is unavailable for a fast density rebuild")
    ion_densities = np.asarray(state.thermal_ion_densities_m3, dtype=float).copy()
    if str(state.collision_scope).startswith("pure_fast_") and str(state.collision_scope).endswith("_benchmark"):
        if ion_densities.size != 1:
            raise ValueError("pure fast benchmark collision states require one ion density entry")
        ion_densities = np.asarray([float(fast_ion_density_m3)], dtype=float)
    elif _nbi_supported_collision_scope(state.collision_scope):
        masses = np.asarray(state.thermal_ion_masses_kg, dtype=float)
        tolerance = 128.0 * np.finfo(float).eps * max(state.fast_ion_species.mass_kg, float(np.max(masses)))
        matches = np.flatnonzero(np.abs(masses - state.fast_ion_species.mass_kg) <= tolerance)
        if matches.size != 1:
            raise ValueError("NBI supported collision state requires one scalar field entry for the test species")
        ion_densities[matches[0]] = float(fast_ion_density_m3)
   
    return build_fbis_collision_parameter_state(electron_density_m3=state.coulomb_logs.electron_density_m3, electron_temperature_J=state.coulomb_logs.electron_temperature_J, ion_densities_m3=ion_densities, ion_charge_numbers=state.thermal_ion_charge_numbers, ion_masses_kg=state.thermal_ion_masses_kg, ion_temperature_J=state.coulomb_logs.ion_temperature_J, fast_ion_species=state.fast_ion_species, background_ion_mass_kg=representative_species.mass_kg, fast_ion_density_m3=float(fast_ion_density_m3), collision_scope=state.collision_scope)

def build_fbis_collision_parameter_state(*, electron_density_m3: float, electron_temperature_J: float, ion_densities_m3: ArrayLike, ion_charge_numbers: ArrayLike, ion_masses_kg: ArrayLike, ion_temperature_J: float, fast_ion_species: IonSpecies, background_ion_mass_kg: float | None = None, fast_ion_density_m3: float = 0.0, collision_scope: str = "single_thermal_background_species") -> FBISCollisionParameterState:
    """Build the staged backend collision state from the fixed plasma closure values"""
    thermal_densities = _require_nonnegative_finite_array(ion_densities_m3, "ion_densities_m3")
    thermal_masses = _require_positive_finite_array(ion_masses_kg, "ion_masses_kg")
    thermal_charges = np.asarray(ion_charge_numbers, dtype=float)
    if not (thermal_densities.shape == thermal_masses.shape == thermal_charges.shape):
        raise ValueError("thermal ion composition arrays must have matching shapes")
    fast_density = float(fast_ion_density_m3)
    if not np.isfinite(fast_density) or fast_density < 0.0:
        raise ValueError("fast_ion_density_m3 must be finite and nonnegative")
    active = thermal_densities > 0.0
    active_masses = thermal_masses[active]
    if active_masses.size == 0:
        raise ValueError("at least one field ion density must be positive")
    if background_ion_mass_kg is None:
        representative_mass = float(active_masses[0])
    else:
        representative_mass = float(background_ion_mass_kg)
    representative_charge = float(abs(thermal_charges[np.flatnonzero(active)[0]]))
    applicability = collision_model_applicability(ion_densities_m3=thermal_densities, ion_masses_kg=thermal_masses, fast_ion_density_present=fast_density > 0.0, collision_scope=collision_scope)
    logs = build_representative_coulomb_log_state(electron_density_m3=electron_density_m3, electron_temperature_J=electron_temperature_J, ion_densities_m3=thermal_densities, ion_charge_numbers=thermal_charges, ion_masses_kg=thermal_masses, ion_temperature_J=ion_temperature_J, fast_ion_species=fast_ion_species, representative_field_ion_mass_kg=representative_mass, representative_field_ion_charge_number=representative_charge)
    tau = spitzer_slowing_down_time_s(electron_temperature_J=electron_temperature_J, electron_density_m3=electron_density_m3, fast_ion_mass_kg=fast_ion_species.mass_kg, fast_ion_charge_number=abs(fast_ion_species.charge_number), coulomb_log_electron=logs.electron_fast_ion)
    vc = critical_velocity_m_s(electron_temperature_J=electron_temperature_J, electron_density_m3=electron_density_m3, ion_densities_m3=ion_densities_m3, ion_charge_numbers=ion_charge_numbers, ion_masses_kg=ion_masses_kg, coulomb_log_electron=logs.electron_fast_ion, coulomb_log_ion=logs.ion_fast_ion)
    zeff = effective_charge(ion_densities_m3=ion_densities_m3, ion_charge_numbers=ion_charge_numbers, electron_density_m3=electron_density_m3)
    bg_mass = representative_mass
    beta_value = float(np.asarray(beta_m(zeff, bg_mass, fast_ion_species.mass_kg)))
    tau_value = float(np.asarray(tau))
    critical_velocity_value = float(np.asarray(vc))
    electron_electron_frequency = float(np.asarray(electron_electron_collision_frequency_s(electron_density_m3, electron_temperature_J, logs.electron_electron)))
    pairwise_state = _build_pairwise_collision_state(fast_ion_species=fast_ion_species, electron_density_m3=float(electron_density_m3), electron_temperature_J=float(electron_temperature_J), ion_densities_m3=thermal_densities, ion_charge_numbers=thermal_charges, ion_masses_kg=thermal_masses, ion_temperature_J=float(ion_temperature_J), fast_ion_density_m3=fast_density, collision_scope=str(collision_scope), logs=logs, spitzer_slowing_down_time_s_value=tau_value, critical_velocity_m_s_value=critical_velocity_value, beta_m_value=beta_value, representative_field_ion_mass_kg=representative_mass, applicability=applicability)

    return FBISCollisionParameterState(
        fast_ion_species=fast_ion_species,
        coulomb_logs=logs,
        spitzer_slowing_down_time_s=tau_value,
        critical_velocity_m_s=critical_velocity_value,
        critical_energy_J=float(np.asarray(critical_energy_J(vc, fast_ion_species.mass_kg))),
        beta_m=beta_value,
        electron_electron_collision_frequency_s=electron_electron_frequency,
        pairwise_collision_state=pairwise_state,
        thermal_ion_densities_m3=thermal_densities.copy(),
        thermal_ion_masses_kg=thermal_masses.copy(),
        thermal_ion_charge_numbers=thermal_charges.copy(),
        fast_ion_density_m3=fast_density,
        multispecies_beta_m_applicable=(applicability.multispecies_beta_m_applicable),
        fast_ion_collision_decomposition_applicable=(applicability.fast_ion_collision_decomposition_applicable),
        thermal_screening_composition_applicable=(applicability.thermal_screening_composition_applicable),
        ion_fast_coulomb_log_applicable=(applicability.ion_fast_coulomb_log_applicable),
        collision_physics_applicable=bool(applicability.collision_physics_applicable),
        collision_scope=str(collision_scope),
    )

def collision_parameter_metadata(state: FBISCollisionParameterState) -> dict[str, object]:
    """JSON friendly metadata for the backend collision state"""
    logs = state.coulomb_logs
    beta_qualified = bool(state.multispecies_beta_m_applicable)
    qualified_beta_m = state.beta_m if beta_qualified else None
    provisional_beta_m = None if beta_qualified else state.beta_m
    ion_log_qualified = bool(state.ion_fast_coulomb_log_applicable)
    qualified_ion_fast_log = logs.ion_fast_ion if ion_log_qualified else None
    provisional_ion_fast_log = None if ion_log_qualified else logs.ion_fast_ion
    metadata = {
        "collision_parameter_backend_model": "fbis_collision_parameter_state",
        "collision_fast_ion_species_id": state.fast_ion_species.species_id,
        "collision_fast_ion_species_name": state.fast_ion_species.name,
        "collision_fast_ion_mass_kg": state.fast_ion_species.mass_kg,
        "collision_fast_ion_charge_number": state.fast_ion_species.charge_number,
        "collision_state_scope": ("single_fast_deuteron_representative_background" if state.fast_ion_species.species_id == "deuterium" else "single_fast_ion_representative_background"),
        "collision_scope": state.collision_scope,
        "collision_thermal_ion_densities_m3": (None if state.thermal_ion_densities_m3 is None else [float(value) for value in state.thermal_ion_densities_m3]),
        "collision_thermal_ion_masses_kg": (None if state.thermal_ion_masses_kg is None else [float(value) for value in state.thermal_ion_masses_kg]),
        "collision_ion_density_array_role": "kinetic_D_T_iteration_reference_for_scalar_Coulomb_logs" if _nbi_supported_collision_scope(state.collision_scope) else "maxwellian_reference_D_T_composition",
        "collision_fast_ion_density_m3": state.fast_ion_density_m3,
        "multispecies_beta_m_applicability": state.multispecies_beta_m_applicable,
        "fast_ion_collision_decomposition_applicability": state.fast_ion_collision_decomposition_applicable,
        "thermal_screening_composition_applicability": state.thermal_screening_composition_applicable,
        "ion_fast_coulomb_log_applicability": state.ion_fast_coulomb_log_applicable,
        "collision_physics_applicability": state.collision_physics_applicable,
        "collision_fast_ion_excluded_from_thermal_ion_temperature_array": not (str(state.collision_scope).startswith("pure_fast_") and str(state.collision_scope).endswith("_benchmark")),
        "collision_state_temperature_closure": logs.temperature_closure_model,
        "collision_state_recompute_policy": "rebuilt_from_each_outer_closure_state_and_each_fast_self_density_iterate",
        "collision_state_is_self_consistent_electron_temperature": False,
        "plasma_temperature_closure_status": logs.temperature_closure_model,
        "coulomb_log_model": logs.model,
        "coulomb_log_electron_fast_ion": logs.electron_fast_ion,
        "coulomb_log_ion_fast_ion": qualified_ion_fast_log,
        "coulomb_log_ion_fast_ion_qualified": ion_log_qualified,
        "coulomb_log_ion_fast_ion_provisional_diagnostic": provisional_ion_fast_log,
        "coulomb_log_ion_fast_ion_status": ("legacy_shared_ion_log_qualified_for_single_species_reference" if ion_log_qualified else "legacy_shared_ion_log_unqualified_for_mixed_D_T_and_retained_only_for_fast_self_and_Eq14_diagnostics"),
        "coulomb_log_ion_fast_ion_representative_assumption": ("not_used_single_thermal_species" if ion_log_qualified else "first_active_thermal_ion_species_mass_and_charge"),
        "coulomb_log_electron_electron": logs.electron_electron,
        "coulomb_log_ion_self": logs.ion_self,
        "coulomb_log_debye_length_m": logs.debye_length_m,
        "coulomb_log_electron_fast_ion_b_min_m": logs.electron_fast_ion_b_min_m,
        "coulomb_log_ion_fast_ion_b_min_m": logs.ion_fast_ion_b_min_m,
        "coulomb_log_electron_electron_b_min_m": logs.electron_electron_b_min_m,
        "coulomb_log_ion_self_b_min_m": logs.ion_self_b_min_m,
        "collision_state_electron_density_m3": logs.electron_density_m3,
        "collision_state_ion_density_m3": logs.ion_density_m3,
        "collision_state_electron_temperature_keV": energy_keV_from_J(logs.electron_temperature_J).item(),
        "collision_state_ion_temperature_keV": energy_keV_from_J(logs.ion_temperature_J).item(),
        "spitzer_slowing_down_time_s": state.spitzer_slowing_down_time_s,
        "critical_velocity_m_s": state.critical_velocity_m_s,
        "critical_energy_keV": energy_keV_from_J(state.critical_energy_J).item(),
        "beta_m": qualified_beta_m,
        "beta_m_qualified": beta_qualified,
        "beta_m_provisional_diagnostic": provisional_beta_m,
        "collision_qualified_aggregate_beta_m": qualified_beta_m,
        "collision_representative_aggregate_beta_m_diagnostic": provisional_beta_m,
        "collision_representative_aggregate_status": ("qualified_single_species_beta_m" if beta_qualified else "provisional_diagnostic_unqualified_mixed_species_representative"),
        "collision_representative_aggregate_is_provisional": not beta_qualified,
        "collision_representative_aggregate_assumption": ("not_used_single_species_beta_m" if beta_qualified else "first_active_thermal_ion_species_mass_and_charge"),
        "electron_electron_collision_frequency_s": state.electron_electron_collision_frequency_s,
        "fast_ion_self_collision_coulomb_log_source": "legacy_shared_ion_log_without_non_Maxwellian_relative_energy",
        **pairwise_collision_metadata(state.pairwise_collision_state),
    }
    if state.fast_ion_species.species_id == "deuterium":
        metadata.update({"collision_fast_D_density_m3": state.fast_ion_density_m3, "collision_fast_D_excluded_from_thermal_ion_temperature_array": not (str(state.collision_scope).startswith("pure_fast_") and str(state.collision_scope).endswith("_benchmark"))})

    return metadata
