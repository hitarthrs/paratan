"""
Convert geometry linked neutral beam attenuation into FBIS speed and Lambda source arrays

The selected axial attenuation rate density supplies q of z
Beam energy fixes the source speed and the local magnetic field with injection pitch gives Lambda birth equal to sin squared pitch divided by B_tilde
Each axial source is normalized with the FBIS velocity cell measure so integrating over speed and Lambda recovers q of z
"""
from __future__ import annotations
from dataclasses import dataclass
import numpy as np
from numpy.typing import ArrayLike
from source_model_revamp.beam.beam_attenuation import MultiEnergyNeutralBeamAttenuation, NeutralBeamAttenuation, attenuation_channel_rate_density, shine_through_power_W, total_charge_exchange_rate_s, total_deuterium_charge_exchange_rate_s, total_deuterium_ionization_rate_s, total_electron_ionization_rate_s, total_ionization_rate_s, total_net_fueling_rate_s, total_tritium_charge_exchange_rate_s, total_tritium_ionization_rate_s
from source_model_revamp.fbis.beam_source_definition import beam_lambda_birth, beam_speed_from_energy_m_s
from source_model_revamp.fbis.beam_source_distribution import  integrate_lambda_source_density, monoenergetic_lambda_source_density
from source_model_revamp.fbis.mirror_invariant_grid import LambdaGrid
from source_model_revamp.fbis.velocity_space_grid import SpeedGrid

@dataclass(frozen=True)
class BeamSpatialComponentSource:
    """
    Spatial and invariant velocity source for one attenuated beam energy component
    
    Axial profiles have one value per attenuation cell
    Lambda_birth_profile follows the representative injection pitch through the local B_tilde profile
    Lambda_birth_sample_profiles and pitch_angle_weights represent an optional finite pitch distribution
    source_v_lambda_per_s has shape axial cell by speed cell by Lambda cell when a Lambda grid is requested
    """
    energy_J: float
    speed_m_s: float
    power_W: float
    birth_rate_s: float
    pitch_angle_rad: float
    pitch_angle_samples_rad: np.ndarray
    pitch_angle_weights: np.ndarray
    B_tilde_birth_profile: np.ndarray
    Lambda_birth_profile: np.ndarray
    Lambda_birth_sample_profiles: np.ndarray
    axial_birth_rate_density_m3_s: np.ndarray
    axial_power_density_W_m3: np.ndarray
    source_v_lambda_per_s: np.ndarray | None

    @property
    def birth_pitch_angle_profile_rad(self) -> np.ndarray:
        """Return the representative injection pitch repeated on the axial birth grid"""
        return np.full(self.B_tilde_birth_profile.shape, self.pitch_angle_rad, dtype=float)

@dataclass(frozen=True)
class AttenuatedBeamComponentSource:
    """
    Selected kinetic source channel for one attenuated beam energy component
    
    birth_rate_s and deposited_birth_power_W refer to the selected channel stored in channel
    The default production channel is fast ion birth
    Separate attenuation totals retain net fueling ionization and target resolved charge exchange bookkeeping
    """
    attenuation: NeutralBeamAttenuation
    spatial_source: BeamSpatialComponentSource
    channel: str
    input_neutral_rate_s: float
    input_power_W: float
    birth_rate_s: float
    deposited_birth_power_W: float
    net_fueling_rate_s: float
    charge_exchange_rate_s: float
    ionization_rate_s: float
    electron_ionization_rate_s: float
    deuterium_ionization_rate_s: float
    tritium_ionization_rate_s: float
    deuterium_charge_exchange_rate_s: float
    tritium_charge_exchange_rate_s: float

@dataclass(frozen=True)
class AttenuatedMultiEnergyBeamSource:
    """
    Combined selected source channel for all energy components of one beam or compatible beam group
    
    Component sources retain their individual beam energies
    Total axial profiles are sums over components and total_source_v_lambda_per_s is present only when a Lambda grid was built
    """
    component_sources: tuple[AttenuatedBeamComponentSource, ...]
    channel: str
    total_input_power_W: float
    total_input_neutral_rate_s: float
    total_birth_rate_s: float
    total_deposited_birth_power_W: float
    total_net_fueling_rate_s: float
    total_charge_exchange_rate_s: float
    total_ionization_rate_s: float
    total_electron_ionization_rate_s: float
    total_deuterium_ionization_rate_s: float
    total_tritium_ionization_rate_s: float
    total_deuterium_charge_exchange_rate_s: float
    total_tritium_charge_exchange_rate_s: float
    shine_through_power_W: float
    total_axial_birth_rate_density_m3_s: np.ndarray
    total_axial_birth_power_density_W_m3: np.ndarray
    total_source_v_lambda_per_s: np.ndarray | None

def _finite_1d_array(name: str, values: ArrayLike) -> np.ndarray:
    arr = np.asarray(values, dtype=float)
    if arr.ndim != 1:
        raise ValueError(f"{name} must be 1D")
    if not np.all(np.isfinite(arr)):
        raise ValueError(f"{name} must contain only finite values")
    
    return arr

def _positive_1d_array(name: str, values: ArrayLike) -> np.ndarray:
    arr = _finite_1d_array(name, values)
    if np.any(arr <= 0.0):
        raise ValueError(f"{name} must contain only positive values")
    
    return arr

def channel_rate_density_from_attenuation(attenuation: NeutralBeamAttenuation, channel: str = "fast_ion_birth") -> np.ndarray:
    """Select a supported attenuation rate density channel"""
    return attenuation_channel_rate_density(attenuation, channel=channel)

def integrate_rate_density_over_attenuation_volume_s(rate_density_m3_s: ArrayLike, attenuation: NeutralBeamAttenuation) -> float:
    """Integral q(z) dV using the attenuation cell volumes"""
    density = np.asarray(rate_density_m3_s, dtype=float)
    if density.shape != attenuation.cell_volumes_m3.shape:
        raise ValueError("rate density and attenuation cell volumes must have the same shape")
    
    return float(np.sum(density * attenuation.cell_volumes_m3))

def pitch_angle_profile_from_representative_pitch(pitch_angle_rad: float, B_tilde_birth_profile: ArrayLike) -> np.ndarray:
    """Repeat a representative on axis pitch on the local magnetic field grid"""
    B_birth = _positive_1d_array("B_tilde_birth_profile", B_tilde_birth_profile)
    pitch = float(pitch_angle_rad)
    if not np.isfinite(pitch):
        raise ValueError("pitch_angle_rad must be finite")

    return np.full(B_birth.shape, pitch, dtype=float)

def lambda_birth_profile_from_pitch_angle(pitch_angle_rad: float, B_tilde_birth_profile: ArrayLike) -> np.ndarray:
    """Λ_b(z) = sin^2(θ_b) / B_tilde_birth(z)"""
    B_birth = _positive_1d_array("B_tilde_birth_profile", B_tilde_birth_profile)
    pitch_profile = pitch_angle_profile_from_representative_pitch(pitch_angle_rad, B_birth)
    Lambda = np.asarray(beam_lambda_birth(pitch_angle_rad=pitch_profile, B_tilde_birth=B_birth), dtype=float,)
    if np.any(Lambda < 0.0) or np.any(Lambda > 1.0):
        raise ValueError("beam birth Lambda profile must stay within [0, 1]")
    
    return Lambda

def spatial_lambda_source_from_profile(speed_grid: SpeedGrid, lambda_grid: LambdaGrid, axial_birth_rate_density_m3_s: ArrayLike, birth_speed_m_s: float, lambda_birth_profile: ArrayLike, B_tilde_birth_profile: ArrayLike, *, lambda_birth_sample_profiles: ArrayLike | None = None, pitch_angle_weights: ArrayLike | None = None) -> np.ndarray:
    """Build S(z, v, Λ) from an attenuation birth rate density profile"""
    q_profile = _finite_1d_array("axial_birth_rate_density_m3_s", axial_birth_rate_density_m3_s)
    Lambda_birth = _finite_1d_array("lambda_birth_profile", lambda_birth_profile)
    B_birth = _positive_1d_array("B_tilde_birth_profile", B_tilde_birth_profile)
    if not (q_profile.shape == Lambda_birth.shape == B_birth.shape):
        raise ValueError("q_profile, lambda_birth_profile, and B_tilde_birth_profile must have the same shape")
    if np.any(q_profile < 0.0):
        raise ValueError("axial_birth_rate_density_m3_s must be nonnegative")
    if lambda_birth_sample_profiles is None:
        sample_profiles = Lambda_birth[None, :]
        sample_weights = np.ones(1, dtype=float)
    else:
        sample_profiles = np.asarray(lambda_birth_sample_profiles, dtype=float)
        sample_weights = np.asarray(pitch_angle_weights, dtype=float)
        if sample_profiles.ndim != 2 or sample_profiles.shape[1] != q_profile.size:
            raise ValueError("lambda_birth_sample_profiles must have shape (n_pitch, n_z)")
        if sample_weights.shape != (sample_profiles.shape[0],):
            raise ValueError("pitch_angle_weights must have one value per pitch sample")
        if np.any(~np.isfinite(sample_profiles)) or np.any((sample_profiles < 0.0) | (sample_profiles > 1.0)):
            raise ValueError("sampled Lambda birth profiles must be finite and inside [0, 1]")
        if np.any(~np.isfinite(sample_weights)) or np.any(sample_weights < 0.0):
            raise ValueError("pitch angle weights must be finite and nonnegative")
        weight_sum = float(np.sum(sample_weights))
        if weight_sum <= 0.0:
            raise ValueError("pitch angle weights must have positive sum")
        sample_weights = sample_weights / weight_sum
    source = np.zeros((q_profile.size, speed_grid.centers_m_s.size, lambda_grid.centers.size), dtype=float)
    for i_z in range(q_profile.size):
        for sample_index, sample_weight in enumerate(sample_weights):
            source[i_z] += float(sample_weight) * monoenergetic_lambda_source_density(
                speed_grid=speed_grid,
                lambda_grid=lambda_grid,
                birth_rate_density_m3_s=float(q_profile[i_z]),
                birth_speed_m_s=float(birth_speed_m_s),
                lambda_birth=float(sample_profiles[sample_index, i_z]),
                B_tilde_birth=float(B_birth[i_z]),
            )

    return source

def integrate_spatial_lambda_source_with_attenuation_s(speed_grid: SpeedGrid, lambda_grid: LambdaGrid, spatial_source_z_v_lambda_per_s: ArrayLike, attenuation: NeutralBeamAttenuation, B_tilde_birth_profile: ArrayLike) -> float:
    """
    Integrate an axial speed Lambda source over velocity space and attenuation cell volumes
    
    Each axial cell uses its local B_tilde birth value for the Lambda pitch measure
    """
    source = np.asarray(spatial_source_z_v_lambda_per_s, dtype=float)
    B_birth = _positive_1d_array("B_tilde_birth_profile", B_tilde_birth_profile)
    if source.shape[0] != attenuation.cell_volumes_m3.size or B_birth.shape != attenuation.cell_volumes_m3.shape:
        raise ValueError("spatial source, cell volumes, and B_tilde_birth_profile have inconsistent axial sizes")
    total = 0.0
    for i_z in range(source.shape[0]):
        local_rate_density = integrate_lambda_source_density(speed_grid=speed_grid, lambda_grid=lambda_grid, source_v_lambda_per_s=source[i_z], B_tilde_birth=float(B_birth[i_z]))
        total += local_rate_density * attenuation.cell_volumes_m3[i_z]

    return float(total)

def attenuated_component_source(attenuation: NeutralBeamAttenuation, beam_energy_J: float, particle_mass_kg: float, pitch_angle_rad: float, speed_grid: SpeedGrid, lambda_grid: LambdaGrid | None = None, channel: str = "fast_ion_birth", *, pitch_angle_samples_rad: ArrayLike | None = None, pitch_angle_weights: ArrayLike | None = None) -> AttenuatedBeamComponentSource:
    """
    Build one FBIS source component from an attenuation result
    
    The selected attenuation channel sets the represented axial rate density
    Beam energy sets source speed and represented power
    The representative pitch and optional finite pitch quadrature are converted to local Lambda birth profiles
    When lambda_grid is None the spatial profiles are built without a velocity space source array
    """
    q_profile = np.asarray(channel_rate_density_from_attenuation(attenuation, channel=channel), dtype=float)
    birth_rate = integrate_rate_density_over_attenuation_volume_s(q_profile, attenuation)
    birth_power = birth_rate * float(beam_energy_J)
    speed = float(beam_speed_from_energy_m_s(beam_energy_J=beam_energy_J, fast_ion_mass_kg=particle_mass_kg))
    power_density = q_profile * float(beam_energy_J)
    B_birth = np.asarray(attenuation.B_tilde_centers, dtype=float)
    Lambda_birth = lambda_birth_profile_from_pitch_angle(pitch_angle_rad, B_birth)
    if pitch_angle_samples_rad is None:
        pitch_samples = np.asarray([pitch_angle_rad], dtype=float)
        pitch_weights = np.ones(1, dtype=float)
    else:
        pitch_samples = _finite_1d_array("pitch_angle_samples_rad", pitch_angle_samples_rad)
        pitch_weights = _finite_1d_array("pitch_angle_weights", pitch_angle_weights)
        if pitch_samples.shape != pitch_weights.shape:
            raise ValueError("pitch angle samples and weights must have the same shape")
        if np.any((pitch_samples < 0.0) | (pitch_samples > 0.5 * np.pi)):
            raise ValueError("pitch angle samples must lie inside [0, pi/2]")
        if np.any(pitch_weights < 0.0) or float(np.sum(pitch_weights)) <= 0.0:
            raise ValueError("pitch angle weights must be nonnegative with positive sum")
        pitch_weights = pitch_weights / np.sum(pitch_weights)
    Lambda_samples = np.asarray([lambda_birth_profile_from_pitch_angle(angle, B_birth) for angle in pitch_samples], dtype=float)
    source_lambda = None
    if lambda_grid is not None:
        source_lambda = spatial_lambda_source_from_profile(
            speed_grid=speed_grid,
            lambda_grid=lambda_grid,
            axial_birth_rate_density_m3_s=q_profile,
            birth_speed_m_s=speed,
            lambda_birth_profile=Lambda_birth,
            B_tilde_birth_profile=B_birth,
            lambda_birth_sample_profiles=Lambda_samples,
            pitch_angle_weights=pitch_weights,
        )
    spatial_source = BeamSpatialComponentSource(
        energy_J=float(beam_energy_J),
        speed_m_s=speed,
        power_W=float(birth_power),
        birth_rate_s=float(birth_rate),
        pitch_angle_rad=float(pitch_angle_rad),
        pitch_angle_samples_rad=pitch_samples,
        pitch_angle_weights=pitch_weights,
        B_tilde_birth_profile=B_birth,
        Lambda_birth_profile=Lambda_birth,
        Lambda_birth_sample_profiles=Lambda_samples,
        axial_birth_rate_density_m3_s=q_profile,
        axial_power_density_W_m3=power_density,
        source_v_lambda_per_s=source_lambda,
    )

    return AttenuatedBeamComponentSource(
        attenuation=attenuation,
        spatial_source=spatial_source,
        channel=str(channel).strip().lower(),
        input_neutral_rate_s=float(attenuation.input_neutral_rate_s),
        input_power_W=float(attenuation.input_neutral_rate_s * beam_energy_J),
        birth_rate_s=float(birth_rate),
        deposited_birth_power_W=float(birth_power),
        net_fueling_rate_s=total_net_fueling_rate_s(attenuation),
        charge_exchange_rate_s=total_charge_exchange_rate_s(attenuation),
        ionization_rate_s=total_ionization_rate_s(attenuation),
        electron_ionization_rate_s=total_electron_ionization_rate_s(attenuation),
        deuterium_ionization_rate_s=total_deuterium_ionization_rate_s(attenuation),
        tritium_ionization_rate_s=total_tritium_ionization_rate_s(attenuation),
        deuterium_charge_exchange_rate_s=total_deuterium_charge_exchange_rate_s(attenuation),
        tritium_charge_exchange_rate_s=total_tritium_charge_exchange_rate_s(attenuation),
    )

def attenuated_multi_energy_source(attenuation: MultiEnergyNeutralBeamAttenuation, particle_mass_kg: float, pitch_angle_rad: float, speed_grid: SpeedGrid, lambda_grid: LambdaGrid | None = None, channel: str = "fast_ion_birth", *, pitch_angle_samples_rad: ArrayLike | None = None, pitch_angle_weights: ArrayLike | None = None) -> AttenuatedMultiEnergyBeamSource:
    """Build and sum FBIS source representations for every energy component of one attenuated beam"""
    components = tuple(attenuated_component_source(
            attenuation=component,
            beam_energy_J=float(attenuation.component_energies_J[i]),
            particle_mass_kg=particle_mass_kg,
            pitch_angle_rad=pitch_angle_rad,
            speed_grid=speed_grid,
            lambda_grid=lambda_grid,
            channel=channel,
            pitch_angle_samples_rad=pitch_angle_samples_rad,
            pitch_angle_weights=pitch_angle_weights,
        ) for i, component in enumerate(attenuation.component_results)
    )
    total_q = np.sum([component.spatial_source.axial_birth_rate_density_m3_s for component in components], axis=0)
    total_power_density = np.sum([component.spatial_source.axial_power_density_W_m3 for component in components],axis=0,)
    total_lambda = None
    if lambda_grid is not None:
        total_lambda = np.sum([component.spatial_source.source_v_lambda_per_s for component in components], axis=0,)

    return AttenuatedMultiEnergyBeamSource(
        component_sources=components,
        channel=str(channel).strip().lower(),
        total_input_power_W=float(attenuation.total_input_power_W),
        total_input_neutral_rate_s=float(attenuation.total_input_neutral_rate_s),
        total_birth_rate_s=float(sum(component.birth_rate_s for component in components)),
        total_deposited_birth_power_W=float(sum(component.deposited_birth_power_W for component in components)),
        total_net_fueling_rate_s=float(sum(component.net_fueling_rate_s for component in components)),
        total_charge_exchange_rate_s=float(sum(component.charge_exchange_rate_s for component in components)),
        total_ionization_rate_s=float(sum(component.ionization_rate_s for component in components)),
        total_electron_ionization_rate_s=float(sum(component.electron_ionization_rate_s for component in components)),
        total_deuterium_ionization_rate_s=float(sum(component.deuterium_ionization_rate_s for component in components)),
        total_tritium_ionization_rate_s=float(sum(component.tritium_ionization_rate_s for component in components)),
        total_deuterium_charge_exchange_rate_s=float(sum(component.deuterium_charge_exchange_rate_s for component in components)),
        total_tritium_charge_exchange_rate_s=float(sum(component.tritium_charge_exchange_rate_s for component in components)),
        shine_through_power_W=shine_through_power_W(attenuation),
        total_axial_birth_rate_density_m3_s=total_q,
        total_axial_birth_power_density_W_m3=total_power_density,
        total_source_v_lambda_per_s=total_lambda,
    )

def combine_attenuated_multi_energy_sources(sources: tuple[AttenuatedMultiEnergyBeamSource, ...]) -> AttenuatedMultiEnergyBeamSource:
    """
    Combine compatible independently attenuated beam sources
    
    All inputs must use the same selected channel and axial profile shape
    Lambda sources must either be present for every input or absent for every input
    Component source tuples are concatenated while scalar totals and axial arrays are summed
    """
    items = tuple(sources)
    if not items:
        raise ValueError("sources must contain at least one attenuated beam")
    channels = {source.channel for source in items}
    if len(channels) != 1:
        raise ValueError("all attenuated beams must use the same source channel")
    rate_profiles = [np.asarray(source.total_axial_birth_rate_density_m3_s, dtype=float) for source in items]
    power_profiles = [np.asarray(source.total_axial_birth_power_density_W_m3, dtype=float) for source in items]
    if any(profile.shape != rate_profiles[0].shape for profile in rate_profiles[1:]) or any(profile.shape != power_profiles[0].shape for profile in power_profiles[1:]):
        raise ValueError("attenuated beam axial profiles must use the same grid")
    source_arrays = [source.total_source_v_lambda_per_s for source in items]
    if any(array is None for array in source_arrays) and not all(array is None for array in source_arrays):
        raise ValueError("attenuated beams must either all provide Lambda sources or all omit them")
    total_lambda = None if source_arrays[0] is None else np.sum([np.asarray(array, dtype=float) for array in source_arrays], axis=0)

    return AttenuatedMultiEnergyBeamSource(
        component_sources=tuple(component for source in items for component in source.component_sources),
        channel=items[0].channel,
        total_input_power_W=float(sum(source.total_input_power_W for source in items)),
        total_input_neutral_rate_s=float(sum(source.total_input_neutral_rate_s for source in items)),
        total_birth_rate_s=float(sum(source.total_birth_rate_s for source in items)),
        total_deposited_birth_power_W=float(sum(source.total_deposited_birth_power_W for source in items)),
        total_net_fueling_rate_s=float(sum(source.total_net_fueling_rate_s for source in items)),
        total_charge_exchange_rate_s=float(sum(source.total_charge_exchange_rate_s for source in items)),
        total_ionization_rate_s=float(sum(source.total_ionization_rate_s for source in items)),
        total_electron_ionization_rate_s=float(sum(source.total_electron_ionization_rate_s for source in items)),
        total_deuterium_ionization_rate_s=float(sum(source.total_deuterium_ionization_rate_s for source in items)),
        total_tritium_ionization_rate_s=float(sum(source.total_tritium_ionization_rate_s for source in items)),
        total_deuterium_charge_exchange_rate_s=float(sum(source.total_deuterium_charge_exchange_rate_s for source in items)),
        total_tritium_charge_exchange_rate_s=float(sum(source.total_tritium_charge_exchange_rate_s for source in items)),
        shine_through_power_W=float(sum(source.shine_through_power_W for source in items)),
        total_axial_birth_rate_density_m3_s=np.sum(rate_profiles, axis=0),
        total_axial_birth_power_density_W_m3=np.sum(power_profiles, axis=0),
        total_source_v_lambda_per_s=total_lambda,
    )

def total_fast_birth_rate_from_attenuated_components_s(source: AttenuatedMultiEnergyBeamSource) -> float:
    """Return the selected component source rates summed over all stored components"""
    return float(sum(component.birth_rate_s for component in source.component_sources))

def total_deposited_birth_power_from_components_W(source: AttenuatedMultiEnergyBeamSource) -> float:
    """Return selected source power summed over all stored components"""
    return float(sum(component.deposited_birth_power_W for component in source.component_sources))