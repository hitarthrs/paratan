"""
Neutral beam attenuation and source assembly for the pipeline

The stage resolves beam geometry, target density authority, atomic stopping rates, attenuated component deposition, and species grouped modal sources
"""
from __future__ import annotations
from typing import TYPE_CHECKING
import numpy as np
from source_model_revamp.beam.beam_attenuation import multi_energy_attenuation_from_power_fractions_on_path, shine_through_power_fraction, total_fast_ion_birth_rate_multi_s
from source_model_revamp.beam.beam_path_geometry import beam_path_on_axisymmetric_grid
from source_model_revamp.beam.hydrogenic_dt_atomic_data import COLLISIONDB_D107356_CROSS_SECTION_UNITS, COLLISIONDB_D107356_DATA_SHA256, COLLISIONDB_D107356_DOI, COLLISIONDB_D107356_ENERGY_FRAME, COLLISIONDB_D107356_ENERGY_UNITS, COLLISIONDB_D107356_MAX_KEV_PER_U, COLLISIONDB_D107356_METHOD, COLLISIONDB_D107356_MIN_KEV_PER_U, COLLISIONDB_D107356_QID, COLLISIONDB_D107356_SOURCE_PDF_SHA256, HYDROGENIC_DT_GROUND_STATE_REFERENCE, IONIZATION_HIGH_ENERGY_TRUNCATION_MAX_RELATIVE_RATE, IONIZATION_LOW_ENERGY_TRUNCATION_MAX_RELATIVE_RATE, heavy_particle_maxwellian_target_energy_moment_J, hydrogenic_dt_atomic_rates_per_m
from source_model_revamp.fbis.species import DEUTERON, TRITON, ion_species
from source_model_revamp.beam.kinetic_target_rates import kinetic_target_atomic_rates_per_m
from source_model_revamp.beam.charge_exchange_redistribution import build_charge_exchange_target_sink_states
from source_model_revamp.beam.source_from_attenuation import AttenuatedMultiEnergyBeamSource, attenuated_multi_energy_source, combine_attenuated_multi_energy_sources
from source_model_revamp.fbis.mirror_invariant_grid import LambdaGrid, uniform_lambda_grid
from source_model_revamp.fbis.velocity_space_grid import PitchGrid, SpeedGrid, uniform_pitch_grid
from source_model_revamp.fusion.reactions import KEV_TO_J
from source_model_revamp.integration.pipeline_common import _speed_grid_for_species, _volume_average_spatial_lambda_source
from source_model_revamp.integration.pipeline_types import BeamDepositionResult, BeamEnsembleResult, GeometryStageResult
from source_model_revamp.integration.config.beam import BeamConfig
from source_model_revamp.integration.config.root import SourceModelRunConfig
from source_model_revamp.integration.plasma_closure import uses_startup_seed_only
from source_model_revamp.numerical_quadrature import gauss_legendre_rule

if TYPE_CHECKING:
    from source_model_revamp.integration.modal_stage.types import OperatingPointDensityState
    from source_model_revamp.integration.pipeline_types import KineticStageResult

def beamline_axis_pitch_angle_deg(direction_unit: np.ndarray) -> float:
    """
    Return the unsigned angle between the beamline and the axial magnetic field in degrees
    
    The sign of the beamline z component does not change the pitch magnitude
    """
    direction = np.asarray(direction_unit, dtype=float)
    if direction.shape != (3,) or not np.all(np.isfinite(direction)):
        raise ValueError("beamline direction must be a finite length 3 vector in x, y, z order")
    norm = float(np.linalg.norm(direction))
    if norm <= 0.0:
        raise ValueError("beamline direction must have nonzero length")
    u_z_magnitude = float(np.clip(abs(direction[2] / norm), 0.0, 1.0))

    return float(np.rad2deg(np.arccos(u_z_magnitude)))

def _resolve_beam_pitch_deg(beam_config: BeamConfig, direction_unit: np.ndarray) -> tuple[str, float, float, float, bool, bool]:
    """
    Resolve the beam pitch from the configured model and endpoint geometry
    
    The result retains the configured and endpoint derived values together with consistency and override flags
    """
    configured_deg = float(beam_config.injection_angle_deg)
    endpoint_deg = beamline_axis_pitch_angle_deg(direction_unit)
    mismatch_deg = abs(endpoint_deg - configured_deg)
    tolerance_deg = float(beam_config.injection_angle_tolerance_deg)
    if not np.isfinite(tolerance_deg) or tolerance_deg < 0.0:
        raise ValueError("beam.injection_angle_tolerance_deg must be finite and nonnegative")
    pitch_source = str(beam_config.injection_pitch_source).strip().lower().replace("-", "_")
    geometry_consistent = bool(mismatch_deg <= tolerance_deg)
    if pitch_source == "beamline_endpoints":
        if not geometry_consistent:
            raise ValueError("beam.injection_angle_deg disagrees with the pitch derived from beamline endpoints " f"by {mismatch_deg:.12g} degrees, exceeding beam.injection_angle_tolerance_deg=" f"{tolerance_deg:.12g}; use injection_pitch_source: configured_override only for a deliberate pitch override")
        effective_deg = endpoint_deg
        override_selected = False
    elif pitch_source == "configured_override":
        effective_deg = configured_deg
        override_selected = True
    else:
        raise ValueError("beam.injection_pitch_source must be beamline_endpoints or configured_override")

    return pitch_source, endpoint_deg, configured_deg, effective_deg, geometry_consistent, override_selected

def _beam_pitch_quadrature(beam_config: BeamConfig, *, representative_pitch_rad: float, beamline_length_m: float) -> tuple[np.ndarray, np.ndarray, float, str]:
    """
    Build the finite beam pitch quadrature in radians
    
    The angular half width includes the configured divergence and the geometric footprint angle from beam radius divided by beamline length
    """
    model = str(beam_config.pitch_distribution_model).strip().lower()
    representative = float(representative_pitch_rad)
    if model == "mono_pitch_benchmark":
        return (np.asarray([representative], dtype=float), np.ones(1, dtype=float), 0.0, "explicit_mono_pitch_benchmark")
    if model != "geometry_linked_finite_width":
        raise ValueError(f"unsupported beam pitch distribution model {model!r}")
    divergence = float(np.deg2rad(beam_config.beamline.divergence_half_angle_deg))
    footprint_angle = float(np.arctan2(beam_config.beamline.radius_m, max(float(beamline_length_m), np.finfo(float).tiny)))
    half_width = max(divergence, footprint_angle)
    if half_width <= 0.0:
        raise ValueError("geometry_linked_finite_width requires nonzero beamline divergence or footprint radius")
    count = int(beam_config.pitch_quadrature_points)
    if count < 2:
        raise ValueError("finite width beam pitch integration requires at least two quadrature points")
    lower = max(0.0, representative - half_width)
    upper = min(0.5 * np.pi, representative + half_width)
    if upper <= lower:
        raise ValueError("beam pitch width interval has no physical support")
    nodes, weights = gauss_legendre_rule(count)
    samples = 0.5 * (lower + upper) + 0.5 * (upper - lower) * nodes
    normalized_weights = 0.5 * weights
    width_source = ("configured_divergence_half_angle" if divergence >= footprint_angle and divergence > 0.0 else "beam_footprint_radius_over_beamline_length")
    
    return samples, normalized_weights, half_width, width_source

def _seed_reference_support_mask(geometry: GeometryStageResult) -> np.ndarray:
    """Return confined cells where the configured startup or reference plasma profile is physically supported"""
    z = np.asarray(geometry.z_centers_m, dtype=float)
    background = geometry.background_profiles
    if background is None:
        return np.ones(z.shape, dtype=bool)

    support = np.asarray(background.support_domain.contains(z), dtype=bool)
    if support.shape != z.shape:
        raise ValueError("configured plasma support mask must match the beam axial grid")

    return support

def _seed_reference_positive_charge_profile(config: SourceModelRunConfig, geometry: GeometryStageResult) -> np.ndarray:
    """Return the configured startup or reference D plus T positive charge number density on confined cells in m⁻³"""
    z = np.asarray(geometry.z_centers_m, dtype=float)
    background = geometry.background_profiles
    if background is None:
        value = (config.plasma_closure.background_deuterium_midplane_density_m3 + config.plasma_closure.background_tritium_midplane_density_m3)
        return np.full(z.shape, float(value), dtype=float)
    profile = np.asarray(background.positive_charge_density_m3_at(z), dtype=float)
    if profile.shape != z.shape:
        raise ValueError("configured plasma positive charge must match the beam axial grid")

    return profile

def _seed_reference_species_profiles(config: SourceModelRunConfig, geometry: GeometryStageResult) -> tuple[np.ndarray, np.ndarray]:
    """Return configured startup or reference deuterium and tritium number density profiles on confined cells in m⁻³"""
    z = np.asarray(geometry.z_centers_m, dtype=float)
    background = geometry.background_profiles
    if background is None:
        support = _seed_reference_support_mask(geometry)
        deuterium = np.where(support, float(config.plasma_closure.background_deuterium_midplane_density_m3), 0.0)
        tritium = np.where(support, float(config.plasma_closure.background_tritium_midplane_density_m3), 0.0)
    else:
        deuterium = np.asarray(background.deuterium_density_m3_at(z), dtype=float)
        tritium = np.asarray(background.tritium_density_m3_at(z), dtype=float)
    if deuterium.shape != z.shape or tritium.shape != z.shape:
        raise ValueError("configured ion profiles must match the beam axial grid")
    if np.any(~np.isfinite(deuterium)) or np.any(~np.isfinite(tritium)) or np.any(deuterium < 0.0) or np.any(tritium < 0.0):
        raise ValueError("configured ion profiles must be finite and nonnegative")
    return deuterium, tritium

def _speed_grids_match(left: SpeedGrid, right: SpeedGrid) -> bool:
    """Return whether two speed grids have the same bounds, cell count, centers, edges, and widths"""
    return all(np.array_equal(getattr(left, name), getattr(right, name)) for name in ("centers_m_s", "faces_m_s", "widths_m_s", "shell_volumes_m3_s3", "face_areas_m2_s2"))

def _defined_mask(values: object, expected_shape: tuple[int, ...]) -> np.ndarray:
    """Return the finite definition mask for one explicit target density array with the expected shape"""
    raw = np.asarray(values)
    if raw.shape != expected_shape:
        raise ValueError("beam target density defined mask must match the beam axial grid")
    if np.issubdtype(raw.dtype, np.bool_):
        return raw.astype(bool, copy=True)
    try:
        numeric = np.asarray(values, dtype=float)
    except (TypeError, ValueError) as exc:
        raise ValueError("beam target density defined mask must contain booleans") from exc
    if np.any(~np.isfinite(numeric)) or np.any((numeric != 0.0) & (numeric != 1.0)):
        raise ValueError("beam target density defined mask must contain only boolean values")

    return numeric.astype(bool)

def _beam_target_density_profile(config: SourceModelRunConfig, geometry: GeometryStageResult, *, effective_path_lengths_m: np.ndarray, target_density_profile_m3: np.ndarray | None, target_density_defined_mask: np.ndarray | None, operating_point_density_state: OperatingPointDensityState | None) -> tuple[np.ndarray, str, np.ndarray, np.ndarray]:
    """
    Select the stopping target density used for beam attenuation on confined cells
    
    Benchmark closure keeps its prescribed density authority
    Startup seed closure uses the reference profile before a kinetic state exists and then uses the solved electron density from OperatingPointDensityState
    """
    expected_shape = np.asarray(geometry.z_centers_m).shape
    if target_density_profile_m3 is not None and operating_point_density_state is not None:
        raise ValueError("provide either target_density_profile_m3 or operating_point_density_state, not both")
    if operating_point_density_state is not None and target_density_defined_mask is not None:
        raise ValueError("operating_point_density_state already owns its target density definition mask")
    if target_density_profile_m3 is None and operating_point_density_state is None:
        if target_density_defined_mask is not None:
            raise ValueError("target_density_defined_mask requires an explicit target_density_profile_m3")
        benchmark_closure = str(config.kinetic_electrostatic.modal_density_closure_model).strip().lower() == "egedal_beam_plasma_quasineutral"
        if benchmark_closure:
            benchmark_density = config.plasma_closure.benchmark_prescribed_electron_density_m3
            if benchmark_density is None:
                raise ValueError("benchmark prescribed electron density is unavailable")
            density = np.full(expected_shape, float(benchmark_density), dtype=float)
            source = "benchmark_prescribed_actual_electron_density"
        else:
            density = _seed_reference_positive_charge_profile(config, geometry)
            source = "startup_or_reference_positive_charge_initial_target"
        defined = np.ones(expected_shape, dtype=bool)
    elif operating_point_density_state is not None:
        density = np.asarray(operating_point_density_state.electron_cell_density_m3, dtype=float)
        defined = _defined_mask(operating_point_density_state.profile_support_mask, expected_shape)
        state_zeta = np.asarray(operating_point_density_state.zeta_cells, dtype=float)
        geometry_zeta = np.asarray(geometry.zeta_centers, dtype=float)
        if (state_zeta.shape != geometry_zeta.shape or np.any(~np.isfinite(state_zeta)) or not np.allclose(state_zeta, geometry_zeta, rtol=0.0, atol=1.0e-12)):
            raise ValueError("operating point density state does not match the beam geometry axial grid")
        source = "solved_confined_eq70_electron_density"
    else:
        density = np.asarray(target_density_profile_m3, dtype=float)
        defined = (np.ones(expected_shape, dtype=bool) if target_density_defined_mask is None else _defined_mask(target_density_defined_mask, expected_shape))
        source = "quasineutral_total_target_density"

    if density.shape != expected_shape:
        raise ValueError("beam target density profile must match the beam axial grid")
    if np.any(~np.isfinite(density)) or np.any(density < 0.0):
        raise ValueError("beam target density profile must be finite and nonnegative")

    if uses_startup_seed_only(config.plasma_closure) and (target_density_profile_m3 is not None or operating_point_density_state is not None):
        support = np.asarray(defined, dtype=bool).copy()
    else:
        support = _seed_reference_support_mask(geometry)
    density = density.copy()
    density[~support] = 0.0
    defined = np.asarray(defined, dtype=bool)
    defined[~support] = True
    effective_path = np.asarray(effective_path_lengths_m, dtype=float)
    if effective_path.shape != expected_shape:
        raise ValueError("beam effective path lengths must match the beam axial grid")
    required = support & (effective_path > 0.0)
    undefined_required = required & ~defined
    if np.any(undefined_required):
        indices = np.flatnonzero(undefined_required).tolist()
        raise ValueError("beam target density is undefined on required beam-intersected target " f"cells with indices {indices}")

    return density, source, defined, support

def _build_beam_deposition(config: SourceModelRunConfig, geometry: GeometryStageResult, beam_config: BeamConfig, speed_grid: SpeedGrid, lambda_grid: LambdaGrid, pitch_grid: PitchGrid, *, electron_temperature_J: float, target_density_profile_m3: np.ndarray | None = None, target_density_defined_mask: np.ndarray | None = None, operating_point_density_state: OperatingPointDensityState | None = None, kinetic_target_state: KineticStageResult | None = None) -> BeamDepositionResult:
    """
    Build attenuation and fast ion births for one enabled neutral beam
    
    The beam path is intersected with the physical plasma domain, component energies are formed from `energy_keV * energy_multiplier`, atomic stopping rates are evaluated on the active target state, and attenuated births are mapped to the species speed and Λ source grid
    """
    projectile_species = ion_species(beam_config.species)
    beamline = beam_config.beamline
    path = beam_path_on_axisymmetric_grid(geometry.z_edges_m, geometry.B_tilde_centers, start_point_m=beamline.start_m, end_point_m=beamline.end_m, reference_area_m2=geometry.midplane_area_m2, plasma_radius_m=geometry.plasma_radius_m, beam_radius_m=beamline.radius_m, divergence_half_angle_rad=np.deg2rad(beamline.divergence_half_angle_deg), radial_overlap_model=beamline.radial_overlap_model)
    pitch_source, beamline_axis_angle_deg, configured_injection_angle_deg, effective_injection_angle_deg, injection_geometry_consistent, override_selected = _resolve_beam_pitch_deg(beam_config, path.direction_unit)
    injection_angle_mismatch_deg = abs(beamline_axis_angle_deg - configured_injection_angle_deg)
    physical_intersection_count = int(np.count_nonzero(path.physical_path_lengths_m > 0.0))
    effective_intersection_count = int(np.count_nonzero(path.effective_path_lengths_m > 0.0))
    if physical_intersection_count == 0:
        raise ValueError("beamline segment does not intersect the configured axial source domain")
    if effective_intersection_count == 0:
        raise ValueError("beamline segment has zero overlap with the configured plasma source volume")
    component_power_fractions = np.asarray([component.power_fraction for component in beam_config.components], dtype=float)
    component_energy_multipliers = np.asarray([component.energy_multiplier for component in beam_config.components], dtype=float)
    component_energies_J = beam_config.energy_keV * KEV_TO_J * component_energy_multipliers
    (target_density, target_density_profile_source, target_density_profile_defined_mask, target_density_profile_support_mask) = _beam_target_density_profile(config, geometry, effective_path_lengths_m=path.effective_path_lengths_m, target_density_profile_m3=target_density_profile_m3, target_density_defined_mask=target_density_defined_mask, operating_point_density_state=operating_point_density_state)
    startup_seed_only = uses_startup_seed_only(config.plasma_closure)
    # Switch heavy particle stopping from the startup seed to the solved kinetic target once it exists
    kinetic_target_active = bool(startup_seed_only and kinetic_target_state is not None)
    if kinetic_target_active:
        kinetic_rates = kinetic_target_atomic_rates_per_m(
            component_energies_J=component_energies_J,
            projectile_species=projectile_species,
            electron_density_profile_m3=target_density,
            electron_temperature_J=electron_temperature_J,
            kinetic_target_state=kinetic_target_state,
            beam_direction_unit=path.direction_unit,
            gyroangle_points=int(config.fusion.num_gyroangle_points),
        )
        rates = kinetic_rates.rates
        deuterium_target_density = kinetic_rates.deuterium_density_profile_m3
        tritium_target_density = kinetic_rates.tritium_density_profile_m3
        deuterium_target_mean_energy_by_component_z = kinetic_rates.deuterium_charge_exchange_target_mean_energy_J
        tritium_target_mean_energy_by_component_z = kinetic_rates.tritium_charge_exchange_target_mean_energy_J
        heavy_target_model = kinetic_rates.target_model
        heavy_target_gyroangle_points = kinetic_rates.gyroangle_points
    else:
        deuterium_target_density, tritium_target_density = _seed_reference_species_profiles(config, geometry)
        rates = hydrogenic_dt_atomic_rates_per_m(
            component_energies_J=component_energies_J,
            projectile_species=projectile_species,
            electron_density_profile_m3=target_density,
            deuterium_density_profile_m3=deuterium_target_density,
            tritium_density_profile_m3=tritium_target_density,
            electron_temperature_J=electron_temperature_J,
            ion_temperature_J=float(config.plasma_closure.background_ion_temperature_keV * KEV_TO_J),
            model_name=beam_config.atomic_data_model,
        )
        deuterium_target_mean_energy_by_component_z = None
        tritium_target_mean_energy_by_component_z = None
        heavy_target_model = "startup_seed_Maxwellian" if startup_seed_only else "startup_or_reference_Maxwellian"
        heavy_target_gyroangle_points = None
    attenuation = multi_energy_attenuation_from_power_fractions_on_path(total_power_W=beam_config.power_W, component_energies_J=component_energies_J, power_fractions=component_power_fractions, path_geometry=path, ionization_rate_per_m=rates.ionization_rate_per_m, charge_exchange_rate_per_m=rates.charge_exchange_rate_per_m, other_loss_rate_per_m=rates.other_loss_rate_per_m, electron_ionization_rate_per_m=rates.electron_ionization_rate_per_m, deuterium_ionization_rate_per_m=rates.deuterium_ionization_rate_per_m, tritium_ionization_rate_per_m=rates.tritium_ionization_rate_per_m, deuterium_charge_exchange_rate_per_m=rates.deuterium_charge_exchange_rate_per_m, tritium_charge_exchange_rate_per_m=rates.tritium_charge_exchange_rate_per_m)
    effective_injection_angle_rad = float(np.deg2rad(effective_injection_angle_deg))
    pitch_samples, pitch_weights, pitch_half_width_rad, pitch_width_source = _beam_pitch_quadrature(beam_config, representative_pitch_rad=effective_injection_angle_rad, beamline_length_m=path.total_length_m)
    attenuated_source = attenuated_multi_energy_source(attenuation=attenuation, particle_mass_kg=projectile_species.mass_kg, pitch_angle_rad=effective_injection_angle_rad, speed_grid=speed_grid, lambda_grid=lambda_grid, channel=beam_config.channel, pitch_angle_samples_rad=pitch_samples, pitch_angle_weights=pitch_weights)
    if attenuated_source.total_source_v_lambda_per_s is None:
        raise RuntimeError("attenuated source did not build Lambda source arrays")
    volume_averaged_source = _volume_average_spatial_lambda_source(attenuated_source.total_source_v_lambda_per_s, geometry.cell_volumes_m3)
    axial_birth_density = np.asarray(attenuated_source.total_axial_birth_rate_density_m3_s, dtype=float)
    axial_power_density = np.asarray(attenuated_source.total_axial_birth_power_density_W_m3, dtype=float)
    axial_birth_rate = axial_birth_density * geometry.cell_volumes_m3
    axial_birth_power = axial_power_density * geometry.cell_volumes_m3
    birth_profile_support = (np.asarray(path.effective_path_lengths_m, dtype=float) > 0.0 ) & np.asarray(target_density_profile_support_mask, dtype=bool)
    component_deposited_power = np.asarray([ float(component.deposited_birth_power_W) for component in attenuated_source.component_sources], dtype=float,)
    component_support = component_power_fractions > 0.0
    ion_temperature_J = float(config.plasma_closure.background_ion_temperature_keV * KEV_TO_J)
    component_deuterium_target_pumpout_rate_s = np.asarray([component.deuterium_charge_exchange_rate_s for component in attenuated_source.component_sources], dtype=float)
    component_tritium_target_pumpout_rate_s = np.asarray([component.tritium_charge_exchange_rate_s for component in attenuated_source.component_sources], dtype=float)
    if kinetic_target_active:
        component_deuterium_target_mean_energy_J = np.asarray([np.sum(component.cell_deuterium_charge_exchange_rates_s * deuterium_target_mean_energy_by_component_z[index]) / rate if rate > 0.0 else 0.0 for index, (component, rate) in enumerate(zip(attenuation.component_results, component_deuterium_target_pumpout_rate_s, strict=True))], dtype=float)
        component_tritium_target_mean_energy_J = np.asarray([np.sum(component.cell_tritium_charge_exchange_rates_s * tritium_target_mean_energy_by_component_z[index]) / rate if rate > 0.0 else 0.0 for index, (component, rate) in enumerate(zip(attenuation.component_results, component_tritium_target_pumpout_rate_s, strict=True))], dtype=float)
    else:
        component_deuterium_target_mean_energy_J = np.asarray([heavy_particle_maxwellian_target_energy_moment_J(float(energy_J), projectile_species.mass_kg, DEUTERON.mass_kg, ion_temperature_J, "charge_exchange") for energy_J in component_energies_J], dtype=float)
        component_tritium_target_mean_energy_J = np.asarray([heavy_particle_maxwellian_target_energy_moment_J(float(energy_J), projectile_species.mass_kg, TRITON.mass_kg, ion_temperature_J, "charge_exchange") for energy_J in component_energies_J], dtype=float)
    if component_deuterium_target_pumpout_rate_s.shape != component_energies_J.shape or component_tritium_target_pumpout_rate_s.shape != component_energies_J.shape:
        raise RuntimeError("beam component charge exchange rates must match component energies")
    component_deuterium_target_pumpout_energy_W = component_deuterium_target_pumpout_rate_s * component_deuterium_target_mean_energy_J
    component_tritium_target_pumpout_energy_W = component_tritium_target_pumpout_rate_s * component_tritium_target_mean_energy_J
    deuterium_target_pumpout_energy_W = float(np.sum(component_deuterium_target_pumpout_energy_W))
    tritium_target_pumpout_energy_W = float(np.sum(component_tritium_target_pumpout_energy_W))
    representative_component = attenuated_source.component_sources[0].spatial_source
    birth_pitch_profile_deg = np.rad2deg(representative_component.birth_pitch_angle_profile_rad)
    birth_lambda_profile = np.asarray(representative_component.Lambda_birth_profile, dtype=float)
    metadata = {
        "beam_source_model": "single_species_resolved_neutral_beam",
        "beam_id": beam_config.id,
        "beam_enabled": beam_config.enabled,
        "beam_species": beam_config.species,
        "beam_attenuation_model": beam_config.model,
        "beam_atomic_data_model": beam_config.atomic_data_model,
        "beam_atomic_data_reference": HYDROGENIC_DT_GROUND_STATE_REFERENCE,
        "beam_atomic_data_collisiondb_qid": COLLISIONDB_D107356_QID,
        "beam_atomic_data_collisiondb_doi": COLLISIONDB_D107356_DOI,
        "beam_atomic_data_collisiondb_method": COLLISIONDB_D107356_METHOD,
        "beam_atomic_data_collisiondb_energy_frame": COLLISIONDB_D107356_ENERGY_FRAME,
        "beam_atomic_data_collisiondb_energy_units": COLLISIONDB_D107356_ENERGY_UNITS,
        "beam_atomic_data_collisiondb_cross_section_units": COLLISIONDB_D107356_CROSS_SECTION_UNITS,
        "beam_atomic_data_collisiondb_data_sha256": COLLISIONDB_D107356_DATA_SHA256,
        "beam_atomic_data_collisiondb_source_pdf_sha256": COLLISIONDB_D107356_SOURCE_PDF_SHA256,
        "beam_atomic_data_collisiondb_min_energy_keV_per_u": COLLISIONDB_D107356_MIN_KEV_PER_U,
        "beam_atomic_data_collisiondb_max_energy_keV_per_u": COLLISIONDB_D107356_MAX_KEV_PER_U,
        "beam_atomic_data_ionization_low_energy_truncation_max_relative_rate": IONIZATION_LOW_ENERGY_TRUNCATION_MAX_RELATIVE_RATE,
        "beam_atomic_data_ionization_high_energy_truncation_max_relative_rate": IONIZATION_HIGH_ENERGY_TRUNCATION_MAX_RELATIVE_RATE,
        "beam_target_density_profile_source": target_density_profile_source,
        "beam_target_density_role": "electron_impact_ionization_target_density",
        "beam_deuterium_target_density_profile_m3": deuterium_target_density,
        "beam_tritium_target_density_profile_m3": tritium_target_density,
        "beam_thermal_deuterium_target_density_profile_m3": deuterium_target_density,
        "beam_thermal_tritium_target_density_profile_m3": tritium_target_density,
        "beam_target_density_profile_m3": target_density,
        "beam_target_density_profile_defined_mask": target_density_profile_defined_mask,
        "beam_target_density_profile_support_mask": target_density_profile_support_mask,
        "beam_target_density_profile_valid": True,
        "beam_stopping_density_model": "electron_and_solved_kinetic_D_T_profiles" if kinetic_target_active else "electron_and_configured_startup_or_reference_D_T_profiles",
        "beam_heavy_particle_target_model": heavy_target_model,
        "beam_heavy_particle_target_uses_solved_kinetic_state": kinetic_target_active,
        "beam_heavy_particle_target_gyroangle_points": heavy_target_gyroangle_points,
        "beam_stopping_species_resolved": True,
        "beam_stopping_target_species": ["electrons", "deuterium", "tritium"],
        "beam_stopping_model_limitation": "ground_state_hydrogenic_equal_relative_speed_isotope_mapping",
        "beam_electron_temperature_keV": float(electron_temperature_J / KEV_TO_J),
        "beam_ion_temperature_keV": None if kinetic_target_active else float(config.plasma_closure.background_ion_temperature_keV),
        "beam_power_W": beam_config.power_W,
        "beam_energy_keV": beam_config.energy_keV,
        "beam_injection_angle_deg": configured_injection_angle_deg,
        "beam_pitch_source": pitch_source,
        "beam_pitch_model": "representative_pitch_from_beamline_relative_to_on_axis_z_directed_magnetic_field" if not override_selected else "explicit_configured_representative_pitch_override",
        "beam_pitch_spatial_model": "finite_uniform_pitch_angle_distribution_on_1d_on_axis_Bz_geometry",
        "beam_pitch_distribution_model": beam_config.pitch_distribution_model,
        "beam_pitch_half_width_deg": float(np.rad2deg(pitch_half_width_rad)),
        "beam_pitch_width_source": pitch_width_source,
        "beam_pitch_quadrature_points": int(pitch_samples.size),
        "beam_pitch_quadrature_samples_deg": np.rad2deg(pitch_samples),
        "beam_pitch_quadrature_weights": pitch_weights,
        "beam_pitch_distribution_is_finite_width": bool(pitch_half_width_rad > 0.0),
        "beam_pitch_endpoint_reversal_invariant": True,
        "beam_coordinate_order": "x_y_z",
        "beam_magnetic_axis_coordinate": "z",
        "beamline_direction_unit": path.direction_unit,
        "beamline_axis_angle_deg": beamline_axis_angle_deg,
        "beam_configured_injection_angle_deg": configured_injection_angle_deg,
        "beam_effective_injection_angle_deg": effective_injection_angle_deg,
        "beam_injection_angle_mismatch_deg": injection_angle_mismatch_deg,
        "beam_injection_geometry_consistent": injection_geometry_consistent,
        "beam_injection_angle_tolerance_deg": float(getattr(beam_config, "injection_angle_tolerance_deg", 1.0e-6)),
        "beam_pitch_override_selected": override_selected,
        "beam_pitch_selection_valid": bool(injection_geometry_consistent or override_selected),
        "beam_birth_pitch_angle_profile_deg": birth_pitch_profile_deg,
        "beam_birth_pitch_profile_scope": "all_axial_cells_with_birth_rate_profile_identifying_active_cells",
        "beam_birth_lambda_profile": birth_lambda_profile,
        "beam_birth_eta_profile": None,
        "beam_birth_eta_profile_status": "deferred_to_geometry_dependent_modal_eta_lambda_map",
        "beam_component_ids": [component.id for component in beam_config.components],
        "beam_power_fractions": component_power_fractions.tolist(),
        "beam_energy_fractions": component_energy_multipliers.tolist(),
        "beam_fast_birth_rate_s": float(attenuated_source.total_birth_rate_s),
        "beam_ionization_birth_rate_s": float(attenuated_source.total_ionization_rate_s),
        "beam_electron_ionization_birth_rate_s": float(attenuated_source.total_electron_ionization_rate_s),
        "beam_deuterium_ionization_birth_rate_s": float(attenuated_source.total_deuterium_ionization_rate_s),
        "beam_tritium_ionization_birth_rate_s": float(attenuated_source.total_tritium_ionization_rate_s),
        "beam_charge_exchange_fast_birth_rate_s": float(attenuated_source.total_charge_exchange_rate_s),
        "beam_deuterium_charge_exchange_fast_birth_rate_s": float(attenuated_source.total_deuterium_charge_exchange_rate_s),
        "beam_tritium_charge_exchange_fast_birth_rate_s": float(attenuated_source.total_tritium_charge_exchange_rate_s),
        "beam_charge_exchange_target_ion_pumpout_rate_s": float(attenuated_source.total_charge_exchange_rate_s),
        "beam_deuterium_target_pumpout_rate_s": float(attenuated_source.total_deuterium_charge_exchange_rate_s),
        "beam_tritium_target_pumpout_rate_s": float(attenuated_source.total_tritium_charge_exchange_rate_s),
        "beam_deuterium_target_pumpout_energy_W": deuterium_target_pumpout_energy_W,
        "beam_tritium_target_pumpout_energy_W": tritium_target_pumpout_energy_W,
        "beam_charge_exchange_target_pumpout_energy_W": deuterium_target_pumpout_energy_W + tritium_target_pumpout_energy_W,
        "beam_component_deuterium_target_pumpout_rate_s": component_deuterium_target_pumpout_rate_s,
        "beam_component_tritium_target_pumpout_rate_s": component_tritium_target_pumpout_rate_s,
        "beam_component_deuterium_target_mean_energy_J": component_deuterium_target_mean_energy_J,
        "beam_component_tritium_target_mean_energy_J": component_tritium_target_mean_energy_J,
        "beam_component_deuterium_target_pumpout_energy_W": component_deuterium_target_pumpout_energy_W,
        "beam_component_tritium_target_pumpout_energy_W": component_tritium_target_pumpout_energy_W,
        "beam_net_plasma_fueling_rate_s": float(attenuated_source.total_net_fueling_rate_s),
        "beam_particle_bookkeeping_scope": "fast_population_birth_equals_ionization_plus_charge_exchange_net_plasma_fueling_equals_ionization",
        "beam_charge_exchange_target_pumpout_kinetic_model": "solved_kinetic_target_reaction_weighted_energy_moment" if kinetic_target_active else "target_resolved_rate_and_Anderson_relative_velocity_weighted_target_kinetic_energy_primary_target_neutral_escapes",
        "beam_deposited_birth_power_W": float(attenuated_source.total_deposited_birth_power_W),
        "beam_shine_through_power_fraction": shine_through_power_fraction(attenuation),
        "beam_path_physical_intersected_cell_count": physical_intersection_count,
        "beam_path_intersected_cell_count": effective_intersection_count,
        "beam_volume_averaged_source_nonzero_cells": int(np.count_nonzero(volume_averaged_source > 0.0)),
        "beam_attenuation_beam_start_m": path.start_point_m,
        "beam_attenuation_beam_end_m": path.end_point_m,
        "beam_attenuation_beam_direction_unit": path.direction_unit,
        "beam_attenuation_beam_radius_m": float(path.beam_radius_m),
        "beam_attenuation_radial_overlap_model": path.radial_overlap_model,
        "beam_path_center_points_m": path.beam_center_points_m,
        "beam_path_beam_radius_profile_m": path.beam_radius_profile_m,
        "beam_path_physical_path_lengths_m": path.physical_path_lengths_m,
        "beam_path_effective_path_lengths_m": path.effective_path_lengths_m,
        "beam_path_radial_overlap_fractions": path.radial_overlap_fractions,
        "beam_path_plasma_radius_profile_m": path.plasma_radius_m,
        "beam_source_axial_physical_birth_rate_density_m3_s": axial_birth_density,
        "beam_source_axial_physical_birth_rate_s": axial_birth_rate,
        "beam_birth_profile_support_mask": birth_profile_support,
        "beam_birth_profile_physical_cell_volumes_m3": np.asarray(geometry.cell_volumes_m3, dtype=float),
        "beam_birth_profile_quantity": "volumetric_birth_rate_density",
        "beam_birth_profile_units": "m^-3 s^-1",
        "beam_birth_profile_weight_model": "physical_source_cell_volume",
        "beam_integrated_birth_profile_weight_model": "integrated_per_axial_interval_unit_weights",
        "beam_component_deposited_birth_power_W": component_deposited_power,
        "beam_component_configured_power_fraction": component_power_fractions,
        "beam_component_support_mask": component_support,
        "beam_component_metric_quantity": "deposited_birth_power",
        "beam_component_metric_units": "W",
        "beam_source_axial_physical_birth_power_density_W_m3": axial_power_density,
        "beam_source_axial_physical_birth_power_W": axial_birth_power,
        "speed_centers_m_s": speed_grid.centers_m_s,
        "speed_faces_m_s": speed_grid.faces_m_s,
        "lambda_centers": lambda_grid.centers,
        "lambda_faces": lambda_grid.faces,
        "pitch_centers": pitch_grid.centers,
        "pitch_faces": pitch_grid.faces,
    }

    return BeamDepositionResult(beam_id=beam_config.id, projectile_species=beam_config.species, path_geometry=path, attenuation=attenuation, attenuated_source=attenuated_source, component_energies_J=component_energies_J, source_z_v_lambda_per_s=attenuated_source.total_source_v_lambda_per_s, volume_averaged_source_v_lambda_per_s=volume_averaged_source, total_fast_birth_rate_s=total_fast_ion_birth_rate_multi_s(attenuation), total_deposited_birth_power_W=float(attenuated_source.total_deposited_birth_power_W), metadata=metadata)

def _beam_deposition_output_metadata(deposition: BeamDepositionResult) -> dict[str, object]:
    """Return one beam diagnostic record without repeating shared ensemble grids"""
    metadata = deposition.metadata

    return {
        "beam_id": deposition.beam_id,
        "enabled": metadata["beam_enabled"],
        "projectile_species": deposition.projectile_species,
        "power_W": metadata["beam_power_W"],
        "energy_keV": metadata["beam_energy_keV"],
        "configured_injection_angle_deg": metadata["beam_configured_injection_angle_deg"],
        "effective_injection_angle_deg": metadata["beam_effective_injection_angle_deg"],
        "beamline_axis_angle_deg": metadata["beamline_axis_angle_deg"],
        "injection_geometry_consistent": metadata["beam_injection_geometry_consistent"],
        "component_ids": metadata["beam_component_ids"],
        "component_energies_J": deposition.component_energies_J,
        "component_configured_power_fraction": metadata["beam_component_configured_power_fraction"],
        "component_deposited_birth_power_W": metadata["beam_component_deposited_birth_power_W"],
        "fast_birth_rate_s": deposition.total_fast_birth_rate_s,
        "electron_ionization_birth_rate_s": metadata["beam_electron_ionization_birth_rate_s"],
        "deuterium_ionization_birth_rate_s": metadata["beam_deuterium_ionization_birth_rate_s"],
        "tritium_ionization_birth_rate_s": metadata["beam_tritium_ionization_birth_rate_s"],
        "ionization_birth_rate_s": metadata["beam_ionization_birth_rate_s"],
        "deuterium_charge_exchange_fast_birth_rate_s": metadata["beam_deuterium_charge_exchange_fast_birth_rate_s"],
        "tritium_charge_exchange_fast_birth_rate_s": metadata["beam_tritium_charge_exchange_fast_birth_rate_s"],
        "charge_exchange_fast_birth_rate_s": metadata["beam_charge_exchange_fast_birth_rate_s"],
        "deuterium_target_pumpout_rate_s": metadata["beam_deuterium_target_pumpout_rate_s"],
        "tritium_target_pumpout_rate_s": metadata["beam_tritium_target_pumpout_rate_s"],
        "deuterium_target_pumpout_energy_W": metadata["beam_deuterium_target_pumpout_energy_W"],
        "tritium_target_pumpout_energy_W": metadata["beam_tritium_target_pumpout_energy_W"],
        "component_deuterium_target_mean_energy_J": metadata["beam_component_deuterium_target_mean_energy_J"],
        "component_tritium_target_mean_energy_J": metadata["beam_component_tritium_target_mean_energy_J"],
        "net_plasma_fueling_rate_s": metadata["beam_net_plasma_fueling_rate_s"],
        "deposited_birth_power_W": deposition.total_deposited_birth_power_W,
        "shine_through_power_W": deposition.attenuated_source.shine_through_power_W,
        "start_m": metadata["beam_attenuation_beam_start_m"],
        "end_m": metadata["beam_attenuation_beam_end_m"],
        "direction_unit": metadata["beam_attenuation_beam_direction_unit"],
        "beam_radius_m": metadata["beam_attenuation_beam_radius_m"],
        "path_center_points_m": metadata["beam_path_center_points_m"],
        "beam_radius_profile_m": metadata["beam_path_beam_radius_profile_m"],
        "physical_path_lengths_m": metadata["beam_path_physical_path_lengths_m"],
        "effective_path_lengths_m": metadata["beam_path_effective_path_lengths_m"],
        "radial_overlap_fractions": metadata["beam_path_radial_overlap_fractions"],
        "plasma_radius_profile_m": metadata["beam_path_plasma_radius_profile_m"],
        "birth_pitch_angle_profile_deg": metadata["beam_birth_pitch_angle_profile_deg"],
        "birth_lambda_profile": metadata["beam_birth_lambda_profile"],
        "axial_birth_rate_density_m3_s": metadata["beam_source_axial_physical_birth_rate_density_m3_s"],
        "axial_birth_rate_s": metadata["beam_source_axial_physical_birth_rate_s"],
        "axial_birth_profile_support_mask": metadata["beam_birth_profile_support_mask"],
        "axial_birth_power_density_W_m3": metadata["beam_source_axial_physical_birth_power_density_W_m3"],
        "axial_birth_power_W": metadata["beam_source_axial_physical_birth_power_W"],
        "source_z_v_lambda_per_s": deposition.source_z_v_lambda_per_s,
        "volume_averaged_source_v_lambda_per_s": deposition.volume_averaged_source_v_lambda_per_s,
    }

def _multi_beam_metadata(config: SourceModelRunConfig, geometry: GeometryStageResult, enabled_beams: tuple[BeamConfig, ...], depositions: tuple[BeamDepositionResult, ...], canonical_depositions: tuple[BeamDepositionResult, ...], attenuated_source_by_species: dict[str, AttenuatedMultiEnergyBeamSource], combined_volume_source_by_species: dict[str, np.ndarray], speed_grid_by_species: dict[str, SpeedGrid], primary_species: str, lambda_grid: LambdaGrid, pitch_grid: PitchGrid) -> dict[str, object]:
    """Build combined and per beam attenuation, deposition, source, and atomic data metadata for the beam ensemble"""
    first = canonical_depositions[0].metadata
    for deposition in canonical_depositions[1:]:
        for key in (
            "beam_target_density_profile_m3",
            "beam_target_density_profile_defined_mask",
            "beam_target_density_profile_support_mask",
            "beam_thermal_deuterium_target_density_profile_m3",
            "beam_thermal_tritium_target_density_profile_m3",
        ):
            if not np.array_equal(np.asarray(deposition.metadata[key]), np.asarray(first[key])):
                raise RuntimeError("all beams in one iteration must use the same target density state")
    beam_by_id = {beam.id: beam for beam in enabled_beams}
    component_energies_J = np.concatenate([deposition.component_energies_J for deposition in canonical_depositions])
    component_deposited_power_W = np.concatenate([np.asarray(deposition.metadata["beam_component_deposited_birth_power_W"], dtype=float) for deposition in canonical_depositions])
    total_input_power_W = float(sum(beam.power_W for beam in enabled_beams))
    component_input_power_W = np.concatenate([beam_by_id[deposition.beam_id].power_W * np.asarray([component.power_fraction for component in beam_by_id[deposition.beam_id].components], dtype=float) for deposition in canonical_depositions])
    component_configured_power_fraction = component_input_power_W / total_input_power_W
    component_ids = [f"{deposition.beam_id}:{component.id}" for deposition in canonical_depositions for component in beam_by_id[deposition.beam_id].components]
    component_beam_ids = [deposition.beam_id for deposition in canonical_depositions for component in beam_by_id[deposition.beam_id].components]
    axial_birth_rate_density = np.sum([np.asarray(deposition.metadata["beam_source_axial_physical_birth_rate_density_m3_s"], dtype=float) for deposition in canonical_depositions], axis=0)
    axial_birth_rate = np.sum([np.asarray(deposition.metadata["beam_source_axial_physical_birth_rate_s"], dtype=float) for deposition in canonical_depositions], axis=0)
    axial_birth_power_density = np.sum([np.asarray(deposition.metadata["beam_source_axial_physical_birth_power_density_W_m3"], dtype=float) for deposition in canonical_depositions], axis=0)
    axial_birth_power = np.sum([np.asarray(deposition.metadata["beam_source_axial_physical_birth_power_W"], dtype=float) for deposition in canonical_depositions], axis=0)
    birth_support = np.any([np.asarray(deposition.metadata["beam_birth_profile_support_mask"], dtype=bool) for deposition in canonical_depositions], axis=0)
    physical_path_support = np.any([np.asarray(deposition.metadata["beam_path_physical_path_lengths_m"], dtype=float) > 0.0 for deposition in canonical_depositions], axis=0)
    effective_path_support = np.any([np.asarray(deposition.metadata["beam_path_effective_path_lengths_m"], dtype=float) > 0.0 for deposition in canonical_depositions], axis=0)
    pitch_weighted_sum = np.sum([np.asarray(deposition.metadata["beam_birth_pitch_angle_profile_deg"], dtype=float) * np.asarray(deposition.metadata["beam_source_axial_physical_birth_rate_s"], dtype=float) for deposition in canonical_depositions], axis=0)
    lambda_weighted_sum = np.sum([np.asarray(deposition.metadata["beam_birth_lambda_profile"], dtype=float) * np.asarray(deposition.metadata["beam_source_axial_physical_birth_rate_s"], dtype=float) for deposition in canonical_depositions], axis=0)
    birth_pitch_profile_deg = np.divide(pitch_weighted_sum, axial_birth_rate, out=np.zeros_like(pitch_weighted_sum), where=axial_birth_rate > 0.0)
    birth_lambda_profile = np.divide(lambda_weighted_sum, axial_birth_rate, out=np.zeros_like(lambda_weighted_sum), where=axial_birth_rate > 0.0)
    scalar_weights = np.asarray([deposition.total_fast_birth_rate_s for deposition in canonical_depositions], dtype=float)
    if float(np.sum(scalar_weights)) <= 0.0:
        scalar_weights = np.asarray([beam_by_id[deposition.beam_id].power_W for deposition in canonical_depositions], dtype=float)
    effective_injection_angle_deg = float(np.average([deposition.metadata["beam_effective_injection_angle_deg"] for deposition in canonical_depositions], weights=scalar_weights))
    configured_injection_angle_deg = float(np.average([deposition.metadata["beam_configured_injection_angle_deg"] for deposition in canonical_depositions], weights=scalar_weights))
    total_birth_rate_s = float(sum(source.total_birth_rate_s for source in attenuated_source_by_species.values()))
    total_deposited_power_W = float(sum(source.total_deposited_birth_power_W for source in attenuated_source_by_species.values()))
    total_shine_through_power_W = float(sum(source.shine_through_power_W for source in attenuated_source_by_species.values()))
    mean_birth_energy_keV = 0.0 if total_birth_rate_s <= 0.0 else total_deposited_power_W / total_birth_rate_s / KEV_TO_J
    active_species = tuple(sorted(attenuated_source_by_species))
    axial_birth_rate_density_by_species = {species: np.asarray(attenuated_source_by_species[species].total_axial_birth_rate_density_m3_s, dtype=float) for species in active_species}
    axial_birth_rate_by_species = {species: axial_birth_rate_density_by_species[species] * np.asarray(geometry.cell_volumes_m3, dtype=float) for species in active_species}
    axial_birth_power_density_by_species = {species: np.asarray(attenuated_source_by_species[species].total_axial_birth_power_density_W_m3, dtype=float) for species in active_species}
    axial_birth_power_by_species = {species: axial_birth_power_density_by_species[species] * np.asarray(geometry.cell_volumes_m3, dtype=float) for species in active_species}
    species_birth_rate_density_sum = np.sum(tuple(axial_birth_rate_density_by_species.values()), axis=0)
    species_birth_power_density_sum = np.sum(tuple(axial_birth_power_density_by_species.values()), axis=0)
    birth_rate_identity_error = float(np.max(np.abs(species_birth_rate_density_sum - axial_birth_rate_density)) / max(float(np.max(np.abs(axial_birth_rate_density))), np.finfo(float).tiny))
    birth_power_identity_error = float(np.max(np.abs(species_birth_power_density_sum - axial_birth_power_density)) / max(float(np.max(np.abs(axial_birth_power_density))), np.finfo(float).tiny))
    nominal_energies_keV = np.asarray([beam.energy_keV for beam in enabled_beams], dtype=float)
    configured_angles_deg = np.asarray([beam.injection_angle_deg for beam in enabled_beams], dtype=float)
    common_nominal_energy_keV = float(nominal_energies_keV[0]) if np.allclose(nominal_energies_keV, nominal_energies_keV[0], rtol=0.0, atol=0.0) else None
    common_configured_angle_deg = float(configured_angles_deg[0]) if np.allclose(configured_angles_deg, configured_angles_deg[0], rtol=0.0, atol=0.0) else None
    attenuation_models = {beam.model for beam in enabled_beams}
    atomic_models = {beam.atomic_data_model for beam in enabled_beams}
    pitch_models = {beam.pitch_distribution_model for beam in enabled_beams}
    primary_speed_grid = speed_grid_by_species[primary_species]

    metadata = {
        "beam_source_model": "species_resolved_neutral_beam_ensemble",
        "beam_attenuation_model": next(iter(attenuation_models)) if len(attenuation_models) == 1 else "mixed_beam_attenuation_models",
        "beam_atomic_data_model": next(iter(atomic_models)) if len(atomic_models) == 1 else "mixed_atomic_data_models",
        "beam_atomic_data_reference": HYDROGENIC_DT_GROUND_STATE_REFERENCE,
        "beam_atomic_data_collisiondb_qid": COLLISIONDB_D107356_QID,
        "beam_atomic_data_collisiondb_doi": COLLISIONDB_D107356_DOI,
        "beam_atomic_data_collisiondb_method": COLLISIONDB_D107356_METHOD,
        "beam_atomic_data_collisiondb_energy_frame": COLLISIONDB_D107356_ENERGY_FRAME,
        "beam_atomic_data_collisiondb_energy_units": COLLISIONDB_D107356_ENERGY_UNITS,
        "beam_atomic_data_collisiondb_cross_section_units": COLLISIONDB_D107356_CROSS_SECTION_UNITS,
        "beam_atomic_data_collisiondb_data_sha256": COLLISIONDB_D107356_DATA_SHA256,
        "beam_atomic_data_collisiondb_source_pdf_sha256": COLLISIONDB_D107356_SOURCE_PDF_SHA256,
        "beam_atomic_data_collisiondb_min_energy_keV_per_u": COLLISIONDB_D107356_MIN_KEV_PER_U,
        "beam_atomic_data_collisiondb_max_energy_keV_per_u": COLLISIONDB_D107356_MAX_KEV_PER_U,
        "beam_atomic_data_ionization_low_energy_truncation_max_relative_rate": IONIZATION_LOW_ENERGY_TRUNCATION_MAX_RELATIVE_RATE,
        "beam_atomic_data_ionization_high_energy_truncation_max_relative_rate": IONIZATION_HIGH_ENERGY_TRUNCATION_MAX_RELATIVE_RATE,
        "beam_target_density_profile_source": first["beam_target_density_profile_source"],
        "beam_target_density_role": first["beam_target_density_role"],
        "beam_target_density_profile_m3": first["beam_target_density_profile_m3"],
        "beam_target_density_profile_defined_mask": first["beam_target_density_profile_defined_mask"],
        "beam_target_density_profile_support_mask": first["beam_target_density_profile_support_mask"],
        "beam_thermal_deuterium_target_density_profile_m3": first["beam_thermal_deuterium_target_density_profile_m3"],
        "beam_thermal_tritium_target_density_profile_m3": first["beam_thermal_tritium_target_density_profile_m3"],
        "beam_target_density_profile_valid": all(deposition.metadata.get("beam_target_density_profile_valid") is True for deposition in canonical_depositions),
        "beam_stopping_density_model": first["beam_stopping_density_model"],
        "beam_stopping_species_resolved": True,
        "beam_stopping_target_species": ["electrons", "deuterium", "tritium"],
        "beam_stopping_model_limitation": first["beam_stopping_model_limitation"],
        "beam_electron_temperature_keV": first["beam_electron_temperature_keV"],
        "beam_ion_temperature_keV": first["beam_ion_temperature_keV"],
        "beam_power_W": total_input_power_W,
        "beam_effective_birth_energy_keV": float(mean_birth_energy_keV),
        "beam_representative_injection_angle_deg": configured_injection_angle_deg,
        "beam_effective_injection_angle_deg": effective_injection_angle_deg,
        "beam_pitch_distribution_model": next(iter(pitch_models)) if len(pitch_models) == 1 else "mixed_beam_pitch_models",
        "beam_pitch_selection_valid": all(deposition.metadata.get("beam_pitch_selection_valid") is True for deposition in canonical_depositions),
        "beam_injection_geometry_consistent": all(deposition.metadata.get("beam_injection_geometry_consistent") is True for deposition in canonical_depositions),
        "beam_pitch_override_selected": any(deposition.metadata.get("beam_pitch_override_selected") is True for deposition in canonical_depositions),
        "beam_birth_pitch_angle_profile_deg": birth_pitch_profile_deg,
        "beam_birth_pitch_profile_scope": "birth_rate_weighted_combined_D_T_beam_profile",
        "beam_birth_lambda_profile": birth_lambda_profile,
        "beam_birth_eta_profile": None,
        "beam_birth_eta_profile_status": "deferred_to_geometry_dependent_modal_eta_lambda_map",
        "beam_component_ids": component_ids,
        "beam_component_beam_ids": component_beam_ids,
        "beam_power_fractions": component_configured_power_fraction,
        "beam_energy_fractions": np.concatenate([np.asarray([component.energy_multiplier for component in beam_by_id[deposition.beam_id].components], dtype=float) for deposition in canonical_depositions]),
        "beam_fast_birth_rate_s": total_birth_rate_s,
        "beam_electron_ionization_birth_rate_s": float(sum(source.total_electron_ionization_rate_s for source in attenuated_source_by_species.values())),
        "beam_deuterium_ionization_birth_rate_s": float(sum(source.total_deuterium_ionization_rate_s for source in attenuated_source_by_species.values())),
        "beam_tritium_ionization_birth_rate_s": float(sum(source.total_tritium_ionization_rate_s for source in attenuated_source_by_species.values())),
        "beam_ionization_birth_rate_s": float(sum(source.total_ionization_rate_s for source in attenuated_source_by_species.values())),
        "beam_deuterium_charge_exchange_fast_birth_rate_s": float(sum(source.total_deuterium_charge_exchange_rate_s for source in attenuated_source_by_species.values())),
        "beam_tritium_charge_exchange_fast_birth_rate_s": float(sum(source.total_tritium_charge_exchange_rate_s for source in attenuated_source_by_species.values())),
        "beam_charge_exchange_fast_birth_rate_s": float(sum(source.total_charge_exchange_rate_s for source in attenuated_source_by_species.values())),
        "beam_charge_exchange_target_ion_pumpout_rate_s": float(sum(source.total_charge_exchange_rate_s for source in attenuated_source_by_species.values())),
        "beam_deuterium_target_pumpout_rate_s": float(sum(source.total_deuterium_charge_exchange_rate_s for source in attenuated_source_by_species.values())),
        "beam_tritium_target_pumpout_rate_s": float(sum(source.total_tritium_charge_exchange_rate_s for source in attenuated_source_by_species.values())),
        "beam_deuterium_target_pumpout_energy_W": float(sum(float(deposition.metadata["beam_deuterium_target_pumpout_energy_W"]) for deposition in canonical_depositions)),
        "beam_tritium_target_pumpout_energy_W": float(sum(float(deposition.metadata["beam_tritium_target_pumpout_energy_W"]) for deposition in canonical_depositions)),
        "beam_charge_exchange_target_pumpout_energy_W": float(sum(float(deposition.metadata["beam_charge_exchange_target_pumpout_energy_W"]) for deposition in canonical_depositions)),
        "beam_net_plasma_fueling_rate_s": float(sum(source.total_net_fueling_rate_s for source in attenuated_source_by_species.values())),
        "beam_particle_bookkeeping_scope": "sum_of_independently_attenuated_D_T_beam_birth_fueling_and_target_pumpout_rates",
        "beam_charge_exchange_target_pumpout_kinetic_model": first["beam_charge_exchange_target_pumpout_kinetic_model"],
        "beam_deposited_birth_power_W": total_deposited_power_W,
        "beam_shine_through_power_fraction": total_shine_through_power_W / total_input_power_W,
        "beam_path_physical_intersected_cell_count": int(np.count_nonzero(physical_path_support)),
        "beam_path_intersected_cell_count": int(np.count_nonzero(effective_path_support)),
        "beam_volume_averaged_source_nonzero_cells": int(sum(np.count_nonzero(source > 0.0) for source in combined_volume_source_by_species.values())),
        "beam_source_axial_physical_birth_rate_density_m3_s": axial_birth_rate_density,
        "beam_source_axial_physical_birth_rate_s": axial_birth_rate,
        "beam_axial_birth_rate_density_m3_s_by_species": axial_birth_rate_density_by_species,
        "beam_axial_birth_rate_s_by_species": axial_birth_rate_by_species,
        "beam_birth_profile_support_mask": birth_support,
        "beam_birth_profile_physical_cell_volumes_m3": np.asarray(geometry.cell_volumes_m3, dtype=float),
        "beam_birth_profile_quantity": "volumetric_birth_rate_density",
        "beam_birth_profile_units": "m^-3 s^-1",
        "beam_birth_profile_weight_model": "physical_source_cell_volume",
        "beam_integrated_birth_profile_weight_model": "integrated_per_axial_interval_unit_weights",
        "beam_component_deposited_birth_power_W": component_deposited_power_W,
        "beam_component_configured_power_fraction": component_configured_power_fraction,
        "beam_component_support_mask": component_input_power_W > 0.0,
        "beam_component_metric_quantity": "deposited_birth_power",
        "beam_component_metric_units": "W",
        "beam_source_axial_physical_birth_power_density_W_m3": axial_birth_power_density,
        "beam_source_axial_physical_birth_power_W": axial_birth_power,
        "beam_axial_birth_power_density_W_m3_by_species": axial_birth_power_density_by_species,
        "beam_axial_birth_power_W_by_species": axial_birth_power_by_species,
        "beam_species_birth_rate_profile_identity_relative_error": birth_rate_identity_error,
        "beam_species_birth_power_profile_identity_relative_error": birth_power_identity_error,
        "beam_species_profile_identity_relative_tolerance": 1.0e-12,
        "beam_species_birth_rate_profile_identity_check_passed": bool(birth_rate_identity_error <= 1.0e-12),
        "beam_species_birth_power_profile_identity_check_passed": bool(birth_power_identity_error <= 1.0e-12),
        "speed_centers_m_s": primary_speed_grid.centers_m_s,
        "speed_faces_m_s": primary_speed_grid.faces_m_s,
        "lambda_centers": lambda_grid.centers,
        "lambda_faces": lambda_grid.faces,
        "pitch_centers": pitch_grid.centers,
        "pitch_faces": pitch_grid.faces,
        "beam_collection_model": "independent_beam_attenuation_with_species_source_aggregation",
        "beam_count": len(config.beams),
        "enabled_beam_count": len(enabled_beams),
        "beam_order": [beam.id for beam in enabled_beams],
        "beam_canonical_order": [deposition.beam_id for deposition in canonical_depositions],
        "beam_ids": [beam.id for beam in config.beams],
        "enabled_beam_ids": [beam.id for beam in enabled_beams],
        "beam_species_by_id": {beam.id: beam.species for beam in config.beams},
        "beam_power_W_by_id": {beam.id: beam.power_W for beam in config.beams},
        "beam_energy_keV_by_id": {beam.id: beam.energy_keV for beam in config.beams},
        "beam_component_ids_by_id": {beam.id: [component.id for component in beam.components] for beam in config.beams},
        "beam_fast_birth_rate_s_by_id": {deposition.beam_id: deposition.total_fast_birth_rate_s for deposition in depositions},
        "beam_deposited_birth_power_W_by_id": {deposition.beam_id: deposition.total_deposited_birth_power_W for deposition in depositions},
        "beam_shine_through_power_W_by_id": {deposition.beam_id: deposition.attenuated_source.shine_through_power_W for deposition in depositions},
        "beam_deposition_by_id": {deposition.beam_id: _beam_deposition_output_metadata(deposition) for deposition in depositions},
        "beam_source_aggregation": "sum_by_projectile_species_before_nonlinear_fbis",
        "beam_source_aggregation_order": "canonical_beam_id_order",
        "beam_active_projectile_species": list(active_species),
        "beam_primary_projectile_species": primary_species,
        "beam_total_fast_birth_rate_s_by_species": {species: float(source.total_birth_rate_s) for species, source in attenuated_source_by_species.items()},
        "beam_total_deposited_birth_power_W_by_species": {species: float(source.total_deposited_birth_power_W) for species, source in attenuated_source_by_species.items()},
        "beam_source_construction_status": "operating_point_and_fusion_compatible_fast_tritium" if TRITON.species_id in active_species else "operating_point_compatible_deuterium_only",
    }
    if common_nominal_energy_keV is not None:
        metadata["beam_energy_keV"] = common_nominal_energy_keV
    if common_configured_angle_deg is not None:
        metadata["beam_injection_angle_deg"] = common_configured_angle_deg

    return metadata

def build_beam_stage(config: SourceModelRunConfig, geometry: GeometryStageResult, *, electron_temperature_J: float | None = None, target_density_profile_m3: np.ndarray | None = None, target_density_defined_mask: np.ndarray | None = None, operating_point_density_state: OperatingPointDensityState | None = None, kinetic_target_state: KineticStageResult | None = None) -> BeamEnsembleResult:
    """
    Build all enabled beams independently and combine their sources by fast ion species
    
    Each beam keeps its own attenuation and deposited power accounting while species sources share one compatible speed grid before entering the modal kinetic stage
    """
    enabled_beams = config.enabled_beams
    temperature_J = float(config.plasma_closure.electron_temperature_initial_guess_keV * KEV_TO_J if electron_temperature_J is None else electron_temperature_J)
    if not np.isfinite(temperature_J) or temperature_J <= 0.0:
        raise ValueError("electron_temperature_J must be positive and finite")
    lambda_grid = uniform_lambda_grid(config.kinetic_electrostatic.lambda_bins, 0.0, 1.0)
    pitch_grid = uniform_pitch_grid(config.kinetic_electrostatic.pitch_bins, -1.0, 1.0)
    startup_seed_only = uses_startup_seed_only(config.plasma_closure)
    beams_by_species: dict[str, tuple[BeamConfig, ...]] = {}
    for species_id in sorted({beam.species for beam in enabled_beams}):
        beams_by_species[species_id] = tuple(beam for beam in enabled_beams if beam.species == species_id)
    speed_grid_by_species: dict[str, SpeedGrid] = {}
    previous_speed_grids = None if kinetic_target_state is None else kinetic_target_state.speed_grid_by_species
    for species_id, species_beams in beams_by_species.items():
        energies_J = np.concatenate([beam.energy_keV * np.asarray([component.energy_multiplier for component in beam.components], dtype=float) * KEV_TO_J for beam in species_beams])
        speed_grid = _speed_grid_for_species(config, ion_species(species_id), energies_J)
        previous_speed_grid = None if previous_speed_grids is None else previous_speed_grids.get(species_id)
        if startup_seed_only and previous_speed_grid is not None:
            if not _speed_grids_match(speed_grid, previous_speed_grid):
                raise ValueError(f"kinetic target {species_id} speed grid does not match the current beam speed grid")
            speed_grid = previous_speed_grid
        speed_grid_by_species[species_id] = speed_grid
    depositions = tuple(
        _build_beam_deposition(
            config,
            geometry,
            beam,
            speed_grid_by_species[beam.species],
            lambda_grid,
            pitch_grid,
            electron_temperature_J=temperature_J,
            target_density_profile_m3=target_density_profile_m3,
            target_density_defined_mask=target_density_defined_mask,
            operating_point_density_state=operating_point_density_state,
            kinetic_target_state=kinetic_target_state,
        ) for beam in enabled_beams
    )
    canonical_depositions = tuple(sorted(depositions, key=lambda deposition: deposition.beam_id))
    attenuated_source_by_species: dict[str, AttenuatedMultiEnergyBeamSource] = {}
    combined_source_by_species: dict[str, np.ndarray] = {}
    combined_volume_source_by_species: dict[str, np.ndarray] = {}
    total_fast_birth_rate_by_species: dict[str, float] = {}
    for species_id in sorted(beams_by_species):
        species_depositions = tuple(deposition for deposition in canonical_depositions if deposition.projectile_species == species_id)
        species_source = species_depositions[0].attenuated_source if len(species_depositions) == 1 else combine_attenuated_multi_energy_sources(tuple(deposition.attenuated_source for deposition in species_depositions))
        if species_source.total_source_v_lambda_per_s is None:
            raise RuntimeError("combined attenuated beam source did not build Lambda source arrays")
        attenuated_source_by_species[species_id] = species_source
        combined_source_by_species[species_id] = np.asarray(species_source.total_source_v_lambda_per_s, dtype=float)
        combined_volume_source_by_species[species_id] = np.sum([deposition.volume_averaged_source_v_lambda_per_s for deposition in species_depositions], axis=0)
        total_fast_birth_rate_by_species[species_id] = float(species_source.total_birth_rate_s)
    primary_species = DEUTERON.species_id if DEUTERON.species_id in attenuated_source_by_species else sorted(attenuated_source_by_species)[0]
    primary_source = attenuated_source_by_species[primary_species]
    primary_speed_grid = speed_grid_by_species[primary_species]
    metadata = _multi_beam_metadata(config, geometry, enabled_beams, depositions, canonical_depositions, attenuated_source_by_species, combined_volume_source_by_species, speed_grid_by_species, primary_species, lambda_grid, pitch_grid)
    charge_exchange_sinks = None
    if startup_seed_only and kinetic_target_state is not None:
        charge_exchange_sinks = build_charge_exchange_target_sink_states(depositions=canonical_depositions, kinetic_target_state=kinetic_target_state, geometry=geometry, gyroangle_points=int(config.fusion.num_gyroangle_points))
    metadata.update({
        "beam_heavy_particle_target_role": "solved_kinetic_D_T" if startup_seed_only and kinetic_target_state is not None else "startup_seed_first_iteration" if startup_seed_only else "prescribed_reference",
        "beam_heavy_particle_target_is_final_plasma_authority": bool(startup_seed_only and kinetic_target_state is not None) or not startup_seed_only,
        "beam_heavy_particle_target_transition": "complete" if startup_seed_only and kinetic_target_state is not None else "awaiting_first_kinetic_state" if startup_seed_only else "not_applicable",
        "beam_charge_exchange_target_sink_model": "inactive_first_seed_iteration" if startup_seed_only and kinetic_target_state is None else "reaction_weighted_speed_resolved_pitch_averaged" if charge_exchange_sinks is not None else "not_modeled_reference_target",
        "beam_charge_exchange_target_sink_rate_s_by_species": {} if charge_exchange_sinks is None else {key: value.reference_event_rate_s for key, value in charge_exchange_sinks.items()},
        "beam_charge_exchange_target_sink_energy_W_by_species": {} if charge_exchange_sinks is None else {key: value.reference_energy_removal_W for key, value in charge_exchange_sinks.items()},
        "beam_charge_exchange_target_sink_reference_identity_error_by_species": {} if charge_exchange_sinks is None else {key: value.reference_rate_identity_relative_error for key, value in charge_exchange_sinks.items()},
        "beam_charge_exchange_target_sink_qualified_by_species": {} if charge_exchange_sinks is None else {key: value.qualified for key, value in charge_exchange_sinks.items()},
    })

    return BeamEnsembleResult(
        beam_results_by_id={deposition.beam_id: deposition for deposition in depositions},
        beam_order=tuple(deposition.beam_id for deposition in depositions),
        attenuated_source=primary_source,
        attenuated_source_by_species=attenuated_source_by_species,
        speed_grid=primary_speed_grid,
        speed_grid_by_species=speed_grid_by_species,
        lambda_grid=lambda_grid,
        pitch_grid=pitch_grid,
        component_energies_J=np.concatenate([deposition.component_energies_J for deposition in canonical_depositions if deposition.projectile_species == primary_species]),
        source_z_v_lambda_per_s=combined_source_by_species[primary_species],
        volume_averaged_source_v_lambda_per_s=combined_volume_source_by_species[primary_species],
        combined_source_by_species=combined_source_by_species,
        combined_volume_averaged_source_by_species=combined_volume_source_by_species,
        total_fast_birth_rate_by_species=total_fast_birth_rate_by_species,
        charge_exchange_sink_by_target_species=charge_exchange_sinks,
        total_fast_birth_rate_s=float(sum(total_fast_birth_rate_by_species.values())),
        total_deposited_birth_power_W=float(sum(source.total_deposited_birth_power_W for source in attenuated_source_by_species.values())),
        total_injected_power_W=float(sum(beam.power_W for beam in enabled_beams)),
        metadata=metadata,
    )
