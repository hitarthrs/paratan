"""Dataclasses for Eq 59 nonlinear convergence warm starts and speed convergence"""
from __future__ import annotations
from dataclasses import dataclass
import numpy as np

@dataclass(frozen=True)
class Eq59ConvergenceDiagnostics:
    """Record nonlinear convergence evidence for the retained Eq 59 modes
    
    Histories track modal updates residuals density mean energy fast self scales and collision operator changes
    """
    assessed: bool
    converged: bool
    iterations: int
    status: str
    failure_reason: str
    relative_change_history: tuple[float, ...]
    absolute_change_history: tuple[float, ...]
    residual_history: tuple[float, ...]
    absolute_residual_history: tuple[float, ...]
    inventory_history: tuple[float, ...]
    inventory_relative_change_history: tuple[float, ...]
    effective_energy_history_J: tuple[float, ...]
    effective_energy_relative_change_history: tuple[float, ...]
    fast_self_density_history_m3: tuple[float, ...]
    fast_self_density_relative_change_history: tuple[float, ...]
    fast_self_g_scale_history_m3_s3: tuple[float, ...]
    fast_self_net_drag_scale_history_m3_s3: tuple[float, ...]
    ion_drag_operator_relative_change_history: tuple[float, ...]
    ion_energy_diffusion_operator_relative_change_history: tuple[float, ...]
    ion_pitch_scattering_operator_relative_change_history: tuple[float, ...]
    stagnated: bool
    oscillatory: bool
    diverged: bool
    fixed_point_iteration_converged: bool
    nonlinear_residual_converged: bool
    inventory_converged: bool
    effective_energy_converged: bool
    fast_self_density_converged: bool
    collision_operator_converged: bool
    speed_resolution_assessed: bool
    speed_resolution_converged: bool
    speed_domain_assessed: bool
    speed_domain_converged: bool
    warm_start_used: bool
    warm_start_remapped: bool
    warm_start_status: str
    warm_start_rejection_reason: str
    solve_method: str
    collision_operator_model: str
    warm_start_fallback_to_pairwise_seed: bool = False
    discarded_warm_start_iterations: int = 0
    discarded_warm_start_status: str = ""
    discarded_warm_start_remapped: bool = False

@dataclass(frozen=True)
class Eq59WarmStartState:
    """Store the Eq 59 state needed to validate and remap a later warm start
    
    modal_distribution_f_j_v and active_eigenvalues_by_speed have shape (n_mode, n_speed)
    """
    speed_faces_m_s: np.ndarray
    modal_distribution_f_j_v: np.ndarray
    magnetic_eigenvalues: np.ndarray
    active_eigenvalues_by_speed: np.ndarray
    source_coefficients_by_component_j: np.ndarray
    source_speeds_m_s: np.ndarray
    spitzer_slowing_down_time_s: float
    collision_operator_model: str
    collision_operator_structure_fingerprint: tuple[object, ...]
    collision_operator_closure_fingerprint: tuple[object, ...]
    particle_mass_kg: float | None
    charge_exchange_loss_frequency_s: np.ndarray | None = None
    species_label: str = "D"
    model: str = "egedal_eq59_hot_ion_rosenbluth"
    nonlinear_converged: bool = True
    mode_eta_integrals: np.ndarray | None = None
    physical_eigenfunctions_eta: np.ndarray | None = None
    eta_grid: np.ndarray | None = None

@dataclass(frozen=True)
class _Eq59WarmStartDecision:
    """Store the validated warm start candidate reuse status remap status and rejection reason"""
    distribution: np.ndarray | None
    used: bool
    remapped: bool
    status: str
    rejection_reason: str

@dataclass(frozen=True)
class Eq59SpeedConvergenceDiagnostics:
    """Record nested speed resolution and upper domain convergence evidence
    
    Histories include physical distribution changes moments losses electron heating upper tails and source to loss ratios
    """
    assessed: bool
    converged: bool
    speed_cell_count_history: np.ndarray
    speed_domain_upper_m_s_history: np.ndarray
    runtime_s_history: np.ndarray
    nonlinear_converged_history: np.ndarray
    common_domain_distribution_relative_change_history: np.ndarray
    inventory_history: np.ndarray
    inventory_relative_change_history: np.ndarray
    effective_energy_history_J: np.ndarray
    effective_energy_relative_change_history: np.ndarray
    modal_particle_loss_rate_history: np.ndarray
    modal_particle_loss_relative_change_history: np.ndarray
    modal_power_loss_history_W_per_m3: np.ndarray
    modal_power_loss_relative_change_history: np.ndarray
    electron_heating_power_history_W_per_m3: np.ndarray
    electron_heating_power_relative_change_history: np.ndarray
    high_speed_tail_population_fraction_history: np.ndarray
    high_speed_tail_energy_fraction_history: np.ndarray
    source_to_modal_loss_ratio_history: np.ndarray
    source_to_modal_loss_ratio_relative_change_history: np.ndarray
    speed_resolution_assessed: bool
    speed_resolution_converged: bool
    speed_domain_assessed: bool
    speed_domain_converged: bool
    electron_heating_speed_resolution_assessed: bool
    electron_heating_speed_resolution_converged: bool
    electron_heating_speed_domain_assessed: bool
    electron_heating_speed_domain_converged: bool
    production_grid_index: int
    production_grid_matched: bool
    selected_speed_cell_count: int
    selected_speed_domain_upper_m_s: float
    failure_reason: str

@dataclass(frozen=True)
class _Eq59IterationMetrics:
    """Store convergence metrics evaluated from one accepted nonlinear iterate"""
    relative_change: float
    absolute_change: float
    solution_scale: float
    relative_residual: float
    absolute_residual: float
    inventory: float | None
    inventory_relative_change: float | None
    effective_energy_J: float | None
    effective_energy_relative_change: float | None
    fast_self_density_m3: float
    fast_self_density_relative_change: float | None
    fast_self_g_scale_m3_s3: float
    fast_self_net_drag_scale_m3_s3: float
    ion_drag_operator_relative_change: float | None
    ion_energy_diffusion_operator_relative_change: float | None
    ion_pitch_scattering_operator_relative_change: float | None
