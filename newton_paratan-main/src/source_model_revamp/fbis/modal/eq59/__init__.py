"""Egedal Eq 59 velocity space solver"""
from source_model_revamp.fbis.modal.eq59.types import Eq59ConvergenceDiagnostics
from source_model_revamp.fbis.modal.eq59.types import Eq59WarmStartState
from source_model_revamp.fbis.modal.eq59.types import Eq59SpeedConvergenceDiagnostics
from source_model_revamp.fbis.modal.eq59.collision_operator_state import Eq59PairOperatorContribution, Eq59CollisionOperatorState, build_eq59_collision_operator_state, eq59_collision_operator_metadata
from source_model_revamp.fbis.modal.eq59.maxwellian_pairs import maxwellian_rosenbluth_coefficients
from source_model_revamp.fbis.modal.eq59.energy_moments import Eq59PairEnergyMoment, FastIonElectronHeatingState, build_fast_ion_electron_heating_state, fast_ion_electron_heating_metadata
from source_model_revamp.fbis.modal.eq59.operator import _scharfetter_gummel_face_coefficients
from source_model_revamp.fbis.modal.eq59.operator import _assemble_eq59_mode_tridiagonal
from source_model_revamp.fbis.modal.eq59.operator import _solve_eq59_mode_finite_difference
from source_model_revamp.fbis.modal.eq59.state import build_eq59_warm_start_state
from source_model_revamp.fbis.modal.eq59.solve import hot_ion_rosenbluth_modal_distribution_eq59
from source_model_revamp.fbis.modal.eq59.convergence import build_eq59_speed_convergence_grid_sequence
from source_model_revamp.fbis.modal.eq59.convergence import assess_eq59_speed_convergence

__all__ = [
    'Eq59ConvergenceDiagnostics',
    'Eq59WarmStartState',
    'Eq59SpeedConvergenceDiagnostics',
    'Eq59PairOperatorContribution',
    'Eq59CollisionOperatorState',
    'build_eq59_collision_operator_state',
    'eq59_collision_operator_metadata',
    'maxwellian_rosenbluth_coefficients',
    'Eq59PairEnergyMoment',
    'FastIonElectronHeatingState',
    'build_fast_ion_electron_heating_state',
    'fast_ion_electron_heating_metadata',
    '_scharfetter_gummel_face_coefficients',
    '_assemble_eq59_mode_tridiagonal',
    '_solve_eq59_mode_finite_difference',
    'build_eq59_warm_start_state',
    'hot_ion_rosenbluth_modal_distribution_eq59',
    'build_eq59_speed_convergence_grid_sequence',
    'assess_eq59_speed_convergence',
]
