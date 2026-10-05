"""Canonical configuration model identifiers and alias tables"""
PRODUCTION_WORKFLOW_MODEL = "egedal_modal_fbis_openmc_source"
MODEL_ALIASES = {PRODUCTION_WORKFLOW_MODEL: PRODUCTION_WORKFLOW_MODEL, "single_beam_egedal_modal_fbis_openmc_source": PRODUCTION_WORKFLOW_MODEL}
GEOMETRY_MODEL_ALIASES = {"paratan_coil_field": "paratan_coil_field", "egedal_analytic": "egedal_analytic"}
BEAM_MODEL_ALIASES = {"geometry_linked_attenuation": "geometry_linked_attenuation"}
NBI_SUPPORTED_STATIONARY_PLASMA_CLOSURE_MODEL = "nbi_supported_stationary"
PLASMA_CLOSURE_MODEL_ALIASES = {NBI_SUPPORTED_STATIONARY_PLASMA_CLOSURE_MODEL: NBI_SUPPORTED_STATIONARY_PLASMA_CLOSURE_MODEL, 'egedal_beam_plasma_quasineutral': 'egedal_beam_plasma_quasineutral'}
QUASINEUTRAL_ELECTRON_DENSITY_MODE = "quasineutral_from_total_positive_charge"
ELECTRON_DENSITY_MODE_ALIASES = {QUASINEUTRAL_ELECTRON_DENSITY_MODE: QUASINEUTRAL_ELECTRON_DENSITY_MODE}
MODAL_DENSITY_CLOSURE_ALIASES = {NBI_SUPPORTED_STATIONARY_PLASMA_CLOSURE_MODEL: NBI_SUPPORTED_STATIONARY_PLASMA_CLOSURE_MODEL, 'egedal_beam_plasma_quasineutral': 'egedal_beam_plasma_quasineutral'}
POWER_BALANCE_MODE_ALIASES = {"fixed_closure": "fixed_closure", "self_consistent_electron_energy": "self_consistent_electron_energy"}
NEUTRON_MODEL_ALIASES = {"distribution_kinematics_spectrum": "distribution_kinematics_spectrum"}
EXPORT_MODEL_ALIASES = {"axisymmetric_volume_bins": "axisymmetric_volume_bins"}
FIELD_STRENGTH_FIT_MODEL = "fit_lf_hf_scaling_to_midplane_and_mirror_fields"
OPENMC_EXPORT_SPATIAL_MODEL = "axisymmetric_cylindrical_volume_bins"
