"""Top level orchestration for the magnetic mirror source model"""
from __future__ import annotations
from dataclasses import replace
from time import perf_counter
import numpy as np
from source_model_revamp.fbis.modal.types import ModalPhysicalDistributionError
from source_model_revamp.integration.pipeline_reporting import represented_population_scope
from source_model_revamp.integration.pipeline_stages.fusion_stage import build_fusion_stage
from source_model_revamp.integration.pipeline_stages.geometry_stage import build_geometry_stage
from source_model_revamp.integration.pipeline_stages.metadata_stage import build_pipeline_metadata
from source_model_revamp.integration.pipeline_stages.neutron_stage import build_neutron_stage
from source_model_revamp.integration.pipeline_stages.openmc_export_stage import build_openmc_export_stage
from source_model_revamp.integration.pipeline_stages.temperature_closure_stage import build_temperature_resolved_operating_point
from source_model_revamp.integration.assessment import build_run_assessment
from source_model_revamp.integration.validation import NEUTRON_OUT_OF_RANGE_FRACTION_TOLERANCE, SourceConstructionError, validate_correlated_event_state, validate_current_identities, validate_density_state, validate_electrostatic_state, validate_expander_state, validate_fusion_neutron_state, validate_geometry_state, validate_kinetic_state, validate_neutron_domain, validate_openmc_file_source, validate_source_constructibility
from source_model_revamp.integration.pipeline_types import SourceModelPipelineResult
from source_model_revamp.integration.config.root import SourceModelRunConfig
from source_model_revamp.pipeline_errors import InvalidPipelineStateError, NeutronEnergyDomainError
from source_model_revamp.runtime_profile import finalize_runtime_profile, new_runtime_profile, record_runtime

def run_source_model_pipeline(config: SourceModelRunConfig, *, write_openmc_source: bool = True, progress_callback=None) -> SourceModelPipelineResult:
    """
    Run geometry, temperature resolved operating point, fusion, neutron, correlated source validation, OpenMC export, assessment, population scope reporting, and metadata assembly
    
    Direct state validators raise stage specific InvalidPipelineStateError before the compact RunAssessment is formed
    Runtime for each major stage is recorded in the final metadata
    """
    pipeline_runtime_started = perf_counter()
    runtime_profile = new_runtime_profile()

    stage_started = perf_counter()
    geometry = build_geometry_stage(config)
    record_runtime(runtime_profile, "geometry", perf_counter() - stage_started)

    stage_started = perf_counter()
    # Reject invalid calculated states immediately while numerical qualification remains an assessment concern
    try:
        validate_geometry_state(geometry)
    except (AttributeError, TypeError, ValueError) as exc:
        raise InvalidPipelineStateError(str(exc), stage="geometry") from exc
    record_runtime(runtime_profile, "geometry_validation", perf_counter() - stage_started)

    stage_started = perf_counter()
    try:
        operating_point = build_temperature_resolved_operating_point(config, geometry, progress_callback=progress_callback)
    except ModalPhysicalDistributionError as exc:
        raise InvalidPipelineStateError(str(exc), stage="operating_point") from exc
    record_runtime(runtime_profile, "operating_point", perf_counter() - stage_started)

    beam = operating_point.beam
    kinetic = operating_point.kinetic
    expander = operating_point.expander
    if expander is None or expander.system_state is None:
        raise InvalidPipelineStateError("temperature resolved operating point did not return the coupled expander state", stage="operating_point")

    stage_started = perf_counter()
    try:
        validate_kinetic_state(kinetic)
        validate_density_state(geometry, kinetic, beam, identity_tolerance=float(kinetic.metadata["density_identity_relative_tolerance"]))
        validate_electrostatic_state(kinetic)
        validate_current_identities({**kinetic.metadata, **expander.metadata})
        validate_expander_state(expander)
    except (AttributeError, KeyError, TypeError, ValueError) as exc:
        raise InvalidPipelineStateError(str(exc), stage="operating_point") from exc
    record_runtime(runtime_profile, "operating_point_validation", perf_counter() - stage_started)

    stage_started = perf_counter()
    fusion = build_fusion_stage(config, geometry, kinetic, expander)
    record_runtime(runtime_profile, "fusion", perf_counter() - stage_started)

    stage_started = perf_counter()

    stage_started = perf_counter()
    try:
        neutrons = build_neutron_stage(config, geometry, kinetic, fusion)
    except NeutronEnergyDomainError as exc:
        raise InvalidPipelineStateError(str(exc), stage="neutrons") from exc
    record_runtime(runtime_profile, "neutrons", perf_counter() - stage_started)

    stage_started = perf_counter()
    event_bank = neutrons.correlated_event_bank
    if event_bank is None and (not np.isfinite(neutrons.total_neutron_rate_s) or neutrons.total_neutron_rate_s <= 0.0):
        raise SourceConstructionError("invalid_total_rate", "must be positive and finite", quantity="total_neutron_rate_s", value=neutrons.total_neutron_rate_s)
    validate_source_constructibility(event_bank)
    try:
        validate_fusion_neutron_state(fusion, neutrons, relative_tolerance=float(neutrons.metadata["neutron_fusion_to_spectrum_conservation_tolerance"]))
        validate_neutron_domain(neutrons.metadata, fraction_tolerance=NEUTRON_OUT_OF_RANGE_FRACTION_TOLERANCE)
        validate_correlated_event_state(event_bank, expected_physical_rate_s=float(neutrons.metadata["neutron_scalar_table_iv_fusion_rate_s"]), relative_tolerance=5.0e-13, mass_shell_tolerance=1.0e-10)
    except (AttributeError, KeyError, TypeError, ValueError) as exc:
        raise InvalidPipelineStateError(str(exc), stage="neutrons") from exc
    record_runtime(runtime_profile, "neutron_validation", perf_counter() - stage_started)

    stage_started = perf_counter()
    export = build_openmc_export_stage(config, neutrons, write_source=write_openmc_source)
    record_runtime(runtime_profile, "openmc_export", perf_counter() - stage_started)
    if write_openmc_source and export.file_source_bundle is None:
        raise SourceConstructionError("missing_correlated_file_source", "required correlated FileSource representation was not created", quantity="openmc_file_source_bundle")

    stage_started = perf_counter()
    if export.file_source_bundle is not None:
        try:
            validate_openmc_file_source(export.file_source_bundle.metadata, probability_tolerance=2.0e-12)
        except (AttributeError, TypeError, ValueError) as exc:
            raise InvalidPipelineStateError(str(exc), stage="openmc_export") from exc
    record_runtime(runtime_profile, "openmc_export_validation", perf_counter() - stage_started)

    stage_started = perf_counter()
    # Build convergence and applicability status only after direct physical state validation succeeds
    assessment = build_run_assessment(operating_point=operating_point, kinetic=kinetic, expander=expander, fusion=fusion)
    record_runtime(runtime_profile, "assessment", perf_counter() - stage_started)

    stage_started = perf_counter()
    included, excluded, limitations = represented_population_scope(config, kinetic)
    record_runtime(runtime_profile, "population_scope", perf_counter() - stage_started)

    partial = SourceModelPipelineResult(config=config, geometry=geometry, beam=beam, kinetic=kinetic, expander=expander, fusion=fusion, neutrons=neutrons, openmc_export=export, metadata={}, operating_point=operating_point, assessment=assessment)
    stage_started = perf_counter()
    metadata = build_pipeline_metadata(partial)
    metadata.update({
        "represented_populations": list(included),
        "excluded_populations": list(excluded),
        "limitations": list(limitations),
        "result_status": assessment.result_status,
        "calculation_valid": True,
        "numerically_converged": assessment.numerically_converged,
        "failed_convergence_checks": list(assessment.failed_convergence_checks),
        "model_applicable": assessment.model_applicable,
        "applicability_limitations": list(assessment.applicability_limitations),
        "source_constructible": True,
        "source_generated": export.file_source_bundle is not None,
        "source_file_valid": True if export.file_source_bundle is not None else None,
        "convergence_checks": [record.to_metadata() for record in assessment.convergence],
        "applicability_checks": [record.to_metadata() for record in assessment.applicability],
    })
    record_runtime(runtime_profile, "metadata_assembly", perf_counter() - stage_started)
    metadata["runtime_pipeline_profile"] = finalize_runtime_profile(runtime_profile, total_s=perf_counter() - pipeline_runtime_started)
    return replace(partial, metadata=metadata)

__all__ = ["run_source_model_pipeline"]
