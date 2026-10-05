"""Exceptions shared by pipeline stages and low level exporters"""
from __future__ import annotations

class PipelineError(RuntimeError):
    """Base class for typed integration pipeline failures"""

class InvalidPipelineStateError(PipelineError):
    """A hard mathematical or data validity condition failed"""
    def __init__(self, message: str, *, stage: str) -> None:
        super().__init__(message)
        self.stage = stage

class PipelineConservationError(PipelineError):
    """A stage detected a mathematically invalid conservation defect"""
    def __init__(self, message: str, *, stage: str, check_name: str) -> None:
        super().__init__(message)
        self.stage = stage
        self.check_name = check_name

class NeutronEnergyDomainError(ValueError):
    """A neutron histogram domain excludes weighted events"""

class OpenMCExportDependencyError(PipelineError):
    """An optional dependency required only for HDF5 export is unavailable"""
