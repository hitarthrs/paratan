"""
Correlated neutron event bank NPZ sidecar input and output

The sidecar retains correlated position direction energy probability and reaction audit arrays together with the separate physical neutron rate
"""
from __future__ import annotations
import json
import os
from pathlib import Path
import tempfile
from typing import Any, Mapping
import numpy as np
from source_model_revamp.neutrons.events import CorrelatedNeutronEventBank

FORMAT_NAME = "source_model_revamp.correlated_neutron_event_bank"
FORMAT_VERSION = 1
DEFAULT_FILENAME = "correlated_neutron_events.npz"
_REQUIRED_ARRAYS = ("positions_m", "directions", "energies_J", "normalized_weights", "reaction_keys", "component_labels", "axial_cell_indices", "center_of_mass_energy_J", "equivalent_deuteron_lab_energy_eV", "cm_emission_mu", "reactant_a_velocity_m_s", "reactant_b_velocity_m_s")

def _json_safe(value: Any) -> Any:
    """
    Convert metadata values to JSON compatible containers and scalar types
    
    Nonfinite floating values are represented by JSON null while numerical event arrays remain stored separately in the NPZ archive
    """
    if isinstance(value, Mapping):
        return {str(key): _json_safe(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_safe(item) for item in value]
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, float) and not np.isfinite(value):
        return None
    return value

def _header(event_bank: CorrelatedNeutronEventBank) -> dict[str, Any]:
    """
    Build the versioned sidecar header for one correlated neutron event bank
    
    The header records units coordinate order event count physical rate probability role and JSON compatible event metadata
    """
    return {
        "format": FORMAT_NAME,
        "format_version": FORMAT_VERSION,
        "physical_total_rate_s": float(event_bank.physical_total_rate_s),
        "event_count": int(event_bank.event_count),
        "coordinate_order": "x_y_z",
        "position_units": "m",
        "direction_units": "unit_vector",
        "energy_units": "J",
        "probability_weight_role": "normalized_event_probability",
        "metadata": _json_safe(event_bank.metadata),
    }

def write_correlated_neutron_event_bank_npz(event_bank: CorrelatedNeutronEventBank, path: str | Path) -> Path:
    """
    Write one correlated event bank to a compressed NPZ sidecar and return its path
    
    The archive is written to a temporary file read back through the companion loader and then atomically replaces the requested destination
    """
    output_path = Path(path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(dir=output_path.parent, prefix=f".{output_path.stem}.", suffix=".npz", delete=False) as temporary:
        temporary_path = Path(temporary.name)
    try:
        np.savez_compressed(
            temporary_path,
            positions_m=np.asarray(event_bank.positions_m, dtype=float),
            directions=np.asarray(event_bank.directions, dtype=float),
            energies_J=np.asarray(event_bank.energies_J, dtype=float),
            normalized_weights=np.asarray(event_bank.normalized_weights, dtype=float),
            reaction_keys=np.asarray(event_bank.reaction_keys, dtype="U8"),
            component_labels=np.asarray(event_bank.component_labels, dtype="U96"),
            axial_cell_indices=np.asarray(event_bank.axial_cell_indices, dtype=np.int64),
            center_of_mass_energy_J=np.asarray(event_bank.center_of_mass_energy_J, dtype=float),
            equivalent_deuteron_lab_energy_eV=np.asarray(event_bank.equivalent_deuteron_lab_energy_eV, dtype=float),
            cm_emission_mu=np.asarray(event_bank.cm_emission_mu, dtype=float),
            reactant_a_velocity_m_s=np.asarray(event_bank.reactant_a_velocity_m_s, dtype=float),
            reactant_b_velocity_m_s=np.asarray(event_bank.reactant_b_velocity_m_s, dtype=float),
            header_json=json.dumps(_header(event_bank), sort_keys=True, separators=(",", ":")),
        )
        read_correlated_neutron_event_bank_npz(temporary_path)
        os.replace(temporary_path, output_path)
    finally:
        if temporary_path.exists():
            temporary_path.unlink()
    return output_path

def read_correlated_neutron_event_bank_npz(path: str | Path,) -> CorrelatedNeutronEventBank:
    """
    Read and validate one correlated neutron event bank NPZ sidecar
    
    The loader checks format identity version units required arrays metadata type and event count before returning CorrelatedNeutronEventBank
    """
    source_path = Path(path)
    with np.load(source_path, allow_pickle=False) as archive:
        if "header_json" not in archive.files:
            raise ValueError("correlated event bank sidecar is missing header_json")
        missing = [name for name in _REQUIRED_ARRAYS if name not in archive.files]
        if missing:
            raise ValueError("correlated event bank sidecar is missing " + ", ".join(missing))
        header = json.loads(str(archive["header_json"].item()))
        if header.get("format") != FORMAT_NAME:
            raise ValueError("correlated event bank sidecar format identity is invalid")
        if int(header.get("format_version", -1)) != FORMAT_VERSION:
            raise ValueError("correlated event bank sidecar version is unsupported")
        if header.get("coordinate_order") != "x_y_z":
            raise ValueError("correlated event bank coordinate order is unsupported")
        if header.get("position_units") != "m":
            raise ValueError("correlated event bank position units are unsupported")
        if header.get("direction_units") != "unit_vector":
            raise ValueError("correlated event bank direction units are unsupported")
        if header.get("energy_units") != "J":
            raise ValueError("correlated event bank energy units are unsupported")
        if header.get("probability_weight_role") != "normalized_event_probability":
            raise ValueError("correlated event bank probability weight role is unsupported")
        metadata = header.get("metadata", {})
        if not isinstance(metadata, Mapping):
            raise ValueError("correlated event bank metadata must be a mapping")
        bank = CorrelatedNeutronEventBank(
            positions_m=np.asarray(archive["positions_m"], dtype=float),
            directions=np.asarray(archive["directions"], dtype=float),
            energies_J=np.asarray(archive["energies_J"], dtype=float),
            normalized_weights=np.asarray(archive["normalized_weights"], dtype=float),
            physical_total_rate_s=float(header["physical_total_rate_s"]),
            reaction_keys=np.asarray(archive["reaction_keys"]).astype("U8"),
            component_labels=np.asarray(archive["component_labels"]).astype("U96"),
            axial_cell_indices=np.asarray(archive["axial_cell_indices"], dtype=np.int64),
            center_of_mass_energy_J=np.asarray(archive["center_of_mass_energy_J"], dtype=float),
            equivalent_deuteron_lab_energy_eV=np.asarray(archive["equivalent_deuteron_lab_energy_eV"], dtype=float),
            cm_emission_mu=np.asarray(archive["cm_emission_mu"], dtype=float),
            reactant_a_velocity_m_s=np.asarray(archive["reactant_a_velocity_m_s"], dtype=float),
            reactant_b_velocity_m_s=np.asarray(archive["reactant_b_velocity_m_s"], dtype=float),
            metadata=dict(metadata),
        )
    if int(header.get("event_count", -1)) != bank.event_count:
        raise ValueError("correlated event bank event count does not match its header")
    
    return bank