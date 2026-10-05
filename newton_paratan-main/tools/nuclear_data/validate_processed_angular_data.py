"""Validate the processed ENDF B VIII.1 fusion angular data used at runtime"""
from __future__ import annotations
import argparse
import hashlib
import json
from importlib import resources
from pathlib import Path
import sys
import numpy as np

REPOSITORY_ROOT = Path(__file__).resolve().parents[2]
SOURCE_ROOT = REPOSITORY_ROOT / "src"
if str(SOURCE_ROOT) not in sys.path:
    sys.path.insert(0, str(SOURCE_ROOT))

from source_model_revamp.nuclear_data.endf_b_viii1 import load_evaluated_fusion_angular_distribution

PACKAGE = "source_model_revamp.nuclear_data.endf_b_viii1"
PROVENANCE_FILENAME = "angular_data_provenance.json"
REACTION_FILES = {"dd_n": "ddn_angular_data.npz", "dt_n": "dtn_angular_data.npz"}

def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()

def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--normalization-tolerance", type=float, default=1.0e-5)
    arguments = parser.parse_args()
    provenance_resource = resources.files(PACKAGE).joinpath(PROVENANCE_FILENAME)
    provenance = None
    if provenance_resource.is_file():
        with resources.as_file(provenance_resource) as path:
            provenance = json.loads(path.read_text(encoding="utf-8"))
    for reaction_key, filename in REACTION_FILES.items():
        distribution = load_evaluated_fusion_angular_distribution(reaction_key)
        validation = distribution.runtime_validation
        if validation.get("angular_law_validation_passed") is not True:
            raise ValueError(f"{reaction_key} angular law validation did not pass")
        energy_min, energy_max = distribution.incident_energy_bounds_eV
        energies = np.asarray([energy_min, 0.5 * (energy_min + energy_max), energy_max], dtype=float)
        normalization_errors = [abs(distribution.normalization_integral(float(energy)) - 1.0) for energy in energies]
        if max(normalization_errors) > arguments.normalization_tolerance:
            raise ValueError(f"{reaction_key} interpolated angular normalization exceeds tolerance")
        runtime_resource = resources.files(PACKAGE).joinpath(filename)
        with resources.as_file(runtime_resource) as path:
            runtime_sha256 = _sha256(path)
        if provenance is not None:
            expected = provenance.get("reactions", {}).get(reaction_key, {}).get("runtime_sha256")
            if expected is not None and runtime_sha256 != expected:
                raise ValueError(f"{reaction_key} runtime file does not match angular_data_provenance.json")
        print(f"{reaction_key}: MAT={distribution.mat} MF={distribution.mf} MT={distribution.mt} LCT={distribution.lct} LAW={distribution.law} LANG={distribution.lang}")
        print(f"{reaction_key}: energy coverage {energy_min:.12g} to {energy_max:.12g} eV")
        print(f"{reaction_key}: maximum checked normalization error {max(normalization_errors):.12g}")
        print(f"{reaction_key}: runtime sha256 {runtime_sha256}")
    print("processed ENDF angular data validation passed")

if __name__ == "__main__":
    main()
