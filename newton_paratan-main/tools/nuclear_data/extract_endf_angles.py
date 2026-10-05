"""Extract DD and DT neutron angular laws from official ENDF B VIII.1 files"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import sys
import zipfile
import numpy as np
from numpy.polynomial.legendre import legval

REPOSITORY_ROOT = Path(__file__).resolve().parents[2]
SOURCE_ROOT = REPOSITORY_ROOT / "src"
if str(SOURCE_ROOT) not in sys.path:
    sys.path.insert(0, str(SOURCE_ROOT))

from endf6 import EndfLaw2Knot, parse_material_header, parse_mf3_section, parse_mf6_section

REACTIONS = {
    "dd_n": {"reaction_label": "D(d,n)3He", "target": "H-2", "mat": 128, "dat_filename": "d_001-H-2_0128.dat", "zip_filename": "d_001-H-2_0128.zip", "expected_lang": 12, "expected_residual_zap": 2003, "runtime_filename": "ddn_angular_data.npz"},
    "dt_n": {"reaction_label": "T(d,n)4He", "target": "H-3", "mat": 131, "dat_filename": "d_001-H-3_0131.dat", "zip_filename": "d_001-H-3_0131.zip", "expected_lang": 0, "expected_residual_zap": 2004, "runtime_filename": "dtn_angular_data.npz"},
}

def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()

def _verify_archive(dat_path: Path, zip_path: Path) -> dict[str, object]:
    if not dat_path.is_file():
        raise FileNotFoundError(dat_path)
    if not zip_path.is_file():
        raise FileNotFoundError(zip_path)
    with zipfile.ZipFile(zip_path) as archive:
        names = archive.namelist()
        if names != [dat_path.name]:
            raise ValueError(f"{zip_path.name} must contain only {dat_path.name}, found {names}")
        archived_bytes = archive.read(dat_path.name)
    raw_bytes = dat_path.read_bytes()
    if archived_bytes != raw_bytes:
        raise ValueError(f"{zip_path.name} does not contain the supplied DAT bytes")
    return {
        "dat_filename": dat_path.name,
        "dat_size_bytes": dat_path.stat().st_size,
        "zip_filename": zip_path.name,
        "zip_size_bytes": zip_path.stat().st_size,
        "archive_member": dat_path.name,
        "archive_matches_dat": True,
        "dat_sha256": _sha256(dat_path),
        "zip_sha256": _sha256(zip_path),
    }

def _validate_common(reaction_key: str, specification: dict[str, object], raw_dir: Path) -> tuple[dict[str, object], dict[str, np.ndarray]]:
    dat_path = raw_dir / str(specification["dat_filename"])
    zip_path = raw_dir / str(specification["zip_filename"])
    archive_info = _verify_archive(dat_path, zip_path)
    mat = int(specification["mat"])
    material = parse_material_header(dat_path, mat=mat)
    mf3 = parse_mf3_section(dat_path, mat=mat, mt=50)
    mf6 = parse_mf6_section(dat_path, mat=mat, mt=50)
    if mf6.reference_frame != 2:
        raise ValueError(f"{reaction_key} MF=6 MT=50 must use LCT=2")
    if len(mf6.products) != 2:
        raise ValueError(f"{reaction_key} MF=6 MT=50 must contain two products")
    neutron = mf6.products[0]
    residual = mf6.products[1]
    incident_energy = np.asarray([knot.incident_energy_eV for knot in neutron.angular_knots], dtype=float)
    langs = {knot.lang for knot in neutron.angular_knots}
    if langs != {int(specification["expected_lang"])}:
        raise ValueError(f"{reaction_key} unexpected LANG values {sorted(langs)}")
    common_arrays = {
        "incident_energy_eV": incident_energy,
        "energy_interpolation_breakpoints": np.asarray(neutron.angular_interpolation.breakpoints, dtype=np.int64),
        "energy_interpolation_laws": np.asarray(neutron.angular_interpolation.laws, dtype=np.int64),
        "yield_incident_energy_eV": np.asarray(neutron.yield_incident_energy_eV, dtype=float),
        "yield_values": np.asarray(neutron.yield_values, dtype=float),
        "yield_interpolation_breakpoints": np.asarray(neutron.yield_interpolation.breakpoints, dtype=np.int64),
        "yield_interpolation_laws": np.asarray(neutron.yield_interpolation.laws, dtype=np.int64),
    }
    audit = {
        "reaction_key": reaction_key,
        "reaction_label": specification["reaction_label"],
        "library": "ENDF/B-VIII.1",
        "release_date": "2024-08-30",
        "sublibrary": "incident deuteron",
        "target": specification["target"],
        "projectile": "deuteron",
        "mat": mat,
        "mf": 6,
        "mt": 50,
        "reference_frame": "center_of_mass",
        "incident_energy_frame": "deuteron_lab_on_stationary_target",
        "angular_cosine_definition": "neutron_direction_relative_to_incident_deuteron_direction_in_center_of_mass",
        "azimuthal_probability_model": "uniform_about_incident_deuteron_axis",
        "lct": mf6.reference_frame,
        "law": neutron.law,
        "lang": int(specification["expected_lang"]),
        "neutron_zap": neutron.zap,
        "neutron_awp": neutron.awp,
        "residual_zap": residual.zap,
        "residual_awp": residual.awp,
        "residual_law": residual.law,
        "product_yield": 1.0,
        "incident_energy_min_eV": float(incident_energy[0]),
        "incident_energy_max_eV": float(incident_energy[-1]),
        "incident_energy_point_count": int(incident_energy.size),
        "incident_energy_interpolation_breakpoints": list(neutron.angular_interpolation.breakpoints),
        "incident_energy_interpolation_laws": list(neutron.angular_interpolation.laws),
        "mf3_mass_difference_Q_eV": mf3.mass_difference_Q_eV,
        "mf3_reaction_Q_eV": mf3.reaction_Q_eV,
        "mf3_breakup_flag": mf3.breakup_flag,
        "mf3_incident_energy_min_eV": mf3.incident_energy_eV[0],
        "mf3_incident_energy_max_eV": mf3.incident_energy_eV[-1],
        "mf3_point_count": len(mf3.incident_energy_eV),
        "mf3_interpolation_breakpoints": list(mf3.interpolation.breakpoints),
        "mf3_interpolation_laws": list(mf3.interpolation.laws),
        "material_za": material.za,
        "material_awr": material.awr,
        "projectile_awr": material.projectile_awr,
        "material_header_emax_eV": material.material_emax_eV,
        "reaction_range_exceeds_material_header_emax": bool(incident_energy[-1] > material.material_emax_eV),
        "material_header_comments": list(material.comments),
        "total_cross_section_role": "audit_only",
        "production_total_cross_section_model": "Bosch_Hale_Table_IV",
        "angular_data_role": "normalized_conditional_probability_in_mu_cm",
        **archive_info,
    }
    return audit, common_arrays

def _extract_tabulated(knots: tuple[EndfLaw2Knot, ...]) -> tuple[dict[str, np.ndarray], dict[str, object]]:
    offsets = [0]
    mu_values: list[float] = []
    probability_values: list[float] = []
    normalization_errors: list[float] = []
    minima: list[float] = []
    point_counts: list[int] = []
    for knot in knots:
        if len(knot.values) != 2 * knot.item_count:
            raise ValueError("LANG=12 LIST length does not match NL")
        mu = np.asarray(knot.values[0::2], dtype=float)
        probability = np.asarray(knot.values[1::2], dtype=float)
        integral = float(np.trapezoid(probability, mu))
        normalization_errors.append(abs(integral - 1.0))
        minima.append(float(np.min(probability)))
        point_counts.append(int(mu.size))
        mu_values.extend(mu.tolist())
        probability_values.extend(probability.tolist())
        offsets.append(len(mu_values))
    arrays = {
        "representation_code": np.asarray(12, dtype=np.int64),
        "angular_offsets": np.asarray(offsets, dtype=np.int64),
        "angular_mu": np.asarray(mu_values, dtype=float),
        "angular_probability_density_per_mu": np.asarray(probability_values, dtype=float),
        "angular_item_count": np.asarray(point_counts, dtype=np.int64),
    }
    audit = {
        "representation": "tabulated_probability_linear_in_mu",
        "angular_point_count_minimum": min(point_counts),
        "angular_point_count_maximum": max(point_counts),
        "maximum_knot_normalization_error": max(normalization_errors),
        "minimum_knot_probability_density_per_mu": min(minima),
    }
    return arrays, audit

def _extract_legendre(knots: tuple[EndfLaw2Knot, ...]) -> tuple[dict[str, np.ndarray], dict[str, object]]:
    order_count = np.asarray([knot.item_count for knot in knots], dtype=np.int64)
    maximum_order = int(np.max(order_count))
    coefficients = np.zeros((len(knots), maximum_order), dtype=float)
    minimum_probability = np.inf
    maximum_normalization_error = 0.0
    mu_check = np.linspace(-1.0, 1.0, 20001)
    for index, knot in enumerate(knots):
        if len(knot.values) != knot.item_count:
            raise ValueError("LANG=0 LIST length does not match NL")
        coefficients[index, : knot.item_count] = knot.values
        polynomial_coefficients = np.zeros(knot.item_count + 1, dtype=float)
        polynomial_coefficients[0] = 0.5
        orders = np.arange(1, knot.item_count + 1, dtype=float)
        polynomial_coefficients[1:] = 0.5 * (2.0 * orders + 1.0) * np.asarray(knot.values, dtype=float)
        probability = legval(mu_check, polynomial_coefficients)
        minimum_probability = min(minimum_probability, float(np.min(probability)))
        maximum_normalization_error = max(maximum_normalization_error, abs(float(np.trapezoid(probability, mu_check)) - 1.0))
    arrays = { "representation_code": np.asarray(0, dtype=np.int64), "legendre_coefficients": coefficients, "legendre_order_count": order_count}
    audit = {
        "representation": "legendre_coefficients",
        "legendre_order_minimum": int(np.min(order_count)),
        "legendre_order_maximum": maximum_order,
        "maximum_knot_normalization_error_on_dense_grid": maximum_normalization_error,
        "minimum_knot_probability_density_per_mu_on_dense_grid": minimum_probability,
    }
    return arrays, audit

def extract_reaction(reaction_key: str, specification: dict[str, object], raw_dir: Path, output_dir: Path) -> dict[str, object]:
    audit, arrays = _validate_common(reaction_key, specification, raw_dir)
    dat_path = raw_dir / str(specification["dat_filename"])
    mf6 = parse_mf6_section(dat_path, mat=int(specification["mat"]), mt=50)
    knots = mf6.products[0].angular_knots
    if int(specification["expected_lang"]) == 12:
        angular_arrays, angular_audit = _extract_tabulated(knots)
    else:
        angular_arrays, angular_audit = _extract_legendre(knots)
    arrays.update(angular_arrays)
    arrays.update(
        {
            "mat": np.asarray(int(specification["mat"]), dtype=np.int64),
            "mf": np.asarray(6, dtype=np.int64),
            "mt": np.asarray(50, dtype=np.int64),
            "lct": np.asarray(2, dtype=np.int64),
            "law": np.asarray(2, dtype=np.int64),
            "lang": np.asarray(int(specification["expected_lang"]), dtype=np.int64),
            "q_value_eV": np.asarray(audit["mf3_reaction_Q_eV"], dtype=float),
            "target_awr": np.asarray(audit["material_awr"], dtype=float),
            "projectile_awr": np.asarray(audit["projectile_awr"], dtype=float),
            "neutron_awp": np.asarray(audit["neutron_awp"], dtype=float),
            "residual_awp": np.asarray(audit["residual_awp"], dtype=float),
            "residual_zap": np.asarray(audit["residual_zap"], dtype=np.int64),
        }
    )
    output_dir.mkdir(parents=True, exist_ok=True)
    runtime_path = output_dir / str(specification["runtime_filename"])
    np.savez_compressed(runtime_path, **arrays)
    audit["runtime_filename"] = runtime_path.name
    audit["runtime_size_bytes"] = runtime_path.stat().st_size
    audit["runtime_sha256"] = _sha256(runtime_path)

    return audit

def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--raw-dir", type=Path, default=REPOSITORY_ROOT / "tools/nuclear_data/raw/endf_b_viii1_deuterons")
    parser.add_argument("--output-dir", type=Path, default=REPOSITORY_ROOT / "src/source_model_revamp/nuclear_data/endf_b_viii1")
    parser.add_argument("--provenance-file", type=Path, default=None)
    arguments = parser.parse_args()
    audits = {reaction_key: extract_reaction(reaction_key, specification, arguments.raw_dir, arguments.output_dir) for reaction_key, specification in REACTIONS.items()}
    provenance_path = arguments.provenance_file or arguments.output_dir / "angular_data_provenance.json"
    provenance_path.parent.mkdir(parents=True, exist_ok=True)
    with provenance_path.open("w", encoding="utf-8") as handle:
        json.dump({"library": "ENDF/B-VIII.1", "reactions": audits}, handle, indent=2, sort_keys=True)
        handle.write("\n")
    for reaction_key, audit in audits.items():
        print(f"{reaction_key}: {audit['runtime_filename']} {audit['runtime_sha256']}")
    print(f"provenance: {provenance_path}")

if __name__ == "__main__":
    main()
