"""Command line interface for source model diagnostic plotting"""
from __future__ import annotations
import argparse
from source_model_revamp.plotting.common import DEFAULT_LOCAL_B_OVER_B0
from source_model_revamp.plotting.registry import run_plots

def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-dir", required=True, help="Run directory containing source_model_metadata.json")
    parser.add_argument("--output-dir", default=None, help="Plot output directory, default is <run-dir>/plots")
    parser.add_argument("--groups", nargs="+", default=["all"], help="Plot groups or names to generate")
    parser.add_argument("--seed", type=int, default=12345, help="Random seed for sampled plotting")
    parser.add_argument("--local-b-over-b0", nargs="+", type=float, default=list(DEFAULT_LOCAL_B_OVER_B0), help="Target B over B0 values for local velocity plots")
    parser.add_argument("--max-modes", type=int, default=12, help="Maximum number of modal eigenfunctions to plot")
    parser.add_argument("--keep-going", action="store_true", default=True, help="Continue after individual plot failures")
    parser.add_argument("--better-cloud-elevation", type=float, default=20.0, help="3D source cloud elevation angle in degrees")
    parser.add_argument("--better-cloud-azimuth", type=float, default=-62.0, help="3D source cloud azimuth angle in degrees")
    parser.add_argument("--better-cloud-samples", type=int, default=100000, help="Probability weighted points shown in the source cloud")
    parser.add_argument("--better-cloud-radial-exaggeration", type=float, default=4.0, help="Display aspect multiplier for source cloud transverse axes")
    parser.add_argument("--better-cloud-color", choices=("energy", "rate_density"), default="energy", help="Quantity used to color the source cloud")
    parser.add_argument("--better-energy-bins", type=int, default=1000, help="Energy bins used in neutron spectra")
    parser.add_argument("--direction-reference", choices=("machine", "beam"), default="machine", help="Reference axis for neutron direction plots")
    parser.add_argument("--direction-beam-id", default=None, help="Beam ID used when direction reference is beam")
    parser.add_argument("--dt-spectrum-bins", type=int, default=220, help="Energy bins used in the component resolved DT spectrum")
    parser.add_argument("--dt-angle-bins", type=int, default=90, help="Angular bins used in the DT PDF and CDF")
    return parser

def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    manifest = run_plots(args)
    made = sum(1 for item in manifest if item.get("status") == "made")
    skipped = sum(1 for item in manifest if item.get("status") == "skipped")
    failed = sum(1 for item in manifest if item.get("status") == "failed")
    print(f"plots made: {made}, skipped: {skipped}, failed: {failed}")
    
    return 1 if failed else 0
