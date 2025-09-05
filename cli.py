#!/usr/bin/env python3
import argparse
from src.paratan.runner import run_tandem_model, run_simple_model_modular, build_simple_model_only, build_tandem_model_only

def main():
    parser = argparse.ArgumentParser(description="Run Paratan model.")
    parser.add_argument("--mode", choices=["simple", "tandem", "simple-build-only", "tandem-build-only"], required=True)
    parser.add_argument("--input", "-i", type=str, required=True, help="Input YAML file path")
    parser.add_argument("--output-dir", "-o", type=str, default="output", help="Directory to write output files")

    args = parser.parse_args()

    if args.mode == "simple":
        run_simple_model_modular(args.input, args.output_dir)
    elif args.mode == "tandem":
        run_tandem_model(args.input, args.output_dir)
    elif args.mode == "simple-build-only":
        build_simple_model_only(args.input, args.output_dir)
        print("✅ Model built successfully! XML files and plots created in:", args.output_dir)
        print("💡 To run the simulation, use: python cli.py --mode simple --input", args.input, "--output-dir", args.output_dir)
    elif args.mode == "tandem-build-only":
        build_tandem_model_only(args.input, args.output_dir)
        print("✅ Model built successfully! XML files and plots created in:", args.output_dir)
        print("💡 To run the simulation, use: python cli.py --mode tandem --input", args.input, "--output-dir", args.output_dir)

if __name__ == "__main__":
    main()
