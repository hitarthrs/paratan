#!/usr/bin/env python3
import argparse
from src.paratan.runner import run_tandem_model

def main():
    parser = argparse.ArgumentParser(description="Run Paratan model.")
    parser.add_argument("--mode", choices=["simple", "tandem"], required=True)
    parser.add_argument("--input", type=str, required=True)
    parser.add_argument("--output-dir", type=str, default="output", help="Directory to write output files")

    args = parser.parse_args()

    if args.mode == "simple":
        # run_simple_model(args.input, args.output_dir)
        print(f"Luls")
    elif args.mode == "tandem":
        run_tandem_model(args.input, args.output_dir)

if __name__ == "__main__":
    main()
