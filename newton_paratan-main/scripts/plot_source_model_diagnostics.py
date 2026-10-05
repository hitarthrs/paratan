from source_model_revamp.plotting.cli import build_parser, main
from source_model_revamp.plotting.registry import load_metadata, run_plots

__all__ = ["build_parser", "load_metadata", "run_plots", "main"]

if __name__ == "__main__":
    main()
