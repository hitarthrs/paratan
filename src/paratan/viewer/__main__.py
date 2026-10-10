"""Paratan interactive geometry viewer.

Usage
-----
    python -m src.paratan.viewer --input input_files/simple_parametric_input_new.yaml

Or from the repo root after setting PYTHONPATH=.
"""
from __future__ import annotations

import argparse
from pathlib import Path


def main(argv: list[str] | None = None) -> int:
    # GL fix BEFORE any pyvista/VTK import (app.py pulls those in).
    from src.paratan.viewer.gl_env import apply_gl_fix_and_reexec

    apply_gl_fix_and_reexec(argv)

    repo_root = Path(__file__).resolve().parents[3]
    default_input = repo_root / "input_files" / "simple_parametric_input_new.yaml"

    parser = argparse.ArgumentParser(
        description="Paratan simple-mirror interactive viewer (Trame + PyVista)."
    )
    parser.add_argument(
        "--input",
        "-i",
        type=Path,
        default=default_input,
        help=f"Simple mirror YAML (default: {default_input})",
    )
    parser.add_argument("--host", default="127.0.0.1")
    parser.add_argument("--port", type=int, default=8080)
    parser.add_argument(
        "--n-theta",
        type=int,
        default=96,
        help="Angular tessellation resolution for a full revolution (lower = faster)",
    )
    parser.add_argument(
        "--no-browser",
        action="store_true",
        help="Do not auto-open a browser tab",
    )
    parser.add_argument(
        "--render",
        choices=("client", "server", "trame"),
        default="trame",
        help=(
            "trame=local WebGL while orbiting (default, snappiest); "
            "server=JPEG from Python; client=WebGL only"
        ),
    )
    parser.add_argument(
        "--timeout",
        type=int,
        default=2,
        help=(
            "Seconds to keep the server after the last browser tab disconnects "
            "(default: 2). Use 0 to leave the process running until Ctrl-C."
        ),
    )
    args = parser.parse_args(argv)

    if not args.input.is_file():
        parser.error(f"Input YAML not found: {args.input}")
    if args.timeout < 0:
        parser.error("--timeout must be >= 0")

    from src.paratan.viewer.app import run_app

    print(f"[paratan-viewer] input  : {args.input.resolve()}")
    print(f"[paratan-viewer] open   : http://{args.host}:{args.port}/")
    print(f"[paratan-viewer] render : {args.render}")
    if args.timeout:
        print(f"[paratan-viewer] timeout: {args.timeout}s after last browser disconnect")
    else:
        print("[paratan-viewer] timeout: disabled (Ctrl-C to stop)")
    print("[paratan-viewer] tip    : orbit should feel local; if blank use --render server")
    print("[paratan-viewer] if blank: hard-refresh the tab (Ctrl+Shift+R)")
    run_app(
        args.input,
        host=args.host,
        port=args.port,
        n_theta=args.n_theta,
        open_browser=not args.no_browser,
        render_mode=args.render,
        timeout=args.timeout,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
