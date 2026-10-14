#!/usr/bin/env python3
"""Launch the ParaTAN interactive geometry viewer from the repository root.

Examples
--------
    ./view_simple_mirror.py
    ./view_simple_mirror.py -i input_files/tandem_parametric_input.yaml --port 8081
    ./view_simple_mirror.py --render server --no-browser
"""
from __future__ import annotations

import sys
from pathlib import Path

# Ensure repo-root imports work when invoked as ./view_simple_mirror.py
_ROOT = Path(__file__).resolve().parent
if str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))

# Must run before importing anything that pulls in VTK/PyVista.
from src.paratan.viewer.gl_env import apply_gl_fix_and_reexec

apply_gl_fix_and_reexec()

from src.paratan.viewer.__main__ import main

if __name__ == "__main__":
    raise SystemExit(main())
