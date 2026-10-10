#!/usr/bin/env python3
"""Launch the Paratan simple-mirror interactive viewer."""
from __future__ import annotations

# Must run before importing anything that pulls in VTK/PyVista.
from src.paratan.viewer.gl_env import apply_gl_fix_and_reexec

apply_gl_fix_and_reexec()

from src.paratan.viewer.__main__ import main

if __name__ == "__main__":
    raise SystemExit(main())
