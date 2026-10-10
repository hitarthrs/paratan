"""System Mesa / libstdc++ bootstrap for conda VTK.

Must run before importing pyvista/VTK.
"""
from __future__ import annotations

import os
import sys
from pathlib import Path


def apply_gl_fix_and_reexec(argv: list[str] | None = None) -> None:
    """Re-exec with system libstdc++ + DRI path so VTK can open a GL context."""
    if os.environ.get("_PARATAN_VIEWER_GL_FIXED") == "1":
        return

    sys_lib = Path("/usr/lib/x86_64-linux-gnu/libstdc++.so.6")
    dri = Path("/usr/lib/x86_64-linux-gnu/dri")
    env = os.environ.copy()
    env["_PARATAN_VIEWER_GL_FIXED"] = "1"

    if dri.is_dir():
        env["LIBGL_DRIVERS_PATH"] = str(dri)
        # Prefer software GL when hardware DRI is flaky under conda.
        env.setdefault("LIBGL_ALWAYS_SOFTWARE", "1")
        env.setdefault("GALLIUM_DRIVER", "llvmpipe")

    if sys_lib.is_file():
        prev = env.get("LD_PRELOAD", "")
        parts = [str(sys_lib)] + [
            p for p in prev.split(":") if p and p != str(sys_lib)
        ]
        env["LD_PRELOAD"] = ":".join(parts)

    # Off-screen by default: Trame client mode still builds a VTK scene server-side.
    env.setdefault("PYVISTA_OFF_SCREEN", "true")

    # Preserve -m launch semantics; sys.argv[0] alone becomes a file launch
    # and loses the repository import path after re-exec.
    launch_args = sys.orig_argv[1:] if argv is None else ["-m", "src.paratan.viewer", *argv]
    os.execve(sys.executable, [sys.executable, *launch_args], env)
