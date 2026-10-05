"""Run the Python or C++ entry point cross-platform.

Usage (same style as ``scripts/bootstrap.py``)::

    python scripts/run.py py [target]  # default target: main.py
    python scripts/run.py cpp [args...]

- ``py``: run a Python file against the compiled ``pps`` module.
- ``cpp``: run the compiled C++ demo (``build/pps``, ``build/pps.exe``
  on Windows). Extra args are forwarded to the executable.
"""

import os
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent


def run_py(target: str) -> int:
    # Fallback only: mise.toml [env] already provides PYTHONPATH when run
    # via mise. This covers direct `python scripts/run.py` without mise.
    build = str(ROOT / "build")
    paths = os.environ.get("PYTHONPATH", "").split(os.pathsep)
    if build not in paths:
        paths = [build, *(p for p in paths if p)]
        os.environ["PYTHONPATH"] = os.pathsep.join(paths)

    completed = subprocess.run([sys.executable, target])
    return completed.returncode


def run_cpp(args: list[str]) -> int:
    exe = ROOT / "build" / ("pps.exe" if os.name == "nt" else "pps")
    if not exe.is_file():
        print(
            f"error: C++ executable not found: {exe}\n"
            "  build it first with: mise run cpp::build",
            file=sys.stderr,
        )
        return 1
    completed = subprocess.run([str(exe), *args])
    return completed.returncode


def main(argv: list[str]) -> int:
    mode = argv[1] if len(argv) > 1 else "py"
    if mode == "cpp":
        return run_cpp(argv[2:])
    rest = argv[1:]
    if rest and rest[0] == "py":
        rest = rest[1:]
    target = rest[0] if rest else "main.py"
    return run_py(target)


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
