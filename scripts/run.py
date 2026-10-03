import os
import subprocess
import sys
from pathlib import Path

# Fallback only: mise.toml [env] already provides PYTHONPATH when run
# via mise. This covers direct `python scripts/run.py` without mise.
build = str(Path("build").resolve())
paths = os.environ.get("PYTHONPATH", "").split(os.pathsep)
if build not in paths:
    paths = [build, *(p for p in paths if p)]
    os.environ["PYTHONPATH"] = os.pathsep.join(paths)

target = sys.argv[1] if len(sys.argv) > 1 else "main.py"

_ = subprocess.run([sys.executable, target], check=True)
