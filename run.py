import subprocess
import sys
from pathlib import Path

executable = Path("build") / ("pps.exe" if sys.platform == "win32" else "pps")

subprocess.run([executable], check=True)
