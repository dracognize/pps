"""Xóa CMakeCache.txt + CMakeFiles (giữ _deps) rồi configure lại.

Usage:
    python scripts/reconfigure.py
"""

import shutil
import subprocess
import sys
from pathlib import Path

ROOT  = Path(__file__).resolve().parent.parent
BUILD = ROOT / "build"

def main() -> int:
    # Xóa cache nhưng GIỮ LẠI _deps (Eigen, fmt, ftxui đã download)
    for name in ["CMakeCache.txt", "CMakeFiles"]:
        target = BUILD / name
        if target.is_file():
            target.unlink()
            print(f"Removed file: {target}")
        elif target.is_dir():
            shutil.rmtree(target)
            print(f"Removed dir:  {target}")

    print("\n--- cmake configure ---")
    result = subprocess.run(
        ["cmake", "-S", str(ROOT), "-B", str(BUILD), "-G", "Ninja"],
        check=False,
    )
    return result.returncode

if __name__ == "__main__":
    raise SystemExit(main())
