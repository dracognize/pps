import shutil
import subprocess

TOOLS = ["cmake", "ninja"]


for tool in TOOLS:
    if shutil.which(tool):
        print(f"{tool} found: {shutil.which(tool)}")
    else:
        print(f"{tool} not found, installing with Mise...")
        subprocess.run(
            ["mise", "use", "--global", f"{tool}@latest"],
            check=True,
        )
