"""Ensure external toolchains are available, installing via Mise if needed.

Groups (normal users only need ``build`` — ``docs`` pulls mdbook):

- ``build``: cmake, ninja (default; what configure/build/run need)
- ``docs``:  mdbook (what ``mise run docs`` needs)
- ``all``:   everything (CI / fresh machine setup)
"""

import shutil
import subprocess
import sys

BUILD_TOOLS = ["cmake", "ninja"]
DOCS_TOOLS = ["mdbook"]

GROUPS = {
    "build": BUILD_TOOLS,
    "docs": DOCS_TOOLS,
    "all": BUILD_TOOLS + DOCS_TOOLS,
}


def ensure(tool: str) -> None:
    path = shutil.which(tool)
    if path:
        print(f"{tool} found: {path}")
        return
    print(f"{tool} not found, installing with Mise...")
    _ = subprocess.run(
        ["mise", "use", "--global", f"{tool}@latest"],
        check=True,
    )


def main(argv: list[str]) -> int:
    group = argv[1] if len(argv) > 1 else "build"
    tools = GROUPS.get(group)
    if tools is None:
        print(
            f"unknown bootstrap group: {group!r} "
            f"(expected one of: {', '.join(sorted(GROUPS))})",
            file=sys.stderr,
        )
        return 2
    for tool in tools:
        ensure(tool)
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
