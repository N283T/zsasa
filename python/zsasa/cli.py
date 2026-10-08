"""CLI entry point that delegates to the native zsasa binary."""

from __future__ import annotations

import os
import sys
from pathlib import Path

from zsasa._ffi import _environment_directories


def _is_script(path: Path) -> bool:
    """Return whether a file is a script (starts with ``#!``), not a native executable.

    A file that cannot be read counts as a script: it cannot be run either.
    """
    try:
        with path.open("rb") as f:
            return f.read(2) == b"#!"
    except OSError:
        return True


def _find_binary() -> str:
    """Locate the native zsasa binary.

    Order: the binary bundled in the package (a wheel), then the environment of the
    running interpreter, for a pure-Python installation whose binary a package manager
    ships on its own (see ``_environment_directories``).

    The ``zsasa`` console script of this package is installed as ``<prefix>/bin/zsasa``,
    the very path the native binary has in such an environment. Running it from here
    would start this function again, without end, so a candidate that is a script is
    never returned. On Windows only ``<prefix>/Library/bin`` is searched: console script
    launchers are executables too, but they are written to ``<prefix>/Scripts``.
    """
    package_dir = Path(__file__).parent
    name = "zsasa.exe" if sys.platform == "win32" else "zsasa"
    bundled = package_dir / name
    if bundled.exists():
        return str(bundled)

    directories = _environment_directories("bin")
    for directory in directories:
        candidate = directory / name
        if candidate.is_file() and not _is_script(candidate):
            return str(candidate)

    msg = (
        f"zsasa binary not found: it is not bundled in the zsasa package ({bundled}) "
        f"and there is no native {name} in {', '.join(str(d) for d in directories)}. "
        "Reinstall with 'pip install --force-reinstall zsasa' (a wheel bundles the binary), "
        "or install the zsasa command-line program into this environment."
    )
    raise FileNotFoundError(msg)


def main() -> None:
    """Run the native zsasa CLI binary."""
    try:
        binary = _find_binary()
    except FileNotFoundError as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)
    _exec_msg = (
        f"Error: failed to execute zsasa binary at {binary}: {{err}}\n"
        "The binary may be corrupted or built for a different platform.\n"
        "Try reinstalling: pip install --force-reinstall zsasa"
    )
    if sys.platform == "win32":
        # Windows doesn't support execvp reliably, use subprocess
        import subprocess

        try:
            sys.exit(subprocess.call([binary, *sys.argv[1:]]))
        except KeyboardInterrupt:
            sys.exit(130)
        except OSError as e:
            print(_exec_msg.format(err=e), file=sys.stderr)
            sys.exit(1)
    else:
        try:
            os.execvp(binary, [binary, *sys.argv[1:]])
        except OSError as e:
            print(_exec_msg.format(err=e), file=sys.stderr)
            sys.exit(1)


if __name__ == "__main__":
    main()
