"""CLI entry point that delegates to the bundled zsasa binary."""

from __future__ import annotations

import os
import sys
from pathlib import Path


def _find_binary() -> str:
    """Locate the bundled zsasa binary."""
    package_dir = Path(__file__).parent
    name = "zsasa.exe" if sys.platform == "win32" else "zsasa"
    binary = package_dir / name
    if not binary.exists():
        msg = (
            f"zsasa binary not found at {binary}. "
            "This may indicate a broken installation. "
            "Try reinstalling: pip install --force-reinstall zsasa"
        )
        raise FileNotFoundError(msg)
    return str(binary)


def main() -> None:
    """Run the bundled zsasa CLI binary."""
    binary = _find_binary()
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
