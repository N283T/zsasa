"""Custom build hook to compile Zig library during pip install."""

from __future__ import annotations

import shutil
import subprocess
import sys
from pathlib import Path
from typing import Any

from hatchling.builders.hooks.plugin.interface import BuildHookInterface


class ZigBuildHook(BuildHookInterface):
    """Build hook that compiles the Zig library before packaging."""

    PLUGIN_NAME = "zig-build"

    def initialize(self, version: str, build_data: dict[str, Any]) -> None:
        """Build the Zig library and copy it to the package directory."""
        # Determine paths
        root_dir = Path(self.root).parent  # Go up from python/ to project root
        python_dir = Path(self.root)
        package_dir = python_dir / "freesasa_zig"

        # Platform-specific library name
        if sys.platform == "darwin":
            lib_name = "libfreesasa_zig.dylib"
        elif sys.platform == "win32":
            lib_name = "freesasa_zig.dll"
        else:
            lib_name = "libfreesasa_zig.so"

        lib_src = root_dir / "zig-out" / "lib" / lib_name
        lib_dst = package_dir / lib_name

        # Check if library already exists and is newer than source
        if lib_dst.exists():
            # For development, always rebuild if zig-out doesn't exist
            if not lib_src.exists():
                self._build_zig(root_dir)
                shutil.copy2(lib_src, lib_dst)
            elif lib_src.stat().st_mtime > lib_dst.stat().st_mtime:
                self._build_zig(root_dir)
                shutil.copy2(lib_src, lib_dst)
        else:
            # Build if not exists
            if not lib_src.exists():
                self._build_zig(root_dir)
            shutil.copy2(lib_src, lib_dst)

        # Include the library in the wheel
        build_data["force_include"][str(lib_dst)] = f"freesasa_zig/{lib_name}"

    def _build_zig(self, root_dir: Path) -> None:
        """Run zig build command."""
        self.app.display_info("Building Zig library...")

        # Check if zig is available
        if shutil.which("zig") is None:
            msg = (
                "Zig compiler not found. Please install Zig 0.15.2+ from https://ziglang.org/download/ "
                "or set FREESASA_ZIG_LIB to point to a pre-built library."
            )
            raise RuntimeError(msg)

        try:
            subprocess.run(
                ["zig", "build", "-Doptimize=ReleaseFast"],
                cwd=root_dir,
                check=True,
                capture_output=True,
                text=True,
                timeout=600,  # 10 minute timeout
            )
            self.app.display_success("Zig library built successfully")
        except subprocess.CalledProcessError as e:
            self.app.display_error(f"Zig build failed:\n{e.stderr}")
            raise RuntimeError("Failed to build Zig library") from e
