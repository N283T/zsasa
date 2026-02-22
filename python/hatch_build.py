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
        # Mark as platform-specific wheel (required for .so/.dylib bundling)
        build_data["infer_tag"] = True

        # Determine paths
        root_dir = Path(self.root).parent  # Go up from python/ to project root
        python_dir = Path(self.root)
        package_dir = python_dir / "zsasa"

        # Platform-specific library name
        if sys.platform == "darwin":
            lib_name = "libzsasa.dylib"
        elif sys.platform == "win32":
            lib_name = "zsasa.dll"
        else:
            lib_name = "libzsasa.so"

        # Zig outputs DLLs to zig-out/bin/ on Windows, .so/.dylib to zig-out/lib/
        if sys.platform == "win32":
            lib_src = root_dir / "zig-out" / "bin" / lib_name
        else:
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
        build_data["force_include"][str(lib_dst)] = f"zsasa/{lib_name}"

    def _find_zig(self) -> list[str]:
        """Find the Zig compiler command."""
        # Check PATH first (system install or setup-zig action)
        if shutil.which("zig"):
            return ["zig"]
        # Check python -m ziglang (PyPI ziglang package)
        try:
            subprocess.run(
                [sys.executable, "-m", "ziglang", "version"],
                capture_output=True,
                check=True,
                timeout=10,
            )
            return [sys.executable, "-m", "ziglang"]
        except (subprocess.CalledProcessError, FileNotFoundError, subprocess.TimeoutExpired):
            return []

    def _build_zig(self, root_dir: Path) -> None:
        """Run zig build command."""
        self.app.display_info("Building Zig library...")

        zig_cmd = self._find_zig()
        if not zig_cmd:
            msg = (
                "Zig compiler not found. Install Zig 0.15.2+ from "
                "https://ziglang.org/download/ or run: pip install ziglang"
            )
            raise RuntimeError(msg)

        try:
            subprocess.run(
                [*zig_cmd, "build", "-Doptimize=ReleaseFast"],
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
