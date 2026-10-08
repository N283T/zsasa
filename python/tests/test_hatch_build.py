"""Tests for the wheel build hook (``python/hatch_build.py``)."""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path
from typing import Any

import pytest

pytest.importorskip("hatchling")

HOOK_FILE = Path(__file__).parent.parent / "hatch_build.py"


def _load_hook_module() -> Any:
    spec = importlib.util.spec_from_file_location("zsasa_hatch_build", HOOK_FILE)
    assert spec is not None
    assert spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


hatch_build = _load_hook_module()


class _App:
    """Collects what the hook reports."""

    def __init__(self) -> None:
        self.messages: list[str] = []

    def display_info(self, message: str) -> None:
        self.messages.append(message)

    display_success = display_info
    display_error = display_info


@pytest.fixture
def project(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> Path:
    """A project directory as the hook sees it: ``<root>/zsasa`` is the package.

    There are no Zig sources next to it and no ``zig`` on ``PATH``, so a hook that
    tried to build or to locate the sources would fail.
    """
    (tmp_path / "python" / "zsasa").mkdir(parents=True)
    monkeypatch.setenv("PATH", str(tmp_path / "empty-path"))
    return tmp_path / "python"


def _hook(project: Path) -> tuple[Any, _App]:
    app = _App()
    hook = hatch_build.ZigBuildHook(
        str(project), {}, None, None, str(project / "dist"), "wheel", app=app
    )
    return hook, app


def test_no_bundle_builds_a_pure_python_wheel(
    monkeypatch: pytest.MonkeyPatch, project: Path
) -> None:
    """Nothing is compiled, nothing is added, and the tag is left to hatchling."""
    monkeypatch.setenv("ZSASA_NO_BUNDLE", "1")
    hook, app = _hook(project)
    build_data: dict[str, Any] = {"force_include": {}, "infer_tag": False, "pure_python": True}

    hook.initialize("standard", build_data)

    assert build_data == {"force_include": {}, "infer_tag": False, "pure_python": True}
    assert app.messages == ["ZSASA_NO_BUNDLE=1: not bundling the Zig library and binary"]
    assert sorted(p.name for p in (project / "zsasa").iterdir()) == []


@pytest.mark.parametrize("name", hatch_build.NATIVE_FILE_NAMES)
def test_no_bundle_refuses_native_files_left_in_the_package(
    monkeypatch: pytest.MonkeyPatch, project: Path, name: str
) -> None:
    """A copy from an earlier build would end up in the ``py3-none-any`` wheel."""
    monkeypatch.setenv("ZSASA_NO_BUNDLE", "1")
    (project / "zsasa" / name).write_bytes(b"stale")
    hook, _ = _hook(project)

    with pytest.raises(RuntimeError, match="ZSASA_NO_BUNDLE=1") as excinfo:
        hook.initialize("standard", {"force_include": {}})

    assert name in str(excinfo.value)


@pytest.mark.parametrize("value", [None, "", "0", "true"])
def test_only_the_value_1_turns_bundling_off(
    monkeypatch: pytest.MonkeyPatch, project: Path, value: str | None
) -> None:
    """Without the switch the hook goes on to bundle and marks the wheel as platform-specific."""
    if value is None:
        monkeypatch.delenv("ZSASA_NO_BUNDLE", raising=False)
    else:
        monkeypatch.setenv("ZSASA_NO_BUNDLE", value)
    hook, _ = _hook(project)
    build_data: dict[str, Any] = {"force_include": {}}

    # No Zig sources in this project: the hook stops where it looks for them.
    with pytest.raises(FileNotFoundError, match="Zig source files were not found"):
        hook.initialize("standard", build_data)

    assert build_data["infer_tag"] is True


def test_bundling_copies_a_finished_build_into_the_wheel(
    monkeypatch: pytest.MonkeyPatch, project: Path
) -> None:
    """The default: ``zig-out`` of the checkout is copied and force-included."""
    monkeypatch.delenv("ZSASA_NO_BUNDLE", raising=False)
    root = project.parent
    (root / "build.zig").write_text("")
    (root / "src").mkdir()
    lib_name = {"darwin": "libzsasa.dylib", "win32": "zsasa.dll"}.get(sys.platform, "libzsasa.so")
    exe_name = "zsasa.exe" if sys.platform == "win32" else "zsasa"
    lib_dir = "bin" if sys.platform == "win32" else "lib"
    (root / "zig-out" / lib_dir).mkdir(parents=True, exist_ok=True)
    (root / "zig-out" / "bin").mkdir(parents=True, exist_ok=True)
    (root / "zig-out" / lib_dir / lib_name).write_bytes(b"library")
    (root / "zig-out" / "bin" / exe_name).write_bytes(b"binary")
    hook, _ = _hook(project)
    build_data: dict[str, Any] = {"force_include": {}}

    hook.initialize("standard", build_data)

    package = project / "zsasa"
    assert build_data["infer_tag"] is True
    assert build_data["force_include"] == {
        str(package / lib_name): f"zsasa/{lib_name}",
        str(package / exe_name): f"zsasa/{exe_name}",
    }
    assert (package / lib_name).read_bytes() == b"library"
    assert (package / exe_name).read_bytes() == b"binary"
