"""Tests for FFI library discovery and the ABI version check."""

from __future__ import annotations

import ctypes.util
import os
import sys
import warnings
from pathlib import Path

import pytest

from zsasa import _ffi


@pytest.fixture
def checkout(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> Path:
    """A fake checkout whose ``python/zsasa/_ffi.py`` is the module under test.

    Returns the checkout root. The library file name is the one of this platform,
    ``ZSASA_LIB`` is unset and the working directory is an empty directory.
    """
    package_dir = tmp_path / "python" / "zsasa"
    package_dir.mkdir(parents=True)
    work_dir = tmp_path / "cwd"
    work_dir.mkdir()
    monkeypatch.delenv("ZSASA_LIB", raising=False)
    monkeypatch.setattr(_ffi, "__file__", str(package_dir / "_ffi.py"))
    monkeypatch.chdir(work_dir)
    return tmp_path


def _write_library(path: Path, content: bytes, mtime: float) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(content)
    os.utime(path, (mtime, mtime))
    return path


def _lib_name() -> str:
    return _ffi._library_file_name()


def test_windows_editable_install_finds_zig_out_bin(
    monkeypatch: pytest.MonkeyPatch,
    checkout: Path,
) -> None:
    """Windows editable installs should discover Zig DLLs under zig-out/bin."""
    dll_path = checkout / "zig-out" / "bin" / "zsasa.dll"
    dll_path.parent.mkdir(parents=True)
    dll_path.write_bytes(b"")

    monkeypatch.setattr(_ffi.sys, "platform", "win32")

    assert _ffi._find_library() == dll_path


def test_zsasa_lib_environment_variable_has_priority(
    monkeypatch: pytest.MonkeyPatch,
    checkout: Path,
) -> None:
    _write_library(checkout / "python" / "zsasa" / _lib_name(), b"bundled", 1000)
    _write_library(checkout / "zig-out" / "lib" / _lib_name(), b"built", 2000)
    chosen = checkout / "custom" / "libzsasa-custom.so"
    monkeypatch.setenv("ZSASA_LIB", str(chosen))

    with warnings.catch_warnings():
        warnings.simplefilter("error")
        assert _ffi._find_library() == chosen


def test_newer_zig_out_build_beats_stale_bundled_library(checkout: Path) -> None:
    """The failure of the report: an old bundled copy shadowed a fresh build."""
    bundled = _write_library(checkout / "python" / "zsasa" / _lib_name(), b"old", 1000)
    built = _write_library(checkout / "zig-out" / "lib" / _lib_name(), b"new", 2000)

    with pytest.warns(UserWarning) as record:
        assert _ffi._find_library() == built

    message = str(record[0].message)
    assert f"Loading {built}" in message
    assert f"skipping {bundled}" in message
    assert "newer" in message


def test_newer_bundled_library_beats_older_zig_out_build(checkout: Path) -> None:
    bundled = _write_library(checkout / "python" / "zsasa" / _lib_name(), b"new", 2000)
    built = _write_library(checkout / "zig-out" / "lib" / _lib_name(), b"old", 1000)

    with pytest.warns(UserWarning, match="skipping") as record:
        assert _ffi._find_library() == bundled

    assert f"skipping {built}" in str(record[0].message)


def test_identical_bundled_copy_of_the_build_is_not_a_conflict(checkout: Path) -> None:
    """The build hook copies zig-out to the package directory: no warning for that."""
    bundled = _write_library(checkout / "python" / "zsasa" / _lib_name(), b"same", 1000)
    _write_library(checkout / "zig-out" / "lib" / _lib_name(), b"same", 2000)

    with warnings.catch_warnings():
        warnings.simplefilter("error")
        assert _ffi._find_library() == bundled


def test_single_library_is_found_without_warning(checkout: Path) -> None:
    built = _write_library(checkout / "zig-out" / "lib" / _lib_name(), b"built", 1000)

    with warnings.catch_warnings():
        warnings.simplefilter("error")
        assert _ffi._find_library() == built

    built.unlink()
    bundled = _write_library(checkout / "python" / "zsasa" / _lib_name(), b"bundled", 1000)
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        assert _ffi._find_library() == bundled


def test_current_directory_is_not_searched(
    monkeypatch: pytest.MonkeyPatch,
    checkout: Path,
) -> None:
    """A library in the working directory must not be loaded."""
    work_dir = checkout / "cwd"
    _write_library(work_dir / _lib_name(), b"stray", 1000)
    _write_library(work_dir / "zig-out" / "lib" / _lib_name(), b"stray", 1000)
    _write_library(work_dir / "zig-out" / "bin" / _lib_name(), b"stray", 1000)

    # Keep a library in the real system directories from making this test depend on the host.
    real_exists = Path.exists
    monkeypatch.setattr(
        Path,
        "exists",
        lambda self: (
            False if str(self).startswith(("/usr/lib", "/usr/local/lib")) else real_exists(self)
        ),
    )

    with pytest.raises(FileNotFoundError, match="Could not find"):
        _ffi._find_library()


def test_loaded_library_has_the_expected_abi_version() -> None:
    _, lib = _ffi._get_lib()
    assert lib.zsasa_abi_version() == _ffi._EXPECTED_ABI_VERSION


class _FakeLib:
    """A library object whose zsasa_abi_version behaves as configured."""

    def __init__(self, version: int | None) -> None:
        self._version = version

    def __getattr__(self, name: str):  # noqa: ANN204
        if name == "zsasa_abi_version" and self._version is not None:
            return lambda: self._version
        raise AttributeError(name)


def test_matching_abi_version_is_accepted() -> None:
    _ffi._check_abi_version(_FakeLib(_ffi._EXPECTED_ABI_VERSION), Path("/lib/libzsasa.so"))


def test_different_abi_version_is_rejected_with_the_library_path() -> None:
    lib = _FakeLib(_ffi._EXPECTED_ABI_VERSION + 1)

    with pytest.raises(ImportError) as excinfo:
        _ffi._check_abi_version(lib, Path("/some/where/libzsasa.so"))

    message = str(excinfo.value)
    assert "/some/where/libzsasa.so" in message
    assert f"ABI version {_ffi._EXPECTED_ABI_VERSION + 1}" in message
    assert f"expects ABI version {_ffi._EXPECTED_ABI_VERSION}" in message


def test_missing_abi_version_symbol_is_rejected_with_the_library_path() -> None:
    with pytest.raises(ImportError, match="does not export zsasa_abi_version") as excinfo:
        _ffi._check_abi_version(_FakeLib(None), Path("/some/where/libzsasa.so"))

    assert "/some/where/libzsasa.so" in str(excinfo.value)


def test_loading_a_library_without_the_symbol_fails_with_import_error(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """A real shared library that is not zsasa (so it has no zsasa_abi_version)."""
    if sys.platform == "win32":
        pytest.skip("no system library with a known path")
    other = ctypes.util.find_library("c") or ctypes.util.find_library("m")
    if other is None:
        pytest.skip("no system C library found")
    # find_library can return a bare name; dlopen resolves it the same way.
    monkeypatch.setenv("ZSASA_LIB", other)

    with pytest.raises(ImportError, match="zsasa_abi_version") as excinfo:
        _ffi._load_library()

    assert other in str(excinfo.value)
