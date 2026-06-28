"""Tests for FFI library discovery."""

from __future__ import annotations

from pathlib import Path

from zsasa import _ffi


def test_windows_editable_install_finds_zig_out_bin(
    monkeypatch,
    tmp_path: Path,
) -> None:
    """Windows editable installs should discover Zig DLLs under zig-out/bin."""
    dll_path = tmp_path.joinpath("zig-out", "bin", "zsasa.dll")
    dll_path.parent.mkdir(parents=True)
    dll_path.write_bytes(b"")

    monkeypatch.delenv("ZSASA_LIB", raising=False)
    monkeypatch.setattr(_ffi.sys, "platform", "win32")
    monkeypatch.chdir(tmp_path)

    assert _ffi._find_library() == dll_path
