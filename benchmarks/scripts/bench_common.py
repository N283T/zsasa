"""Shared benchmark utilities.

Pure library module imported by bench.py and bench_lr.py.
Contains binary path resolution, parsing helpers, system info, and I/O utilities.
"""

from __future__ import annotations

import json
import os
import platform
import subprocess
from pathlib import Path


def parse_threads(threads_str: str) -> list[int]:
    """Parse thread specification like '1-10' or '1,4,8'."""
    result = []
    for part in threads_str.split(","):
        if "-" in part:
            start, end = part.split("-", 1)
            result.extend(range(int(start), int(end) + 1))
        else:
            result.append(int(part))
    return sorted(set(result))


def get_binary_path(tool: str) -> Path:
    """Get binary path for a tool.

    Tools:
        zig      -> zig-out/bin/zsasa
        freesasa -> benchmarks/external/freesasa/src/freesasa
        rust     -> benchmarks/external/rustsasa/target/release/rust-sasa
        lahuta   -> benchmarks/external/lahuta/build/cli/lahuta
    """
    root = Path(__file__).parent.parent.parent

    if tool == "zig":
        return root.joinpath("zig-out", "bin", "zsasa")
    elif tool == "freesasa":
        return root.joinpath("benchmarks", "external", "freesasa", "src", "freesasa")
    elif tool == "rust":
        return root.joinpath(
            "benchmarks", "external", "rustsasa", "target", "release", "rust-sasa"
        )
    elif tool == "lahuta":
        return root.joinpath(
            "benchmarks", "external", "lahuta", "build", "cli", "lahuta"
        )
    else:
        raise ValueError(f"Unknown tool: {tool}")


def get_system_info() -> dict:
    """Get system information."""
    info = {
        "os": platform.system(),
        "os_version": platform.release(),
        "arch": platform.machine(),
        "cpu_cores": os.cpu_count() or 1,
    }

    if platform.system() == "Darwin":
        try:
            result = subprocess.run(
                ["sysctl", "-n", "machdep.cpu.brand_string"],
                capture_output=True,
                text=True,
            )
            if result.returncode == 0:
                info["cpu_model"] = result.stdout.strip()
            result = subprocess.run(
                ["sysctl", "-n", "hw.memsize"], capture_output=True, text=True
            )
            if result.returncode == 0:
                info["memory_gb"] = int(result.stdout.strip()) // (1024**3)
        except Exception:
            pass

    return info


def get_n_atoms_from_pdb(pdb_path: Path) -> int:
    """Count ATOM records in a PDB file."""
    count = 0
    try:
        with open(pdb_path) as f:
            for line in f:
                if line.startswith("ATOM  "):
                    count += 1
    except Exception:
        pass
    return count


def scan_input_directory(input_dir: Path) -> list[tuple[str, int]]:
    """Scan directory for .pdb files and return (id, n_atoms=0) list.

    Uses os.scandir for fast scanning. n_atoms is resolved lazily during run.
    """
    entries = []
    with os.scandir(input_dir) as it:
        for entry in it:
            if entry.is_file() and entry.name.endswith(".pdb"):
                entries.append(entry.name)

    entries.sort()

    structures = []
    for filename in entries:
        pdb_id = filename[:-4]  # strip .pdb
        structures.append((pdb_id, 0))

    return structures


def load_sample_file(sample_path: Path) -> list[str]:
    """Load sample file and return list of IDs.

    Supports both formats:
    - v1: {"samples": ["id1", "id2", ...]}
    - v2: {"samples": {"bin": [{"id": "...", ...}, ...], ...}}
    """
    with open(sample_path) as f:
        data = json.load(f)

    if "samples" not in data:
        raise ValueError(f"Invalid sample file: missing 'samples' key in {sample_path}")

    samples = data["samples"]

    # v1: flat list
    if isinstance(samples, list):
        return samples

    # v2: dict of bins
    ids = []
    for entries in samples.values():
        for entry in entries:
            ids.append(entry["id"])
    return ids
