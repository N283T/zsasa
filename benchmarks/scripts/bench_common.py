"""Shared benchmark utilities.

Pure library module imported by bench.py and bench_lr.py.
Contains tool definitions, parsing, system info, run functions, and summary display.
"""

from __future__ import annotations

import csv
import json
import math
import os
import platform
import re
import subprocess
import tempfile
from pathlib import Path

from rich.console import Console
from rich.table import Table

console = Console()

TOOLS = ["zig_f64", "zig_f32", "freesasa", "rust"]
TOOL_ALIASES = {"zig": "zig_f64"}
ALGORITHMS = ["sr", "lr"]


def parse_tool(tool: str) -> tuple[str, str, str]:
    """Parse tool name into (canonical, base, precision).

    Examples:
        "zig_f64" -> ("zig_f64", "zig", "f64")
        "zig_f32" -> ("zig_f32", "zig", "f32")
        "zig"     -> ("zig_f64", "zig", "f64")  # alias
        "freesasa" -> ("freesasa", "freesasa", "f64")
        "rust"     -> ("rust", "rust", "f64")
    """
    tool = TOOL_ALIASES.get(tool, tool)
    if tool.startswith("zig_f"):
        return tool, "zig", tool.split("_")[1]
    return tool, tool, "f64"


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
    """Get binary path for a tool."""
    root = Path(__file__).parent.parent.parent

    if tool == "zig":
        return root.joinpath("zig-out", "bin", "zsasa")
    elif tool == "freesasa":
        return root.joinpath(
            "benchmarks", "external", "freesasa-bench", "src", "freesasa"
        )
    elif tool == "rust":
        return root.joinpath(
            "benchmarks",
            "external",
            "rustsasa-bench",
            "target",
            "release",
            "rust-sasa",
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


def run_zig(
    pdb_path: Path,
    algorithm: str,
    n_threads: int,
    precision: str = "f64",
    n_points: int = 100,
    n_slices: int = 20,
    use_bitmask: bool = False,
) -> tuple[float, float]:
    """Run Zig benchmark. Returns (sasa_time_ms, total_sasa)."""
    binary = get_binary_path("zig")
    if not binary.exists():
        raise FileNotFoundError(f"Zig binary not found: {binary}")

    output_path: Path | None = None

    try:
        with tempfile.NamedTemporaryFile(suffix=".json", delete=False) as f:
            output_path = Path(f.name)

        cmd = [
            str(binary),
            "calc",
            "--timing",
            f"--algorithm={algorithm}",
            f"--threads={n_threads}",
            f"--precision={precision}",
        ]
        if algorithm == "lr":
            cmd.append(f"--n-slices={n_slices}")
        else:
            cmd.append(f"--n-points={n_points}")
        if use_bitmask:
            cmd.append("--use-bitmask")
        cmd.extend([str(pdb_path), str(output_path)])

        result = subprocess.run(cmd, capture_output=True, text=True, timeout=600)

        if result.returncode != 0:
            raise RuntimeError(f"Zig failed: {result.stderr}")

        sasa_time_ms = 0.0
        for line in result.stderr.split("\n"):
            match = re.search(r"SASA calculation:\s*([\d.]+)\s*ms", line)
            if match:
                sasa_time_ms = float(match.group(1))
                break

        with open(output_path) as f:
            data = json.load(f)

        return sasa_time_ms, data["total_area"]

    finally:
        if output_path is not None:
            output_path.unlink(missing_ok=True)


def run_freesasa(
    pdb_path: Path,
    algorithm: str,
    n_threads: int,
    n_points: int = 100,
    n_slices: int = 20,
) -> tuple[float, float]:
    """Run FreeSASA C benchmark. Returns (sasa_time_ms, total_sasa)."""
    binary = get_binary_path("freesasa")
    if not binary.exists():
        raise FileNotFoundError(f"FreeSASA binary not found: {binary}")

    cmd = [
        str(binary),
        str(pdb_path),
        f"--n-threads={n_threads}",
    ]

    if algorithm == "sr":
        cmd.extend(["--shrake-rupley", f"--resolution={n_points}"])
    else:
        cmd.extend(["--lee-richards", f"--resolution={n_slices}"])

    result = subprocess.run(cmd, capture_output=True, text=True, timeout=600)

    if result.returncode != 0:
        raise RuntimeError(f"FreeSASA failed: {result.stderr}")

    total_sasa = 0.0
    for line in result.stdout.split("\n"):
        if line.startswith("Total"):
            match = re.search(r"[\d.]+", line)
            if match:
                total_sasa = float(match.group())
                break

    sasa_time_ms = 0.0
    for line in result.stderr.split("\n"):
        match = re.search(r"SASA calculation time:\s*([\d.]+)\s*ms", line)
        if match:
            sasa_time_ms = float(match.group(1))
            break

    return sasa_time_ms, total_sasa


def run_rust(
    pdb_path: Path,
    n_threads: int,
    n_points: int = 100,
) -> tuple[float, float]:
    """Run RustSASA benchmark. Returns (sasa_time_ms, total_sasa)."""
    binary = get_binary_path("rust")
    if not binary.exists():
        raise FileNotFoundError(f"RustSASA binary not found: {binary}")

    cmd = [
        str(binary),
        str(pdb_path),
        "-n",
        str(n_points),
        "-t",
        str(n_threads),
    ]

    result = subprocess.run(cmd, capture_output=True, text=True, timeout=600)

    if result.returncode != 0:
        raise RuntimeError(f"RustSASA failed: {result.stderr}")

    sasa_time_ms = 0.0
    for line in result.stderr.split("\n"):
        if "SASA_TIME_US:" in line:
            try:
                us = int(line.split(":")[1].strip())
                sasa_time_ms = us / 1000.0
            except (ValueError, IndexError):
                pass
            break

    total_sasa = 0.0
    for line in result.stdout.split("\n"):
        match = re.search(r"Total SASA:\s*([\d.]+)", line)
        if match:
            total_sasa = float(match.group(1))
            break

    return sasa_time_ms, total_sasa


def run_benchmark(
    tool: str,
    pdb_path: Path,
    algorithm: str,
    n_threads: int,
    precision: str = "f64",
    n_points: int = 100,
    n_slices: int = 20,
    use_bitmask: bool = False,
) -> tuple[float, float]:
    """Run benchmark for a specific tool. Returns (sasa_time_ms, total_sasa)."""
    if tool == "zig":
        return run_zig(
            pdb_path, algorithm, n_threads, precision, n_points, n_slices, use_bitmask
        )
    elif tool == "freesasa":
        return run_freesasa(pdb_path, algorithm, n_threads, n_points, n_slices)
    elif tool == "rust":
        if algorithm != "sr":
            raise ValueError("RustSASA only supports SR algorithm")
        return run_rust(pdb_path, n_threads, n_points)
    else:
        raise ValueError(f"Unknown tool: {tool}")


def print_summary(csv_path: Path, warmup: int, runs: int) -> None:
    """Print summary statistics from results CSV."""
    from collections import defaultdict

    by_threads: dict[int, list[float]] = defaultdict(list)

    with open(csv_path) as f:
        reader = csv.DictReader(f)
        for row in reader:
            t = int(row["threads"])
            by_threads[t].append(float(row["sasa_time_ms"]))

    if not by_threads:
        return

    table = Table(title="Summary (SASA time per structure)")
    table.add_column("Threads", style="cyan", justify="right")
    table.add_column("N", justify="right")
    table.add_column("Mean (ms)", justify="right")
    table.add_column("Std (ms)", justify="right")
    table.add_column("Min (ms)", justify="right")
    table.add_column("Max (ms)", justify="right")
    table.add_column("Total (s)", justify="right")

    for t in sorted(by_threads):
        times = by_threads[t]
        n = len(times)
        mean = sum(times) / n
        variance = sum((x - mean) ** 2 for x in times) / (n - 1) if n > 1 else 0.0
        std = math.sqrt(variance)
        total_s = sum(times) / 1000.0
        table.add_row(
            str(t),
            f"{n:,}",
            f"{mean:.3f}",
            f"{std:.3f}",
            f"{min(times):.3f}",
            f"{max(times):.3f}",
            f"{total_s:.1f}",
        )

    console.print()
    console.print(table)
