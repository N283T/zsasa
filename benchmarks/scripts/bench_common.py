"""Shared benchmark utilities.

Pure library module imported by bench.py, bench_lr.py, and bench_batch.py.
Contains binary path resolution, parsing helpers, system info, hyperfine runner,
and I/O utilities.
"""

from __future__ import annotations

import csv
import json
import os
import platform
import shlex
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

from rich.console import Console
from rich.table import Table

console = Console()

TOOL_ALIASES = {"zig": "zig_f64"}


def parse_tool(tool: str) -> tuple[str, str, str]:
    """Parse tool name into (canonical, base, precision).

    Examples:
        "zig_f64"  -> ("zig_f64", "zig", "f64")
        "zig_f32"  -> ("zig_f32", "zig", "f32")
        "zig"      -> ("zig_f64", "zig", "f64")  # alias
        "freesasa" -> ("freesasa", "freesasa", "f64")
        "rust"     -> ("rust", "rust", "f64")
        "lahuta"   -> ("lahuta", "lahuta", "f64")
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


def _bin_dir() -> Path:
    """Get benchmarks/external/bin/ directory (populated by setup.sh)."""
    return Path(__file__).parent.parent.joinpath("external", "bin")


# tool name -> binary filename in bin/
_TOOL_BINARIES = {
    "zig": "zsasa",
    "freesasa": "freesasa",
    "rust": "rust-sasa",
    "lahuta": "lahuta",
}


def get_binary_path(tool: str) -> Path:
    """Get binary path for a tool (all symlinked in external/bin/ by setup.sh).

    Tools:
        zig      -> bin/zsasa
        freesasa -> bin/freesasa
        rust     -> bin/rust-sasa
        lahuta   -> bin/lahuta
    """
    if tool not in _TOOL_BINARIES:
        raise ValueError(f"Unknown tool: {tool}")
    return _bin_dir().joinpath(_TOOL_BINARIES[tool])


def get_system_info() -> dict:
    """Get system information for benchmark config."""
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
                timeout=5,
            )
            if result.returncode == 0:
                info["cpu_model"] = result.stdout.strip()
        except (subprocess.SubprocessError, OSError) as e:
            info["cpu_model_error"] = str(e)

        try:
            result = subprocess.run(
                ["sysctl", "-n", "hw.memsize"],
                capture_output=True,
                text=True,
                timeout=5,
            )
            if result.returncode == 0:
                info["memory_gb"] = int(result.stdout.strip()) // (1024**3)
        except (subprocess.SubprocessError, OSError, ValueError) as e:
            info["memory_gb_error"] = str(e)

    return info


def get_n_atoms_from_pdb(pdb_path: Path) -> int:
    """Count ATOM records in a PDB file (excludes HETATM).

    Returns 0 and logs a warning to stderr on read failure.
    """
    count = 0
    try:
        with open(pdb_path) as f:
            for line in f:
                if line.startswith("ATOM  "):
                    count += 1
    except OSError as e:
        print(f"Warning: could not read {pdb_path}: {e}", file=sys.stderr)
    return count


def scan_input_directory(input_dir: Path) -> list[tuple[str, int]]:
    """Scan directory for .pdb files and return (id, 0) pairs.

    Uses os.scandir for fast scanning. n_atoms is set to 0 as a placeholder;
    callers resolve it lazily via get_n_atoms_from_pdb.
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


def quote_path(path: Path) -> str:
    """Shell-quote a path for safe embedding in command strings."""
    return shlex.quote(str(path))


def run_hyperfine(cmd: str, warmup: int, runs: int, json_path: Path) -> dict | None:
    """Run hyperfine and return results dict, or None on failure.

    Times out after 600 seconds. Logs diagnostic details on failure.
    """
    hyperfine_cmd = [
        "hyperfine",
        "--warmup",
        str(warmup),
        "--runs",
        str(runs),
        "--export-json",
        str(json_path),
        cmd,
    ]
    try:
        subprocess.run(
            hyperfine_cmd, check=True, capture_output=True, text=True, timeout=600
        )
    except subprocess.TimeoutExpired:
        console.print(f"[red]  Timeout: command exceeded 600s: {cmd}[/red]")
        return None
    except subprocess.CalledProcessError as e:
        stderr_snippet = (e.stderr or "").strip()[-500:]
        console.print(f"[red]  hyperfine exited with code {e.returncode}[/red]")
        if stderr_snippet:
            console.print(f"[dim]  stderr: {stderr_snippet}[/dim]")
        return None

    try:
        if not json_path.exists():
            console.print(f"[red]  No JSON output at {json_path}[/red]")
            return None
        with open(json_path) as f:
            data = json.load(f)
        if not data.get("results"):
            console.print("[red]  Empty results in hyperfine JSON[/red]")
            return None
        return data["results"][0]
    except (json.JSONDecodeError, KeyError) as e:
        console.print(f"[red]  Failed to parse hyperfine JSON: {e}[/red]")
        return None


def print_hyperfine_summary(csv_path: Path) -> None:
    """Print summary table from hyperfine results CSV (mean_s column)."""
    by_threads: dict[int, list[float]] = defaultdict(list)

    with open(csv_path) as f:
        reader = csv.DictReader(f)
        for row in reader:
            by_threads[int(row["threads"])].append(float(row["mean_s"]))

    if not by_threads:
        return

    table = Table(title="Summary (wall-clock per structure)")
    table.add_column("Threads", style="cyan", justify="right")
    table.add_column("N", justify="right")
    table.add_column("Mean (s)", justify="right")
    table.add_column("Min (s)", justify="right")
    table.add_column("Max (s)", justify="right")

    for t in sorted(by_threads):
        times = by_threads[t]
        n = len(times)
        mean = sum(times) / n
        table.add_row(
            str(t),
            f"{n:,}",
            f"{mean:.4f}",
            f"{min(times):.4f}",
            f"{max(times):.4f}",
        )

    console.print()
    console.print(table)
