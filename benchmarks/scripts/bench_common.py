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
from collections import defaultdict
from pathlib import Path

from rich.console import Console
from rich.table import Table

console = Console()

TOOL_ALIASES = {"zig": "zig_f64", "zig_bitmask": "zig_f64_bitmask"}

LAHUTA_BITMASK_POINTS = {64, 128, 256}

BITMASK_CAPABLE_BASES = {"zig", "lahuta"}


def parse_tool(tool: str) -> tuple[str, str, str, bool]:
    """Parse tool name into (canonical, base, precision, use_bitmask).

    Raises:
        ValueError: If a _bitmask suffix is used with a tool that does not
            support bitmask mode.

    Examples:
        "zig_f64"          -> ("zig_f64", "zig", "f64", False)
        "zig_f32"          -> ("zig_f32", "zig", "f32", False)
        "zig"              -> ("zig_f64", "zig", "f64", False)  # alias
        "zig_f64_bitmask"  -> ("zig_f64_bitmask", "zig", "f64", True)
        "zig_f32_bitmask"  -> ("zig_f32_bitmask", "zig", "f32", True)
        "zig_bitmask"      -> ("zig_f64_bitmask", "zig", "f64", True)  # alias
        "freesasa"         -> ("freesasa", "freesasa", "f64", False)
        "rust"             -> ("rust", "rust", "f64", False)
        "lahuta"           -> ("lahuta", "lahuta", "f64", False)
        "lahuta_bitmask"   -> ("lahuta_bitmask", "lahuta", "f64", True)
    """
    tool = TOOL_ALIASES.get(tool, tool)

    # Strip _bitmask suffix
    use_bitmask = tool.endswith("_bitmask")
    base_tool = tool.removesuffix("_bitmask") if use_bitmask else tool

    if base_tool.startswith("zig_f"):
        effective_base = "zig"
        precision = base_tool.split("_")[1]
    else:
        effective_base = base_tool
        precision = "f64"

    if use_bitmask and effective_base not in BITMASK_CAPABLE_BASES:
        raise ValueError(
            f"Tool '{tool}' does not support bitmask mode. "
            f"Bitmask is only supported for: {sorted(BITMASK_CAPABLE_BASES)}"
        )

    return tool, effective_base, precision, use_bitmask


def parse_threads(threads_str: str) -> list[int]:
    """Parse thread specification like '1-10' or '1,4,8'.

    Raises:
        ValueError: If the specification contains invalid or non-positive values.
    """
    result = []
    for part in threads_str.split(","):
        part = part.strip()
        if not part:
            continue
        try:
            if "-" in part and not part.startswith("-"):
                start, end = part.split("-", 1)
                result.extend(range(int(start), int(end) + 1))
            else:
                result.append(int(part))
        except ValueError:
            raise ValueError(
                f"Invalid thread specification '{part}' in '{threads_str}'. "
                f"Expected format: '1,4,8' or '1-10'"
            ) from None
    if not result:
        raise ValueError(f"No valid thread counts in '{threads_str}'")
    if any(t <= 0 for t in result):
        raise ValueError(f"Thread counts must be positive, got: {sorted(set(result))}")
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


def get_n_atoms_from_pdb(pdb_path: Path) -> int | None:
    """Count ATOM records in a PDB file (excludes HETATM).

    Returns None and logs a warning on read failure.
    """
    count = 0
    try:
        with open(pdb_path) as f:
            for line in f:
                if line.startswith("ATOM  "):
                    count += 1
    except OSError as e:
        console.print(f"[yellow]Warning: could not read {pdb_path}: {e}[/yellow]")
        return None
    return count


def scan_input_directory(input_dir: Path) -> list[tuple[str, int]]:
    """Scan directory for .pdb files and return (id, 0) pairs.

    Uses os.scandir for fast scanning. n_atoms is set to 0 as a placeholder;
    callers resolve it lazily via get_n_atoms_from_pdb.
    """
    entries = []
    try:
        with os.scandir(input_dir) as it:
            for entry in it:
                if entry.is_file() and entry.name.endswith(".pdb"):
                    entries.append(entry.name)
    except OSError as e:
        raise RuntimeError(f"Cannot scan directory {input_dir}: {e}") from e

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
    for bin_name, entries in samples.items():
        for i, entry in enumerate(entries):
            if "id" not in entry:
                raise ValueError(
                    f"Sample entry {i} in bin '{bin_name}' missing 'id' key "
                    f"in {sample_path}"
                )
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
        result = data["results"][0]
        required_keys = {"mean", "stddev", "min", "max", "median"}
        missing = required_keys - result.keys()
        if missing:
            console.print(f"[red]  Hyperfine result missing keys: {missing}[/red]")
            return None
        return result
    except (json.JSONDecodeError, KeyError) as e:
        console.print(f"[red]  Failed to parse hyperfine JSON: {e}[/red]")
        return None


def print_hyperfine_summary(csv_path: Path) -> None:
    """Print summary table from hyperfine results CSV.

    Groups rows by thread count and summarizes the per-structure mean_s
    values (mean, min, max) for each thread count.
    """
    by_threads: dict[int, list[float]] = defaultdict(list)

    try:
        with open(csv_path) as f:
            reader = csv.DictReader(f)
            for row in reader:
                by_threads[int(row["threads"])].append(float(row["mean_s"]))
    except (OSError, KeyError, ValueError) as e:
        console.print(f"[yellow]Warning: could not parse summary CSV: {e}[/yellow]")
        return

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
