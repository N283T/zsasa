#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.12"
# dependencies = ["rich>=13.0", "typer>=0.9.0", "numpy>=1.26", "gemmi>=0.7", "cffi>=1.16"]
# ///
"""Benchmark: Python bindings vs Zig CLI overhead comparison.

Demonstrates that Python bindings have minimal overhead compared to
direct Zig CLI usage, since both use the same Zig SASA engine.

Compares four metrics:
  1. CLI SASA-only   : Internal SASA calculation time (from --timing)
  2. Python SASA-only: calculate_sasa() with pre-loaded data
  3. CLI total       : Full wall-clock (file I/O + parse + SASA + output)
  4. Python total    : gemmi load + classify + calculate_sasa()

Structures are downloaded from AlphaFold DB on first run.

Usage:
    ./benchmarks/scripts/bench_python.py
    ./benchmarks/scripts/bench_python.py --runs 20 --threads 1
    ./benchmarks/scripts/bench_python.py --threads 0 --output benchmarks/results/python

Output:
    benchmarks/results/python/
    ├── config.json
    └── results.json
"""

from __future__ import annotations

import json
import os
import platform
import re
import subprocess
import sys
import tempfile
import time
import urllib.request
from datetime import datetime
from pathlib import Path
from typing import Annotated

import numpy as np
import typer
from rich.console import Console
from rich.table import Table

app = typer.Typer(help="Python bindings vs Zig CLI benchmark")
console = Console()

# AlphaFold structures for benchmarking (various sizes)
STRUCTURES = [
    ("AF-P68871-F1", "P68871", "Hemoglobin β"),  # ~147 aa, ~1.1k atoms
    ("AF-P04637-F1", "P04637", "p53"),  # ~393 aa, ~3k atoms
    ("AF-P02768-F1", "P02768", "Albumin (BSA)"),  # ~609 aa, ~4.8k atoms
    ("AF-P00533-F1", "P00533", "EGFR"),  # ~1210 aa, ~9.4k atoms
]

ALPHAFOLD_VERSION = 6


def get_root_dir() -> Path:
    """Get project root directory."""
    return Path(__file__).resolve().parent.parent.parent


def get_pdb_cache_dir() -> Path:
    """Get cache directory for AlphaFold PDB files."""
    d = get_root_dir().joinpath("benchmarks", "alphafold_pdb")
    d.mkdir(parents=True, exist_ok=True)
    return d


def download_structure(uniprot_id: str, model_id: str) -> Path:
    """Download AlphaFold structure if not cached."""
    cache_dir = get_pdb_cache_dir()
    pdb_path = cache_dir.joinpath(f"{model_id}.pdb")

    if pdb_path.exists():
        return pdb_path

    url = f"https://alphafold.ebi.ac.uk/files/{model_id}-model_v{ALPHAFOLD_VERSION}.pdb"
    console.print(f"  Downloading {model_id} from AlphaFold DB...")
    urllib.request.urlretrieve(url, pdb_path)  # noqa: S310

    # Verify it's a valid PDB
    with open(pdb_path) as f:
        first_line = f.readline()
    if not first_line.startswith(("HEADER", "REMARK", "TITLE")):
        pdb_path.unlink()
        msg = f"Downloaded file is not a valid PDB: {url}"
        raise RuntimeError(msg)

    return pdb_path


def count_atoms_pdb(pdb_path: Path) -> int:
    """Count ATOM records in PDB file (excluding H)."""
    count = 0
    with open(pdb_path) as f:
        for line in f:
            if line.startswith("ATOM"):
                element = line[76:78].strip()
                if element not in ("H", "D"):
                    count += 1
    return count


def bench_cli_timing(
    pdb_path: Path,
    n_runs: int,
    threads: int,
    n_points: int = 100,
    use_bitmask: bool = False,
) -> dict:
    """Benchmark Zig CLI, extracting SASA-only time from --timing."""
    zsasa_bin = get_root_dir().joinpath("zig-out", "bin", "zsasa")
    if not zsasa_bin.exists():
        console.print(
            "[red]zsasa binary not found. Run: zig build -Doptimize=ReleaseFast[/red]"
        )
        raise typer.Exit(1)

    thread_arg = f"--threads={threads}" if threads > 0 else "--threads=0"

    sasa_times: list[float] = []
    total_times: list[float] = []

    for _ in range(n_runs):
        with tempfile.NamedTemporaryFile(suffix=".json") as tmp:
            cmd = [
                str(zsasa_bin),
                "calc",
                "--timing",
                "-q",
                thread_arg,
                f"--n-points={n_points}",
                str(pdb_path),
                tmp.name,
            ]
            if use_bitmask:
                cmd.insert(-2, "--use-bitmask")
            start = time.perf_counter()
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
            )
            wall = time.perf_counter() - start

            if result.returncode != 0:
                console.print(f"[red]CLI error: {result.stderr}[/red]")
                continue

            total_times.append(wall)

            # Parse SASA calculation time from --timing output
            for line in result.stderr.splitlines():
                m = re.search(r"SASA calculation:\s+([\d.]+)\s+ms", line)
                if m:
                    sasa_times.append(float(m.group(1)) / 1000)  # ms → s
                    break

    return {
        "sasa_times": sasa_times,
        "total_times": total_times,
    }


def bench_python_sasa(
    coords: np.ndarray,
    radii: np.ndarray,
    n_runs: int,
    threads: int,
    n_points: int = 100,
    use_bitmask: bool = False,
) -> list[float]:
    """Benchmark Python calculate_sasa() with pre-loaded data."""
    # Import here to avoid top-level dependency on installed zsasa
    sys.path.insert(0, str(get_root_dir().joinpath("python")))
    from zsasa import calculate_sasa

    # Warmup
    calculate_sasa(
        coords, radii, n_threads=threads, n_points=n_points, use_bitmask=use_bitmask
    )

    times = []
    for _ in range(n_runs):
        start = time.perf_counter()
        calculate_sasa(
            coords, radii, n_threads=threads, n_points=n_points, use_bitmask=use_bitmask
        )
        elapsed = time.perf_counter() - start
        times.append(elapsed)
    return times


def bench_python_e2e(
    pdb_path: Path,
    n_runs: int,
    threads: int,
    n_points: int = 100,
) -> list[float]:
    """Benchmark Python end-to-end: gemmi load + classify + SASA."""
    sys.path.insert(0, str(get_root_dir().joinpath("python")))
    from zsasa.core import ClassifierType
    from zsasa.integrations.gemmi import calculate_sasa_from_structure

    # Warmup
    calculate_sasa_from_structure(
        str(pdb_path),
        classifier=ClassifierType.PROTOR,
        n_threads=threads,
        n_points=n_points,
    )

    times = []
    for _ in range(n_runs):
        start = time.perf_counter()
        calculate_sasa_from_structure(
            str(pdb_path),
            classifier=ClassifierType.PROTOR,
            n_threads=threads,
            n_points=n_points,
        )
        elapsed = time.perf_counter() - start
        times.append(elapsed)
    return times


def median(values: list[float]) -> float:
    """Calculate median of a list."""
    s = sorted(values)
    n = len(s)
    if n == 0:
        return 0.0
    return s[n // 2]


def prepare_python_data(pdb_path: Path, threads: int) -> tuple[np.ndarray, np.ndarray]:
    """Load structure with gemmi and classify atoms to get coords + radii."""
    sys.path.insert(0, str(get_root_dir().joinpath("python")))
    from zsasa.core import ClassifierType, classify_atoms
    from zsasa.integrations.gemmi import extract_atoms_from_model

    import gemmi

    structure = gemmi.read_structure(str(pdb_path))
    atom_data = extract_atoms_from_model(structure[0], include_hydrogens=False)
    classification = classify_atoms(
        atom_data.residue_names,
        atom_data.atom_names,
        ClassifierType.PROTOR,
    )
    return atom_data.coords, classification.radii


@app.callback(invoke_without_command=True)
def main(
    runs: Annotated[int, typer.Option(help="Number of benchmark runs")] = 10,
    warmup: Annotated[int, typer.Option(help="Number of warmup runs (CLI only)")] = 2,
    threads: Annotated[int, typer.Option(help="Thread count (0 = auto)")] = 0,
    output: Annotated[
        str | None, typer.Option(help="Output directory for results")
    ] = None,
    n_points: Annotated[
        int,
        typer.Option(
            "--n-points",
            "-N",
            help="Number of sphere test points per atom (default: 100)",
        ),
    ] = 100,
    use_bitmask: Annotated[
        bool,
        typer.Option(
            "--use-bitmask",
            help="Use bitmask neighbor list for zsasa",
        ),
    ] = False,
) -> None:
    """Compare Python bindings vs Zig CLI performance."""
    console.print("[bold]Python Bindings vs Zig CLI Benchmark[/bold]\n")

    # Download structures
    console.print("[bold]Downloading structures...[/bold]")
    pdb_files: list[tuple[str, str, Path, int]] = []
    for model_id, uniprot_id, name in STRUCTURES:
        pdb_path = download_structure(uniprot_id, model_id)
        n_atoms = count_atoms_pdb(pdb_path)
        pdb_files.append((model_id, name, pdb_path, n_atoms))
        console.print(f"  {name}: {n_atoms:,} atoms")

    console.print()

    # CLI warmup
    if warmup > 0:
        console.print(f"[dim]Warming up CLI ({warmup} runs)...[/dim]")
        zsasa_bin = get_root_dir().joinpath("zig-out", "bin", "zsasa")
        thread_arg = f"--threads={threads}" if threads > 0 else "--threads=0"
        for _ in range(warmup):
            for _, _, pdb_path, _ in pdb_files:
                with tempfile.NamedTemporaryFile(suffix=".json") as tmp:
                    warmup_cmd = [
                        str(zsasa_bin),
                        "calc",
                        "-q",
                        thread_arg,
                        f"--n-points={n_points}",
                        str(pdb_path),
                        tmp.name,
                    ]
                    if use_bitmask:
                        warmup_cmd.insert(-2, "--use-bitmask")
                    subprocess.run(
                        warmup_cmd,
                        capture_output=True,
                    )

    # Run benchmarks
    results = []

    for model_id, name, pdb_path, n_atoms in pdb_files:
        console.print(f"[bold]Benchmarking {name} ({n_atoms:,} atoms)...[/bold]")

        # CLI benchmark
        cli = bench_cli_timing(pdb_path, runs, threads, n_points, use_bitmask)

        # Python SASA-only
        coords, radii = prepare_python_data(pdb_path, threads)
        py_sasa_times = bench_python_sasa(
            coords, radii, runs, threads, n_points, use_bitmask
        )

        # Python end-to-end
        py_e2e_times = bench_python_e2e(pdb_path, runs, threads, n_points)

        results.append(
            {
                "model_id": model_id,
                "name": name,
                "n_atoms": n_atoms,
                "cli_sasa_ms": median(cli["sasa_times"]) * 1000,
                "cli_total_ms": median(cli["total_times"]) * 1000,
                "python_sasa_ms": median(py_sasa_times) * 1000,
                "python_e2e_ms": median(py_e2e_times) * 1000,
                "raw": {
                    "cli_sasa": cli["sasa_times"],
                    "cli_total": cli["total_times"],
                    "python_sasa": py_sasa_times,
                    "python_e2e": py_e2e_times,
                },
            }
        )

    # Display results
    console.print()

    # Table 1: SASA-only comparison
    t1 = Table(title="SASA Calculation Only (same Zig engine)")
    t1.add_column("Structure", style="cyan")
    t1.add_column("Atoms", justify="right")
    t1.add_column("CLI (ms)", justify="right", style="green")
    t1.add_column("Python (ms)", justify="right", style="blue")
    t1.add_column("Overhead", justify="right")

    for r in results:
        overhead = (
            ((r["python_sasa_ms"] / r["cli_sasa_ms"]) - 1) * 100
            if r["cli_sasa_ms"] > 0
            else 0
        )
        t1.add_row(
            r["name"],
            f"{r['n_atoms']:,}",
            f"{r['cli_sasa_ms']:.2f}",
            f"{r['python_sasa_ms']:.2f}",
            f"{overhead:+.1f}%",
        )

    console.print(t1)
    console.print()

    # Table 2: End-to-end comparison
    t2 = Table(title="End-to-End (file load + classify + SASA)")
    t2.add_column("Structure", style="cyan")
    t2.add_column("Atoms", justify="right")
    t2.add_column("CLI total (ms)", justify="right", style="green")
    t2.add_column("Python total (ms)", justify="right", style="blue")
    t2.add_column("Ratio", justify="right")

    for r in results:
        ratio = r["python_e2e_ms"] / r["cli_total_ms"] if r["cli_total_ms"] > 0 else 0
        t2.add_row(
            r["name"],
            f"{r['n_atoms']:,}",
            f"{r['cli_total_ms']:.1f}",
            f"{r['python_e2e_ms']:.1f}",
            f"{ratio:.2f}x",
        )

    console.print(t2)
    console.print()
    console.print("[dim]Both CLI and Python use ProtOr classifier.[/dim]")
    console.print(
        "[dim]SASA-only: CLI from --timing output, Python from calculate_sasa() call.[/dim]"
    )
    console.print(
        f"[dim]Threads: {'auto' if threads == 0 else threads}, Runs: {runs}[/dim]"
    )

    # Save results
    if output is None:
        output = str(get_root_dir().joinpath("benchmarks", "results", "python"))

    out_dir = Path(output)
    out_dir.mkdir(parents=True, exist_ok=True)

    config = {
        "timestamp": datetime.now().strftime("%Y-%m-%d_%H%M%S"),
        "system": {
            "os": platform.system(),
            "os_version": platform.release(),
            "arch": platform.machine(),
            "cpu_cores": os.cpu_count(),
            "cpu_model": platform.processor() or "unknown",
        },
        "parameters": {
            "runs": runs,
            "warmup": warmup,
            "threads": threads,
            "n_points": n_points,
            "use_bitmask": use_bitmask,
            "probe_radius": 1.4,
            "algorithm": "sr",
            "structures": [
                {"model_id": m, "name": n, "n_atoms": a} for m, n, _, a in pdb_files
            ],
        },
    }

    with open(out_dir.joinpath("config.json"), "w") as f:
        json.dump(config, f, indent=2)

    # Save results (without raw times for readability)
    summary = [{k: v for k, v in r.items() if k != "raw"} for r in results]
    with open(out_dir.joinpath("results.json"), "w") as f:
        json.dump(summary, f, indent=2)

    console.print(f"\n[green]Results saved to {out_dir}[/green]")


if __name__ == "__main__":
    app()
