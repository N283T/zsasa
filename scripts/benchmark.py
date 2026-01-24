#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.11"
# dependencies = [
#     "freesasa>=2.0",
#     "numpy>=1.24",
# ]
# ///
"""Unified benchmark comparing Zig implementation with FreeSASA across all structures.

Measures performance across multiple structure sizes with both algorithms.
Compares implementations:
- Zig CLI (subprocess)
- Zig Python bindings (library)
- FreeSASA Python (library)
- FreeSASA C (optional, subprocess with multi-threading support)

Usage:
    ./scripts/benchmark.py [--runs=N] [--threads=N] [--structure=PDB]
    ./scripts/benchmark.py --use-c [--fs-c-path=PATH] [--threads=N]

Examples:
    ./scripts/benchmark.py                      # Run all benchmarks
    ./scripts/benchmark.py --runs=5             # 5 runs per benchmark
    ./scripts/benchmark.py --structure=1a0q     # Single structure only
    ./scripts/benchmark.py --use-c              # Include FreeSASA C benchmark
    ./scripts/benchmark.py --use-c --threads=4  # FreeSASA C with 4 threads

Note on benchmark fairness:
    Execution order is: Zig CLI -> Zig Python -> FreeSASA Python

    CPU cache warming effect: FreeSASA runs after Zig, benefiting from cached data
    (coordinates, NumPy internals, math functions, allocator state). Testing shows
    FreeSASA is ~20% faster when run after Zig vs. cold start. This means the
    benchmark slightly favors FreeSASA, making Zig's speedup numbers conservative.

    Memory interference: Tested and confirmed no negative impact. Zig's allocations
    are properly freed via defer, and Python GC has no significant effect.
"""

from __future__ import annotations

import json
import re
import subprocess
import sys
import tempfile
import time
from dataclasses import dataclass
from pathlib import Path

import freesasa
import numpy as np

# Add the python package to path for Zig Python bindings
sys.path.insert(0, str(Path(__file__).parent.parent / "python"))
try:
    from freesasa_zig import calculate_sasa as zig_calculate_sasa

    ZIG_PYTHON_AVAILABLE = True
except ImportError:
    ZIG_PYTHON_AVAILABLE = False


@dataclass
class BenchmarkResult:
    """Result of a single benchmark run."""

    pdb_id: str
    n_atoms: int
    algorithm: str
    tool: str
    time_ms: float
    total_area: float
    sasa_only_ms: float | None = (
        None  # SASA calculation time only (for fair comparison)
    )


def run_zig_benchmark(
    input_path: Path,
    algorithm: str = "sr",
    n_threads: int = 0,
    classifier: str = "protor",
) -> tuple[float, float, dict[str, float]]:
    """Run Zig benchmark with timing. Returns (time_ms, total_area, timing_breakdown)."""
    zig_binary = Path(__file__).parent.parent / "zig-out" / "bin" / "freesasa_zig"

    with tempfile.NamedTemporaryFile(suffix=".json", delete=False) as f:
        output_path = Path(f.name)

    try:
        cmd = [
            str(zig_binary),
            "--timing",
            f"--algorithm={algorithm}",
            f"--classifier={classifier}",
            str(input_path),
            str(output_path),
        ]
        if n_threads > 0:
            cmd.insert(3, f"--threads={n_threads}")

        result = subprocess.run(cmd, capture_output=True, text=True, timeout=600)

        if result.returncode != 0:
            raise RuntimeError(f"Zig failed: {result.stderr}")

        # Parse timing from stderr
        timing = {}
        for line in result.stderr.split("\n"):
            if ":" in line and "ms" in line:
                match = re.match(r"\s*(.+?):\s*([\d.]+)\s*ms", line)
                if match:
                    key = match.group(1).strip().lower().replace(" ", "_")
                    timing[key] = float(match.group(2))

        # Get total time
        total_time = timing.get("total", 0.0)

        # Parse output
        with open(output_path) as f:
            data = json.load(f)

        return total_time, data["total_area"], timing

    finally:
        output_path.unlink(missing_ok=True)


def run_zig_python_benchmark(
    input_path: Path,
    algorithm: str = "sr",
    n_points: int = 100,
    n_slices: int = 20,
    n_threads: int = 0,
) -> tuple[float, float]:
    """Run Zig Python bindings benchmark. Returns (time_ms, total_area).

    Measures only the SASA calculation time for fair comparison.
    """
    if not ZIG_PYTHON_AVAILABLE:
        raise RuntimeError("Zig Python bindings not available")

    # Load input data
    with open(input_path) as f:
        data = json.load(f)

    coords = np.column_stack([data["x"], data["y"], data["z"]])
    radii = np.array(data["r"])

    # Time only the SASA calculation
    start = time.perf_counter()
    if algorithm == "sr":
        result = zig_calculate_sasa(
            coords, radii, algorithm="sr", n_points=n_points, n_threads=n_threads
        )
    else:
        result = zig_calculate_sasa(
            coords, radii, algorithm="lr", n_slices=n_slices, n_threads=n_threads
        )
    elapsed = time.perf_counter() - start

    return elapsed * 1000, result.total_area


def run_freesasa_python_benchmark(
    input_path: Path,
    algorithm: str = "sr",
    n_points: int = 100,
    n_slices: int = 20,
) -> tuple[float, float]:
    """Run FreeSASA Python benchmark. Returns (time_ms, total_area).

    Measures only the freesasa.calc() time for fair comparison.
    """
    # Load pre-computed ProtOr radii input
    with open(input_path) as f:
        data = json.load(f)

    n_atoms = len(data["x"])
    coords = []
    for i in range(n_atoms):
        coords.extend([data["x"][i], data["y"][i], data["z"][i]])

    radii = data["r"]

    # Set algorithm parameters
    if algorithm == "sr":
        params = freesasa.Parameters(
            {"algorithm": freesasa.ShrakeRupley, "n-points": n_points}
        )
    else:
        params = freesasa.Parameters(
            {"algorithm": freesasa.LeeRichards, "n-slices": n_slices}
        )

    # Time only the SASA calculation
    start = time.perf_counter()
    result = freesasa.calcCoord(coords, radii, params)
    elapsed = time.perf_counter() - start

    return elapsed * 1000, result.totalArea()


def run_freesasa_c_benchmark(
    cif_path: Path,
    algorithm: str = "sr",
    n_points: int = 100,
    n_slices: int = 20,
    n_threads: int = 1,
    fs_c_binary: Path | None = None,
) -> tuple[float, float]:
    """Run FreeSASA C benchmark. Returns (time_ms, total_area).

    Requires FreeSASA C binary compiled with thread support.
    """
    if fs_c_binary is None:
        # Default location: freesasa-c/src/freesasa in project root
        fs_c_binary = Path(__file__).parent.parent / "freesasa-c" / "src" / "freesasa"

    if not fs_c_binary.exists():
        raise FileNotFoundError(f"FreeSASA C binary not found: {fs_c_binary}")

    if not cif_path.exists():
        raise FileNotFoundError(f"CIF file not found: {cif_path}")

    # Build command
    cmd = [
        str(fs_c_binary),
        "--cif",
        "--radii=protor",
        f"--n-threads={n_threads}",
    ]

    if algorithm == "sr":
        cmd.extend(["--shrake-rupley", f"--resolution={n_points}"])
    else:
        cmd.extend(["--lee-richards", f"--resolution={n_slices}"])

    cmd.append(str(cif_path))

    # Run and time
    start = time.perf_counter()
    result = subprocess.run(cmd, capture_output=True, text=True, timeout=600)
    elapsed = time.perf_counter() - start

    if result.returncode != 0:
        raise RuntimeError(f"FreeSASA C failed: {result.stderr}")

    # Parse total area from output
    total_area = 0.0
    for line in result.stdout.split("\n"):
        if line.startswith("Total"):
            match = re.search(r"[\d.]+", line)
            if match:
                total_area = float(match.group())
                break

    return elapsed * 1000, total_area


def download_cif_if_needed(pdb_id: str, cif_dir: Path) -> Path:
    """Download CIF file from RCSB if not present."""
    cif_path = cif_dir / f"{pdb_id}.cif"
    if cif_path.exists():
        return cif_path

    import urllib.error
    import urllib.request

    url = f"https://files.rcsb.org/download/{pdb_id.upper()}.cif"
    print(f"  Downloading {pdb_id}.cif...")
    cif_dir.mkdir(parents=True, exist_ok=True)

    try:
        urllib.request.urlretrieve(url, cif_path)
    except (urllib.error.URLError, urllib.error.HTTPError) as e:
        # Clean up partial download
        cif_path.unlink(missing_ok=True)
        raise FileNotFoundError(f"Failed to download {pdb_id}.cif: {e}") from e

    return cif_path


def run_benchmarks(
    structures: list[tuple[str, str, str]],
    base_dir: Path,
    n_runs: int = 3,
    n_threads: int = 0,
    use_c: bool = False,
    fs_c_path: Path | None = None,
) -> list[BenchmarkResult]:
    """Run all benchmarks and return results."""
    results = []
    algorithms = ["sr", "lr"]

    for pdb_id, category, description in structures:
        # Use ProtOr-radii inputs for both Zig and FreeSASA
        input_path = base_dir / "inputs_protor" / f"{pdb_id}.json"

        if not input_path.exists():
            print(f"  {pdb_id}: Skipping (input not found)")
            continue

        # Get atom count
        with open(input_path) as f:
            data = json.load(f)
            n_atoms = len(data.get("x", []))

        print(f"\n{pdb_id.upper()} ({n_atoms} atoms):")

        for algo in algorithms:
            # Zig benchmark
            zig_times = []
            zig_sasa_times = []
            zig_area = 0.0
            for _ in range(n_runs):
                try:
                    t, area, timing = run_zig_benchmark(
                        input_path, algo, n_threads, "protor"
                    )
                    zig_times.append(t)
                    zig_area = area
                    # Extract SASA-only time from timing breakdown
                    if "sasa_calculation" in timing:
                        zig_sasa_times.append(timing["sasa_calculation"])
                except Exception as e:
                    print(f"    Zig {algo.upper()}: ERROR - {e}")
                    break

            if zig_times:
                avg_time = sum(zig_times) / len(zig_times)
                avg_sasa_time = (
                    sum(zig_sasa_times) / len(zig_sasa_times)
                    if zig_sasa_times
                    else None
                )
                results.append(
                    BenchmarkResult(
                        pdb_id,
                        n_atoms,
                        algo,
                        "zig-cli",
                        avg_time,
                        zig_area,
                        avg_sasa_time,
                    )
                )
                sasa_str = f", SASA: {avg_sasa_time:.2f}" if avg_sasa_time else ""
                print(
                    f"  Zig CLI {algo.upper():2s}: {avg_time:8.2f} ms{sasa_str} (area: {zig_area:.2f})"
                )

            # Zig Python bindings benchmark
            if not ZIG_PYTHON_AVAILABLE:
                print(f"  Zig Py  {algo.upper():2s}: SKIPPED (bindings not available)")
            elif ZIG_PYTHON_AVAILABLE:
                zig_py_times = []
                zig_py_area = 0.0
                for _ in range(n_runs):
                    try:
                        t, area = run_zig_python_benchmark(
                            input_path, algo, n_threads=n_threads
                        )
                        zig_py_times.append(t)
                        zig_py_area = area
                    except Exception as e:
                        print(f"    Zig Py {algo.upper()}: ERROR - {e}")
                        break

                if zig_py_times:
                    avg_time = sum(zig_py_times) / len(zig_py_times)
                    results.append(
                        BenchmarkResult(
                            pdb_id, n_atoms, algo, "zig-py", avg_time, zig_py_area
                        )
                    )
                    print(
                        f"  Zig Py  {algo.upper():2s}: {avg_time:8.2f} ms (area: {zig_py_area:.2f})"
                    )

            # FreeSASA Python benchmark
            fs_times = []
            fs_area = 0.0
            for _ in range(n_runs):
                try:
                    t, area = run_freesasa_python_benchmark(input_path, algo)
                    fs_times.append(t)
                    fs_area = area
                except Exception as e:
                    print(f"    FreeSASA {algo.upper()}: ERROR - {e}")
                    break

            if fs_times:
                avg_time = sum(fs_times) / len(fs_times)
                results.append(
                    BenchmarkResult(
                        pdb_id, n_atoms, algo, "freesasa", avg_time, fs_area
                    )
                )
                print(
                    f"  FS      {algo.upper():2s}: {avg_time:8.2f} ms (area: {fs_area:.2f})"
                )

            # FreeSASA C benchmark (optional)
            if use_c:
                cif_dir = base_dir / "cif"
                try:
                    cif_path = download_cif_if_needed(pdb_id, cif_dir)
                    fsc_times = []
                    fsc_area = 0.0
                    for _ in range(n_runs):
                        try:
                            t, area = run_freesasa_c_benchmark(
                                cif_path,
                                algo,
                                n_threads=n_threads if n_threads > 0 else 1,
                                fs_c_binary=fs_c_path,
                            )
                            fsc_times.append(t)
                            fsc_area = area
                        except Exception as e:
                            print(f"    FS-C {algo.upper()}: ERROR - {e}")
                            break

                    if fsc_times:
                        avg_time = sum(fsc_times) / len(fsc_times)
                        results.append(
                            BenchmarkResult(
                                pdb_id, n_atoms, algo, "freesasa-c", avg_time, fsc_area
                            )
                        )
                        print(
                            f"  FS-C    {algo.upper():2s}: {avg_time:8.2f} ms (area: {fsc_area:.2f})"
                        )
                except FileNotFoundError as e:
                    print(f"    FS-C {algo.upper()}: SKIP - {e}")

    return results


def print_summary(results: list[BenchmarkResult]) -> None:
    """Print benchmark summary table."""
    print("\n" + "=" * 100)
    print("BENCHMARK SUMMARY (SASA-only for fair comparison)")
    print("=" * 100)

    # Group by PDB
    pdbs = sorted(
        set(r.pdb_id for r in results),
        key=lambda p: next(r.n_atoms for r in results if r.pdb_id == p),
    )

    # Check if we have Zig Python results
    has_zig_py = any(r.tool == "zig-py" for r in results)

    if has_zig_py:
        print(
            f"\n{'PDB':<6} {'Atoms':>7} {'Algo':<3} "
            f"{'Zig CLI':>10} {'Zig Py':>10} {'FreeSASA':>10} "
            f"{'CLI vs FS':>10} {'Py vs FS':>10}"
        )
        print("-" * 90)
    else:
        print(
            f"\n{'PDB':<8} {'Atoms':>8} {'Algo':<4} "
            f"{'Zig SASA':>10} {'FS (ms)':>10} {'Speedup':>10}"
        )
        print("-" * 60)

    for pdb in pdbs:
        pdb_results = [r for r in results if r.pdb_id == pdb]
        n_atoms = pdb_results[0].n_atoms if pdb_results else 0

        for algo in ["sr", "lr"]:
            zig_cli_r = next(
                (r for r in pdb_results if r.algorithm == algo and r.tool == "zig-cli"),
                None,
            )
            zig_py_r = next(
                (r for r in pdb_results if r.algorithm == algo and r.tool == "zig-py"),
                None,
            )
            fs_r = next(
                (
                    r
                    for r in pdb_results
                    if r.algorithm == algo and r.tool == "freesasa"
                ),
                None,
            )

            if has_zig_py and zig_cli_r and zig_py_r and fs_r:
                # Use SASA-only time for CLI (or total if not available)
                cli_time = (
                    zig_cli_r.sasa_only_ms
                    if zig_cli_r.sasa_only_ms
                    else zig_cli_r.time_ms
                )
                py_time = zig_py_r.time_ms
                fs_time = fs_r.time_ms
                cli_speedup = fs_time / cli_time if cli_time > 0 else 0
                py_speedup = fs_time / py_time if py_time > 0 else 0
                print(
                    f"{pdb:<6} {n_atoms:>7} {algo.upper():<3} "
                    f"{cli_time:>9.2f}ms {py_time:>9.2f}ms {fs_time:>9.2f}ms "
                    f"{cli_speedup:>9.2f}x {py_speedup:>9.2f}x"
                )
            elif zig_cli_r and fs_r:
                # Fallback: no Zig Python results
                zig_time = (
                    zig_cli_r.sasa_only_ms
                    if zig_cli_r.sasa_only_ms
                    else zig_cli_r.time_ms
                )
                speedup = fs_r.time_ms / zig_time if zig_time > 0 else 0
                print(
                    f"{pdb:<8} {n_atoms:>8} {algo.upper():<4} "
                    f"{zig_time:>10.2f} {fs_r.time_ms:>10.2f} {speedup:>9.2f}x"
                )


def main() -> int:
    # Parse arguments
    n_runs = 3
    n_threads = 0
    structure_filter = None
    use_c = False
    fs_c_path = None

    for arg in sys.argv[1:]:
        if arg.startswith("--runs="):
            n_runs = int(arg.split("=")[1])
        elif arg.startswith("--threads="):
            n_threads = int(arg.split("=")[1])
        elif arg.startswith("--structure="):
            structure_filter = arg.split("=")[1].lower()
        elif arg == "--use-c":
            use_c = True
        elif arg.startswith("--fs-c-path="):
            fs_c_path = Path(arg.split("=")[1])
            use_c = True  # Implicitly enable C benchmark

    # Available structures
    structures = [
        ("1crn", "tiny", "Crambin"),
        ("1ubq", "small", "Ubiquitin"),
        ("1a0q", "medium", "Lipid transfer protein"),
        ("3hhb", "medium", "Hemoglobin"),
        ("1aon", "large", "GroEL-GroES"),
        ("4v6x", "xlarge", "Ribosome"),
    ]

    if structure_filter:
        structures = [(p, c, d) for p, c, d in structures if p == structure_filter]
        if not structures:
            print(f"Error: Unknown structure: {structure_filter}")
            return 1

    base_dir = Path(__file__).parent.parent / "benchmarks"

    print("=" * 80)
    print(
        f"SASA Benchmark (runs={n_runs}, threads={n_threads or 'auto'}, use_c={use_c})"
    )
    print("=" * 80)

    results = run_benchmarks(structures, base_dir, n_runs, n_threads, use_c, fs_c_path)
    print_summary(results)

    return 0


if __name__ == "__main__":
    sys.exit(main())
