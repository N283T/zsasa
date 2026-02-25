#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.12"
# dependencies = ["polars>=1.0", "matplotlib>=3.8", "numpy>=1.26", "rich>=13.0", "typer>=0.9.0"]
# ///
"""Adaptive SR benchmark: accuracy + timing.

Phase 1 (Validation): Run each config, collect per-file SASA, compare vs reference.
Phase 2 (Timing): Measure wall time with hyperfine.

Reference: f64, no bitmask, same n_points (gold standard).

Usage:
    ./benchmarks/scripts/bench_adaptive.py \
        -i benchmarks/UP000005640_9606_HUMAN_v6/pdb \
        -n human --n-points 128

Output:
    benchmarks/results/adaptive/<name>/
    ├── config.json
    ├── validation.csv       # Per-file SASA for each config
    ├── accuracy.png         # Scatter plots vs reference
    ├── stats.json           # Accuracy statistics
    └── bench_*.json         # Hyperfine timing results
"""

from __future__ import annotations

import json
import os
import platform
import shutil
import subprocess
import tempfile
from datetime import datetime
from pathlib import Path
from typing import Annotated

import typer
from rich.console import Console
from rich.table import Table

app = typer.Typer(help="Adaptive SR benchmark: accuracy + timing")
console = Console()


# ---------------------------------------------------------------------------
# Config
# ---------------------------------------------------------------------------

ADAPTIVE_CONFIGS = [
    # (label, extra_flags)
    ("adaptive_64_010_090", "--adaptive --coarse-points 64 --adaptive-low 0.10 --adaptive-high 0.90"),
    ("adaptive_64_020_080", "--adaptive --coarse-points 64 --adaptive-low 0.20 --adaptive-high 0.80"),
    ("adaptive_64_030_070", "--adaptive --coarse-points 64 --adaptive-low 0.30 --adaptive-high 0.70"),
    ("adaptive_64_040_060", "--adaptive --coarse-points 64 --adaptive-low 0.40 --adaptive-high 0.60"),
]


def get_root_dir() -> Path:
    return Path(__file__).parent.parent.parent


def get_zsasa_path() -> Path:
    return get_root_dir().joinpath("zig-out", "bin", "zsasa")


def get_system_info() -> dict:
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
                capture_output=True, text=True,
            )
            if result.returncode == 0:
                info["cpu_model"] = result.stdout.strip()
            result = subprocess.run(
                ["sysctl", "-n", "hw.memsize"],
                capture_output=True, text=True,
            )
            if result.returncode == 0:
                info["memory_gb"] = int(result.stdout.strip()) // (1024**3)
        except Exception:
            pass
    return info


# ---------------------------------------------------------------------------
# Phase 1: Validation
# ---------------------------------------------------------------------------


def run_zsasa_collect(
    zsasa: Path,
    input_dir: Path,
    n_points: int,
    threads: int,
    extra_flags: str = "",
) -> dict[str, float]:
    """Run zsasa batch, return {stem: total_sasa}."""
    with tempfile.TemporaryDirectory(prefix="adaptive_val_") as tmp:
        out_dir = Path(tmp)
        cmd = (
            f"{zsasa} batch {input_dir} {out_dir}"
            f" --n-points={n_points} --threads={threads} --precision=f64"
        )
        if extra_flags:
            cmd += f" {extra_flags}"

        console.print(f"  [dim]$ {cmd}[/]")
        proc = subprocess.run(cmd, shell=True, capture_output=True, text=True, timeout=3600)
        if proc.returncode != 0:
            console.print(f"[red]Failed: {proc.stderr[:300]}[/]")
            return {}

        results: dict[str, float] = {}
        for json_file in sorted(out_dir.glob("*.json")):
            try:
                data = json.loads(json_file.read_text())
                results[json_file.stem] = float(data["total_area"])
            except (json.JSONDecodeError, KeyError, ValueError):
                continue
        return results


def run_validation(
    zsasa: Path,
    input_dir: Path,
    n_points: int,
    threads: int,
    results_dir: Path,
) -> Path:
    """Run all configs, collect SASA, write validation.csv."""
    import polars as pl

    configs = [
        ("ref_f64", ""),
        ("bitmask_f64", "--use-bitmask"),
    ]
    configs.extend(
        (label, f"--use-bitmask {flags}")
        for label, flags in ADAPTIVE_CONFIGS
    )

    all_results: dict[str, dict[str, float]] = {}

    for label, flags in configs:
        console.print(f"[bold cyan]{label}[/]")
        all_results[label] = run_zsasa_collect(zsasa, input_dir, n_points, threads, flags)
        console.print(f"  → {len(all_results[label])} files")

    # Merge into DataFrame
    all_stems = set()
    for r in all_results.values():
        all_stems.update(r.keys())

    rows = []
    for stem in sorted(all_stems):
        row: dict = {"structure": stem}
        for label in [c[0] for c in configs]:
            row[label] = all_results[label].get(stem)
        rows.append(row)

    columns = ["structure"] + [c[0] for c in configs]
    df = pl.DataFrame(rows).select(columns)

    csv_path = results_dir.joinpath("validation.csv")
    df.write_csv(csv_path)
    console.print(f"[green]Saved:[/] {csv_path} ({df.height} structures)")
    return csv_path


# ---------------------------------------------------------------------------
# Accuracy analysis
# ---------------------------------------------------------------------------


def analyze_accuracy(csv_path: Path, results_dir: Path) -> dict:
    """Compute accuracy stats and generate plots."""
    import matplotlib.pyplot as plt
    import numpy as np
    import polars as pl

    df = pl.read_csv(csv_path)
    reference = "ref_f64"

    if reference not in df.columns:
        console.print(f"[red]Reference column '{reference}' not found[/]")
        return {}

    skip_cols = {"structure", reference}
    compare_cols = [c for c in df.columns if c not in skip_cols]

    all_stats: dict[str, dict] = {}

    # Stats table
    table = Table(title=f"Accuracy vs {reference}")
    table.add_column("Config", style="cyan")
    table.add_column("N", justify="right")
    table.add_column("R²", justify="right")
    table.add_column("Mean Err %", justify="right")
    table.add_column("Median Err %", justify="right")
    table.add_column("p95 Err %", justify="right")
    table.add_column("Max Err %", justify="right")

    for col in compare_cols:
        pair = df.select([reference, col]).drop_nulls()
        if pair.height < 2:
            continue

        ref_arr = pair[reference].to_numpy()
        comp_arr = pair[col].to_numpy()

        mask = ref_arr > 0
        ref_arr = ref_arr[mask]
        comp_arr = comp_arr[mask]

        if len(ref_arr) < 2:
            continue

        rel_errors = np.abs(comp_arr - ref_arr) / ref_arr * 100
        correlation = np.corrcoef(comp_arr, ref_arr)[0, 1]

        stats = {
            "n": int(len(ref_arr)),
            "r_squared": float(correlation ** 2),
            "mean_rel_error": float(np.mean(rel_errors)),
            "median_rel_error": float(np.median(rel_errors)),
            "p95_rel_error": float(np.percentile(rel_errors, 95)),
            "max_rel_error": float(np.max(rel_errors)),
        }
        all_stats[col] = stats

        table.add_row(
            col,
            str(stats["n"]),
            f"{stats['r_squared']:.8f}",
            f"{stats['mean_rel_error']:.4f}",
            f"{stats['median_rel_error']:.4f}",
            f"{stats['p95_rel_error']:.4f}",
            f"{stats['max_rel_error']:.4f}",
        )

    console.print()
    console.print(table)

    # Save stats
    stats_path = results_dir.joinpath("stats.json")
    stats_path.write_text(json.dumps(all_stats, indent=2))
    console.print(f"[green]Saved:[/] {stats_path}")

    # Scatter plots
    n_plots = len(compare_cols)
    if n_plots == 0:
        return all_stats

    cols_per_row = min(3, n_plots)
    n_rows = (n_plots + cols_per_row - 1) // cols_per_row
    fig, axes = plt.subplots(
        n_rows, cols_per_row,
        figsize=(7 * cols_per_row, 6 * n_rows),
        squeeze=False,
    )

    for idx, col in enumerate(compare_cols):
        ax = axes[idx // cols_per_row][idx % cols_per_row]
        pair = df.select([reference, col]).drop_nulls()
        if pair.height < 2:
            ax.set_title(f"{col}: no data")
            continue

        ref_arr = pair[reference].to_numpy()
        comp_arr = pair[col].to_numpy()
        mask = ref_arr > 0
        ref_arr = ref_arr[mask]
        comp_arr = comp_arr[mask]

        ax.scatter(comp_arr, ref_arr, alpha=0.15, s=4, color="#3498db")
        max_val = max(float(np.max(ref_arr)), float(np.max(comp_arr)))
        ax.plot([0, max_val], [0, max_val], "r--", linewidth=1, label="y = x")

        s = all_stats.get(col, {})
        stats_text = (
            f"R² = {s.get('r_squared', 0):.8f}\n"
            f"Mean = {s.get('mean_rel_error', 0):.4f}%\n"
            f"p95 = {s.get('p95_rel_error', 0):.4f}%\n"
            f"Max = {s.get('max_rel_error', 0):.4f}%"
        )
        ax.text(
            0.05, 0.95, stats_text,
            transform=ax.transAxes, fontsize=9,
            verticalalignment="top",
            bbox=dict(boxstyle="round", facecolor="wheat", alpha=0.5),
        )
        ax.set_xlabel(f"{col} SASA (Å²)")
        ax.set_ylabel(f"{reference} SASA (Å²)")
        ax.set_title(col)
        ax.set_aspect("equal")

    # Hide empty subplots
    for idx in range(n_plots, n_rows * cols_per_row):
        axes[idx // cols_per_row][idx % cols_per_row].set_visible(False)

    fig.tight_layout()
    plot_path = results_dir.joinpath("accuracy.png")
    fig.savefig(plot_path, dpi=150)
    plt.close(fig)
    console.print(f"[green]Saved:[/] {plot_path}")

    return all_stats


# ---------------------------------------------------------------------------
# Phase 2: Timing
# ---------------------------------------------------------------------------


def run_timing(
    zsasa: Path,
    input_dir: Path,
    n_points: int,
    threads: int,
    warmup: int,
    runs: int,
    results_dir: Path,
) -> None:
    """Run hyperfine timing for all configs."""
    if not shutil.which("hyperfine"):
        console.print("[red]hyperfine not found, skipping timing[/]")
        return

    configs = [
        ("bitmask_f64", f"--use-bitmask --n-points={n_points} --threads={threads} --precision=f64"),
    ]
    for label, flags in ADAPTIVE_CONFIGS:
        configs.append(
            (label, f"--use-bitmask --n-points={n_points} --threads={threads} --precision=f64 {flags}")
        )

    timing_results: list[dict] = []

    for label, flags in configs:
        json_out = results_dir.joinpath(f"bench_{label}.json")
        cmd = f"{zsasa} batch -q {input_dir} {flags}"

        console.print(f"[bold cyan]Timing: {label}[/]")
        console.print(f"  [dim]$ {cmd}[/]")

        hyperfine_cmd = [
            "hyperfine",
            "--warmup", str(warmup),
            "--runs", str(runs),
            "--export-json", str(json_out),
            cmd,
        ]

        try:
            subprocess.run(hyperfine_cmd, check=True, timeout=3600)
        except subprocess.CalledProcessError as e:
            console.print(f"[red]Failed: {e}[/]")
            continue

        if json_out.exists():
            data = json.loads(json_out.read_text())
            if data.get("results"):
                r = data["results"][0]
                timing_results.append({"name": label, **r})

    # Summary table
    if not timing_results:
        return

    table = Table(title="Timing Results")
    table.add_column("Config", style="cyan")
    table.add_column("Mean (s)", justify="right")
    table.add_column("σ (s)", justify="right")
    table.add_column("Min (s)", justify="right")
    table.add_column("Speedup", justify="right")

    baseline = next((r for r in timing_results if r["name"] == "bitmask_f64"), None)
    baseline_mean = baseline["mean"] if baseline else 1.0

    for r in timing_results:
        speedup = baseline_mean / r["mean"]
        style = "green" if speedup > 1.02 else ("red" if speedup < 0.98 else "")
        table.add_row(
            r["name"],
            f"{r['mean']:.3f}",
            f"±{r['stddev']:.3f}",
            f"{r['min']:.3f}",
            f"[{style}]{speedup:.3f}x[/{style}]" if style else f"{speedup:.3f}x",
        )

    console.print()
    console.print(table)


# ---------------------------------------------------------------------------
# Main command
# ---------------------------------------------------------------------------


@app.command()
def main(
    input_dir: Annotated[
        Path,
        typer.Option("--input", "-i", help="PDB directory", exists=True, dir_okay=True, file_okay=False),
    ],
    name: Annotated[
        str,
        typer.Option("--name", "-n", help="Benchmark name"),
    ],
    n_points: Annotated[
        int,
        typer.Option("--n-points", "-N", help="Fine n_points (default: 128)"),
    ] = 128,
    threads: Annotated[
        int,
        typer.Option("--threads", "-T", help="Thread count"),
    ] = 8,
    warmup: Annotated[
        int,
        typer.Option("--warmup", "-w", help="Hyperfine warmup runs"),
    ] = 2,
    runs: Annotated[
        int,
        typer.Option("--runs", "-r", help="Hyperfine timing runs"),
    ] = 5,
    output_dir: Annotated[
        Path | None,
        typer.Option("--output", "-o", help="Output directory"),
    ] = None,
    skip_validation: Annotated[
        bool,
        typer.Option("--skip-validation", help="Skip validation phase"),
    ] = False,
    skip_timing: Annotated[
        bool,
        typer.Option("--skip-timing", help="Skip timing phase"),
    ] = False,
) -> None:
    """Benchmark adaptive SR: accuracy + timing."""
    zsasa = get_zsasa_path()
    if not zsasa.exists():
        console.print(f"[red]zsasa not found: {zsasa}[/]")
        console.print("Run: zig build -Doptimize=ReleaseFast")
        raise typer.Exit(1)

    root = get_root_dir()
    if output_dir is None:
        results_dir = root.joinpath("benchmarks", "results", "adaptive", name)
    else:
        results_dir = output_dir
    results_dir.mkdir(parents=True, exist_ok=True)

    n_files = sum(1 for f in input_dir.iterdir() if f.is_file())

    # Save config
    config = {
        "timestamp": datetime.now().strftime("%Y-%m-%d_%H%M%S"),
        "system": get_system_info(),
        "parameters": {
            "name": name,
            "input_dir": str(input_dir),
            "n_files": n_files,
            "n_points": n_points,
            "threads": threads,
            "warmup": warmup,
            "runs": runs,
            "adaptive_configs": [
                {"label": label, "flags": flags}
                for label, flags in ADAPTIVE_CONFIGS
            ],
        },
    }
    config_path = results_dir.joinpath("config.json")
    config_path.write_text(json.dumps(config, indent=2))

    console.print(f"[bold]=== Adaptive SR Benchmark: {name} ===[/]")
    console.print(f"Input: {input_dir} ({n_files} files)")
    console.print(f"n_points: {n_points}, threads: {threads}")
    console.print(f"Output: {results_dir}")
    console.print()

    # Phase 1: Validation
    if not skip_validation:
        console.print("[bold]=== Phase 1: Validation ===[/]")
        csv_path = run_validation(zsasa, input_dir, n_points, threads, results_dir)
        analyze_accuracy(csv_path, results_dir)
        console.print()

    # Phase 2: Timing
    if not skip_timing:
        console.print("[bold]=== Phase 2: Timing ===[/]")
        run_timing(zsasa, input_dir, n_points, threads, warmup, runs, results_dir)
        console.print()

    console.print(f"[bold green]=== Done! Results: {results_dir} ===[/]")


if __name__ == "__main__":
    app()
