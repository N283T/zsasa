#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.12"
# dependencies = [
#     "polars",
#     "matplotlib",
#     "typer",
#     "rich",
# ]
# ///
"""Analyze batch benchmark results (hyperfine JSON format).

Usage:
    ./benchmarks/scripts/analyze_batch.py summary           # Show summary table
    ./benchmarks/scripts/analyze_batch.py summary -n ecoli  # Specific benchmark
    ./benchmarks/scripts/analyze_batch.py plot -n ecoli     # Generate charts
    ./benchmarks/scripts/analyze_batch.py all -n ecoli      # Both
"""

from __future__ import annotations

import json
import re
from pathlib import Path
from typing import Annotated

import matplotlib.pyplot as plt
import polars as pl
import typer
from rich import print as rprint
from rich.table import Table

app = typer.Typer(help="Analyze batch benchmark results")

RESULTS_DIR = Path(__file__).parent.parent.joinpath("results", "batch")
PLOTS_DIR = Path(__file__).parent.parent.joinpath("results", "plots", "batch")

TOOL_ORDER = {"freesasa": 0, "rustsasa": 1, "zsasa": 2}

COLOR_MAP = {
    ("rustsasa", None): "#e74c3c",
    ("freesasa", None): "#3498db",
    ("zsasa", "f32"): "#2ecc71",
    ("zsasa", "f64"): "#27ae60",
}


def parse_benchmark_name(filename: str) -> dict:
    """Parse benchmark filename to extract tool, precision, threads.

    Examples:
        bench_zsasa_f64_8t.json -> {tool: zsasa, precision: f64, threads: 8}
        bench_freesasa_1t.json -> {tool: freesasa, precision: None, threads: 1}
        bench_rustsasa_8t.json -> {tool: rustsasa, precision: None, threads: 8}
    """
    name = filename.replace("bench_", "").replace(".json", "")

    # Extract threads (e.g., "8t" -> 8)
    threads_match = re.search(r"(\d+)t$", name)
    threads = int(threads_match.group(1)) if threads_match else 1
    name = re.sub(r"_?\d+t$", "", name)

    # Extract precision (f32 or f64)
    precision = None
    if "_f32" in name:
        precision = "f32"
        name = name.replace("_f32", "")
    elif "_f64" in name:
        precision = "f64"
        name = name.replace("_f64", "")

    return {"tool": name, "precision": precision, "threads": threads}


def load_batch_data(benchmark_name: str | None = None) -> pl.DataFrame:
    """Load batch results from hyperfine JSON files."""
    if benchmark_name:
        search_dirs = [RESULTS_DIR.joinpath(benchmark_name)]
    else:
        search_dirs = list(RESULTS_DIR.iterdir())

    rows = []
    for bench_dir in search_dirs:
        if not bench_dir.is_dir():
            continue

        for json_file in bench_dir.glob("bench_*.json"):
            try:
                with open(json_file) as f:
                    data = json.load(f)

                if not data.get("results"):
                    continue

                result = data["results"][0]
                info = parse_benchmark_name(json_file.name)

                rows.append(
                    {
                        "benchmark": bench_dir.name,
                        "tool": info["tool"],
                        "precision": info["precision"],
                        "threads": info["threads"],
                        "mean_s": result["mean"],
                        "stddev_s": result["stddev"],
                        "min_s": result["min"],
                        "max_s": result["max"],
                        "runs": len(result.get("times", [])),
                    }
                )
            except (json.JSONDecodeError, KeyError) as e:
                rprint(f"[yellow]Warning: Failed to parse {json_file}: {e}[/yellow]")
                continue

    if not rows:
        raise FileNotFoundError(
            f"No batch results found in {RESULTS_DIR}"
            + (f"/{benchmark_name}" if benchmark_name else "")
        )

    return pl.DataFrame(rows)


def setup_style():
    """Set up matplotlib style."""
    plt.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.size": 12,
            "axes.titlesize": 14,
            "axes.labelsize": 12,
            "figure.dpi": 150,
            "savefig.dpi": 150,
            "savefig.bbox": "tight",
            "axes.spines.top": False,
            "axes.spines.right": False,
        }
    )


def _format_ratio(ratio: float, *, is_baseline: bool) -> str:
    """Format a speedup ratio as a Rich-styled string."""
    if is_baseline:
        return "[dim]baseline[/dim]"
    color = "green" if ratio > 1 else "red"
    return f"[{color}]{ratio:.2f}x[/{color}]"


@app.command()
def summary(
    name: Annotated[
        str | None,
        typer.Option("--name", "-n", help="Benchmark name (e.g., ecoli)"),
    ] = None,
):
    """Print summary table with FreeSASA as baseline."""
    df = load_batch_data(name)

    # Get FreeSASA baseline
    freesasa = df.filter(pl.col("tool") == "freesasa")
    if freesasa.height == 0:
        rprint("[yellow]No FreeSASA baseline found, using slowest as baseline[/yellow]")
        baseline_time = df["mean_s"].max()
    else:
        baseline_time = freesasa["mean_s"][0]

    # Sort by tool order then by precision
    df = df.with_columns(
        pl.col("tool").replace_strict(TOOL_ORDER, default=99).alias("_order")
    ).sort(["_order", "precision", "threads"], descending=[False, False, True])

    benchmark = df["benchmark"][0] if df.height > 0 else "unknown"

    table = Table(title=f"Batch Benchmark: {benchmark}")
    table.add_column("Tool", style="cyan")
    table.add_column("Precision")
    table.add_column("Threads", justify="right")
    table.add_column("Mean", justify="right")
    table.add_column("Std Dev", justify="right")
    table.add_column("vs FreeSASA", justify="right")

    for row in df.iter_rows(named=True):
        tool = row["tool"]
        prec = row["precision"] or "-"
        threads = str(row["threads"])
        mean = row["mean_s"]
        stddev = row["stddev_s"]

        is_baseline = tool == "freesasa"
        speedup = _format_ratio(baseline_time / mean, is_baseline=is_baseline)

        table.add_row(
            tool,
            prec,
            threads,
            f"{mean:.3f}s",
            f"±{stddev:.3f}",
            speedup,
        )

    rprint(table)
    rprint("\n[dim]>1x = faster than FreeSASA baseline[/dim]")


@app.command()
def plot(
    name: Annotated[
        str | None,
        typer.Option("--name", "-n", help="Benchmark name (e.g., ecoli)"),
    ] = None,
):
    """Generate comparison bar chart."""
    setup_style()
    df = load_batch_data(name)

    # Filter to multi-threaded results only for cleaner chart
    df_mt = df.filter(pl.col("threads") > 1)
    if df_mt.height == 0:
        df_mt = df

    # Sort by tool order then by precision
    df_mt = df_mt.with_columns(
        pl.col("tool").replace_strict(TOOL_ORDER, default=99).alias("_order")
    ).sort(["_order", "precision"])

    # Get FreeSASA baseline for speedup calculation
    freesasa = df.filter(pl.col("tool") == "freesasa")
    baseline_time = freesasa["mean_s"][0] if freesasa.height > 0 else df["mean_s"].max()

    # Prepare data
    labels = []
    times = []
    colors = []

    for row in df_mt.iter_rows(named=True):
        tool = row["tool"]
        prec = row["precision"]
        threads = row["threads"]

        if prec:
            label = f"{tool}\n({prec}, {threads}t)"
        else:
            label = f"{tool}\n({threads}t)"

        labels.append(label)
        times.append(row["mean_s"])
        colors.append(COLOR_MAP.get((tool, prec), "#95a5a6"))

    # Create figure
    fig, ax = plt.subplots(figsize=(10, 6))

    x = range(len(labels))
    bars = ax.bar(x, times, color=colors, width=0.6, edgecolor="white")

    # Add labels on bars
    for bar, t in zip(bars, times):
        speedup = baseline_time / t
        label = f"{t:.2f}s\n({speedup:.1f}x)"
        ax.text(
            bar.get_x() + bar.get_width() / 2,
            bar.get_height() + 0.5,
            label,
            ha="center",
            va="bottom",
            fontsize=10,
            fontweight="bold",
        )

    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.set_ylabel("Time (seconds)")
    ax.grid(True, alpha=0.3, axis="y")

    benchmark = df["benchmark"][0] if df.height > 0 else "unknown"
    ax.set_title(f"Batch Benchmark: {benchmark} (vs FreeSASA)")

    # Save
    PLOTS_DIR.mkdir(parents=True, exist_ok=True)
    out_path = PLOTS_DIR.joinpath(f"{benchmark}.png")
    fig.savefig(out_path)
    plt.close(fig)
    rprint(f"[green]Saved:[/green] {out_path}")


@app.command("all")
def all_commands(
    name: Annotated[
        str | None,
        typer.Option("--name", "-n", help="Benchmark name (e.g., ecoli)"),
    ] = None,
):
    """Generate summary and plot."""
    summary(name)
    rprint()
    plot(name)


if __name__ == "__main__":
    app()
