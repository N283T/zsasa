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
"""Analyze batch benchmark results (Zig vs Rust).

Usage:
    ./benchmarks/scripts/analyze_batch.py summary  # Show summary table
    ./benchmarks/scripts/analyze_batch.py plot     # Generate comparison chart
    ./benchmarks/scripts/analyze_batch.py all      # Both
"""

from pathlib import Path

import matplotlib.pyplot as plt
import polars as pl
import typer
from rich import print as rprint
from rich.table import Table

app = typer.Typer(help="Analyze batch benchmark results")

RESULTS_DIR = Path(__file__).parent.parent / "results"
PLOTS_DIR = RESULTS_DIR / "plots" / "batch"


def load_batch_data() -> pl.DataFrame:
    """Load batch results (Zig and Rust only)."""
    batch_dirs = list(RESULTS_DIR.glob("batch_*/results.csv"))
    if not batch_dirs:
        raise FileNotFoundError("No batch results found (batch_*/results.csv)")

    dfs = []
    for csv_file in batch_dirs:
        # Skip FreeSASA
        if "freesasa" in csv_file.parent.name:
            continue

        df = pl.read_csv(csv_file)

        # Add precision column if missing
        if "precision" not in df.columns:
            dir_name = csv_file.parent.name
            # Rust uses f32 internally
            if "_f32" in dir_name or "rust" in dir_name:
                df = df.with_columns(pl.lit("f32").alias("precision"))
            else:
                df = df.with_columns(pl.lit("f64").alias("precision"))

        # Normalize column order
        df = df.select(
            [
                "tool",
                "algorithm",
                "precision",
                "threads",
                "run",
                "files",
                "successful",
                "failed",
                "sasa_time_ms",
                "total_time_ms",
            ]
        )
        dfs.append(df)

    return pl.concat(dfs)


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


@app.command()
def summary():
    """Print summary table with Rust as baseline."""
    df = load_batch_data()

    # Aggregate by tool/algorithm/precision (mean across runs)
    stats = (
        df.group_by(["tool", "algorithm", "precision"])
        .agg(
            pl.col("files").first(),
            pl.col("threads").first(),
            pl.col("total_time_ms").mean().alias("total_ms"),
            pl.len().alias("runs"),
        )
        .sort(["algorithm", "tool", "precision"])
    )

    # Get Rust baseline
    rust = stats.filter(pl.col("tool") == "rust")
    if rust.height == 0:
        rprint("[red]No Rust results found[/red]")
        return

    rust_total = rust["total_ms"][0]
    files = rust["files"][0]
    threads = rust["threads"][0]
    rust_throughput = files / (rust_total / 1000)

    table = Table(title=f"Batch Benchmark ({files:,} files, {threads} threads)")
    table.add_column("Tool", style="cyan")
    table.add_column("Precision")
    table.add_column("Time", justify="right")
    table.add_column("vs Rust", justify="right")
    table.add_column("Throughput", justify="right")
    table.add_column("vs Rust", justify="right")

    for row in stats.iter_rows(named=True):
        tool = row["tool"]
        prec = row["precision"]
        total = row["total_ms"]
        throughput = row["files"] / (total / 1000)

        # Calculate speedup vs Rust
        time_ratio = rust_total / total
        tp_ratio = throughput / rust_throughput

        # Format speedup
        if tool == "rust":
            time_vs = "[dim]baseline[/dim]"
            tp_vs = "[dim]baseline[/dim]"
        else:
            color = "green" if time_ratio > 1 else "red"
            time_vs = f"[{color}]{time_ratio:.2f}x[/{color}]"
            color = "green" if tp_ratio > 1 else "red"
            tp_vs = f"[{color}]{tp_ratio:.2f}x[/{color}]"

        table.add_row(
            tool.capitalize(),
            prec,
            f"{total / 1000:.1f}s",
            time_vs,
            f"{throughput:.1f}/s",
            tp_vs,
        )

    rprint(table)
    rprint("\n[dim]>1x = Zig faster than Rust[/dim]")


@app.command()
def plot():
    """Generate comparison bar chart (Rust as baseline)."""
    setup_style()
    PLOTS_DIR.mkdir(parents=True, exist_ok=True)

    df = load_batch_data()

    stats = (
        df.group_by(["tool", "algorithm", "precision"])
        .agg(
            pl.col("files").first(),
            pl.col("threads").first(),
            pl.col("total_time_ms").mean().alias("total_ms"),
        )
        .sort(["tool", "precision"])
    )

    # Get Rust baseline
    rust = stats.filter(pl.col("tool") == "rust")
    if rust.height == 0:
        rprint("[red]No Rust results found[/red]")
        return

    rust_total = rust["total_ms"][0]
    rust_throughput = rust["files"][0] / (rust_total / 1000)
    files = rust["files"][0]
    threads = rust["threads"][0]

    # Prepare data
    labels = []
    total_times = []
    throughputs = []
    colors = []

    color_map = {
        ("rust", "f32"): "#e74c3c",
        ("rust", "f64"): "#e74c3c",
        ("zig", "f32"): "#27ae60",
        ("zig", "f64"): "#2ecc71",
    }

    for row in stats.sort(["tool", "precision"], descending=[True, False]).iter_rows(
        named=True
    ):
        tool = row["tool"]
        prec = row["precision"]
        total_sec = row["total_ms"] / 1000

        label = f"{tool.capitalize()}\n({prec})"
        labels.append(label)
        total_times.append(total_sec)
        throughputs.append(row["files"] / total_sec)
        colors.append(color_map.get((tool, prec), "#95a5a6"))

    # Create figure with two subplots
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    x = range(len(labels))
    width = 0.6

    # Left: Wall-clock Time
    bars1 = axes[0].bar(x, total_times, color=colors, width=width, edgecolor="white")
    for bar, t in zip(bars1, total_times):
        speedup = (rust_total / 1000) / t
        label = f"{t:.0f}s"
        if speedup != 1.0:
            label += f" ({speedup:.2f}x)"
        axes[0].text(
            bar.get_x() + bar.get_width() / 2,
            bar.get_height() - (max(total_times) - min(total_times)) * 0.15,
            label,
            ha="center",
            va="top",
            fontsize=10,
            color="white",
            fontweight="bold",
        )
    axes[0].axhline(y=rust_total / 1000, color="#e74c3c", linestyle="--", alpha=0.7)
    axes[0].set_xticks(x)
    axes[0].set_xticklabels(labels)
    axes[0].set_ylabel("Time (seconds)")
    axes[0].set_title("Wall-clock Time")
    axes[0].grid(True, alpha=0.3, axis="y")
    # Zoom y-axis
    total_min = min(total_times) * 0.95
    total_max = max(total_times) * 1.02
    axes[0].set_ylim(total_min, total_max)

    # Right: Throughput (files/sec)
    bars2 = axes[1].bar(x, throughputs, color=colors, width=width, edgecolor="white")
    for bar, tp in zip(bars2, throughputs):
        speedup = tp / rust_throughput
        label = f"{tp:.1f}"
        if speedup != 1.0:
            label += f" ({speedup:.2f}x)"
        axes[1].text(
            bar.get_x() + bar.get_width() / 2,
            bar.get_height() - (max(throughputs) - min(throughputs)) * 0.15,
            label,
            ha="center",
            va="top",
            fontsize=10,
            color="white",
            fontweight="bold",
        )
    axes[1].axhline(y=rust_throughput, color="#e74c3c", linestyle="--", alpha=0.7)
    axes[1].set_xticks(x)
    axes[1].set_xticklabels(labels)
    axes[1].set_ylabel("Files / second")
    axes[1].set_title("Throughput")
    axes[1].grid(True, alpha=0.3, axis="y")
    # Zoom y-axis
    tp_min = min(throughputs) * 0.95
    tp_max = max(throughputs) * 1.02
    axes[1].set_ylim(tp_min, tp_max)

    fig.suptitle(f"Batch SASA Benchmark ({files:,} files, {threads} threads)")
    fig.tight_layout()
    out_path = PLOTS_DIR / "comparison.png"
    fig.savefig(out_path)
    plt.close(fig)
    rprint(f"[green]Saved:[/green] {out_path}")


@app.command()
def all():
    """Generate summary and plot."""
    summary()
    rprint()
    plot()


if __name__ == "__main__":
    app()
