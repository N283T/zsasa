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
"""Analyze batch benchmark results (Zig vs Rust vs FreeSASA).

Outputs two variants:
- Rust baseline: Rust vs Zig only (FreeSASA excluded for scale)
- All tools: Rust vs Zig vs FreeSASA

Usage:
    ./benchmarks/scripts/analyze_batch.py summary  # Show summary tables
    ./benchmarks/scripts/analyze_batch.py plot     # Generate comparison charts
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

TOOL_ORDER = {"freesasa": 0, "rust": 1, "zig": 2}

COLOR_MAP = {
    ("rust", "f32"): "#e74c3c",
    ("rust", "f64"): "#e74c3c",
    ("freesasa", "f64"): "#3498db",
    ("zig", "f32"): "#f39c12",
    ("zig", "f64"): "#e67e22",
}


def load_batch_data() -> pl.DataFrame:
    """Load batch results (Zig, Rust, and FreeSASA)."""
    batch_dirs = list(RESULTS_DIR.glob("batch_*/results.csv"))
    if not batch_dirs:
        raise FileNotFoundError("No batch results found (batch_*/results.csv)")

    dfs = []
    for csv_file in batch_dirs:
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


def _format_ratio(ratio: float, *, is_baseline: bool) -> str:
    """Format a speedup ratio as a Rich-styled string."""
    if is_baseline:
        return "[dim]baseline[/dim]"
    color = "green" if ratio > 1 else "red"
    return f"[{color}]{ratio:.2f}x[/{color}]"


def _print_summary_table(
    stats: pl.DataFrame,
    *,
    title_suffix: str,
    include_freesasa_baseline: bool = False,
) -> None:
    """Print a summary table with Rust as baseline (optionally also vs FreeSASA)."""
    rust = stats.filter(pl.col("tool") == "rust")
    if rust.height == 0:
        rprint("[red]No Rust results found[/red]")
        return

    rust_total = rust["total_ms"][0]
    files = rust["files"][0]
    threads = rust["threads"][0]
    rust_throughput = files / (rust_total / 1000)

    freesasa_total = None
    freesasa_throughput = None
    if include_freesasa_baseline:
        freesasa = stats.filter(pl.col("tool") == "freesasa")
        if freesasa.height > 0:
            freesasa_total = freesasa["total_ms"][0]
            freesasa_throughput = freesasa["files"][0] / (freesasa_total / 1000)

    table = Table(
        title=f"Batch Benchmark: {title_suffix} ({files:,} files, {threads} threads)"
    )
    table.add_column("Tool", style="cyan")
    table.add_column("Precision")
    table.add_column("Time", justify="right")
    if freesasa_total is not None:
        table.add_column("vs FreeSASA", justify="right")
    table.add_column("vs Rust", justify="right")
    table.add_column("Throughput", justify="right")
    if freesasa_total is not None:
        table.add_column("vs FreeSASA", justify="right")
    table.add_column("vs Rust", justify="right")

    for row in stats.iter_rows(named=True):
        tool = row["tool"]
        prec = row["precision"]
        total = row["total_ms"]
        throughput = row["files"] / (total / 1000)

        time_vs_rust = _format_ratio(rust_total / total, is_baseline=tool == "rust")
        tp_vs_rust = _format_ratio(
            throughput / rust_throughput, is_baseline=tool == "rust"
        )

        cells = [
            tool.capitalize(),
            prec,
            f"{total / 1000:.1f}s",
        ]

        if freesasa_total is not None:
            time_vs_fs = _format_ratio(
                freesasa_total / total, is_baseline=tool == "freesasa"
            )
            cells.append(time_vs_fs)

        cells.append(time_vs_rust)
        cells.append(f"{throughput:.1f}/s")

        if freesasa_total is not None:
            tp_vs_fs = _format_ratio(
                throughput / freesasa_throughput, is_baseline=tool == "freesasa"
            )
            cells.append(tp_vs_fs)

        cells.append(tp_vs_rust)

        table.add_row(*cells)

    rprint(table)


@app.command()
def summary():
    """Print summary tables with Rust as baseline."""
    df = load_batch_data()

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

    # Rust vs Zig only
    rust_zig = stats.filter(pl.col("tool") != "freesasa")
    _print_summary_table(rust_zig, title_suffix="Rust vs Zig")

    # All tools (if FreeSASA data exists)
    if stats.filter(pl.col("tool") == "freesasa").height > 0:
        rprint()
        _print_summary_table(
            stats, title_suffix="All Tools", include_freesasa_baseline=True
        )

    rprint("\n[dim]>1x = faster than baseline[/dim]")


def _plot_chart(
    stats: pl.DataFrame,
    *,
    output_name: str,
    title_suffix: str,
    baseline_tool: str = "rust",
) -> None:
    """Generate a comparison bar chart with the given tool as baseline."""
    baseline = stats.filter(pl.col("tool") == baseline_tool)
    if baseline.height == 0:
        rprint(f"[red]No {baseline_tool} results found[/red]")
        return

    baseline_total = baseline["total_ms"][0]
    baseline_throughput = baseline["files"][0] / (baseline_total / 1000)
    files = baseline["files"][0]
    threads = baseline["threads"][0]
    baseline_color = COLOR_MAP.get((baseline_tool, baseline["precision"][0]), "#95a5a6")

    # Prepare data (sorted: FreeSASA -> Rust -> Zig)
    labels = []
    total_times = []
    throughputs = []
    colors = []

    sorted_stats = stats.with_columns(
        pl.col("tool").replace_strict(TOOL_ORDER, default=99).alias("_order")
    ).sort(["_order", "precision"])

    for row in sorted_stats.iter_rows(named=True):
        tool = row["tool"]
        prec = row["precision"]
        total_sec = row["total_ms"] / 1000

        label = f"{tool.capitalize()}\n({prec})"
        labels.append(label)
        total_times.append(total_sec)
        throughputs.append(row["files"] / total_sec)
        colors.append(COLOR_MAP.get((tool, prec), "#95a5a6"))

    # Create figure with two subplots
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    x = range(len(labels))
    width = 0.6

    # Left: Wall-clock Time
    bars1 = axes[0].bar(x, total_times, color=colors, width=width, edgecolor="white")
    for bar, t in zip(bars1, total_times):
        speedup = (baseline_total / 1000) / t
        label = f"{t:.0f}s"
        if speedup != 1.0:
            label += f" ({speedup:.2f}x)"
        # Place above bar if it's short relative to the tallest
        if t < max(total_times) * 0.3:
            y_pos = bar.get_height() + (max(total_times) - min(total_times)) * 0.01
            va, text_color = "bottom", "black"
        else:
            y_pos = bar.get_height() - (max(total_times) - min(total_times)) * 0.15
            va, text_color = "top", "white"
        axes[0].text(
            bar.get_x() + bar.get_width() / 2,
            y_pos,
            label,
            ha="center",
            va=va,
            fontsize=10,
            color=text_color,
            fontweight="bold",
        )
    axes[0].axhline(
        y=baseline_total / 1000, color=baseline_color, linestyle="--", alpha=0.7
    )
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
        speedup = tp / baseline_throughput
        label = f"{tp:.1f}"
        if speedup != 1.0:
            label += f"\n({speedup:.2f}x)"
        # Place above bar if it's short relative to the tallest
        if tp < max(throughputs) * 0.3:
            y_pos = bar.get_height() + (max(throughputs) - min(throughputs)) * 0.01
            va, text_color = "bottom", "black"
        else:
            y_pos = bar.get_height() - (max(throughputs) - min(throughputs)) * 0.15
            va, text_color = "top", "white"
        axes[1].text(
            bar.get_x() + bar.get_width() / 2,
            y_pos,
            label,
            ha="center",
            va=va,
            fontsize=10,
            color=text_color,
            fontweight="bold",
        )
    axes[1].axhline(
        y=baseline_throughput, color=baseline_color, linestyle="--", alpha=0.7
    )
    axes[1].set_xticks(x)
    axes[1].set_xticklabels(labels)
    axes[1].set_ylabel("Files / second")
    axes[1].set_title("Throughput")
    axes[1].grid(True, alpha=0.3, axis="y")
    # Zoom y-axis
    tp_min = min(throughputs) * 0.95
    tp_max = max(throughputs) * 1.02
    axes[1].set_ylim(tp_min, tp_max)

    baseline_label = baseline_tool.capitalize()
    fig.suptitle(
        f"Batch SASA Benchmark: {title_suffix}"
        f" ({files:,} files, {threads} threads, vs {baseline_label})"
    )
    fig.tight_layout()
    out_path = PLOTS_DIR / output_name
    fig.savefig(out_path)
    plt.close(fig)
    rprint(f"[green]Saved:[/green] {out_path}")


@app.command()
def plot():
    """Generate comparison bar charts (Rust baseline and all tools)."""
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

    # Rust vs Zig only
    rust_zig = stats.filter(pl.col("tool") != "freesasa")
    _plot_chart(rust_zig, output_name="comparison.png", title_suffix="Rust vs Zig")

    # All tools with FreeSASA baseline (if FreeSASA data exists)
    if stats.filter(pl.col("tool") == "freesasa").height > 0:
        _plot_chart(
            stats,
            output_name="comparison_all.png",
            title_suffix="All Tools",
            baseline_tool="freesasa",
        )


@app.command()
def all():
    """Generate summary and plot."""
    summary()
    rprint()
    plot()


if __name__ == "__main__":
    app()
