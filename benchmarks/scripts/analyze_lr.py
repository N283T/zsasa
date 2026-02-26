#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.12"
# dependencies = [
#     "polars>=1.0",
#     "matplotlib>=3.8",
#     "typer>=0.9.0",
#     "rich>=13.0",
# ]
# ///
"""LR benchmark analysis CLI.

Generates summary statistics and plots from bench_lr.py results.

Usage:
    ./benchmarks/scripts/analyze_lr.py summary       # Summary tables
    ./benchmarks/scripts/analyze_lr.py all            # All plots + summary
    ./benchmarks/scripts/analyze_lr.py scatter        # Atoms vs time
    ./benchmarks/scripts/analyze_lr.py threads        # Thread scaling
    ./benchmarks/scripts/analyze_lr.py validation     # SASA validation
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import polars as pl
import typer
from rich import print as rprint
from rich.table import Table

app = typer.Typer(help="LR benchmark analysis CLI")

_BENCHMARKS_DIR = Path(__file__).parent.parent
RESULTS_BASE_LR = _BENCHMARKS_DIR.joinpath("results", "single_lr")
PLOTS_DIR = _BENCHMARKS_DIR.joinpath("results", "plots", "lr")

COLORS = {
    "zsasa": "#f39c12",
    "zsasa_f32": "#e67e22",
    "zsasa_f64": "#f39c12",
    "freesasa": "#3498db",
}

LINESTYLES = {
    "zsasa": "-",
    "zsasa_f32": "--",
    "zsasa_f64": "-",
    "freesasa": "-",
}

DISPLAY_NAMES = {
    "zsasa_f64": "zsasa (f64)",
    "zsasa_f32": "zsasa (f32)",
    "zsasa": "zsasa",
    "freesasa": "FreeSASA",
}


def display_name(tool_label: str) -> str:
    """Get display name for a tool label."""
    return DISPLAY_NAMES.get(tool_label, tool_label)


def setup_style():
    """Set up matplotlib style for clean plots."""
    plt.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.size": 11,
            "axes.titlesize": 13,
            "axes.labelsize": 11,
            "xtick.labelsize": 10,
            "ytick.labelsize": 10,
            "legend.fontsize": 10,
            "figure.titlesize": 14,
            "figure.dpi": 150,
            "savefig.dpi": 150,
            "savefig.bbox": "tight",
            "axes.spines.top": False,
            "axes.spines.right": False,
        }
    )


def load_lr_data(n_slices: int = 20) -> pl.DataFrame:
    """Load LR benchmark results from results/single_lr/<n_slices>/**/results.csv.

    Derives tool_label from the CSV data (tool + precision for zig).
    Aggregates multiple runs into mean values per (structure, threads).
    """
    lr_dir = RESULTS_BASE_LR.joinpath(str(n_slices))
    csv_files = sorted(lr_dir.glob("*/results.csv"))

    if not csv_files:
        raise FileNotFoundError(f"No results.csv files found in {lr_dir}")

    dfs = []
    for f in csv_files:
        df = pl.read_csv(f)
        dfs.append(df)

    df = pl.concat(dfs, how="diagonal")

    # Ensure precision column exists
    if "precision" not in df.columns:
        df = df.with_columns(pl.lit(None).cast(pl.Utf8).alias("precision"))

    # Remap tool names: zig -> zsasa
    df = df.with_columns(
        pl.col("tool").replace({"zig": "zsasa"}).alias("tool"),
    )

    # Create tool_label: zsasa_f64, zsasa_f32, freesasa
    df = df.with_columns(
        pl.when((pl.col("tool") == "zsasa") & pl.col("precision").is_not_null())
        .then(pl.concat_str([pl.col("tool"), pl.lit("_"), pl.col("precision")]))
        .otherwise(pl.col("tool"))
        .alias("tool_label")
    )

    # Aggregate by structure (mean across runs)
    return (
        df.group_by(
            [
                "tool",
                "tool_label",
                "structure",
                "n_atoms",
                "algorithm",
                "precision",
                "threads",
            ]
        )
        .agg(
            pl.col("sasa_time_ms").mean().alias("time_ms"),
            pl.col("sasa_time_ms").std().alias("time_std"),
            pl.col("sasa_time_ms").len().alias("n_runs"),
            pl.col("total_sasa").first(),
        )
        .sort(["tool_label", "n_atoms"])
    )


# === CLI Commands ===


@app.command()
def summary(
    n_slices: int = typer.Option(
        20, "--n-slices", help="Number of slices per atom diameter"
    ),
):
    """Print summary statistics table."""
    df = load_lr_data(n_slices)

    df_t1 = df.filter(pl.col("threads") == 1)

    table = Table(title="LR Single-Threaded Performance Summary")
    table.add_column("Tool", style="green")
    table.add_column("Structures", justify="right")
    table.add_column("Median (ms)", justify="right")
    table.add_column("Mean (ms)", justify="right")
    table.add_column("Std (ms)", justify="right")
    table.add_column("P95 (ms)", justify="right")

    stats = (
        df_t1.group_by("tool_label")
        .agg(
            pl.len().alias("n"),
            pl.col("time_ms").median().alias("median"),
            pl.col("time_ms").mean().alias("mean"),
            pl.col("time_ms").std().alias("std"),
            pl.col("time_ms").quantile(0.95).alias("p95"),
        )
        .sort("tool_label")
    )

    for row in stats.iter_rows(named=True):
        table.add_row(
            row["tool_label"],
            f"{row['n']:,}",
            f"{row['median']:.2f}",
            f"{row['mean']:.2f}",
            f"{row['std']:.2f}",
            f"{row['p95']:.2f}",
        )

    rprint(table)

    # Speedup summary
    rprint("\n[bold]Speedup Ratios (zsasa vs others):[/bold]")

    pivot = (
        df_t1.select(["structure", "tool_label", "time_ms"])
        .pivot(on="tool_label", index="structure", values="time_ms")
        .drop_nulls()
    )

    zsasa_col = "zsasa_f64" if "zsasa_f64" in pivot.columns else "zsasa"

    if zsasa_col in pivot.columns and "freesasa" in pivot.columns:
        speedup = pivot.select(
            (pl.col("freesasa") / pl.col(zsasa_col)).alias("speedup")
        )
        med = speedup["speedup"].median()
        mean = speedup["speedup"].mean()
        rprint(
            f"  LR: zsasa vs FreeSASA = [green]{med:.2f}x[/green] (median), {mean:.2f}x (mean)"
        )

    if "zsasa_f32" in pivot.columns and "zsasa_f64" in pivot.columns:
        speedup = pivot.select(
            (pl.col("zsasa_f64") / pl.col("zsasa_f32")).alias("speedup")
        )
        med = speedup["speedup"].median()
        rprint(
            f"  LR: zsasa(f32) vs zsasa(f64) = [green]{med:.2f}x[/green] (median)"
        )


@app.command()
def scatter(
    n_slices: int = typer.Option(
        20, "--n-slices", help="Number of slices per atom diameter"
    ),
):
    """Generate atoms vs time scatter plot."""
    setup_style()
    df = load_lr_data(n_slices)

    plot_dir = PLOTS_DIR.joinpath("scatter")
    individual_dir = plot_dir.joinpath("individual")
    individual_dir.mkdir(parents=True, exist_ok=True)

    thread_counts = sorted(df["threads"].unique().to_list())

    for threads in thread_counts:
        df_t = df.filter(pl.col("threads") == threads)
        df_sampled = df_t.sample(n=min(5000, df_t.height), seed=42)

        fig, ax = plt.subplots(figsize=(10, 6))
        tool_labels = [
            t
            for t in sorted(df_sampled["tool_label"].unique().to_list())
            if "f32" not in t
        ]
        for tool_label in tool_labels:
            df_tool = df_sampled.filter(pl.col("tool_label") == tool_label)
            ax.scatter(
                df_tool["n_atoms"].to_list(),
                df_tool["time_ms"].to_list(),
                label=display_name(tool_label),
                alpha=0.4,
                s=10,
                color=COLORS.get(tool_label, "#95a5a6"),
            )
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel("Number of Atoms")
        ax.set_ylabel("Execution Time (ms)")
        ax.set_title(f"LR: Time vs Size (threads={threads})")
        ax.legend()
        fig.tight_layout()
        out_path = individual_dir.joinpath(f"t{threads}.png")
        fig.savefig(out_path)
        plt.close(fig)
        rprint(f"[green]Saved:[/green] {out_path}")

    # Grid of all threads
    n_threads = len(thread_counts)
    n_cols = min(3, n_threads)
    n_rows = (n_threads + n_cols - 1) // n_cols

    fig, axes = plt.subplots(n_rows, n_cols, figsize=(6 * n_cols, 5 * n_rows))
    if n_threads == 1:
        axes = [[axes]]
    elif n_rows == 1:
        axes = [axes]

    for idx, threads in enumerate(thread_counts):
        row, col = idx // n_cols, idx % n_cols
        ax = axes[row][col] if n_rows > 1 else axes[0][col]
        df_t = df.filter(pl.col("threads") == threads)
        df_sampled = df_t.sample(n=min(5000, df_t.height), seed=42)
        tool_labels = [
            t
            for t in sorted(df_sampled["tool_label"].unique().to_list())
            if "f32" not in t
        ]
        for tool_label in tool_labels:
            df_tool = df_sampled.filter(pl.col("tool_label") == tool_label)
            ax.scatter(
                df_tool["n_atoms"].to_list(),
                df_tool["time_ms"].to_list(),
                label=display_name(tool_label),
                alpha=0.4,
                s=10,
                color=COLORS.get(tool_label, "#95a5a6"),
            )
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel("Number of Atoms")
        ax.set_ylabel("Execution Time (ms)")
        ax.set_title(f"threads={threads}")
        ax.legend()

    for idx in range(n_threads, n_rows * n_cols):
        row, col = idx // n_cols, idx % n_cols
        axes[row][col].set_visible(False)

    fig.suptitle("LR: Execution Time vs Structure Size", fontsize=14)
    fig.tight_layout()
    out_path = plot_dir.joinpath("grid.png")
    fig.savefig(out_path)
    plt.close(fig)
    rprint(f"[green]Saved:[/green] {out_path}")


@app.command()
def threads(
    n_slices: int = typer.Option(
        20, "--n-slices", help="Number of slices per atom diameter"
    ),
):
    """Generate thread scaling plot."""
    setup_style()
    df = load_lr_data(n_slices)

    plot_dir = PLOTS_DIR.joinpath("thread_scaling")
    plot_dir.mkdir(parents=True, exist_ok=True)

    scaling = (
        df.group_by(["tool_label", "threads"])
        .agg(pl.col("time_ms").median().alias("median_time"))
        .sort(["tool_label", "threads"])
    )

    fig, ax = plt.subplots(figsize=(10, 6))
    tool_labels = [
        t for t in sorted(scaling["tool_label"].unique().to_list()) if "f32" not in t
    ]
    for tool_label in tool_labels:
        df_tool = scaling.filter(pl.col("tool_label") == tool_label)
        ax.plot(
            df_tool["threads"].to_list(),
            df_tool["median_time"].to_list(),
            marker="o",
            label=display_name(tool_label),
            color=COLORS.get(tool_label, "#95a5a6"),
            linestyle=LINESTYLES.get(tool_label, "-"),
            linewidth=2,
        )

    ax.set_xlabel("Thread Count")
    ax.set_ylabel("Median Execution Time (ms)")
    ax.set_title("LR: Thread Scaling")
    ax.legend()
    ax.grid(True, alpha=0.3)
    thread_values = sorted(scaling["threads"].unique().to_list())
    ax.set_xticks(thread_values)
    ax.set_xticklabels([str(t) for t in thread_values])

    fig.tight_layout()
    out_path = plot_dir.joinpath("lr.png")
    fig.savefig(out_path)
    plt.close(fig)
    rprint(f"[green]Saved:[/green] {out_path}")


@app.command()
def validation(
    n_slices: int = typer.Option(
        20, "--n-slices", help="Number of slices per atom diameter"
    ),
):
    """Generate SASA validation scatter plot (zsasa vs FreeSASA)."""
    setup_style()
    df = load_lr_data(n_slices)

    plot_dir = PLOTS_DIR.joinpath("validation")
    plot_dir.mkdir(parents=True, exist_ok=True)

    df_t1 = df.filter(pl.col("threads") == 1)
    if df_t1.height == 0:
        rprint("[yellow]No single-threaded LR data found[/yellow]")
        return

    pivot = (
        df_t1.select(["structure", "tool_label", "total_sasa"])
        .pivot(on="tool_label", index="structure", values="total_sasa")
        .drop_nulls()
    )

    zsasa_col = "zsasa_f64" if "zsasa_f64" in pivot.columns else "zsasa"
    if zsasa_col not in pivot.columns or "freesasa" not in pivot.columns:
        rprint("[yellow]Need both zsasa and freesasa data for validation[/yellow]")
        return

    zsasa_sasa = pivot[zsasa_col].to_list()
    fs_sasa = pivot["freesasa"].to_list()

    fig, ax = plt.subplots(figsize=(8, 8))
    ax.scatter(zsasa_sasa, fs_sasa, alpha=0.3, s=10, color="#3498db")

    max_val = max(max(zsasa_sasa), max(fs_sasa))
    ax.plot([0, max_val], [0, max_val], "r--", linewidth=1.5, label="y = x")

    zsasa_arr = np.array(zsasa_sasa)
    fs_arr = np.array(fs_sasa)
    correlation = np.corrcoef(zsasa_arr, fs_arr)[0, 1]
    r_squared = correlation**2

    rel_errors = np.abs(zsasa_arr - fs_arr) / fs_arr * 100
    max_rel_error = np.max(rel_errors)
    mean_rel_error = np.mean(rel_errors)

    zsasa_label = "zsasa (f64)" if zsasa_col == "zsasa_f64" else "zsasa"
    ax.set_xlabel(f"{zsasa_label} SASA")
    ax.set_ylabel("FreeSASA SASA")
    ax.set_title(f"LR: SASA Validation (n={len(zsasa_sasa):,})")
    ax.legend()

    stats_text = f"R² = {r_squared:.6f}\nMean error = {mean_rel_error:.4f}%\nMax error = {max_rel_error:.4f}%"
    ax.text(
        0.05,
        0.95,
        stats_text,
        transform=ax.transAxes,
        fontsize=10,
        verticalalignment="top",
        bbox=dict(boxstyle="round", facecolor="wheat", alpha=0.5),
    )

    ax.set_aspect("equal")
    fig.tight_layout()

    out_path = plot_dir.joinpath("lr.png")
    fig.savefig(out_path)
    plt.close(fig)
    rprint(f"[green]Saved:[/green] {out_path}")


@app.command(name="all")
def all_plots(
    n_slices: int = typer.Option(
        20, "--n-slices", help="Number of slices per atom diameter"
    ),
):
    """Generate all plots and summary."""
    summary(n_slices)
    rprint("\n[bold]Generating plots...[/bold]\n")
    validation(n_slices)
    scatter(n_slices)
    threads(n_slices)
    rprint(f"\n[bold green]All plots saved to:[/bold green] {PLOTS_DIR}")


if __name__ == "__main__":
    app()
