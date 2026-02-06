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
"""Benchmark analysis CLI.

Generates summary statistics and plots from bench.py results.

Modules:
- analyze_data.py: Constants, data loading, utilities
- analyze_plots.py: All plotting functions

Usage:
    ./benchmarks/scripts/analyze.py summary      # Summary tables
    ./benchmarks/scripts/analyze.py all          # All plots + summary
    ./benchmarks/scripts/analyze.py scatter      # Atoms vs time
    ./benchmarks/scripts/analyze.py threads      # Thread scaling
    ./benchmarks/scripts/analyze.py grid         # Speedup by size/threads
    ./benchmarks/scripts/analyze.py validation   # SASA validation
    ./benchmarks/scripts/analyze.py samples      # Per-bin sample plots
    ./benchmarks/scripts/analyze.py large        # Large structure analysis
    ./benchmarks/scripts/analyze.py efficiency   # Parallel efficiency
    ./benchmarks/scripts/analyze.py speedup      # Best speedup structures
    ./benchmarks/scripts/analyze.py export-csv   # Export to CSV
"""

from pathlib import Path

import polars as pl
import typer
from rich import print as rprint
from rich.table import Table

from analyze_data import (
    BINS,
    PLOTS_DIR,
    RESULTS_DIR,
    add_size_bin,
    compute_speedup_by_bin,
    load_data,
)
from analyze_plots import (
    plot_efficiency,
    plot_grid,
    plot_large,
    plot_samples,
    plot_scatter,
    plot_speedup,
    plot_threads,
    plot_validation,
)

app = typer.Typer(help="Benchmark analysis CLI")


# === CLI Commands ===


@app.command()
def summary():
    """Print summary statistics table."""
    df = load_data()

    # Single-threaded comparison
    df_t1 = df.filter(pl.col("threads") == 1)

    table = Table(title="Single-Threaded Performance Summary")
    table.add_column("Algorithm", style="cyan")
    table.add_column("Tool", style="green")
    table.add_column("Structures", justify="right")
    table.add_column("Median (ms)", justify="right")
    table.add_column("Mean (ms)", justify="right")
    table.add_column("Std (ms)", justify="right")
    table.add_column("P95 (ms)", justify="right")

    stats = (
        df_t1.group_by(["algorithm", "tool_label"])
        .agg(
            pl.len().alias("n"),
            pl.col("time_ms").median().alias("median"),
            pl.col("time_ms").mean().alias("mean"),
            pl.col("time_ms").std().alias("std"),
            pl.col("time_ms").quantile(0.95).alias("p95"),
        )
        .sort(["algorithm", "tool_label"])
    )

    for row in stats.iter_rows(named=True):
        table.add_row(
            row["algorithm"].upper(),
            row["tool_label"],
            f"{row['n']:,}",
            f"{row['median']:.2f}",
            f"{row['mean']:.2f}",
            f"{row['std']:.2f}",
            f"{row['p95']:.2f}",
        )

    rprint(table)

    # Speedup summary
    rprint("\n[bold]Speedup Ratios (Zig vs others):[/bold]")
    for algo in ["sr", "lr"]:
        df_algo = df_t1.filter(pl.col("algorithm") == algo)
        if df_algo.height == 0:
            continue

        pivot = (
            df_algo.select(["structure", "tool_label", "time_ms"])
            .pivot(on="tool_label", index="structure", values="time_ms")
            .drop_nulls()
        )

        zig_col = "zig_f64" if "zig_f64" in pivot.columns else "zig"

        if zig_col in pivot.columns and "freesasa" in pivot.columns:
            speedup = pivot.select(
                (pl.col("freesasa") / pl.col(zig_col)).alias("speedup")
            )
            med = speedup["speedup"].median()
            mean = speedup["speedup"].mean()
            rprint(
                f"  {algo.upper()}: Zig vs FreeSASA = [green]{med:.2f}x[/green] (median), {mean:.2f}x (mean)"
            )

        if "zig_f32" in pivot.columns and "zig_f64" in pivot.columns:
            speedup = pivot.select(
                (pl.col("zig_f64") / pl.col("zig_f32")).alias("speedup")
            )
            med = speedup["speedup"].median()
            rprint(
                f"  {algo.upper()}: Zig(f32) vs Zig(f64) = [green]{med:.2f}x[/green] (median)"
            )

        if algo == "sr" and "rust" in pivot.columns and zig_col in pivot.columns:
            speedup_rust = pivot.select(
                (pl.col("rust") / pl.col(zig_col)).alias("speedup")
            )
            med = speedup_rust["speedup"].median()
            rprint(
                f"  {algo.upper()}: Zig vs Rust = [green]{med:.2f}x[/green] (median)"
            )

    # Speedup by size bin table (SR only)
    rprint("\n")
    df_sr = df_t1.filter(pl.col("algorithm") == "sr")
    if df_sr.height == 0:
        return

    speedup_by_bin = compute_speedup_by_bin(df_sr, threads=1)

    bin_table = Table(title="SR Speedup by Structure Size (threads=1)")
    bin_table.add_column("Size Bin", style="cyan")
    bin_table.add_column("Count", justify="right")
    bin_table.add_column("Zig vs FreeSASA", justify="right")
    bin_table.add_column("Zig vs Rust", justify="right")
    bin_table.add_column("f32 vs f64", justify="right")

    bin_order = [b[2] for b in BINS]
    data_dict = {row["size_bin"]: row for row in speedup_by_bin.iter_rows(named=True)}

    for bin_name in bin_order:
        if bin_name not in data_dict:
            continue
        row = data_dict[bin_name]
        zig_fs = row.get("zig_f64_vs_freesasa")
        zig_rust = row.get("zig_f64_vs_rust")
        f32_f64 = row.get("zig_f32_vs_zig_f64")

        fs_str = f"{zig_fs:.2f}x" if zig_fs else "-"
        rust_str = f"{zig_rust:.2f}x" if zig_rust else "-"
        f32_str = f"{f32_f64:.2f}x" if f32_f64 else "-"

        if zig_fs and zig_fs > 1.0:
            fs_str = f"[green]{fs_str}[/green]"
        elif zig_fs:
            fs_str = f"[red]{fs_str}[/red]"

        if zig_rust and zig_rust > 1.0:
            rust_str = f"[green]{rust_str}[/green]"
        elif zig_rust:
            rust_str = f"[red]{rust_str}[/red]"

        if f32_f64 and f32_f64 > 1.0:
            f32_str = f"[green]{f32_str}[/green]"
        elif f32_f64:
            f32_str = f"[red]{f32_str}[/red]"

        bin_table.add_row(bin_name, f"{row['count']:,}", fs_str, rust_str, f32_str)

    rprint(bin_table)
    rprint("\n[dim]Green = faster, Red = slower[/dim]")


@app.command(name="export-csv")
def export_csv():
    """Export summary tables as CSV files per thread count."""
    df = load_data()
    df = add_size_bin(df)

    csv_dir = RESULTS_DIR.joinpath("csv")
    csv_dir.mkdir(parents=True, exist_ok=True)

    thread_counts = sorted(df["threads"].unique().to_list())

    for threads in thread_counts:
        t_dir = csv_dir.joinpath(f"t{threads}")
        t_dir.mkdir(parents=True, exist_ok=True)

        df_t = df.filter(pl.col("threads") == threads)

        # Performance summary
        perf_stats = (
            df_t.group_by(["algorithm", "tool_label", "precision"])
            .agg(
                pl.len().alias("structures"),
                pl.col("time_ms").median().alias("median_ms"),
                pl.col("time_ms").mean().alias("mean_ms"),
                pl.col("time_ms").std().alias("std_ms"),
                pl.col("time_ms").quantile(0.95).alias("p95_ms"),
            )
            .sort(["algorithm", "tool_label"])
        )
        perf_stats.write_csv(t_dir.joinpath("performance_summary.csv"))

        # Speedup by size bin (SR)
        df_sr = df_t.filter(pl.col("algorithm") == "sr")
        if df_sr.height > 0:
            speedup_by_bin = compute_speedup_by_bin(df_sr, threads=threads)
            bin_order = [b[2] for b in BINS]
            rows = []
            data_dict = {
                row["size_bin"]: row for row in speedup_by_bin.iter_rows(named=True)
            }
            for bin_name in bin_order:
                if bin_name in data_dict:
                    rows.append(data_dict[bin_name])

            if rows:
                speedup_df = pl.DataFrame(rows)
                speedup_df.write_csv(t_dir.joinpath("speedup_by_bin_sr.csv"))

        # Speedup by size bin (LR)
        df_lr = df_t.filter(pl.col("algorithm") == "lr")
        if df_lr.height > 0:
            speedup_by_bin_lr = compute_speedup_by_bin(df_lr, threads=threads)
            rows_lr = []
            data_dict_lr = {
                row["size_bin"]: row for row in speedup_by_bin_lr.iter_rows(named=True)
            }
            for bin_name in bin_order:
                if bin_name in data_dict_lr:
                    rows_lr.append(data_dict_lr[bin_name])

            if rows_lr:
                speedup_df_lr = pl.DataFrame(rows_lr)
                speedup_df_lr.write_csv(t_dir.joinpath("speedup_by_bin_lr.csv"))

        rprint(f"[green]Saved:[/green] {t_dir}/")

    rprint(f"\n[bold]Exported CSV files to {csv_dir}[/bold]")


# === Plot Command Wrappers ===


@app.command()
def scatter():
    """Generate atoms vs time scatter plot."""
    plot_scatter()


@app.command()
def threads():
    """Generate thread scaling plot."""
    plot_threads()


@app.command()
def grid():
    """Generate grid of speedup plots for all thread counts."""
    plot_grid()


@app.command()
def validation():
    """Generate SASA validation scatter plot (Zig vs FreeSASA)."""
    plot_validation()


@app.command()
def samples():
    """Generate thread scaling plots for representative structures per size bin."""
    plot_samples()


@app.command()
def large():
    """Generate speedup bar chart for large structures (50k+ atoms)."""
    plot_large()


@app.command()
def efficiency():
    """Calculate and plot parallel efficiency."""
    plot_efficiency()


@app.command()
def speedup(
    min_atoms: int = typer.Option(50000, help="Minimum atom count for filtering"),
    top_n: int = typer.Option(5, help="Number of top entries to show"),
):
    """Find structures with best Zig speedup at any thread count."""
    plot_speedup(min_atoms=min_atoms, top_n=top_n)


@app.command()
def all():
    """Generate all plots and summary."""
    summary()
    rprint("\n[bold]Generating plots...[/bold]\n")
    plot_validation()
    plot_scatter()
    plot_threads()
    plot_grid()
    plot_samples()
    plot_large()
    plot_efficiency()
    plot_speedup(min_atoms=50000, top_n=5)
    rprint(f"\n[bold green]All plots saved to:[/bold green] {PLOTS_DIR}")


if __name__ == "__main__":
    app()
