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
"""SR benchmark analysis CLI.

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
    RESULTS_BASE,
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

app = typer.Typer(help="SR benchmark analysis CLI")


# === CLI Commands ===


@app.command()
def summary(
    n_points: int = typer.Option(
        100, "--n-points", "-N", help="Number of sphere test points"
    ),
):
    """Print summary statistics table."""
    df = load_data(n_points)

    # Single-threaded comparison
    df_t1 = df.filter(pl.col("threads") == 1)

    table = Table(title="Single-Threaded Performance Summary (SR)")
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
            f"  SR: zsasa vs FreeSASA = [green]{med:.2f}x[/green] (median), {mean:.2f}x (mean)"
        )

    if "zsasa_f32" in pivot.columns and "zsasa_f64" in pivot.columns:
        speedup = pivot.select(
            (pl.col("zsasa_f64") / pl.col("zsasa_f32")).alias("speedup")
        )
        med = speedup["speedup"].median()
        rprint(
            f"  SR: zsasa(f32) vs zsasa(f64) = [green]{med:.2f}x[/green] (median)"
        )

    if "rustsasa" in pivot.columns and zsasa_col in pivot.columns:
        speedup_rust = pivot.select(
            (pl.col("rustsasa") / pl.col(zsasa_col)).alias("speedup")
        )
        med = speedup_rust["speedup"].median()
        rprint(
            f"  SR: zsasa vs RustSASA = [green]{med:.2f}x[/green] (median)"
        )

    # Speedup by size bin table
    rprint("\n")
    speedup_by_bin = compute_speedup_by_bin(df_t1, threads=1)

    bin_table = Table(title="SR Speedup by Structure Size (threads=1)")
    bin_table.add_column("Size Bin", style="cyan")
    bin_table.add_column("Count", justify="right")
    bin_table.add_column("zsasa vs FreeSASA", justify="right")
    bin_table.add_column("zsasa vs RustSASA", justify="right")
    bin_table.add_column("f32 vs f64", justify="right")

    bin_order = [b[2] for b in BINS]
    data_dict = {row["size_bin"]: row for row in speedup_by_bin.iter_rows(named=True)}

    for bin_name in bin_order:
        if bin_name not in data_dict:
            continue
        row = data_dict[bin_name]
        zig_fs = row.get("zsasa_f64_vs_freesasa")
        zig_rust = row.get("zsasa_f64_vs_rustsasa")
        f32_f64 = row.get("zsasa_f32_vs_zsasa_f64")

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
def export_csv(
    n_points: int = typer.Option(
        100, "--n-points", "-N", help="Number of sphere test points"
    ),
):
    """Export summary tables as CSV files per thread count."""
    df = load_data(n_points)
    df = add_size_bin(df)

    csv_dir = RESULTS_BASE.joinpath(str(n_points), "csv")
    csv_dir.mkdir(parents=True, exist_ok=True)

    thread_counts = sorted(df["threads"].unique().to_list())

    for threads in thread_counts:
        t_dir = csv_dir.joinpath(f"t{threads}")
        t_dir.mkdir(parents=True, exist_ok=True)

        df_t = df.filter(pl.col("threads") == threads)

        # Performance summary
        perf_stats = (
            df_t.group_by(["tool_label", "precision"])
            .agg(
                pl.len().alias("structures"),
                pl.col("time_ms").median().alias("median_ms"),
                pl.col("time_ms").mean().alias("mean_ms"),
                pl.col("time_ms").std().alias("std_ms"),
                pl.col("time_ms").quantile(0.95).alias("p95_ms"),
            )
            .sort("tool_label")
        )
        perf_stats.write_csv(t_dir.joinpath("performance_summary.csv"))

        # Speedup by size bin
        if df_t.height > 0:
            speedup_by_bin = compute_speedup_by_bin(df_t, threads=threads)
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

        rprint(f"[green]Saved:[/green] {t_dir}/")

    rprint(f"\n[bold]Exported CSV files to {csv_dir}[/bold]")


# === Plot Command Wrappers ===


@app.command()
def scatter(
    n_points: int = typer.Option(
        100, "--n-points", "-N", help="Number of sphere test points"
    ),
):
    """Generate atoms vs time scatter plot."""
    plot_scatter(n_points)


@app.command()
def threads(
    n_points: int = typer.Option(
        100, "--n-points", "-N", help="Number of sphere test points"
    ),
):
    """Generate thread scaling plot."""
    plot_threads(n_points)


@app.command()
def grid(
    n_points: int = typer.Option(
        100, "--n-points", "-N", help="Number of sphere test points"
    ),
):
    """Generate grid of speedup plots for all thread counts."""
    plot_grid(n_points)


@app.command()
def validation(
    n_points: int = typer.Option(
        100, "--n-points", "-N", help="Number of sphere test points"
    ),
):
    """Generate SASA validation scatter plot (zsasa vs FreeSASA)."""
    plot_validation(n_points)


@app.command()
def samples(
    n_points: int = typer.Option(
        100, "--n-points", "-N", help="Number of sphere test points"
    ),
):
    """Generate thread scaling plots for representative structures per size bin."""
    plot_samples(n_points)


@app.command()
def large(
    n_points: int = typer.Option(
        100, "--n-points", "-N", help="Number of sphere test points"
    ),
):
    """Generate speedup bar chart for large structures (50k+ atoms)."""
    plot_large(n_points)


@app.command()
def efficiency(
    n_points: int = typer.Option(
        100, "--n-points", "-N", help="Number of sphere test points"
    ),
):
    """Calculate and plot parallel efficiency."""
    plot_efficiency(n_points)


@app.command()
def speedup(
    min_atoms: int = typer.Option(50000, help="Minimum atom count for filtering"),
    top_n: int = typer.Option(5, help="Number of top entries to show"),
    n_points: int = typer.Option(
        100, "--n-points", "-N", help="Number of sphere test points"
    ),
):
    """Find structures with best zsasa speedup at any thread count."""
    plot_speedup(min_atoms=min_atoms, top_n=top_n, n_points=n_points)


@app.command()
def all(
    n_points: int = typer.Option(
        100, "--n-points", "-N", help="Number of sphere test points"
    ),
):
    """Generate all plots and summary."""
    summary(n_points)
    rprint("\n[bold]Generating plots...[/bold]\n")
    plot_validation(n_points)
    plot_scatter(n_points)
    plot_threads(n_points)
    plot_grid(n_points)
    plot_samples(n_points)
    plot_large(n_points)
    plot_efficiency(n_points)
    plot_speedup(min_atoms=50000, top_n=5, n_points=n_points)
    rprint(f"\n[bold green]All plots saved to:[/bold green] {PLOTS_DIR}")


if __name__ == "__main__":
    app()
