#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.11"
# dependencies = ["pandas>=2.0", "matplotlib>=3.7", "rich>=13.0", "typer>=0.9.0"]
# ///
"""Analyze benchmark results and generate reports/graphs.

Reads results from benchmarks/results/*/ directories.

Usage:
    ./benchmarks/scripts/analyze.py summary   # Summary tables
    ./benchmarks/scripts/analyze.py validate  # SASA validation
    ./benchmarks/scripts/analyze.py plot      # Generate graphs
    ./benchmarks/scripts/analyze.py all       # All of the above
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Annotated

import matplotlib.pyplot as plt
import pandas as pd
import typer
from rich.console import Console
from rich.table import Table

app = typer.Typer(help="Analyze benchmark results")
console = Console()

RESULTS_DIR = Path(__file__).parent.parent / "results"


def load_all_results() -> pd.DataFrame:
    """Load all benchmark results into a DataFrame."""
    all_data = []

    for result_dir in RESULTS_DIR.iterdir():
        if not result_dir.is_dir():
            continue

        csv_path = result_dir / "results.csv"
        if not csv_path.exists():
            continue

        df = pd.read_csv(csv_path)
        all_data.append(df)

    if not all_data:
        console.print("[yellow]No results found[/yellow]")
        raise typer.Exit(1)

    return pd.concat(all_data, ignore_index=True)


def load_configs() -> dict[str, dict]:
    """Load all config.json files."""
    configs = {}

    for result_dir in RESULTS_DIR.iterdir():
        if not result_dir.is_dir():
            continue

        config_path = result_dir / "config.json"
        if config_path.exists():
            configs[result_dir.name] = json.loads(config_path.read_text())

    return configs


@app.command()
def summary() -> None:
    """Show summary tables of benchmark results."""
    df = load_all_results()

    # Group by tool, structure, algorithm, threads and compute mean
    grouped = (
        df.groupby(["tool", "structure", "algorithm", "threads"])
        .agg({"sasa_time_ms": "mean", "n_atoms": "first", "total_sasa": "mean"})
        .reset_index()
    )

    # For each algorithm, show comparison at max threads
    for algo in grouped["algorithm"].unique():
        algo_df = grouped[grouped["algorithm"] == algo]
        max_threads = algo_df["threads"].max()
        max_df = algo_df[algo_df["threads"] == max_threads]

        table = Table(title=f"{algo.upper()} (threads={max_threads})")
        table.add_column("Structure")
        table.add_column("Atoms", justify="right")

        tools = sorted(max_df["tool"].unique())
        for tool in tools:
            table.add_column(f"{tool} (ms)", justify="right")

        # Sort by n_atoms
        structures = max_df.drop_duplicates("structure").sort_values("n_atoms")

        for _, row in structures.iterrows():
            struct = row["structure"]
            n_atoms = row["n_atoms"]
            values = [struct.upper(), f"{n_atoms:,}"]

            for tool in tools:
                tool_row = max_df[
                    (max_df["structure"] == struct) & (max_df["tool"] == tool)
                ]
                if len(tool_row) > 0:
                    ms = tool_row["sasa_time_ms"].values[0]
                    values.append(f"{ms:.2f}")
                else:
                    values.append("-")

            table.add_row(*values)

        console.print(table)
        console.print()


@app.command()
def validate(
    tolerance: Annotated[
        float,
        typer.Option("--tolerance", "-t", help="SASA difference tolerance (%)"),
    ] = 0.1,
) -> None:
    """Validate SASA values across tools."""
    df = load_all_results()

    # Group by structure, algorithm, threads
    grouped = (
        df.groupby(["structure", "algorithm", "threads", "tool"])
        .agg({"total_sasa": "mean"})
        .reset_index()
    )

    # Compare tools for each structure/algorithm/threads
    mismatches = []

    for (struct, algo, threads), group in grouped.groupby(
        ["structure", "algorithm", "threads"]
    ):
        tools = group["tool"].unique()
        if len(tools) < 2:
            continue

        sasa_values = {row["tool"]: row["total_sasa"] for _, row in group.iterrows()}

        # Compare all pairs
        tool_list = list(sasa_values.keys())
        for i, t1 in enumerate(tool_list):
            for t2 in tool_list[i + 1 :]:
                s1, s2 = sasa_values[t1], sasa_values[t2]
                if s1 == 0 or s2 == 0:
                    continue

                diff_pct = abs(s1 - s2) / max(s1, s2) * 100

                if diff_pct > tolerance:
                    mismatches.append(
                        {
                            "structure": struct,
                            "algorithm": algo,
                            "threads": threads,
                            "tool1": t1,
                            "tool2": t2,
                            "sasa1": s1,
                            "sasa2": s2,
                            "diff_pct": diff_pct,
                        }
                    )

    if mismatches:
        console.print(
            f"[yellow]Found {len(mismatches)} mismatches (>{tolerance}%)[/yellow]\n"
        )

        table = Table(title="SASA Mismatches")
        table.add_column("Structure")
        table.add_column("Algo")
        table.add_column("Threads", justify="right")
        table.add_column("Tool 1")
        table.add_column("SASA 1", justify="right")
        table.add_column("Tool 2")
        table.add_column("SASA 2", justify="right")
        table.add_column("Diff %", justify="right")

        for m in mismatches[:20]:  # Show first 20
            table.add_row(
                m["structure"],
                m["algorithm"],
                str(m["threads"]),
                m["tool1"],
                f"{m['sasa1']:.2f}",
                m["tool2"],
                f"{m['sasa2']:.2f}",
                f"{m['diff_pct']:.3f}",
            )

        console.print(table)
    else:
        console.print(
            f"[green]All SASA values match within {tolerance}% tolerance[/green]"
        )


@app.command()
def plot(
    output_dir: Annotated[
        Path | None,
        typer.Option("--output", "-o", help="Output directory for plots"),
    ] = None,
) -> None:
    """Generate benchmark plots."""
    df = load_all_results()

    if output_dir is None:
        output_dir = RESULTS_DIR / "plots"
    output_dir.mkdir(parents=True, exist_ok=True)

    # Group and average
    grouped = (
        df.groupby(["tool", "structure", "algorithm", "threads", "n_atoms"])
        .agg({"sasa_time_ms": "mean"})
        .reset_index()
    )

    # Plot 1: Execution time vs structure size (log-log) at max threads
    max_threads = grouped["threads"].max()
    max_df = grouped[grouped["threads"] == max_threads]

    for algo in max_df["algorithm"].unique():
        algo_df = max_df[max_df["algorithm"] == algo]

        fig, ax = plt.subplots(figsize=(10, 6))

        for tool in sorted(algo_df["tool"].unique()):
            tool_df = algo_df[algo_df["tool"] == tool].sort_values("n_atoms")
            ax.loglog(
                tool_df["n_atoms"],
                tool_df["sasa_time_ms"],
                "o-",
                label=tool,
                markersize=8,
            )

        ax.set_xlabel("Number of Atoms")
        ax.set_ylabel("Time (ms)")
        ax.set_title(f"{algo.upper()} - Execution Time vs Structure Size")
        ax.legend()
        ax.grid(True, alpha=0.3)

        path = output_dir / f"time_vs_size_{algo}.png"
        fig.savefig(path, dpi=150, bbox_inches="tight")
        plt.close(fig)
        console.print(f"Saved: {path}")

    # Plot 2: Thread scaling for largest structure
    largest = grouped["n_atoms"].max()
    large_df = grouped[grouped["n_atoms"] == largest]

    for algo in large_df["algorithm"].unique():
        algo_df = large_df[large_df["algorithm"] == algo]

        fig, ax = plt.subplots(figsize=(10, 6))

        for tool in sorted(algo_df["tool"].unique()):
            tool_df = algo_df[algo_df["tool"] == tool].sort_values("threads")
            ax.plot(
                tool_df["threads"],
                tool_df["sasa_time_ms"],
                "o-",
                label=tool,
                markersize=8,
            )

        ax.set_xlabel("Threads")
        ax.set_ylabel("Time (ms)")
        struct_name = large_df["structure"].iloc[0].upper()
        ax.set_title(f"{algo.upper()} - Thread Scaling ({struct_name})")
        ax.legend()
        ax.grid(True, alpha=0.3)

        path = output_dir / f"thread_scaling_{algo}.png"
        fig.savefig(path, dpi=150, bbox_inches="tight")
        plt.close(fig)
        console.print(f"Saved: {path}")

    # Plot 3: Speedup ratio (Zig vs others)
    if "zig" in grouped["tool"].unique():
        for algo in grouped["algorithm"].unique():
            algo_df = max_df[max_df["algorithm"] == algo]

            zig_times = algo_df[algo_df["tool"] == "zig"].set_index("structure")[
                "sasa_time_ms"
            ]

            fig, ax = plt.subplots(figsize=(10, 6))

            other_tools = [t for t in algo_df["tool"].unique() if t != "zig"]
            x = range(len(zig_times))
            width = 0.8 / len(other_tools)

            for i, tool in enumerate(other_tools):
                tool_times = algo_df[algo_df["tool"] == tool].set_index("structure")[
                    "sasa_time_ms"
                ]
                speedups = tool_times / zig_times
                speedups = speedups.reindex(zig_times.index)

                ax.bar(
                    [xi + i * width for xi in x],
                    speedups.values,
                    width,
                    label=f"vs {tool}",
                )

            ax.set_xlabel("Structure")
            ax.set_ylabel("Speedup (X times faster)")
            ax.set_title(f"{algo.upper()} - Zig Speedup Ratio")
            ax.set_xticks([xi + width * (len(other_tools) - 1) / 2 for xi in x])
            ax.set_xticklabels([s.upper() for s in zig_times.index])
            ax.axhline(y=1, color="gray", linestyle="--", alpha=0.5)
            ax.legend()
            ax.grid(True, alpha=0.3, axis="y")

            path = output_dir / f"speedup_{algo}.png"
            fig.savefig(path, dpi=150, bbox_inches="tight")
            plt.close(fig)
            console.print(f"Saved: {path}")

    console.print(f"\n[green]Plots saved to: {output_dir}[/green]")


@app.command(name="all")
def run_all() -> None:
    """Run all analyses."""
    console.print("[bold]== Summary ==[/bold]\n")
    summary()

    console.print("\n[bold]== Validation ==[/bold]\n")
    validate()

    console.print("\n[bold]== Plots ==[/bold]\n")
    plot()


if __name__ == "__main__":
    app()
