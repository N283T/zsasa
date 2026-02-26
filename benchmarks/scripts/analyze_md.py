#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.12"
# dependencies = [
#     "matplotlib>=3.8",
#     "typer>=0.9.0",
#     "rich>=13.0",
# ]
# ///
"""MD trajectory benchmark analysis CLI.

Generates summary statistics and plots from bench_md.py results.

Usage:
    ./benchmarks/scripts/analyze_md.py summary --name 6sup_R1
    ./benchmarks/scripts/analyze_md.py bar --name 6sup_R1
    ./benchmarks/scripts/analyze_md.py memory --name 6sup_R1
    ./benchmarks/scripts/analyze_md.py all --name 6sup_R1
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Annotated

import matplotlib.pyplot as plt
import typer
from rich import print as rprint
from rich.table import Table

app = typer.Typer(help="MD trajectory benchmark analysis CLI")

# === Constants ===

RESULTS_BASE = Path(__file__).parent.parent.joinpath("results", "md")


def _results_dir(n_points: int) -> Path:
    return RESULTS_BASE.joinpath(str(n_points))


DISPLAY_NAMES: dict[str, str] = {
    "zig_f32": "zsasa CLI (f32)",
    "zig_f64": "zsasa CLI (f64)",
    "zig_f32_bitmask": "zsasa CLI (f32, bitmask)",
    "zig_f64_bitmask": "zsasa CLI (f64, bitmask)",
    "zsasa_mdtraj": "zsasa.mdtraj",
    "zsasa_mdanalysis": "zsasa.mdanalysis",
    "mdtraj": "MDTraj native",
    "mdsasa_bolt": "mdsasa-bolt",
}

COLORS: dict[str, str] = {
    "zig_f32": "#e67e22",
    "zig_f64": "#f39c12",
    "zig_f32_bitmask": "#d35400",
    "zig_f64_bitmask": "#e08e0b",
    "zsasa_mdtraj": "#27ae60",
    "zsasa_mdanalysis": "#16a085",
    "mdtraj": "#3498db",
    "mdsasa_bolt": "#e74c3c",
}

# Ordered for display (fastest tools first in typical results)
TOOL_ORDER = [
    "zsasa_mdtraj",
    "zsasa_mdanalysis",
    "zig_f32",
    "zig_f64",
    "zig_f32_bitmask",
    "zig_f64_bitmask",
    "mdtraj",
    "mdsasa_bolt",
]


# === Data Model ===


@dataclass(frozen=True)
class BenchResult:
    """Single benchmark result parsed from hyperfine JSON."""

    tool: str
    threads: int | None
    mean: float
    stddev: float | None
    median: float
    min_time: float
    max_time: float
    memory_bytes: int


def _parse_filename(filename: str) -> tuple[str, int | None]:
    """Parse tool name and thread count from benchmark filename.

    Examples:
        bench_zig_f32_1t.json -> ("zig_f32", 1)
        bench_mdsasa_bolt_all.json -> ("mdsasa_bolt", None)
        bench_mdtraj_1t.json -> ("mdtraj", 1)
        bench_zsasa_mdtraj_8t.json -> ("zsasa_mdtraj", 8)
    """
    stem = filename.replace("bench_", "").replace(".json", "")

    # Split off thread suffix
    parts = stem.rsplit("_", 1)
    if len(parts) == 2 and parts[1].endswith("t") and parts[1][:-1].isdigit():
        return parts[0], int(parts[1][:-1])
    if len(parts) == 2 and parts[1] == "all":
        return parts[0], None
    return stem, None


def load_results(name: str, n_points: int = 100) -> list[BenchResult]:
    """Load all benchmark results for a given dataset name."""
    results_dir = _results_dir(n_points).joinpath(name)
    if not results_dir.exists():
        base = _results_dir(n_points)
        available = (
            [d.name for d in base.iterdir() if d.is_dir()] if base.exists() else []
        )
        raise typer.BadParameter(
            f"Dataset '{name}' not found in {base}. Available: {', '.join(sorted(available))}"
        )

    # Resolve "all" threads from config cpu_cores
    config = load_config(name, n_points)
    cpu_cores = (
        config.get("system", {}).get("cpu_cores") if config is not None else None
    )

    results: list[BenchResult] = []
    for json_file in sorted(results_dir.glob("bench_*.json")):
        try:
            data = json.loads(json_file.read_text())
            if not data.get("results"):
                continue
            r = data["results"][0]
            tool, threads = _parse_filename(json_file.name)

            # "all" threads -> cpu_cores from config
            if threads is None and cpu_cores is not None:
                threads = cpu_cores

            mem_list = r.get("memory_usage_byte", [0])
            peak_mem = max(mem_list) if mem_list else 0

            results.append(
                BenchResult(
                    tool=tool,
                    threads=threads,
                    mean=r["mean"],
                    stddev=r.get("stddev"),
                    median=r["median"],
                    min_time=r["min"],
                    max_time=r["max"],
                    memory_bytes=peak_mem,
                )
            )
        except (json.JSONDecodeError, KeyError) as e:
            rprint(f"[yellow]Warning: Skipping {json_file.name}: {e}[/yellow]")
            continue

    return results


def load_config(name: str, n_points: int = 100) -> dict | None:
    """Load config.json for dataset metadata."""
    config_path = _results_dir(n_points).joinpath(name, "config.json")
    if not config_path.exists():
        return None
    try:
        return json.loads(config_path.read_text())
    except json.JSONDecodeError:
        rprint(f"[yellow]Warning: Could not parse {config_path}[/yellow]")
        return None


def get_plots_dir(name: str, n_points: int = 100) -> Path:
    """Get (and create) plots output directory."""
    plots_dir = _results_dir(n_points).joinpath(name, "plots")
    plots_dir.mkdir(parents=True, exist_ok=True)
    return plots_dir


def display_name(tool: str) -> str:
    """Get display name for a tool."""
    return DISPLAY_NAMES.get(tool, tool)


def get_best_per_tool(results: list[BenchResult]) -> dict[str, BenchResult]:
    """Get the result at the highest thread count per tool.

    Picks the result with the most threads (typically 10t) to give a
    consistent comparison point.  For single-threaded tools like MDTraj,
    the only available result (1t) is returned.
    """
    best: dict[str, BenchResult] = {}
    for r in results:
        threads = r.threads if r.threads is not None else 0
        prev = best.get(r.tool)
        prev_threads = (
            prev.threads if prev is not None and prev.threads is not None else 0
        )
        if prev is None or threads > prev_threads:
            best[r.tool] = r
    return best


def get_mdtraj_baseline(results: list[BenchResult]) -> float | None:
    """Get MDTraj native 1t mean time as baseline."""
    for r in results:
        if r.tool == "mdtraj":
            return r.mean
    return None


def setup_style() -> None:
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


def _format_memory(memory_bytes: int) -> str:
    """Format memory size with appropriate unit."""
    mem_gb = memory_bytes / (1024**3)
    if mem_gb >= 1.0:
        return f"{mem_gb:.1f} GB"
    mem_mb = memory_bytes / (1024**2)
    return f"{mem_mb:.0f} MB"


def _threads_label(threads: int | None) -> str:
    """Format thread count for display."""
    if threads is None:
        return "all"
    return str(threads)


def build_subtitle(config: dict | None) -> str:
    """Build subtitle string from config metadata (atom_count, total_frames, stride)."""
    if not config:
        return ""
    params = config.get("parameters", {})
    parts: list[str] = []
    if "atom_count" in params:
        parts.append(f"{params['atom_count']:,} atoms")
    if "total_frames" in params:
        parts.append(f"{params['total_frames']:,} frames")
    stride = params.get("stride", 1)
    if stride and stride > 1:
        parts.append(f"stride={stride}")
    return ", ".join(parts)


# === CLI Commands ===


@app.command()
def summary(
    name: Annotated[str, typer.Option("--name", "-n", help="Dataset name")],
    n_points: Annotated[
        int,
        typer.Option("--n-points", "-N", help="Number of sphere test points"),
    ] = 100,
) -> None:
    """Print summary table of all benchmark results."""
    results = load_results(name, n_points)
    config = load_config(name, n_points)

    if config:
        params = config.get("parameters", {})
        system = config.get("system", {})
        rprint(f"[bold]Dataset:[/bold] {name}")
        rprint(
            f"[dim]System: {system.get('cpu_model', 'N/A')} "
            f"({system.get('cpu_cores', '?')} cores, "
            f"{system.get('memory_gb', '?')} GB RAM)[/dim]"
        )
        traj_parts: list[str] = []
        if "atom_count" in params:
            traj_parts.append(f"{params['atom_count']:,} atoms")
        if "total_frames" in params:
            traj_parts.append(f"{params['total_frames']:,} frames")
        if traj_parts:
            rprint(f"[dim]Trajectory: {', '.join(traj_parts)}[/dim]")
        rprint(
            f"[dim]Params: stride={params.get('stride', '?')}, "
            f"n_points={params.get('n_points', '?')}, "
            f"runs={params.get('runs', '?')}[/dim]"
        )
        rprint()

    mdtraj_baseline = get_mdtraj_baseline(results)

    subtitle = build_subtitle(config)
    table_title = f"MD Benchmark Results: {name}"
    if subtitle:
        table_title += f" ({subtitle})"
    table = Table(title=table_title)
    table.add_column("Tool", style="cyan")
    table.add_column("Threads", justify="right")
    table.add_column("Mean (s)", justify="right")
    table.add_column("Stddev", justify="right")
    table.add_column("Min (s)", justify="right")
    table.add_column("Max (s)", justify="right")
    table.add_column("Memory", justify="right")
    if mdtraj_baseline is not None:
        table.add_column("vs MDTraj", justify="right")

    sorted_results = sorted(results, key=lambda r: r.mean)

    for r in sorted_results:
        stddev_str = f"{r.stddev:.3f}" if r.stddev is not None else "-"

        row = [
            display_name(r.tool),
            _threads_label(r.threads),
            f"{r.mean:.3f}",
            stddev_str,
            f"{r.min_time:.3f}",
            f"{r.max_time:.3f}",
            _format_memory(r.memory_bytes),
        ]

        if mdtraj_baseline is not None and r.mean > 0:
            speedup = mdtraj_baseline / r.mean
            speedup_str = f"{speedup:.1f}x"
            if speedup > 1.0:
                speedup_str = f"[green]{speedup_str}[/green]"
            elif speedup < 1.0:
                speedup_str = f"[red]{speedup_str}[/red]"
            row.append(speedup_str)

        table.add_row(*row)

    rprint(table)


@app.command()
def bar(
    name: Annotated[str, typer.Option("--name", "-n", help="Dataset name")],
    n_points: Annotated[
        int,
        typer.Option("--n-points", "-N", help="Number of sphere test points"),
    ] = 100,
) -> None:
    """Generate vertical bar chart comparing best time per tool."""
    setup_style()
    results = load_results(name, n_points)
    config = load_config(name, n_points)
    best = get_best_per_tool(results)

    ordered = [t for t in TOOL_ORDER if t in best]
    if not ordered:
        rprint("[yellow]No matching tools found for bar chart.[/yellow]")
        return
    ordered.sort(key=lambda t: best[t].mean)

    labels = []
    times = []
    colors = []
    for tool in ordered:
        r = best[tool]
        thread_info = f" ({r.threads}t)" if r.threads is not None else " (all)"
        labels.append(f"{display_name(tool)}{thread_info}")
        times.append(r.mean)
        colors.append(COLORS.get(tool, "#95a5a6"))

    fig, ax = plt.subplots(figsize=(max(6, len(labels) * 1.8), 6))
    x = range(len(labels))
    bars = ax.bar(x, times, color=colors, width=0.6, edgecolor="white")

    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=30, ha="right")
    ax.set_ylabel("Time (seconds)")
    subtitle = build_subtitle(config)
    title = f"MD Trajectory SASA Benchmark: {name}"
    if subtitle:
        title += f"\n({subtitle})"
    ax.set_title(title)
    ax.grid(True, alpha=0.3, axis="y")

    for bar_item, t in zip(bars, times):
        ax.text(
            bar_item.get_x() + bar_item.get_width() / 2,
            bar_item.get_height() + max(times) * 0.02,
            f"{t:.1f}s",
            ha="center",
            va="bottom",
            fontsize=10,
            fontweight="bold",
        )

    ax.set_ylim(0, max(times) * 1.15)

    out_path = get_plots_dir(name, n_points).joinpath("bar.png")
    fig.savefig(out_path)
    plt.close(fig)
    rprint(f"[green]Saved:[/green] {out_path}")


@app.command()
def memory(
    name: Annotated[str, typer.Option("--name", "-n", help="Dataset name")],
    n_points: Annotated[
        int,
        typer.Option("--n-points", "-N", help="Number of sphere test points"),
    ] = 100,
) -> None:
    """Generate memory usage comparison bar chart."""
    setup_style()
    config = load_config(name, n_points)
    results = load_results(name, n_points)
    best = get_best_per_tool(results)

    ordered = [t for t in TOOL_ORDER if t in best]
    if not ordered:
        rprint("[yellow]No matching tools found for memory chart.[/yellow]")
        return
    ordered.sort(key=lambda t: best[t].memory_bytes)

    labels = []
    mem_mbs = []
    colors = []
    for tool in ordered:
        r = best[tool]
        labels.append(display_name(tool))
        mem_mbs.append(r.memory_bytes / (1024**2))
        colors.append(COLORS.get(tool, "#95a5a6"))

    fig, ax = plt.subplots(figsize=(max(6, len(labels) * 1.8), 6))
    x = range(len(labels))
    bars = ax.bar(x, mem_mbs, color=colors, width=0.6, edgecolor="white")

    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=30, ha="right")
    ax.set_ylabel("Peak Memory (MB)")
    subtitle = build_subtitle(config)
    title = f"Memory Usage: {name}"
    if subtitle:
        title += f"\n({subtitle})"
    ax.set_title(title)
    ax.grid(True, alpha=0.3, axis="y")

    for bar_item, mem_mb in zip(bars, mem_mbs):
        ax.text(
            bar_item.get_x() + bar_item.get_width() / 2,
            bar_item.get_height() + max(mem_mbs) * 0.02,
            _format_memory(int(mem_mb * (1024**2))),
            ha="center",
            va="bottom",
            fontsize=10,
            fontweight="bold",
        )

    ax.set_ylim(0, max(mem_mbs) * 1.15)

    out_path = get_plots_dir(name, n_points).joinpath("memory.png")
    fig.savefig(out_path)
    plt.close(fig)
    rprint(f"[green]Saved:[/green] {out_path}")


@app.command(name="all")
def all_cmd(
    name: Annotated[str, typer.Option("--name", "-n", help="Dataset name")],
    n_points: Annotated[
        int,
        typer.Option("--n-points", "-N", help="Number of sphere test points"),
    ] = 100,
) -> None:
    """Run all analysis commands."""
    summary(name=name, n_points=n_points)
    rprint("\n[bold]Generating plots...[/bold]\n")
    bar(name=name, n_points=n_points)
    memory(name=name, n_points=n_points)
    rprint(
        f"\n[bold green]All plots saved to:[/bold green] {get_plots_dir(name, n_points)}"
    )


if __name__ == "__main__":
    app()
