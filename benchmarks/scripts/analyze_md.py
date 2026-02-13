#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.12"
# dependencies = [
#     "matplotlib",
#     "typer>=0.9.0",
#     "rich>=13.0",
# ]
# ///
"""MD trajectory benchmark analysis CLI.

Generates summary statistics and plots from bench_md.py results.

Usage:
    ./benchmarks/scripts/analyze_md.py summary --name 6sup_R1
    ./benchmarks/scripts/analyze_md.py bar --name 6sup_R1
    ./benchmarks/scripts/analyze_md.py threads --name 6sup_R1
    ./benchmarks/scripts/analyze_md.py memory --name 6sup_R1
    ./benchmarks/scripts/analyze_md.py speedup --name 6sup_R1
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

RESULTS_DIR = Path(__file__).parent.parent.joinpath("results", "md")

DISPLAY_NAMES: dict[str, str] = {
    "zig_f32": "zsasa CLI (f32)",
    "zig_f64": "zsasa CLI (f64)",
    "zsasa_mdtraj": "zsasa.mdtraj",
    "zsasa_mdanalysis": "zsasa.mdanalysis",
    "mdtraj": "MDTraj native",
    "mdsasa_bolt": "mdsasa-bolt",
}

COLORS: dict[str, str] = {
    "zig_f32": "#e67e22",
    "zig_f64": "#f39c12",
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


def load_results(name: str) -> list[BenchResult]:
    """Load all benchmark results for a given dataset name."""
    results_dir = RESULTS_DIR.joinpath(name)
    if not results_dir.exists():
        available = [d.name for d in RESULTS_DIR.iterdir() if d.is_dir()]
        raise typer.BadParameter(
            f"Dataset '{name}' not found. Available: {', '.join(sorted(available))}"
        )

    # Resolve "all" threads from config cpu_cores
    config = load_config(name)
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


def load_config(name: str) -> dict | None:
    """Load config.json for dataset metadata."""
    config_path = RESULTS_DIR.joinpath(name, "config.json")
    if not config_path.exists():
        return None
    try:
        return json.loads(config_path.read_text())
    except json.JSONDecodeError:
        rprint(f"[yellow]Warning: Could not parse {config_path}[/yellow]")
        return None


def get_plots_dir(name: str) -> Path:
    """Get (and create) plots output directory."""
    plots_dir = RESULTS_DIR.joinpath(name, "plots")
    plots_dir.mkdir(parents=True, exist_ok=True)
    return plots_dir


def display_name(tool: str) -> str:
    """Get display name for a tool."""
    return DISPLAY_NAMES.get(tool, tool)


def get_best_per_tool(results: list[BenchResult]) -> dict[str, BenchResult]:
    """Get the fastest result per tool (best thread count)."""
    best: dict[str, BenchResult] = {}
    for r in results:
        if r.tool not in best or r.mean < best[r.tool].mean:
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


# === CLI Commands ===


@app.command()
def summary(
    name: Annotated[str, typer.Option("--name", "-n", help="Dataset name")],
) -> None:
    """Print summary table of all benchmark results."""
    results = load_results(name)
    config = load_config(name)

    if config:
        params = config.get("parameters", {})
        system = config.get("system", {})
        rprint(f"[bold]Dataset:[/bold] {name}")
        rprint(
            f"[dim]System: {system.get('cpu_model', 'N/A')} "
            f"({system.get('cpu_cores', '?')} cores, "
            f"{system.get('memory_gb', '?')} GB RAM)[/dim]"
        )
        rprint(
            f"[dim]Params: stride={params.get('stride', '?')}, "
            f"n_points={params.get('n_points', '?')}, "
            f"runs={params.get('runs', '?')}[/dim]"
        )
        rprint()

    mdtraj_baseline = get_mdtraj_baseline(results)

    table = Table(title=f"MD Benchmark Results: {name}")
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
) -> None:
    """Generate horizontal bar chart comparing best time per tool."""
    setup_style()
    results = load_results(name)
    best = get_best_per_tool(results)

    # Sort by time (fastest first from top)
    ordered = [t for t in TOOL_ORDER if t in best]
    if not ordered:
        rprint("[yellow]No matching tools found for bar chart.[/yellow]")
        return
    ordered.sort(key=lambda t: best[t].mean, reverse=True)

    labels = []
    times = []
    colors = []
    for tool in ordered:
        r = best[tool]
        thread_info = f" ({r.threads}t)" if r.threads is not None else " (all)"
        labels.append(f"{display_name(tool)}{thread_info}")
        times.append(r.mean)
        colors.append(COLORS.get(tool, "#95a5a6"))

    fig, ax = plt.subplots(figsize=(10, max(4, len(labels) * 0.8)))
    y_pos = range(len(labels))
    bars = ax.barh(y_pos, times, color=colors, height=0.6, edgecolor="white")

    ax.set_yticks(y_pos)
    ax.set_yticklabels(labels)
    ax.set_xlabel("Time (seconds)")
    ax.set_title(f"MD Trajectory SASA Benchmark: {name}")
    ax.grid(True, alpha=0.3, axis="x")

    for bar_item, t in zip(bars, times):
        ax.text(
            bar_item.get_width() + max(times) * 0.01,
            bar_item.get_y() + bar_item.get_height() / 2,
            f"{t:.1f}s",
            va="center",
            fontsize=10,
        )

    ax.set_xlim(0, max(times) * 1.15)
    fig.tight_layout()

    out_path = get_plots_dir(name).joinpath("bar.png")
    fig.savefig(out_path)
    plt.close(fig)
    rprint(f"[green]Saved:[/green] {out_path}")


@app.command()
def threads(
    name: Annotated[str, typer.Option("--name", "-n", help="Dataset name")],
) -> None:
    """Generate thread scaling plot."""
    setup_style()
    results = load_results(name)

    # Tools that don't scale with threads - always show as reference line
    reference_tools = {"mdtraj", "mdsasa_bolt"}

    # Group by tool
    tool_threads: dict[str, list[tuple[int, float]]] = {}
    ref_lines: dict[str, float] = {}

    for r in results:
        if r.tool in reference_tools:
            ref_lines[r.tool] = r.mean
        elif r.threads is not None:
            tool_threads.setdefault(r.tool, []).append((r.threads, r.mean))

    fig, ax = plt.subplots(figsize=(10, 6))

    # Plot scaling lines
    for tool in TOOL_ORDER:
        if tool not in tool_threads:
            continue
        data = sorted(tool_threads[tool])
        t_vals = [d[0] for d in data]
        time_vals = [d[1] for d in data]
        ax.plot(
            t_vals,
            time_vals,
            marker="o",
            label=display_name(tool),
            color=COLORS.get(tool, "#95a5a6"),
            linewidth=2,
            markersize=6,
        )

    # Plot reference lines
    for tool, time_val in ref_lines.items():
        ax.axhline(
            y=time_val,
            linestyle="--",
            color=COLORS.get(tool, "#95a5a6"),
            linewidth=1.5,
            label=f"{display_name(tool)} ({time_val:.0f}s)",
            alpha=0.7,
        )

    ax.set_xlabel("Thread Count")
    ax.set_ylabel("Time (seconds)")
    ax.set_title(f"Thread Scaling: {name}")
    ax.legend(loc="upper right")
    ax.grid(True, alpha=0.3)

    # Set x ticks to actual thread counts
    all_threads = sorted({t for pts in tool_threads.values() for t, _ in pts})
    if all_threads:
        ax.set_xticks(all_threads)
        ax.set_xticklabels([str(t) for t in all_threads])

    fig.tight_layout()
    out_path = get_plots_dir(name).joinpath("threads.png")
    fig.savefig(out_path)
    plt.close(fig)
    rprint(f"[green]Saved:[/green] {out_path}")


@app.command()
def memory(
    name: Annotated[str, typer.Option("--name", "-n", help="Dataset name")],
) -> None:
    """Generate memory usage comparison bar chart."""
    setup_style()
    results = load_results(name)
    best = get_best_per_tool(results)

    ordered = [t for t in TOOL_ORDER if t in best]
    if not ordered:
        rprint("[yellow]No matching tools found for memory chart.[/yellow]")
        return
    ordered.sort(key=lambda t: best[t].memory_bytes, reverse=True)

    labels = []
    mem_mbs = []
    colors = []
    for tool in ordered:
        r = best[tool]
        labels.append(display_name(tool))
        mem_mbs.append(r.memory_bytes / (1024**2))
        colors.append(COLORS.get(tool, "#95a5a6"))

    fig, ax = plt.subplots(figsize=(10, max(4, len(labels) * 0.8)))
    y_pos = range(len(labels))
    bars = ax.barh(y_pos, mem_mbs, color=colors, height=0.6, edgecolor="white")

    ax.set_yticks(y_pos)
    ax.set_yticklabels(labels)
    ax.set_xlabel("Peak Memory (MB)")
    ax.set_title(f"Memory Usage: {name}")
    ax.grid(True, alpha=0.3, axis="x")

    for bar_item, mem_mb in zip(bars, mem_mbs):
        ax.text(
            bar_item.get_width() + max(mem_mbs) * 0.01,
            bar_item.get_y() + bar_item.get_height() / 2,
            _format_memory(int(mem_mb * (1024**2))),
            va="center",
            fontsize=10,
        )

    ax.set_xlim(0, max(mem_mbs) * 1.15)
    fig.tight_layout()

    out_path = get_plots_dir(name).joinpath("memory.png")
    fig.savefig(out_path)
    plt.close(fig)
    rprint(f"[green]Saved:[/green] {out_path}")


@app.command()
def speedup(
    name: Annotated[str, typer.Option("--name", "-n", help="Dataset name")],
) -> None:
    """Show speedup ratios and generate comparison chart."""
    setup_style()
    results = load_results(name)
    best = get_best_per_tool(results)

    mdtraj_time = get_mdtraj_baseline(results)
    bolt_result = best.get("mdsasa_bolt")
    bolt_time = bolt_result.mean if bolt_result is not None else None

    # Rich table
    table = Table(title=f"Speedup Summary: {name}")
    table.add_column("Tool", style="cyan")
    table.add_column("Threads", justify="right")
    table.add_column("Time (s)", justify="right")
    if mdtraj_time is not None:
        table.add_column("vs MDTraj", justify="right")
    if bolt_time is not None:
        table.add_column("vs mdsasa-bolt", justify="right")

    ordered = [t for t in TOOL_ORDER if t in best]
    if not ordered:
        rprint("[yellow]No matching tools found for speedup analysis.[/yellow]")
        return
    ordered.sort(key=lambda t: best[t].mean)

    for tool in ordered:
        r = best[tool]
        row = [
            display_name(tool),
            _threads_label(r.threads),
            f"{r.mean:.1f}",
        ]

        if mdtraj_time is not None and r.mean > 0:
            s = mdtraj_time / r.mean
            s_str = f"{s:.1f}x"
            if s > 1.0:
                s_str = f"[green]{s_str}[/green]"
            elif s < 1.0:
                s_str = f"[red]{s_str}[/red]"
            row.append(s_str)

        if bolt_time is not None and r.mean > 0:
            s = bolt_time / r.mean
            s_str = f"{s:.1f}x"
            if s > 1.0:
                s_str = f"[green]{s_str}[/green]"
            elif s < 1.0:
                s_str = f"[red]{s_str}[/red]"
            row.append(s_str)

        table.add_row(*row)

    rprint(table)

    # Bar chart: speedup vs MDTraj
    if mdtraj_time is None:
        return

    chart_tools = [t for t in ordered if t != "mdtraj" and best[t].mean > 0]
    if not chart_tools:
        return
    labels = [display_name(t) for t in chart_tools]
    speedups = [mdtraj_time / best[t].mean for t in chart_tools]
    colors = [COLORS.get(t, "#95a5a6") for t in chart_tools]

    fig, ax = plt.subplots(figsize=(10, max(4, len(labels) * 0.8)))
    y_pos = range(len(labels))
    bars = ax.barh(y_pos, speedups, color=colors, height=0.6, edgecolor="white")

    ax.set_yticks(y_pos)
    ax.set_yticklabels(labels)
    ax.set_xlabel("Speedup vs MDTraj native (higher = faster)")
    ax.set_title(f"Speedup vs MDTraj native: {name}")
    ax.axvline(x=1.0, color="gray", linestyle="--", linewidth=1.5)
    ax.grid(True, alpha=0.3, axis="x")

    for bar_item, s in zip(bars, speedups):
        ax.text(
            bar_item.get_width() + max(speedups) * 0.01,
            bar_item.get_y() + bar_item.get_height() / 2,
            f"{s:.1f}x",
            va="center",
            fontsize=10,
        )

    ax.set_xlim(0, max(speedups) * 1.15)
    fig.tight_layout()

    out_path = get_plots_dir(name).joinpath("speedup.png")
    fig.savefig(out_path)
    plt.close(fig)
    rprint(f"[green]Saved:[/green] {out_path}")


@app.command(name="all")
def all_cmd(
    name: Annotated[str, typer.Option("--name", "-n", help="Dataset name")],
) -> None:
    """Run all analysis commands."""
    summary(name=name)
    rprint("\n[bold]Generating plots...[/bold]\n")
    bar(name=name)
    threads(name=name)
    memory(name=name)
    speedup(name=name)
    rprint(f"\n[bold green]All plots saved to:[/bold green] {get_plots_dir(name)}")


if __name__ == "__main__":
    app()
