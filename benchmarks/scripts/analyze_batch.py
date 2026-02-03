#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.12"
# dependencies = [
#     "numpy",
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
    ./benchmarks/scripts/analyze_batch.py plot -n ecoli     # Time comparison chart
    ./benchmarks/scripts/analyze_batch.py memory -n ecoli   # Memory comparison chart
    ./benchmarks/scripts/analyze_batch.py validate -n ecoli # SASA area validation
    ./benchmarks/scripts/analyze_batch.py all               # Everything
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
    ("freesasa", None): "#3498db",  # Blue (same as analyze.py)
    ("rustsasa", None): "#e74c3c",  # Red (same as analyze.py)
    ("zsasa", "f32"): "#e67e22",  # Dark orange (same as analyze.py)
    ("zsasa", "f64"): "#f39c12",  # Light orange (same as analyze.py)
}


def _get_color(tool: str, precision: str | None) -> str:
    return COLOR_MAP.get((tool, precision), "#95a5a6")


def _make_label(tool: str, precision: str | None, threads: int) -> str:
    if precision:
        return f"{tool} ({precision}, {threads}t)"
    return f"{tool} ({threads}t)"


def parse_benchmark_name(filename: str) -> dict:
    """Parse benchmark filename to extract tool, precision, threads.

    Examples:
        bench_zsasa_f64_8t.json -> {tool: zsasa, precision: f64, threads: 8}
        bench_freesasa_1t.json -> {tool: freesasa, precision: None, threads: 1}
        bench_rustsasa_8t.json -> {tool: rustsasa, precision: None, threads: 8}
    """
    name = filename.replace("bench_", "").replace(".json", "")

    threads_match = re.search(r"(\d+)t$", name)
    threads = int(threads_match.group(1)) if threads_match else 1
    name = re.sub(r"_?\d+t$", "", name)

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
        search_dirs = sorted(d for d in RESULTS_DIR.iterdir() if d.is_dir())

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

                memory_bytes = result.get("memory_usage_byte", [])
                max_rss = max(memory_bytes) if memory_bytes else 0

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
                        "max_rss_mb": max_rss / (1024 * 1024),
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


def _sort_df(df: pl.DataFrame) -> pl.DataFrame:
    """Sort by threads ascending, then tool order, then precision."""
    return df.with_columns(
        pl.col("tool").replace_strict(TOOL_ORDER, default=99).alias("_order")
    ).sort(["threads", "_order", "precision"], descending=[False, False, False])


def _get_baseline(df: pl.DataFrame, tool: str) -> float | None:
    """Get baseline time for a tool. Returns None if not found."""
    rows = df.filter(pl.col("tool") == tool)
    if rows.height == 0:
        return None
    return rows["mean_s"][0]


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


def _format_ratio(ratio: float | None, *, is_baseline: bool) -> str:
    """Format a speedup ratio as a Rich-styled string."""
    if is_baseline:
        return "[dim]baseline[/dim]"
    if ratio is None:
        return "[dim]-[/dim]"
    color = "green" if ratio > 1 else "red"
    return f"[{color}]{ratio:.2f}x[/{color}]"


def _get_benchmarks(name: str | None) -> list[str]:
    """Get list of benchmark names to process."""
    if name:
        return [name]
    return sorted(d.name for d in RESULTS_DIR.iterdir() if d.is_dir())


# --- Commands ---


@app.command()
def summary(
    name: Annotated[
        str | None,
        typer.Option("--name", "-n", help="Benchmark name (e.g., ecoli)"),
    ] = None,
):
    """Print summary table with vs FreeSASA and vs RustSASA."""
    df = load_batch_data(name)

    for bench_name in df["benchmark"].unique().sort().to_list():
        df_bench = _sort_df(df.filter(pl.col("benchmark") == bench_name))

        fs_time = _get_baseline(df_bench, "freesasa")

        # Build RustSASA baseline per thread count
        rs_by_threads: dict[int, float] = {}
        for row in df_bench.filter(pl.col("tool") == "rustsasa").iter_rows(named=True):
            rs_by_threads[row["threads"]] = row["mean_s"]

        table = Table(title=f"Batch Benchmark: {bench_name}")
        table.add_column("Tool", style="cyan")
        table.add_column("Precision")
        table.add_column("Threads", justify="right")
        table.add_column("Mean", justify="right")
        table.add_column("Std Dev", justify="right")
        table.add_column("RSS", justify="right")
        if fs_time is not None:
            table.add_column("vs FreeSASA", justify="right")
        if rs_by_threads:
            table.add_column("vs RustSASA", justify="right")

        for row in df_bench.iter_rows(named=True):
            tool = row["tool"]
            prec = row["precision"] or "-"
            threads = row["threads"]
            mean = row["mean_s"]
            stddev = row["stddev_s"]
            rss = row["max_rss_mb"]

            is_fs = tool == "freesasa"
            is_rs = tool == "rustsasa"

            vs_fs = _format_ratio(
                fs_time / mean if fs_time else None, is_baseline=is_fs
            )

            rs_time = rs_by_threads.get(threads)
            vs_rs = _format_ratio(
                rs_time / mean if rs_time else None, is_baseline=is_rs
            )

            rss_str = f"{rss:.0f} MB" if rss > 0 else "-"

            cols = [
                tool,
                prec,
                str(threads),
                f"{mean:.3f}s",
                f"±{stddev:.3f}",
                rss_str,
            ]
            if fs_time is not None:
                cols.append(vs_fs)
            if rs_by_threads:
                cols.append(vs_rs)

            table.add_row(*cols)

        rprint(table)
        rprint()

    rprint("[dim]>1x = faster than baseline[/dim]")


@app.command()
def plot(
    name: Annotated[
        str | None,
        typer.Option("--name", "-n", help="Benchmark name (e.g., ecoli)"),
    ] = None,
):
    """Generate time comparison bar chart."""
    setup_style()

    for bench_name in _get_benchmarks(name):
        df = _sort_df(load_batch_data(bench_name))
        _plot_time(df, bench_name)


def _plot_time(df: pl.DataFrame, bench_name: str):
    """Generate time comparison bar chart for one benchmark."""
    fs_time = _get_baseline(df, "freesasa")

    # RustSASA baseline per thread count
    rs_by_threads: dict[int, float] = {}
    for row in df.filter(pl.col("tool") == "rustsasa").iter_rows(named=True):
        rs_by_threads[row["threads"]] = row["mean_s"]

    labels = []
    times = []
    errors = []
    colors = []
    thread_counts = []

    for row in df.iter_rows(named=True):
        labels.append(_make_label(row["tool"], row["precision"], row["threads"]))
        times.append(row["mean_s"])
        errors.append(row["stddev_s"])
        colors.append(_get_color(row["tool"], row["precision"]))
        thread_counts.append(row["threads"])

    fig, ax = plt.subplots(figsize=(max(8, len(labels) * 1.5), 6))

    x = range(len(labels))
    bars = ax.bar(
        x, times, yerr=errors, color=colors, width=0.6, edgecolor="white", capsize=4
    )

    # Annotate bars with time and speedup (same thread count comparison)
    for bar, t, tc in zip(bars, times, thread_counts):
        rs_time = rs_by_threads.get(tc)
        parts = [f"{t:.1f}s"]
        if fs_time and t != fs_time:
            parts.append(f"vs FS: {fs_time / t:.1f}x")
        if rs_time and t != rs_time:
            parts.append(f"vs RS: {rs_time / t:.1f}x")

        ax.text(
            bar.get_x() + bar.get_width() / 2,
            bar.get_height() + bar.get_height() * 0.03,
            "\n".join(parts),
            ha="center",
            va="bottom",
            fontsize=9,
            fontweight="bold",
        )

    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=30, ha="right")
    ax.set_ylabel("Time (seconds)")
    ax.grid(True, alpha=0.3, axis="y")
    ax.set_title(f"Batch Processing Time: {bench_name}")

    # Add margin at top for labels
    ymax = max(times) * 1.25
    ax.set_ylim(0, ymax)

    PLOTS_DIR.mkdir(parents=True, exist_ok=True)
    out_path = PLOTS_DIR.joinpath(f"{bench_name}_time.png")
    fig.savefig(out_path)
    plt.close(fig)
    rprint(f"[green]Saved:[/green] {out_path}")


@app.command()
def memory(
    name: Annotated[
        str | None,
        typer.Option("--name", "-n", help="Benchmark name (e.g., ecoli)"),
    ] = None,
):
    """Generate memory usage (RSS) bar chart."""
    setup_style()

    for bench_name in _get_benchmarks(name):
        df = load_batch_data(bench_name)
        df = df.filter(pl.col("max_rss_mb") > 0)
        if df.height == 0:
            rprint(f"[yellow]No memory data for {bench_name}[/yellow]")
            continue
        df = _sort_df(df)
        _plot_memory(df, bench_name)


def _plot_memory(df: pl.DataFrame, bench_name: str):
    """Generate memory usage bar chart for one benchmark."""
    labels = []
    rss_values = []
    colors = []

    for row in df.iter_rows(named=True):
        labels.append(_make_label(row["tool"], row["precision"], row["threads"]))
        rss_values.append(row["max_rss_mb"])
        colors.append(_get_color(row["tool"], row["precision"]))

    fig, ax = plt.subplots(figsize=(max(8, len(labels) * 1.5), 6))

    x = range(len(labels))
    bars = ax.bar(x, rss_values, color=colors, width=0.6, edgecolor="white")

    for bar, rss in zip(bars, rss_values):
        ax.text(
            bar.get_x() + bar.get_width() / 2,
            bar.get_height() + max(rss_values) * 0.02,
            f"{rss:.0f} MB",
            ha="center",
            va="bottom",
            fontsize=10,
            fontweight="bold",
        )

    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=30, ha="right")
    ax.set_ylabel("Peak RSS (MB)")
    ax.grid(True, alpha=0.3, axis="y")
    ax.set_title(f"Memory Usage: {bench_name}")
    ax.set_ylim(0, max(rss_values) * 1.15)

    PLOTS_DIR.mkdir(parents=True, exist_ok=True)
    out_path = PLOTS_DIR.joinpath(f"{bench_name}_memory.png")
    fig.savefig(out_path)
    plt.close(fig)
    rprint(f"[green]Saved:[/green] {out_path}")


def _load_sasa_values(
    temp_out: Path,
) -> dict[str, dict[str, float]]:
    """Load total SASA from each tool's output files.

    Returns {tool_name: {structure_name: total_sasa}}.
    """
    tools: dict[str, dict[str, float]] = {}

    for tool_dir in sorted(temp_out.iterdir()):
        if not tool_dir.is_dir():
            continue

        tool_name = tool_dir.name
        sasa_map: dict[str, float] = {}

        # Determine parser based on directory name
        if tool_name.startswith("zig"):
            parser = "zig"
        elif tool_name == "freesasa":
            parser = "freesasa"
        elif tool_name == "rustsasa":
            parser = "rustsasa"
        else:
            continue

        for json_file in sorted(tool_dir.glob("*.json")):
            structure = json_file.stem
            try:
                with open(json_file) as f:
                    data = json.load(f)

                if parser == "zig":
                    total = data.get("total_area", 0.0)
                elif parser == "freesasa":
                    total = 0.0
                    for struct in data.get("results", [{}])[0].get("structure", []):
                        for chain in struct.get("chains", []):
                            total += chain.get("area", {}).get("total", 0.0)
                elif parser == "rustsasa":
                    total = sum(r.get("value", 0.0) for r in data.get("Residue", []))
                else:
                    continue

                sasa_map[structure] = total
            except (json.JSONDecodeError, KeyError, IndexError):
                continue

        if sasa_map:
            tools[tool_name] = sasa_map

    return tools


def _compute_pairwise_stats(
    values_a: list[float], values_b: list[float]
) -> dict[str, float]:
    """Compute pairwise comparison statistics."""
    import numpy as np

    a = np.array(values_a)
    b = np.array(values_b)
    diff = a - b
    abs_diff = np.abs(diff)

    # Relative error (avoid division by zero)
    nonzero = b != 0
    rel_errors = np.zeros_like(a)
    rel_errors[nonzero] = np.abs(diff[nonzero]) / b[nonzero] * 100

    r_squared = np.corrcoef(a, b)[0, 1] ** 2

    return {
        "n": len(a),
        "r_squared": float(r_squared),
        "mean_abs_diff": float(np.mean(abs_diff)),
        "max_abs_diff": float(np.max(abs_diff)),
        "mean_rel_error_pct": float(np.mean(rel_errors)),
        "max_rel_error_pct": float(np.max(rel_errors)),
        "median_rel_error_pct": float(np.median(rel_errors)),
    }


@app.command()
def validate(
    name: Annotated[
        str,
        typer.Option("--name", "-n", help="Benchmark name (e.g., ecoli)"),
    ] = "ecoli",
):
    """Validate SASA areas across tools (zig vs FreeSASA vs RustSASA)."""
    bench_dir = RESULTS_DIR.joinpath(name)
    temp_out = bench_dir.joinpath("temp_out")

    if not temp_out.exists():
        rprint(f"[red]No temp_out directory found in {bench_dir}[/red]")
        raise typer.Exit(1)

    tools = _load_sasa_values(temp_out)
    tool_names = sorted(tools.keys())
    rprint(f"[bold]SASA Validation: {name}[/bold]")
    rprint(f"Tools: {', '.join(tool_names)}")
    for t in tool_names:
        rprint(f"  {t}: {len(tools[t])} structures")
    rprint()

    # Pairwise comparison table
    pairs = []
    for i, t1 in enumerate(tool_names):
        for t2 in tool_names[i + 1 :]:
            pairs.append((t1, t2))

    table = Table(title="Pairwise SASA Comparison")
    table.add_column("Pair", style="cyan")
    table.add_column("N", justify="right")
    table.add_column("R\u00b2", justify="right")
    table.add_column("Mean |diff|", justify="right")
    table.add_column("Max |diff|", justify="right")
    table.add_column("Mean rel err", justify="right")
    table.add_column("Max rel err", justify="right")
    table.add_column("Median rel err", justify="right")

    pair_data: dict[tuple[str, str], tuple[list[str], list[float], list[float]]] = {}
    for t1, t2 in pairs:
        common = sorted(set(tools[t1].keys()) & set(tools[t2].keys()))
        if not common:
            continue

        vals_a = [tools[t1][s] for s in common]
        vals_b = [tools[t2][s] for s in common]
        stats = _compute_pairwise_stats(vals_a, vals_b)
        pair_data[(t1, t2)] = (common, vals_a, vals_b)

        table.add_row(
            f"{t1} vs {t2}",
            str(stats["n"]),
            f"{stats['r_squared']:.8f}",
            f"{stats['mean_abs_diff']:.4f} \u00c5\u00b2",
            f"{stats['max_abs_diff']:.4f} \u00c5\u00b2",
            f"{stats['mean_rel_error_pct']:.4f}%",
            f"{stats['max_rel_error_pct']:.4f}%",
            f"{stats['median_rel_error_pct']:.4f}%",
        )

    rprint(table)
    rprint()

    # Generate scatter plots
    setup_style()
    plot_dir = PLOTS_DIR.joinpath(f"{name}")
    plot_dir.mkdir(parents=True, exist_ok=True)

    for (t1, t2), (common, vals_a, vals_b) in pair_data.items():
        import numpy as np

        a = np.array(vals_a)
        b = np.array(vals_b)
        stats = _compute_pairwise_stats(vals_a, vals_b)

        fig, axes = plt.subplots(1, 2, figsize=(14, 6))

        # Left: scatter plot
        ax = axes[0]
        ax.scatter(a, b, alpha=0.3, s=8, color="#3498db", edgecolors="none")
        max_val = max(a.max(), b.max()) * 1.05
        ax.plot([0, max_val], [0, max_val], "r--", linewidth=1, label="y = x")
        ax.set_xlabel(f"{t1} SASA (\u00c5\u00b2)")
        ax.set_ylabel(f"{t2} SASA (\u00c5\u00b2)")
        ax.set_title(f"{t1} vs {t2} (n={len(common):,})")
        ax.set_aspect("equal")
        ax.legend(loc="lower right")

        stats_text = (
            f"R\u00b2 = {stats['r_squared']:.8f}\n"
            f"Mean err = {stats['mean_rel_error_pct']:.4f}%\n"
            f"Max err = {stats['max_rel_error_pct']:.4f}%"
        )
        ax.text(
            0.05,
            0.95,
            stats_text,
            transform=ax.transAxes,
            fontsize=9,
            verticalalignment="top",
            bbox=dict(boxstyle="round", facecolor="wheat", alpha=0.7),
        )

        # Right: relative error histogram
        ax2 = axes[1]
        nonzero = b != 0
        rel_err = np.abs(a[nonzero] - b[nonzero]) / b[nonzero] * 100
        ax2.hist(rel_err, bins=50, color="#3498db", edgecolor="white", alpha=0.8)
        ax2.axvline(
            np.median(rel_err),
            color="red",
            linestyle="--",
            label=f"median={np.median(rel_err):.4f}%",
        )
        ax2.set_xlabel("Relative Error (%)")
        ax2.set_ylabel("Count")
        ax2.set_title("Relative Error Distribution")
        ax2.legend()

        fig.suptitle(f"SASA Validation: {name}", fontsize=14, fontweight="bold")
        fig.tight_layout()

        out_path = plot_dir.joinpath(f"{name}_validate_{t1}_vs_{t2}.png")
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
    """Generate summary, time chart, memory chart, and validation."""
    summary(name)
    rprint()
    plot(name)
    rprint()
    memory(name)
    if name:
        rprint()
        validate(name)


if __name__ == "__main__":
    app()
