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
"""Analyze batch benchmark results (hyperfine JSON format).

Usage:
    ./benchmarks/scripts/analyze_batch.py summary            # Summary table + files/sec
    ./benchmarks/scripts/analyze_batch.py summary -n ecoli   # Specific benchmark
    ./benchmarks/scripts/analyze_batch.py plot -n ecoli      # Time comparison charts (per thread)
    ./benchmarks/scripts/analyze_batch.py memory -n ecoli    # Memory comparison charts (per thread)
    ./benchmarks/scripts/analyze_batch.py all                # Everything
"""

from __future__ import annotations

import json
import re
import statistics
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

TOOL_ORDER = {"freesasa": 0, "rustsasa": 1, "zsasa": 2, "lahuta": 3}

COLOR_MAP = {
    ("freesasa", None): "#3498db",
    ("rustsasa", None): "#e74c3c",
    ("zsasa", "f32"): "#e67e22",
    ("zsasa", "f64"): "#f39c12",
    ("lahuta", None): "#9b59b6",
}


def _get_color(tool: str, precision: str | None) -> str:
    return COLOR_MAP.get((tool, precision), "#95a5a6")


def _make_label(tool: str, precision: str | None) -> str:
    if precision:
        return f"{tool} ({precision})"
    return tool


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


def load_config(benchmark_name: str) -> dict | None:
    """Load config.json for a benchmark."""
    config_path = RESULTS_DIR.joinpath(benchmark_name, "config.json")
    if not config_path.exists():
        return None
    with open(config_path) as f:
        return json.load(f)


def _filter_outliers(times: list[float]) -> list[float]:
    """Remove outliers using IQR method. Requires 5+ data points.

    Only removes points that deviate >10% from the median to avoid
    over-filtering low-variance data where IQR is extremely tight.
    """
    if len(times) < 5:
        return times
    sorted_t = sorted(times)
    n = len(sorted_t)
    q1 = sorted_t[n // 4]
    q3 = sorted_t[3 * n // 4]
    iqr = q3 - q1
    lower = q1 - 1.5 * iqr
    upper = q3 + 1.5 * iqr
    med = statistics.median(times)
    threshold = med * 0.1
    filtered = [t for t in times if lower <= t <= upper or abs(t - med) <= threshold]
    if filtered and len(filtered) < len(times):
        removed = len(times) - len(filtered)
        rprint(
            f"[dim]  Removed {removed} outlier(s): kept {len(filtered)}/{len(times)} runs[/dim]"
        )
    return filtered if filtered else times


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

                times = result.get("times", [])
                filtered = _filter_outliers(times)
                mean_s = statistics.mean(filtered) if filtered else result["mean"]
                stddev_s = (
                    statistics.stdev(filtered)
                    if len(filtered) >= 2
                    else result["stddev"]
                )

                memory_bytes = result.get("memory_usage_byte", [])
                max_rss = max(memory_bytes) if memory_bytes else 0

                rows.append(
                    {
                        "benchmark": bench_dir.name,
                        "tool": info["tool"],
                        "precision": info["precision"],
                        "threads": info["threads"],
                        "mean_s": mean_s,
                        "stddev_s": stddev_s,
                        "min_s": min(filtered) if filtered else result["min"],
                        "max_s": max(filtered) if filtered else result["max"],
                        "runs": len(filtered),
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
    """Print summary table with speedup ratios and throughput."""
    df = load_batch_data(name)

    for bench_name in df["benchmark"].unique().sort().to_list():
        df_bench = _sort_df(df.filter(pl.col("benchmark") == bench_name))

        # Load config for n_files
        config = load_config(bench_name)
        n_files = config["parameters"]["n_files"] if config else None

        # FreeSASA baseline (always 1t)
        fs_rows = df_bench.filter(pl.col("tool") == "freesasa")
        fs_time = fs_rows["mean_s"][0] if fs_rows.height > 0 else None

        # RustSASA baseline per thread count
        rs_by_threads: dict[int, float] = {}
        for row in df_bench.filter(pl.col("tool") == "rustsasa").iter_rows(named=True):
            rs_by_threads[row["threads"]] = row["mean_s"]

        # Build table
        table = Table(title=f"Batch Benchmark: {bench_name}")
        table.add_column("Tool", style="cyan")
        table.add_column("Precision")
        table.add_column("Threads", justify="right")
        table.add_column("Mean", justify="right")
        table.add_column("Std Dev", justify="right")
        table.add_column("RSS", justify="right")
        if n_files:
            table.add_column("files/sec", justify="right")
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
            fps_str = f"{n_files / mean:,.0f}" if n_files else ""

            cols = [
                tool,
                prec,
                str(threads),
                f"{mean:.3f}s",
                f"±{stddev:.3f}",
                rss_str,
            ]
            if n_files:
                cols.append(fps_str)
            if fs_time is not None:
                cols.append(vs_fs)
            if rs_by_threads:
                cols.append(vs_rs)

            table.add_row(*cols)

        rprint(table)
        if n_files:
            rprint(f"[dim]  {n_files:,} files in dataset[/dim]")
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
    """Generate time comparison bar charts per thread count."""
    PLOTS_DIR.mkdir(parents=True, exist_ok=True)
    config = load_config(bench_name)

    thread_counts = sorted(df["threads"].unique().to_list())

    for threads in thread_counts:
        df_t = _sort_df(df.filter(pl.col("threads") == threads))

        labels = []
        times = []
        errors = []
        colors = []

        for row in df_t.iter_rows(named=True):
            labels.append(_make_label(row["tool"], row["precision"]))
            times.append(row["mean_s"])
            errors.append(row["stddev_s"])
            colors.append(_get_color(row["tool"], row["precision"]))

        fig, ax = plt.subplots(figsize=(max(6, len(labels) * 1.8), 6))

        x = range(len(labels))
        bars = ax.bar(
            x,
            times,
            yerr=errors,
            color=colors,
            width=0.6,
            edgecolor="white",
            capsize=4,
        )

        for bar, t in zip(bars, times):
            ax.text(
                bar.get_x() + bar.get_width() / 2,
                bar.get_height() + max(times) * 0.02,
                f"{t:.2f}s",
                ha="center",
                va="bottom",
                fontsize=10,
                fontweight="bold",
            )

        ax.set_xticks(x)
        ax.set_xticklabels(labels, rotation=30, ha="right")
        ax.set_ylabel("Time (seconds)")
        ax.grid(True, alpha=0.3, axis="y")

        title = f"Batch Processing Time: {bench_name} ({threads}t)"
        if config:
            p = config["parameters"]
            title += f"\n{p['n_files']:,} files, warmup={p['warmup']}, runs={p['runs']}"
        ax.set_title(title)
        ax.set_ylim(0, max(times) * 1.15)

        out_path = PLOTS_DIR.joinpath(f"{bench_name}_time_{threads}t.png")
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
    """Generate memory usage bar charts per thread count."""
    PLOTS_DIR.mkdir(parents=True, exist_ok=True)
    config = load_config(bench_name)

    thread_counts = sorted(df["threads"].unique().to_list())

    for threads in thread_counts:
        df_t = _sort_df(df.filter(pl.col("threads") == threads))

        labels = []
        rss_values = []
        colors = []

        for row in df_t.iter_rows(named=True):
            labels.append(_make_label(row["tool"], row["precision"]))
            rss_values.append(row["max_rss_mb"])
            colors.append(_get_color(row["tool"], row["precision"]))

        if not rss_values or max(rss_values) == 0:
            continue

        fig, ax = plt.subplots(figsize=(max(6, len(labels) * 1.8), 6))

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

        title = f"Memory Usage: {bench_name} ({threads}t)"
        if config:
            p = config["parameters"]
            title += f"\n{p['n_files']:,} files"
        ax.set_title(title)
        ax.set_ylim(0, max(rss_values) * 1.15)

        out_path = PLOTS_DIR.joinpath(f"{bench_name}_memory_{threads}t.png")
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
    """Generate summary, time charts, and memory charts."""
    summary(name)
    rprint()
    plot(name)
    rprint()
    memory(name)


if __name__ == "__main__":
    app()
