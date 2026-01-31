#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.12"
# dependencies = [
#     "polars",
#     "matplotlib",
#     "typer",
#     "rich",
#     "pyyaml",
# ]
# ///
"""Generate benchmark comparison plots."""

from pathlib import Path

import matplotlib.pyplot as plt
import polars as pl
import typer
from rich import print as rprint
from rich.table import Table

app = typer.Typer(help="Generate benchmark plots")

RESULTS_DIR = Path(__file__).parent.parent / "results"
PLOTS_DIR = RESULTS_DIR / "plots"

# Color palette for tools
COLORS = {
    "zig": "#f39c12",  # Orange (legacy, keep for backward compat)
    "zig_f32": "#e67e22",  # Dark orange
    "zig_f64": "#f39c12",  # Light orange
    "freesasa": "#3498db",  # Blue
    "rust": "#e74c3c",  # Red
}

# Line styles for tools (solid by default, dashed for f32)
LINESTYLES = {
    "zig": "-",
    "zig_f32": "--",  # Dashed for f32
    "zig_f64": "-",
    "freesasa": "-",
    "rust": "-",
}


def load_config() -> dict:
    """Load benchmark configuration from YAML file."""
    import yaml

    config_path = RESULTS_DIR / "config.yaml"
    if not config_path.exists():
        raise FileNotFoundError(f"Config file not found: {config_path}")
    with open(config_path) as f:
        return yaml.safe_load(f)


def load_data() -> pl.DataFrame:
    """Load benchmark results specified in config.yaml."""
    config = load_config()

    # Collect all directories from config (sr + lr)
    dirs = config.get("sr", []) + config.get("lr", [])
    csv_files = [RESULTS_DIR / d / "results.csv" for d in dirs]
    csv_files = [f for f in csv_files if f.exists()]

    if not csv_files:
        raise FileNotFoundError("No results.csv files found in configured directories")

    dfs = []
    for f in csv_files:
        df = pl.read_csv(f)
        dir_name = f.parent.name  # e.g., "zig_sr_f32", "zig_sr_f64", "freesasa_sr"

        # Extract precision from directory name (only for zig)
        if "_f32" in dir_name:
            precision = "f32"
        elif "_f64" in dir_name:
            precision = "f64"
        else:
            precision = None  # Non-zig tools

        # Add precision column
        df = df.with_columns(pl.lit(precision).alias("precision"))
        dfs.append(df)

    df = pl.concat(dfs, how="diagonal")  # diagonal handles schema differences

    # Create compound tool identifier for zig variants
    df = df.with_columns(
        pl.when((pl.col("tool") == "zig") & pl.col("precision").is_not_null())
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
            pl.col("total_sasa").first(),
        )
        .sort(["algorithm", "tool_label", "n_atoms"])
    )


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
    table.add_column("P95 (ms)", justify="right")

    stats = (
        df_t1.group_by(["algorithm", "tool_label"])
        .agg(
            pl.len().alias("n"),
            pl.col("time_ms").median().alias("median"),
            pl.col("time_ms").mean().alias("mean"),
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
            f"{row['p95']:.2f}",
        )

    rprint(table)

    # Speedup summary - use tool_label for pivoting
    rprint("\n[bold]Speedup Ratios (Zig vs others):[/bold]")
    for algo in ["sr", "lr"]:
        df_algo = df_t1.filter(pl.col("algorithm") == algo)
        pivot = (
            df_algo.select(["structure", "tool_label", "time_ms"])
            .pivot(on="tool_label", index="structure", values="time_ms")
            .drop_nulls()
        )

        # zig_f64 vs FreeSASA (primary comparison)
        if "zig_f64" in pivot.columns and "freesasa" in pivot.columns:
            speedup = pivot.select(
                (pl.col("freesasa") / pl.col("zig_f64")).alias("speedup")
            )
            med = speedup["speedup"].median()
            mean = speedup["speedup"].mean()
            rprint(
                f"  {algo.upper()}: Zig(f64) vs FreeSASA = [green]{med:.2f}x[/green] (median), {mean:.2f}x (mean)"
            )

        # zig_f32 vs zig_f64 (f32 performance relative to f64)
        if "zig_f32" in pivot.columns and "zig_f64" in pivot.columns:
            speedup = pivot.select(
                (pl.col("zig_f64") / pl.col("zig_f32")).alias("speedup")
            )
            med = speedup["speedup"].median()
            rprint(
                f"  {algo.upper()}: Zig(f32) vs Zig(f64) = [green]{med:.2f}x[/green] (median)"
            )

        if algo == "sr" and "rust" in pivot.columns and "zig_f64" in pivot.columns:
            speedup_rust = pivot.select(
                (pl.col("rust") / pl.col("zig_f64")).alias("speedup")
            )
            med = speedup_rust["speedup"].median()
            rprint(
                f"  {algo.upper()}: Zig(f64) vs Rust = [green]{med:.2f}x[/green] (median)"
            )

    # Speedup by size bin table (SR only)
    rprint("\n")
    df_sr = df_t1.filter(pl.col("algorithm") == "sr")
    speedup_by_bin = compute_speedup_by_bin(df_sr, threads=1)

    bin_table = Table(title="SR Speedup by Structure Size (threads=1)")
    bin_table.add_column("Size Bin", style="cyan")
    bin_table.add_column("Count", justify="right")
    bin_table.add_column("Zig(f64) vs FreeSASA", justify="right")
    bin_table.add_column("Zig(f64) vs Rust", justify="right")
    bin_table.add_column("Zig(f32) vs Zig(f64)", justify="right")

    # Sort by bin order
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

        # Color based on speedup
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


@app.command()
def export_csv():
    """Export summary tables as CSV files per thread count."""
    df = load_data()
    df = add_size_bin(df)

    csv_dir = RESULTS_DIR / "csv"
    csv_dir.mkdir(parents=True, exist_ok=True)

    thread_counts = sorted(df["threads"].unique().to_list())

    for threads in thread_counts:
        t_dir = csv_dir / f"t{threads}"
        t_dir.mkdir(parents=True, exist_ok=True)

        df_t = df.filter(pl.col("threads") == threads)

        # 1. Performance summary (with tool_label for zig variants)
        perf_stats = (
            df_t.group_by(["algorithm", "tool_label", "precision"])
            .agg(
                pl.len().alias("structures"),
                pl.col("time_ms").median().alias("median_ms"),
                pl.col("time_ms").mean().alias("mean_ms"),
                pl.col("time_ms").quantile(0.95).alias("p95_ms"),
            )
            .sort(["algorithm", "tool_label"])
        )
        perf_path = t_dir / "performance_summary.csv"
        perf_stats.write_csv(perf_path)

        # 2. Speedup by size bin (SR)
        df_sr = df_t.filter(pl.col("algorithm") == "sr")
        speedup_by_bin = compute_speedup_by_bin(df_sr, threads=threads)

        # Reorder by bin order
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
            speedup_path = t_dir / "speedup_by_bin_sr.csv"
            speedup_df.write_csv(speedup_path)

        # 3. Speedup by size bin (LR)
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
                speedup_path_lr = t_dir / "speedup_by_bin_lr.csv"
                speedup_df_lr.write_csv(speedup_path_lr)

        rprint(f"[green]Saved:[/green] {t_dir}/")

    rprint(f"\n[bold]Exported CSV files to {csv_dir}[/bold]")


def _plot_scatter(df_algo: pl.DataFrame, algo: str, ax):
    """Plot scatter for a single algorithm."""
    df_sampled = df_algo.sample(n=min(5000, df_algo.height), seed=42)

    # Exclude f32 from main plots (almost no difference from f64)
    tool_labels = [
        t for t in sorted(df_sampled["tool_label"].unique().to_list()) if "f32" not in t
    ]
    for tool_label in tool_labels:
        df_tool = df_sampled.filter(pl.col("tool_label") == tool_label)
        ax.scatter(
            df_tool["n_atoms"].to_list(),
            df_tool["time_ms"].to_list(),
            label=tool_label.replace("_", " ").title(),
            alpha=0.4,
            s=10,
            color=COLORS.get(tool_label, "#95a5a6"),
        )

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Number of Atoms")
    ax.set_ylabel("Execution Time (ms)")
    ax.legend()


@app.command()
def scatter():
    """Generate atoms vs time scatter plot."""
    setup_style()
    df = load_data()

    for algo in ["sr", "lr"]:
        plot_dir = PLOTS_DIR / "scatter" / algo
        individual_dir = plot_dir / "individual"
        individual_dir.mkdir(parents=True, exist_ok=True)

        df_algo = df.filter(pl.col("algorithm") == algo)
        thread_counts = sorted(df_algo["threads"].unique().to_list())

        # Individual plots per thread
        for threads in thread_counts:
            df_t = df_algo.filter(pl.col("threads") == threads)
            fig, ax = plt.subplots(figsize=(10, 6))
            _plot_scatter(df_t, algo, ax)
            ax.set_title(f"{algo.upper()}: Time vs Size (threads={threads})")
            fig.tight_layout()
            out_path = individual_dir / f"t{threads}.png"
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
            df_t = df_algo.filter(pl.col("threads") == threads)
            _plot_scatter(df_t, algo, ax)
            ax.set_title(f"threads={threads}")

        for idx in range(n_threads, n_rows * n_cols):
            row, col = idx // n_cols, idx % n_cols
            axes[row][col].set_visible(False)

        fig.suptitle(f"{algo.upper()}: Execution Time vs Structure Size", fontsize=14)
        fig.tight_layout()
        out_path = plot_dir / "grid.png"
        fig.savefig(out_path)
        plt.close(fig)
        rprint(f"[green]Saved:[/green] {out_path}")


def _plot_threads(df_algo: pl.DataFrame, algo: str, ax):
    """Plot thread scaling for a single algorithm."""
    scaling = (
        df_algo.group_by(["tool_label", "threads"])
        .agg(pl.col("time_ms").median().alias("median_time"))
        .sort(["tool_label", "threads"])
    )

    # Exclude f32 from main plots (almost no difference from f64)
    tool_labels = [
        t for t in sorted(scaling["tool_label"].unique().to_list()) if "f32" not in t
    ]
    for tool_label in tool_labels:
        df_tool = scaling.filter(pl.col("tool_label") == tool_label)
        ax.plot(
            df_tool["threads"].to_list(),
            df_tool["median_time"].to_list(),
            marker="o",
            label=tool_label.replace("_", " ").title(),
            color=COLORS.get(tool_label, "#95a5a6"),
            linestyle=LINESTYLES.get(tool_label, "-"),
            linewidth=2,
        )

    ax.set_xlabel("Thread Count")
    ax.set_ylabel("Median Execution Time (ms)")
    ax.legend()
    ax.grid(True, alpha=0.3)
    # Set x-axis to actual thread values
    thread_values = sorted(scaling["threads"].unique().to_list())
    ax.set_xticks(thread_values)
    ax.set_xticklabels([str(t) for t in thread_values])


@app.command()
def threads():
    """Generate thread scaling plot."""
    setup_style()
    plot_dir = PLOTS_DIR / "thread_scaling"
    individual_dir = plot_dir / "individual"
    individual_dir.mkdir(parents=True, exist_ok=True)

    df = load_data()

    # Individual plots
    for algo in ["sr", "lr"]:
        df_algo = df.filter(pl.col("algorithm") == algo)
        fig, ax = plt.subplots(figsize=(10, 6))
        _plot_threads(df_algo, algo, ax)
        ax.set_title(f"{algo.upper()}: Thread Scaling")
        fig.tight_layout()
        out_path = individual_dir / f"{algo}.png"
        fig.savefig(out_path)
        plt.close(fig)
        rprint(f"[green]Saved:[/green] {out_path}")

    # Combined grid
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    for idx, algo in enumerate(["sr", "lr"]):
        df_algo = df.filter(pl.col("algorithm") == algo)
        _plot_threads(df_algo, algo, axes[idx])
        axes[idx].set_title(f"{algo.upper()}")

    fig.suptitle("Thread Scaling", fontsize=14)
    fig.tight_layout()
    out_path = plot_dir / "grid.png"
    fig.savefig(out_path)
    plt.close(fig)
    rprint(f"[green]Saved:[/green] {out_path}")


# Bin definitions for size-based analysis
BINS = [
    (0, 500, "0-500"),
    (500, 1000, "500-1k"),
    (1000, 2000, "1k-2k"),
    (2000, 5000, "2k-5k"),
    (5000, 10000, "5k-10k"),
    (10000, 20000, "10k-20k"),
    (20000, 50000, "20k-50k"),
    (50000, 100000, "50k-100k"),
    (100000, 200000, "100k-200k"),
    (200000, float("inf"), "200k+"),
]


def add_size_bin(df: pl.DataFrame) -> pl.DataFrame:
    """Add size bin column to DataFrame."""
    expr = pl.lit("200k+")
    for low, high, label in reversed(BINS[:-1]):
        expr = pl.when(pl.col("n_atoms") < high).then(pl.lit(label)).otherwise(expr)
    return df.with_columns(expr.alias("size_bin"))


def compute_speedup_by_bin(df: pl.DataFrame, threads: int = 1) -> pl.DataFrame:
    """Compute speedup ratios by size bin for given thread count."""
    df_t = df.filter(pl.col("threads") == threads)

    pivot = (
        df_t.select(["structure", "tool_label", "n_atoms", "time_ms"])
        .pivot(on="tool_label", index=["structure", "n_atoms"], values="time_ms")
        .drop_nulls()
    )

    # Add speedup columns using tool_label (zig_f64, zig_f32, freesasa, rust)
    cols = pivot.columns

    # zig_f64 vs FreeSASA (primary comparison)
    if "zig_f64" in cols and "freesasa" in cols:
        pivot = pivot.with_columns(
            (pl.col("freesasa") / pl.col("zig_f64")).alias("zig_f64_vs_freesasa")
        )

    # zig_f64 vs Rust
    if "zig_f64" in cols and "rust" in cols:
        pivot = pivot.with_columns(
            (pl.col("rust") / pl.col("zig_f64")).alias("zig_f64_vs_rust")
        )

    # zig_f32 vs zig_f64 (f32 performance relative to f64)
    if "zig_f32" in cols and "zig_f64" in cols:
        pivot = pivot.with_columns(
            (pl.col("zig_f64") / pl.col("zig_f32")).alias("zig_f32_vs_zig_f64")
        )

    # FreeSASA vs Rust
    if "freesasa" in cols and "rust" in cols:
        pivot = pivot.with_columns(
            (pl.col("rust") / pl.col("freesasa")).alias("freesasa_vs_rust")
        )

    pivot = add_size_bin(pivot)

    # Aggregate by bin - median with IQR for error bars
    agg_cols = [pl.len().alias("count")]
    if "zig_f64_vs_freesasa" in pivot.columns:
        agg_cols.extend(
            [
                pl.col("zig_f64_vs_freesasa").median().alias("zig_f64_vs_freesasa"),
                pl.col("zig_f64_vs_freesasa")
                .quantile(0.25)
                .alias("zig_f64_vs_freesasa_q25"),
                pl.col("zig_f64_vs_freesasa")
                .quantile(0.75)
                .alias("zig_f64_vs_freesasa_q75"),
            ]
        )
    if "zig_f64_vs_rust" in pivot.columns:
        agg_cols.extend(
            [
                pl.col("zig_f64_vs_rust").median().alias("zig_f64_vs_rust"),
                pl.col("zig_f64_vs_rust").quantile(0.25).alias("zig_f64_vs_rust_q25"),
                pl.col("zig_f64_vs_rust").quantile(0.75).alias("zig_f64_vs_rust_q75"),
            ]
        )
    if "zig_f32_vs_zig_f64" in pivot.columns:
        agg_cols.extend(
            [
                pl.col("zig_f32_vs_zig_f64").median().alias("zig_f32_vs_zig_f64"),
                pl.col("zig_f32_vs_zig_f64")
                .quantile(0.25)
                .alias("zig_f32_vs_zig_f64_q25"),
                pl.col("zig_f32_vs_zig_f64")
                .quantile(0.75)
                .alias("zig_f32_vs_zig_f64_q75"),
            ]
        )
    if "freesasa_vs_rust" in pivot.columns:
        agg_cols.extend(
            [
                pl.col("freesasa_vs_rust").median().alias("freesasa_vs_rust"),
                pl.col("freesasa_vs_rust").quantile(0.25).alias("freesasa_vs_rust_q25"),
                pl.col("freesasa_vs_rust").quantile(0.75).alias("freesasa_vs_rust_q75"),
            ]
        )

    return pivot.group_by("size_bin").agg(agg_cols).sort("size_bin")


def _plot_speedup_single(
    df_sr: pl.DataFrame, threads: int, ax, show_legend: bool = True
):
    """Plot speedup for a single thread count on given axes."""
    bin_labels = [b[2] for b in BINS]
    speedup_data = compute_speedup_by_bin(df_sr, threads=threads)
    data_dict = {r["size_bin"]: r for r in speedup_data.iter_rows(named=True)}

    x_labels = [b for b in bin_labels if b in data_dict]
    x_pos = list(range(len(x_labels)))

    # zig_f64 vs FreeSASA (primary comparison)
    if "zig_f64_vs_freesasa" in speedup_data.columns:
        y_fs = [data_dict[label]["zig_f64_vs_freesasa"] for label in x_labels]
        y_fs_q25 = [data_dict[label]["zig_f64_vs_freesasa_q25"] for label in x_labels]
        y_fs_q75 = [data_dict[label]["zig_f64_vs_freesasa_q75"] for label in x_labels]
        yerr_fs = [
            [m - q25 for m, q25 in zip(y_fs, y_fs_q25)],
            [q75 - m for m, q75 in zip(y_fs, y_fs_q75)],
        ]
        ax.errorbar(
            x_pos,
            y_fs,
            yerr=yerr_fs,
            marker="o",
            linewidth=1.5,
            markersize=5,
            capsize=3,
            label="Zig(f64) vs FreeSASA",
            color="#3498db",
        )

    # zig_f64 vs Rust
    if "zig_f64_vs_rust" in speedup_data.columns:
        y_rust = [data_dict[label]["zig_f64_vs_rust"] for label in x_labels]
        y_rust_q25 = [data_dict[label]["zig_f64_vs_rust_q25"] for label in x_labels]
        y_rust_q75 = [data_dict[label]["zig_f64_vs_rust_q75"] for label in x_labels]
        yerr_rust = [
            [m - q25 for m, q25 in zip(y_rust, y_rust_q25)],
            [q75 - m for m, q75 in zip(y_rust, y_rust_q75)],
        ]
        ax.errorbar(
            x_pos,
            y_rust,
            yerr=yerr_rust,
            marker="s",
            linewidth=1.5,
            markersize=5,
            capsize=3,
            label="Zig(f64) vs Rust",
            color="#e74c3c",
        )

    # zig_f32 vs zig_f64 (omitted - almost no difference, see Appendix)
    # if "zig_f32_vs_zig_f64" in speedup_data.columns:
    #     y_f32 = [data_dict[label]["zig_f32_vs_zig_f64"] for label in x_labels]
    #     ...

    # FreeSASA vs Rust
    if "freesasa_vs_rust" in speedup_data.columns:
        y_fsr = [data_dict[label]["freesasa_vs_rust"] for label in x_labels]
        y_fsr_q25 = [data_dict[label]["freesasa_vs_rust_q25"] for label in x_labels]
        y_fsr_q75 = [data_dict[label]["freesasa_vs_rust_q75"] for label in x_labels]
        yerr_fsr = [
            [m - q25 for m, q25 in zip(y_fsr, y_fsr_q25)],
            [q75 - m for m, q75 in zip(y_fsr, y_fsr_q75)],
        ]
        ax.errorbar(
            x_pos,
            y_fsr,
            yerr=yerr_fsr,
            marker="D",
            linewidth=1.5,
            markersize=5,
            capsize=3,
            label="FreeSASA vs Rust",
            color="#9b59b6",
        )

    ax.axhline(y=1.0, color="gray", linestyle="--", linewidth=1)
    ax.set_xticks(x_pos)
    ax.set_xticklabels(x_labels, rotation=45, ha="right", fontsize=8)
    ax.set_ylim(0.5, 2.5)
    ax.grid(True, alpha=0.3)

    if show_legend:
        ax.legend(fontsize=8)


@app.command()
def grid():
    """Generate grid of speedup plots for all thread counts."""
    setup_style()
    plot_dir = PLOTS_DIR / "speedup_by_bin"
    individual_dir = plot_dir / "individual"
    individual_dir.mkdir(parents=True, exist_ok=True)

    df = load_data()
    df_sr = df.filter(pl.col("algorithm") == "sr")

    # Get available thread counts
    thread_counts = sorted(df_sr["threads"].unique().to_list())

    # Generate individual plots for each thread count
    for threads in thread_counts:
        fig_single, ax_single = plt.subplots(figsize=(10, 6))
        _plot_speedup_single(df_sr, threads, ax_single, show_legend=True)
        ax_single.set_xlabel("Structure Size (atoms)")
        ax_single.set_ylabel("Speedup Ratio (>1 = Zig faster)")
        ax_single.set_title(f"SR Algorithm: Zig Speedup (threads={threads})")
        fig_single.tight_layout()
        out_path = individual_dir / f"t{threads}.png"
        fig_single.savefig(out_path)
        plt.close(fig_single)
        rprint(f"[green]Saved:[/green] {out_path}")

    # Generate grid plot
    n_threads = len(thread_counts)
    n_cols = min(3, n_threads)
    n_rows = (n_threads + n_cols - 1) // n_cols

    fig, axes = plt.subplots(n_rows, n_cols, figsize=(5 * n_cols, 4 * n_rows))
    if n_threads == 1:
        axes = [[axes]]
    elif n_rows == 1:
        axes = [axes]

    for idx, threads in enumerate(thread_counts):
        row, col = idx // n_cols, idx % n_cols
        ax = axes[row][col] if n_rows > 1 else axes[0][col]
        _plot_speedup_single(df_sr, threads, ax, show_legend=(idx == 0))
        ax.set_title(f"threads={threads}")

    # Hide empty subplots
    for idx in range(n_threads, n_rows * n_cols):
        row, col = idx // n_cols, idx % n_cols
        axes[row][col].set_visible(False)

    fig.suptitle(
        "SR Algorithm: Zig Speedup by Structure Size and Thread Count", fontsize=14
    )
    fig.tight_layout()

    out_path = plot_dir / "grid.png"
    fig.savefig(out_path)
    plt.close(fig)
    rprint(f"[green]Saved:[/green] {out_path}")


@app.command()
def validation():
    """Generate SASA validation scatter plot (Zig f64 vs FreeSASA)."""
    setup_style()
    plot_dir = PLOTS_DIR / "validation"
    plot_dir.mkdir(parents=True, exist_ok=True)

    df = load_data()

    for algo in ["sr", "lr"]:
        df_algo = df.filter((pl.col("algorithm") == algo) & (pl.col("threads") == 1))

        # Pivot to get zig_f64 and freesasa SASA values side by side
        pivot = (
            df_algo.select(["structure", "tool_label", "total_sasa"])
            .pivot(on="tool_label", index="structure", values="total_sasa")
            .drop_nulls()
        )

        if "zig_f64" not in pivot.columns or "freesasa" not in pivot.columns:
            continue

        zig_sasa = pivot["zig_f64"].to_list()
        fs_sasa = pivot["freesasa"].to_list()

        fig, ax = plt.subplots(figsize=(8, 8))

        ax.scatter(zig_sasa, fs_sasa, alpha=0.3, s=10, color="#3498db")

        # y=x line
        max_val = max(max(zig_sasa), max(fs_sasa))
        ax.plot([0, max_val], [0, max_val], "r--", linewidth=1.5, label="y = x")

        # Calculate R² and max relative error
        import numpy as np

        zig_arr = np.array(zig_sasa)
        fs_arr = np.array(fs_sasa)
        correlation = np.corrcoef(zig_arr, fs_arr)[0, 1]
        r_squared = correlation**2

        rel_errors = np.abs(zig_arr - fs_arr) / fs_arr * 100
        max_rel_error = np.max(rel_errors)
        mean_rel_error = np.mean(rel_errors)

        ax.set_xlabel("Zig(f64) SASA (Å²)")
        ax.set_ylabel("FreeSASA SASA (Å²)")
        ax.set_title(f"{algo.upper()}: SASA Validation (n={len(zig_sasa):,})")
        ax.legend()

        # Add stats annotation
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

        out_path = plot_dir / f"{algo}.png"
        fig.savefig(out_path)
        plt.close(fig)
        rprint(f"[green]Saved:[/green] {out_path}")


@app.command()
def samples():
    """Generate thread scaling plots for representative structures per size bin."""
    setup_style()
    df = load_data()
    df_sr = df.filter(pl.col("algorithm") == "sr")
    df_sr = add_size_bin(df_sr)

    plot_dir = PLOTS_DIR / "samples"
    plot_dir.mkdir(parents=True, exist_ok=True)

    bin_order = [b[2] for b in BINS]

    # Get thread counts from data
    thread_counts = sorted(df_sr["threads"].unique().to_list())

    for bin_name in bin_order:
        df_bin = df_sr.filter(pl.col("size_bin") == bin_name)
        if df_bin.height == 0:
            continue

        # Pick up to 3 representative structures (small, medium, large within bin)
        struct_sizes = (
            df_bin.filter(pl.col("threads") == 1)
            .select(["structure", "n_atoms"])
            .unique()
            .sort("n_atoms")
        )
        n = struct_sizes.height
        if n <= 3:
            selected = struct_sizes["structure"].to_list()
        else:
            # Pick first, middle, last
            indices = [0, n // 2, n - 1]
            selected = [struct_sizes["structure"][i] for i in indices]

        if not selected:
            continue

        # Create subplot for this bin
        n_structs = len(selected)
        n_cols = min(3, n_structs)
        n_rows = (n_structs + n_cols - 1) // n_cols

        fig, axes = plt.subplots(
            n_rows, n_cols, figsize=(5 * n_cols, 4 * n_rows), squeeze=False
        )

        for idx, struct in enumerate(selected):
            row, col = idx // n_cols, idx % n_cols
            ax = axes[row][col]

            df_struct = df_bin.filter(pl.col("structure") == struct)
            n_atoms = df_struct["n_atoms"][0]

            for tool_label in sorted(df_struct["tool_label"].unique().to_list()):
                df_tool = df_struct.filter(pl.col("tool_label") == tool_label).sort(
                    "threads"
                )
                ax.plot(
                    df_tool["threads"].to_list(),
                    df_tool["time_ms"].to_list(),
                    marker="o",
                    label=tool_label.replace("_", " ").title(),
                    color=COLORS.get(tool_label, "#95a5a6"),
                    linestyle=LINESTYLES.get(tool_label, "-"),
                    linewidth=1.5,
                )

            ax.set_xlabel("Threads")
            ax.set_ylabel("Time (ms)")
            ax.set_xticks(thread_counts)
            ax.set_title(f"{struct}\n({n_atoms:,} atoms)", fontsize=10)
            ax.grid(True, alpha=0.3)
            if idx == 0:
                ax.legend(fontsize=8)

        # Hide empty subplots
        for idx in range(n_structs, n_rows * n_cols):
            row, col = idx // n_cols, idx % n_cols
            axes[row][col].set_visible(False)

        fig.suptitle(f"SR Thread Scaling: {bin_name} atoms", fontsize=14)
        fig.tight_layout()

        # Sanitize bin name for filename
        safe_name = bin_name.replace("+", "plus")
        out_path = plot_dir / f"{safe_name}.png"
        fig.savefig(out_path)
        plt.close(fig)
        rprint(f"[green]Saved:[/green] {out_path}")

    # Special plot for max size structure
    max_struct_row = (
        df_sr.filter(pl.col("threads") == 1).sort("n_atoms", descending=True).head(1)
    )
    if max_struct_row.height > 0:
        max_struct = max_struct_row["structure"][0]
        max_atoms = max_struct_row["n_atoms"][0]

        df_max = df_sr.filter(pl.col("structure") == max_struct)

        fig, ax = plt.subplots(figsize=(8, 6))

        for tool_label in sorted(df_max["tool_label"].unique().to_list()):
            df_tool = df_max.filter(pl.col("tool_label") == tool_label).sort("threads")
            ax.plot(
                df_tool["threads"].to_list(),
                df_tool["time_ms"].to_list(),
                marker="o",
                label=tool_label.replace("_", " ").title(),
                color=COLORS.get(tool_label, "#95a5a6"),
                linestyle=LINESTYLES.get(tool_label, "-"),
                linewidth=2,
                markersize=8,
            )

        ax.set_xlabel("Threads")
        ax.set_ylabel("Time (ms)")
        ax.set_xticks(thread_counts)
        ax.set_title(f"SR Thread Scaling: {max_struct} ({max_atoms:,} atoms)")
        ax.grid(True, alpha=0.3)
        ax.legend()

        fig.tight_layout()
        out_path = plot_dir / "max_structure.png"
        fig.savefig(out_path)
        plt.close(fig)
        rprint(f"[green]Saved:[/green] {out_path}")


@app.command()
def large():
    """Generate speedup bar chart for large structures (100k+ atoms) at threads=10."""
    setup_style()
    plot_dir = PLOTS_DIR / "large"
    plot_dir.mkdir(parents=True, exist_ok=True)

    df = load_data()
    df = add_size_bin(df)

    # Filter to large structures (100k+) and threads=10
    large_bins = {"100k-200k", "200k+"}
    df_large = df.filter(
        pl.col("size_bin").is_in(large_bins) & (pl.col("threads") == 10)
    )

    if df_large.height == 0:
        rprint("[yellow]No large structures found[/yellow]")
        return

    results = []
    colors_algo = {"sr": "#3498db", "lr": "#e74c3c"}  # SR=blue, LR=red

    # Collect data per algorithm using tool_label
    pivots = {}
    for algo in ["sr", "lr"]:
        df_algo = df_large.filter(pl.col("algorithm") == algo)
        if df_algo.height == 0:
            continue
        pivots[algo] = (
            df_algo.select(["structure", "tool_label", "time_ms"])
            .pivot(on="tool_label", index="structure", values="time_ms")
            .drop_nulls()
        )

    # Build results in reverse order (bottom to top in chart)
    # LR section first (will appear at bottom)
    lr_count = 0

    # For LR, zig_lr directory has no precision suffix, so tool_label is just "zig"
    zig_lr_col = (
        "zig_f64" if "zig_f64" in pivots.get("lr", pl.DataFrame()).columns else "zig"
    )
    if (
        "lr" in pivots
        and zig_lr_col in pivots["lr"].columns
        and "freesasa" in pivots["lr"].columns
    ):
        speedup = pivots["lr"]["freesasa"] / pivots["lr"][zig_lr_col]
        results.append(
            {
                "label": "vs FreeSASA (LR)",
                "speedup": speedup.median(),
                "q25": speedup.quantile(0.25),
                "q75": speedup.quantile(0.75),
                "color": colors_algo["lr"],
            }
        )
        lr_count += 1

    # SR section (will appear at top)
    if (
        "sr" in pivots
        and "zig_f64" in pivots["sr"].columns
        and "rust" in pivots["sr"].columns
    ):
        speedup = pivots["sr"]["rust"] / pivots["sr"]["zig_f64"]
        results.append(
            {
                "label": "vs Rust (SR)",
                "speedup": speedup.median(),
                "q25": speedup.quantile(0.25),
                "q75": speedup.quantile(0.75),
                "color": colors_algo["sr"],
            }
        )

    if (
        "sr" in pivots
        and "zig_f64" in pivots["sr"].columns
        and "freesasa" in pivots["sr"].columns
    ):
        speedup = pivots["sr"]["freesasa"] / pivots["sr"]["zig_f64"]
        results.append(
            {
                "label": "vs FreeSASA (SR)",
                "speedup": speedup.median(),
                "q25": speedup.quantile(0.25),
                "q75": speedup.quantile(0.75),
                "color": colors_algo["sr"],
            }
        )

    if not results:
        rprint("[yellow]No comparison data available[/yellow]")
        return

    # Plot horizontal bar chart
    fig, ax = plt.subplots(figsize=(10, 6))

    labels = [r["label"] for r in results]
    speedups = [r["speedup"] for r in results]
    colors = [r["color"] for r in results]
    xerr = [
        [r["speedup"] - r["q25"] for r in results],
        [r["q75"] - r["speedup"] for r in results],
    ]

    y_pos = range(len(labels))
    ax.barh(y_pos, speedups, xerr=xerr, color=colors, capsize=4, height=0.6)

    # Add separator between LR and SR sections
    if lr_count > 0 and lr_count < len(results):
        ax.axhline(y=lr_count - 0.5, color="gray", linestyle="-", linewidth=0.8)

    ax.set_yticks(y_pos)
    ax.set_yticklabels(labels)
    ax.set_xlabel("Speedup (median, IQR)")
    ax.axvline(x=1.0, color="gray", linestyle="--", linewidth=1.5)

    # Add value labels
    for i, s in enumerate(speedups):
        ax.text(s + 0.05, i, f"{s:.2f}x", va="center", fontsize=10)

    n_structs = (
        df_large.filter(pl.col("algorithm") == "sr").select("structure").unique().height
    )
    ax.set_title(
        f"Zig f64 Speedup on Large Structures (100k+ atoms, n={n_structs}, threads=10)"
    )
    ax.set_xlim(0, max(speedups) * 1.2)
    ax.grid(True, alpha=0.3, axis="x")

    fig.tight_layout()
    out_path = plot_dir / "speedup_bar.png"
    fig.savefig(out_path)
    plt.close(fig)
    rprint(f"[green]Saved:[/green] {out_path}")


@app.command()
def efficiency():
    """Calculate and plot parallel efficiency from existing benchmark data."""
    setup_style()
    df = load_data()
    df = add_size_bin(df)

    plot_dir = PLOTS_DIR / "efficiency"
    plot_dir.mkdir(parents=True, exist_ok=True)

    # Focus on SR algorithm
    df_sr = df.filter(pl.col("algorithm") == "sr")

    # Get threads=1 baseline for each structure/tool_label
    df_t1 = df_sr.filter(pl.col("threads") == 1).select(
        [
            "tool_label",
            "structure",
            "n_atoms",
            "size_bin",
            pl.col("time_ms").alias("t1_ms"),
        ]
    )

    # Join with all thread counts to calculate efficiency
    df_eff = df_sr.join(df_t1, on=["tool_label", "structure", "n_atoms", "size_bin"])
    df_eff = df_eff.with_columns(
        (pl.col("t1_ms") / (pl.col("time_ms") * pl.col("threads"))).alias("efficiency")
    )

    # Summary by size bin and tool_label
    thread_counts = sorted(df_eff["threads"].unique().to_list())
    # Use tool_labels that exist in data (exclude f32 - almost no difference from f64)
    tool_labels = ["zig_f64", "freesasa", "rust"]

    # Table: efficiency by size bin (threads=10)
    t_target = 10
    df_t10 = df_eff.filter(pl.col("threads") == t_target)

    summary_table = Table(title=f"Parallel Efficiency by Size (SR, threads={t_target})")
    summary_table.add_column("Size Bin", style="cyan")
    summary_table.add_column("Count", justify="right")
    summary_table.add_column("Zig(f64)", justify="right")
    summary_table.add_column("FreeSASA", justify="right")
    summary_table.add_column("Rust", justify="right")

    bin_order = [b[2] for b in BINS]

    for bin_name in bin_order:
        df_bin = df_t10.filter(pl.col("size_bin") == bin_name)
        if df_bin.height == 0:
            continue

        count = df_bin.filter(pl.col("tool_label") == "zig_f64").height
        row = [bin_name, f"{count:,}"]

        for tool_label in tool_labels:
            df_tool = df_bin.filter(pl.col("tool_label") == tool_label)
            if df_tool.height > 0:
                eff = df_tool["efficiency"].median()
                if eff >= 0.3:
                    row.append(f"[green]{eff:.2f}[/green]")
                elif eff >= 0.2:
                    row.append(f"[yellow]{eff:.2f}[/yellow]")
                else:
                    row.append(f"[red]{eff:.2f}[/red]")
            else:
                row.append("-")

        summary_table.add_row(*row)

    rprint(summary_table)
    rprint(
        "\n[dim]Efficiency = T1 / (TN × N). Higher = better thread utilization[/dim]"
    )

    # Export CSV for all thread counts
    csv_dir = RESULTS_DIR / "csv" / "efficiency"
    csv_dir.mkdir(parents=True, exist_ok=True)

    for threads in thread_counts:
        df_t = df_eff.filter(pl.col("threads") == threads)
        rows = []
        for bin_name in bin_order:
            df_bin = df_t.filter(pl.col("size_bin") == bin_name)
            if df_bin.height == 0:
                continue
            row = {"size_bin": bin_name}
            for tool_label in tool_labels:
                df_tool = df_bin.filter(pl.col("tool_label") == tool_label)
                if df_tool.height > 0:
                    row[f"{tool_label}_median"] = df_tool["efficiency"].median()
                    row[f"{tool_label}_count"] = df_tool.height
            rows.append(row)
        if rows:
            pl.DataFrame(rows).write_csv(csv_dir / f"t{threads}.csv")

    rprint(f"[green]Saved CSV:[/green] {csv_dir}/")

    # Plot: average efficiency by thread count
    fig, ax = plt.subplots(figsize=(10, 6))

    for tool_label in tool_labels:
        df_tool = df_eff.filter(pl.col("tool_label") == tool_label)
        if df_tool.height == 0:
            continue

        avg_by_thread = (
            df_tool.group_by("threads")
            .agg(pl.col("efficiency").median().alias("eff"))
            .sort("threads")
        )

        ax.plot(
            avg_by_thread["threads"].to_list(),
            avg_by_thread["eff"].to_list(),
            marker="o",
            label=tool_label.replace("_", " ").title(),
            color=COLORS.get(tool_label, "#95a5a6"),
            linestyle=LINESTYLES.get(tool_label, "-"),
            linewidth=2,
            markersize=8,
        )

    ax.axhline(y=1.0, color="gray", linestyle="--", alpha=0.5, label="Ideal")
    ax.set_xlabel("Threads")
    ax.set_ylabel("Parallel Efficiency (median)")
    ax.set_title("Parallel Efficiency Comparison (SR)")
    ax.set_xticks(thread_counts)
    ax.set_ylim(0, 1.1)
    ax.legend()
    ax.grid(True, alpha=0.3)

    fig.tight_layout()
    out_path = plot_dir / "summary.png"
    fig.savefig(out_path)
    plt.close(fig)
    rprint(f"\n[green]Saved:[/green] {out_path}")

    # Plot: efficiency by size bin (comparing tools at threads=10)
    fig, ax = plt.subplots(figsize=(14, 6))

    x_labels = []
    x_pos = []
    pos = 0
    width = 0.2  # Narrower bars to fit more tools

    for bin_name in bin_order:
        df_bin = df_t10.filter(pl.col("size_bin") == bin_name)
        if df_bin.height == 0:
            continue

        x_labels.append(bin_name)
        for i, tool_label in enumerate(tool_labels):
            df_tool = df_bin.filter(pl.col("tool_label") == tool_label)
            if df_tool.height > 0:
                eff = df_tool["efficiency"].median()
                ax.bar(
                    pos + i * width,
                    eff,
                    width,
                    label=tool_label.replace("_", " ").title()
                    if len(x_pos) == 0
                    else "",
                    color=COLORS.get(tool_label, "#95a5a6"),
                )
        x_pos.append(pos + width * 1.5)  # Center of 4 bars
        pos += 1

    ax.set_xticks(x_pos)
    ax.set_xticklabels(x_labels, rotation=45, ha="right")
    ax.set_ylabel("Parallel Efficiency")
    ax.set_title(f"Parallel Efficiency by Structure Size (SR, threads={t_target})")
    ax.legend()
    ax.grid(True, alpha=0.3, axis="y")

    fig.tight_layout()
    out_path = plot_dir / "by_size.png"
    fig.savefig(out_path)
    plt.close(fig)
    rprint(f"[green]Saved:[/green] {out_path}")


@app.command()
def speedup(
    min_atoms: int = typer.Option(50000, help="Minimum atom count for filtering"),
    top_n: int = typer.Option(5, help="Number of top entries to show"),
):
    """Find structures with best Zig speedup at any thread count."""
    setup_style()
    df = load_data()
    df = add_size_bin(df)

    # Filter to SR algorithm and large structures
    df_sr = df.filter((pl.col("algorithm") == "sr") & (pl.col("n_atoms") >= min_atoms))

    # Pivot to get tools side by side
    pivot = (
        df_sr.select(["structure", "n_atoms", "threads", "tool_label", "time_ms"])
        .pivot(
            on="tool_label", index=["structure", "n_atoms", "threads"], values="time_ms"
        )
        .drop_nulls()
    )

    # Calculate speedup
    if "zig_f64" not in pivot.columns:
        rprint("[red]No zig_f64 data found[/red]")
        return

    speedup_cols = []
    if "freesasa" in pivot.columns:
        pivot = pivot.with_columns(
            (pl.col("freesasa") / pl.col("zig_f64")).alias("vs_freesasa")
        )
        speedup_cols.append("vs_freesasa")
    if "rust" in pivot.columns:
        pivot = pivot.with_columns(
            (pl.col("rust") / pl.col("zig_f64")).alias("vs_rust")
        )
        speedup_cols.append("vs_rust")

    plot_dir = PLOTS_DIR / "speedup"
    plot_dir.mkdir(parents=True, exist_ok=True)

    # Print and collect top entries for each comparison
    results = {}

    for col in speedup_cols:
        comparison = col.replace("vs_", "")
        top_entries = pivot.sort(col, descending=True).head(top_n)

        table = Table(
            title=f"Top {top_n} Zig Speedup vs {comparison.title()} ({min_atoms:,}+ atoms)"
        )
        table.add_column("Rank", style="dim")
        table.add_column("Structure", style="cyan")
        table.add_column("Atoms", justify="right")
        table.add_column("Threads", justify="right")
        table.add_column("Speedup", justify="right", style="green")

        entries = []
        for i, row in enumerate(top_entries.iter_rows(named=True), 1):
            table.add_row(
                str(i),
                row["structure"],
                f"{row['n_atoms']:,}",
                str(row["threads"]),
                f"{row[col]:.2f}x",
            )
            entries.append(
                {
                    "structure": row["structure"],
                    "n_atoms": row["n_atoms"],
                    "threads": row["threads"],
                    "speedup": row[col],
                }
            )

        rprint(table)
        rprint()
        results[comparison] = entries

    # Check for Rust outliers (unusually high speedup that might be measurement noise)
    if "vs_rust" in speedup_cols:
        # Find entries where vs_rust is much higher than vs_freesasa
        outlier_check = pivot.filter(
            (pl.col("vs_rust") > 5.0) & (pl.col("n_atoms") < 50000)
        ).sort("vs_rust", descending=True)

        if outlier_check.height > 0:
            rprint(
                "[yellow]Note: Rust has outliers in small structures (likely measurement noise):[/yellow]"
            )
            for row in outlier_check.head(3).iter_rows(named=True):
                rprint(
                    f"  {row['structure']}: {row['vs_rust']:.1f}x @ threads={row['threads']} ({row['n_atoms']:,} atoms)"
                )
            rprint()

    # Generate comparison plot (similar to user's image)
    if len(results) < 2:
        rprint("[yellow]Need both FreeSASA and Rust data for comparison plot[/yellow]")
        return

    # Get best structure for each comparison
    best_fs = results["freesasa"][0]
    best_rust = results["rust"][0]

    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    for idx, (comparison, best) in enumerate(
        [("freesasa", best_fs), ("rust", best_rust)]
    ):
        ax = axes[idx]
        struct = best["structure"]
        n_atoms = best["n_atoms"]

        # Get all thread data for this structure
        df_struct = pivot.filter(pl.col("structure") == struct).sort("threads")

        threads_list = df_struct["threads"].to_list()
        vs_fs = df_struct["vs_freesasa"].to_list()
        vs_rust = df_struct["vs_rust"].to_list()

        ax.plot(
            threads_list,
            vs_fs,
            marker="o",
            label="vs FreeSASA",
            color="#3498db",
            linewidth=2,
            markersize=8,
        )
        ax.plot(
            threads_list,
            vs_rust,
            marker="s",
            label="vs Rust",
            color="#e74c3c",
            linewidth=2,
            markersize=8,
        )

        ax.axhline(y=1.0, color="gray", linestyle="--", linewidth=1)
        ax.set_xlabel("Threads")
        ax.set_ylabel("Speedup (nx)")
        ax.set_xticks(threads_list)
        ax.legend()
        ax.grid(True, alpha=0.3)

        # Find best speedup for this structure
        best_fs_val = max(vs_fs)
        best_fs_t = threads_list[vs_fs.index(best_fs_val)]
        best_rust_val = max(vs_rust)
        best_rust_t = threads_list[vs_rust.index(best_rust_val)]

        # Annotate endpoint values (at threads=10)
        ax.annotate(
            f"{vs_fs[-1]:.2f}x",
            xy=(threads_list[-1], vs_fs[-1]),
            xytext=(5, 0),
            textcoords="offset points",
            color="#3498db",
            fontsize=10,
        )
        ax.annotate(
            f"{vs_rust[-1]:.2f}x",
            xy=(threads_list[-1], vs_rust[-1]),
            xytext=(5, 0),
            textcoords="offset points",
            color="#e74c3c",
            fontsize=10,
        )

        # Title with best info
        if comparison == "freesasa":
            title = f"{struct} ({n_atoms:,} atoms)\nBest vs FreeSASA: {best_fs_val:.2f}x @ threads={best_fs_t}"
        else:
            title = f"{struct} ({n_atoms:,} atoms)\nBest vs Rust (large): {best_rust_val:.2f}x @ threads={best_rust_t}"
        ax.set_title(title)

    fig.suptitle("Zig Speedup vs FreeSASA and Rust (SR Algorithm)", fontsize=14)
    fig.tight_layout()

    out_path = plot_dir / "comparison.png"
    fig.savefig(out_path)
    plt.close(fig)
    rprint(f"[green]Saved:[/green] {out_path}")


@app.command()
def all():
    """Generate all plots and summary."""
    summary()
    rprint("\n[bold]Generating plots...[/bold]\n")
    validation()
    scatter()
    threads()
    grid()
    samples()
    large()
    efficiency()
    speedup()
    rprint(f"\n[bold green]All plots saved to:[/bold green] {PLOTS_DIR}")


if __name__ == "__main__":
    app()
