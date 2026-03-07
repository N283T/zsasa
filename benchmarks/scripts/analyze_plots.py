"""Plotting functions for SR benchmark analysis."""

import matplotlib.pyplot as plt
import polars as pl
from rich import print as rprint

from analyze_data import (
    BINS,
    COLORS,
    LINESTYLES,
    PLOTS_DIR,
    add_size_bin,
    compute_speedup_by_bin,
    display_name,
    load_data,
    metric_label,
    setup_style,
)


# === Internal Plot Helpers ===


def _plot_scatter(df_sr: pl.DataFrame, ax, time_col: str = "time_ms"):
    """Plot scatter for SR algorithm."""
    df_sampled = df_sr.sample(n=min(5000, df_sr.height), seed=42)

    tool_labels = [
        t for t in sorted(df_sampled["tool_label"].unique().to_list()) if "f32" not in t
    ]
    for tool_label in tool_labels:
        df_tool = df_sampled.filter(pl.col("tool_label") == tool_label)
        ax.scatter(
            df_tool["n_atoms"].to_list(),
            df_tool[time_col].to_list(),
            label=display_name(tool_label),
            alpha=0.4,
            s=10,
            color=COLORS.get(tool_label, "#95a5a6"),
        )

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Number of Atoms")
    ax.set_ylabel(f"{metric_label(time_col)} (ms)")
    ax.legend()


def _plot_threads(df_sr: pl.DataFrame, ax, time_col: str = "time_ms"):
    """Plot thread scaling for SR algorithm."""
    scaling = (
        df_sr.group_by(["tool_label", "threads"])
        .agg(pl.col(time_col).median().alias("median_time"))
        .sort(["tool_label", "threads"])
    )

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
    ax.set_ylabel(f"Median {metric_label(time_col)} (ms)")
    ax.legend()
    ax.grid(True, alpha=0.3)
    thread_values = sorted(scaling["threads"].unique().to_list())
    ax.set_xticks(thread_values)
    ax.set_xticklabels([str(t) for t in thread_values])


def _plot_speedup_single(
    df_sr: pl.DataFrame,
    threads: int,
    ax,
    show_legend: bool = True,
    time_col: str = "time_ms",
):
    """Plot speedup for a single thread count on given axes."""
    bin_labels = [b[2] for b in BINS]
    speedup_data = compute_speedup_by_bin(df_sr, threads=threads, time_col=time_col)
    data_dict = {r["size_bin"]: r for r in speedup_data.iter_rows(named=True)}

    x_labels = [b for b in bin_labels if b in data_dict]
    x_pos = list(range(len(x_labels)))

    comparisons = [
        ("zsasa_f64_vs_freesasa", "zsasa(f64) vs FreeSASA", "o", "#3498db"),
        ("zsasa_f64_vs_rustsasa", "zsasa(f64) vs RustSASA", "s", "#e74c3c"),
        ("freesasa_vs_rustsasa", "FreeSASA vs RustSASA", "D", "#9b59b6"),
    ]

    for col_name, label, marker, color in comparisons:
        if col_name not in speedup_data.columns:
            continue
        y = [data_dict[b][col_name] for b in x_labels]
        y_q25 = [data_dict[b][f"{col_name}_q25"] for b in x_labels]
        y_q75 = [data_dict[b][f"{col_name}_q75"] for b in x_labels]
        yerr = [
            [m - q25 for m, q25 in zip(y, y_q25)],
            [q75 - m for m, q75 in zip(y, y_q75)],
        ]
        ax.errorbar(
            x_pos,
            y,
            yerr=yerr,
            marker=marker,
            linewidth=1.5,
            markersize=5,
            capsize=3,
            label=label,
            color=color,
        )

    ax.axhline(y=1.0, color="gray", linestyle="--", linewidth=1)
    ax.set_xticks(x_pos)
    ax.set_xticklabels(x_labels, rotation=45, ha="right", fontsize=7)
    ax.set_ylim(0.5, 3.0)
    ax.grid(True, alpha=0.3)

    if show_legend:
        ax.legend(fontsize=8)


# === Plot Command Implementations ===


def plot_scatter(n_points: int = 100, time_col: str = "time_ms"):
    """Generate atoms vs time scatter plot."""
    setup_style()
    df = load_data(n_points)

    suffix = "_sasa" if time_col == "sasa_time_ms" else ""
    plot_dir = PLOTS_DIR.joinpath("scatter", f"sr{suffix}")
    individual_dir = plot_dir.joinpath("individual")
    individual_dir.mkdir(parents=True, exist_ok=True)

    thread_counts = sorted(df["threads"].unique().to_list())
    mlabel = metric_label(time_col)

    for threads in thread_counts:
        df_t = df.filter(pl.col("threads") == threads)
        fig, ax = plt.subplots(figsize=(10, 6))
        _plot_scatter(df_t, ax, time_col=time_col)
        ax.set_title(f"SR: {mlabel} vs Size (threads={threads})")
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
        _plot_scatter(df_t, ax, time_col=time_col)
        ax.set_title(f"threads={threads}")

    for idx in range(n_threads, n_rows * n_cols):
        row, col = idx // n_cols, idx % n_cols
        axes[row][col].set_visible(False)

    fig.suptitle(f"SR: {mlabel} vs Structure Size", fontsize=14)
    fig.tight_layout()
    out_path = plot_dir.joinpath("grid.png")
    fig.savefig(out_path)
    plt.close(fig)
    rprint(f"[green]Saved:[/green] {out_path}")


def plot_threads(n_points: int = 100, time_col: str = "time_ms"):
    """Generate thread scaling plot."""
    setup_style()
    df = load_data(n_points)

    suffix = "_sasa" if time_col == "sasa_time_ms" else ""
    plot_dir = PLOTS_DIR.joinpath("thread_scaling")
    plot_dir.mkdir(parents=True, exist_ok=True)

    fig, ax = plt.subplots(figsize=(10, 6))
    _plot_threads(df, ax, time_col=time_col)
    ax.set_title(f"SR: Thread Scaling ({metric_label(time_col)})")
    fig.tight_layout()
    out_path = plot_dir.joinpath(f"sr{suffix}.png")
    fig.savefig(out_path)
    plt.close(fig)
    rprint(f"[green]Saved:[/green] {out_path}")


def plot_grid(n_points: int = 100, time_col: str = "time_ms"):
    """Generate grid of speedup plots for all thread counts."""
    setup_style()
    df = load_data(n_points)

    if df.height == 0:
        rprint("[yellow]No SR data found[/yellow]")
        return

    suffix = "_sasa" if time_col == "sasa_time_ms" else ""
    mlabel = metric_label(time_col)
    plot_dir = PLOTS_DIR.joinpath(f"speedup_by_bin{suffix}")
    individual_dir = plot_dir.joinpath("individual")
    individual_dir.mkdir(parents=True, exist_ok=True)

    thread_counts = sorted(df["threads"].unique().to_list())

    for threads in thread_counts:
        fig_single, ax_single = plt.subplots(figsize=(12, 6))
        _plot_speedup_single(
            df, threads, ax_single, show_legend=True, time_col=time_col
        )
        ax_single.set_xlabel("Structure Size (atoms)")
        ax_single.set_ylabel("Speedup Ratio (>1 = zsasa faster)")
        ax_single.set_title(f"SR: zsasa Speedup [{mlabel}] (threads={threads})")
        fig_single.tight_layout()
        out_path = individual_dir.joinpath(f"t{threads}.png")
        fig_single.savefig(out_path)
        plt.close(fig_single)
        rprint(f"[green]Saved:[/green] {out_path}")

    n_threads = len(thread_counts)
    n_cols = min(3, n_threads)
    n_rows = (n_threads + n_cols - 1) // n_cols

    fig, axes = plt.subplots(n_rows, n_cols, figsize=(6 * n_cols, 4.5 * n_rows))
    if n_threads == 1:
        axes = [[axes]]
    elif n_rows == 1:
        axes = [axes]

    for idx, threads in enumerate(thread_counts):
        row, col = idx // n_cols, idx % n_cols
        ax = axes[row][col] if n_rows > 1 else axes[0][col]
        _plot_speedup_single(df, threads, ax, show_legend=(idx == 0), time_col=time_col)
        ax.set_title(f"threads={threads}")

    for idx in range(n_threads, n_rows * n_cols):
        row, col = idx // n_cols, idx % n_cols
        axes[row][col].set_visible(False)

    fig.suptitle(f"SR: zsasa Speedup by Size/Threads [{mlabel}]", fontsize=14)
    fig.tight_layout()

    out_path = plot_dir.joinpath("grid.png")
    fig.savefig(out_path)
    plt.close(fig)
    rprint(f"[green]Saved:[/green] {out_path}")


def plot_validation(n_points: int = 100):
    """Generate SASA validation scatter plot (zsasa f64 vs FreeSASA)."""
    import numpy as np

    setup_style()
    df = load_data(n_points)

    plot_dir = PLOTS_DIR.joinpath("validation")
    plot_dir.mkdir(parents=True, exist_ok=True)

    df_t1 = df.filter(pl.col("threads") == 1)
    if df_t1.height == 0:
        rprint("[yellow]No single-threaded SR data found[/yellow]")
        return

    if "total_sasa" not in df_t1.columns or df_t1["total_sasa"].is_null().all():
        rprint("[yellow]No total_sasa data available for validation[/yellow]")
        return

    pivot = (
        df_t1.select(["structure", "tool_label", "total_sasa"])
        .pivot(on="tool_label", index="structure", values="total_sasa")
        .drop_nulls()
    )

    # Try zig_f64 first, fall back to zig
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
    ax.set_title(f"SR: SASA Validation (n={len(zsasa_sasa):,})")
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

    out_path = plot_dir.joinpath("sr.png")
    fig.savefig(out_path)
    plt.close(fig)
    rprint(f"[green]Saved:[/green] {out_path}")


def plot_samples(n_points: int = 100, time_col: str = "time_ms"):
    """Generate thread scaling plots for representative structures per size bin."""
    setup_style()
    df = load_data(n_points)

    if df.height == 0:
        rprint("[yellow]No SR data found[/yellow]")
        return

    df = add_size_bin(df)

    suffix = "_sasa" if time_col == "sasa_time_ms" else ""
    plot_dir = PLOTS_DIR.joinpath(f"samples{suffix}")
    plot_dir.mkdir(parents=True, exist_ok=True)

    bin_order = [b[2] for b in BINS]
    thread_counts = sorted(df["threads"].unique().to_list())

    for bin_name in bin_order:
        df_bin = df.filter(pl.col("size_bin") == bin_name)
        if df_bin.height == 0:
            continue

        struct_sizes = (
            df_bin.filter(pl.col("threads") == thread_counts[0])
            .select(["structure", "n_atoms"])
            .unique()
            .sort("n_atoms")
        )
        n = struct_sizes.height
        if n <= 3:
            selected = struct_sizes["structure"].to_list()
        else:
            indices = [0, n // 2, n - 1]
            selected = [struct_sizes["structure"][i] for i in indices]

        if not selected:
            continue

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
                    df_tool[time_col].to_list(),
                    marker="o",
                    label=display_name(tool_label),
                    color=COLORS.get(tool_label, "#95a5a6"),
                    linestyle=LINESTYLES.get(tool_label, "-"),
                    linewidth=1.5,
                )

            ax.set_xlabel("Threads")
            ax.set_ylabel(f"{metric_label(time_col)} (ms)")
            ax.set_xticks(thread_counts)
            ax.set_title(f"{struct}\n({n_atoms:,} atoms)", fontsize=10)
            ax.grid(True, alpha=0.3)
            if idx == 0:
                ax.legend(fontsize=8)

        for idx in range(n_structs, n_rows * n_cols):
            row, col = idx // n_cols, idx % n_cols
            axes[row][col].set_visible(False)

        fig.suptitle(f"SR Thread Scaling: {bin_name} atoms", fontsize=14)
        fig.tight_layout()

        safe_name = bin_name.replace("+", "plus")
        out_path = plot_dir.joinpath(f"{safe_name}.png")
        fig.savefig(out_path)
        plt.close(fig)
        rprint(f"[green]Saved:[/green] {out_path}")

    # Max size structure plot
    max_struct_row = (
        df.filter(pl.col("threads") == thread_counts[0])
        .sort("n_atoms", descending=True)
        .head(1)
    )
    if max_struct_row.height > 0:
        max_struct = max_struct_row["structure"][0]
        max_atoms = max_struct_row["n_atoms"][0]
        df_max = df.filter(pl.col("structure") == max_struct)

        fig, ax = plt.subplots(figsize=(8, 6))

        for tool_label in sorted(df_max["tool_label"].unique().to_list()):
            df_tool = df_max.filter(pl.col("tool_label") == tool_label).sort("threads")
            ax.plot(
                df_tool["threads"].to_list(),
                df_tool[time_col].to_list(),
                marker="o",
                label=display_name(tool_label),
                color=COLORS.get(tool_label, "#95a5a6"),
                linestyle=LINESTYLES.get(tool_label, "-"),
                linewidth=2,
                markersize=8,
            )

        ax.set_xlabel("Threads")
        ax.set_ylabel(f"{metric_label(time_col)} (ms)")
        ax.set_xticks(thread_counts)
        ax.set_title(f"SR Thread Scaling: {max_struct} ({max_atoms:,} atoms)")
        ax.grid(True, alpha=0.3)
        ax.legend()

        fig.tight_layout()
        out_path = plot_dir.joinpath("max_structure.png")
        fig.savefig(out_path)
        plt.close(fig)
        rprint(f"[green]Saved:[/green] {out_path}")


def plot_large(n_points: int = 100, time_col: str = "time_ms"):
    """Generate speedup bar chart for large structures (50k+ atoms)."""
    setup_style()
    df = load_data(n_points)
    df = add_size_bin(df)

    suffix = "_sasa" if time_col == "sasa_time_ms" else ""
    plot_dir = PLOTS_DIR.joinpath(f"large{suffix}")
    plot_dir.mkdir(parents=True, exist_ok=True)

    # Large bins (50k+)
    large_bins = {"50k-75k", "75k-100k", "100k-150k", "150k-200k", "200k-500k", "500k+"}

    # Use max thread count available
    max_threads = df["threads"].max()
    df_large = df.filter(
        pl.col("size_bin").is_in(large_bins) & (pl.col("threads") == max_threads)
    )

    if df_large.height == 0:
        rprint("[yellow]No large structures found[/yellow]")
        return

    results = []

    pivot = (
        df_large.select(["structure", "tool_label", time_col])
        .pivot(on="tool_label", index="structure", values=time_col)
        .drop_nulls()
    )

    if "zsasa_f64" in pivot.columns and "rustsasa" in pivot.columns:
        speedup = pivot["rustsasa"] / pivot["zsasa_f64"]
        results.append(
            {
                "label": "vs RustSASA (SR)",
                "speedup": speedup.median(),
                "q25": speedup.quantile(0.25),
                "q75": speedup.quantile(0.75),
                "color": "#3498db",
            }
        )

    if "zsasa_f64" in pivot.columns and "freesasa" in pivot.columns:
        speedup = pivot["freesasa"] / pivot["zsasa_f64"]
        results.append(
            {
                "label": "vs FreeSASA (SR)",
                "speedup": speedup.median(),
                "q25": speedup.quantile(0.25),
                "q75": speedup.quantile(0.75),
                "color": "#3498db",
            }
        )

    if not results:
        rprint("[yellow]No comparison data available[/yellow]")
        return

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

    ax.set_yticks(y_pos)
    ax.set_yticklabels(labels)
    ax.set_xlabel("Speedup (median, IQR)")
    ax.axvline(x=1.0, color="gray", linestyle="--", linewidth=1.5)

    for i, s in enumerate(speedups):
        ax.text(s + 0.05, i, f"{s:.2f}x", va="center", fontsize=10)

    n_structs = df_large.select("structure").unique().height
    ax.set_title(
        f"zsasa Speedup on Large Structures (50k+ atoms, n={n_structs}, threads={max_threads})"
    )
    ax.set_xlim(0, max(speedups) * 1.2)
    ax.grid(True, alpha=0.3, axis="x")

    fig.tight_layout()
    out_path = plot_dir.joinpath("speedup_bar.png")
    fig.savefig(out_path)
    plt.close(fig)
    rprint(f"[green]Saved:[/green] {out_path}")

    # --- Thread scaling for large structures ---
    df_large_all = df.filter(pl.col("size_bin").is_in(large_bins))
    if df_large_all.height > 0:
        fig2, ax2 = plt.subplots(figsize=(10, 6))
        _plot_threads(df_large_all, ax2, time_col=time_col)
        n_structs = df_large_all.select("structure").unique().height
        ax2.set_title(
            f"SR: Thread Scaling on Large Structures (50k+ atoms, n={n_structs})"
        )
        fig2.tight_layout()
        out_path2 = plot_dir.joinpath("speedup_by_threads.png")
        fig2.savefig(out_path2)
        plt.close(fig2)
        rprint(f"[green]Saved:[/green] {out_path2}")


def plot_efficiency(n_points: int = 100, time_col: str = "time_ms"):
    """Calculate and plot parallel efficiency from existing benchmark data."""
    from rich.table import Table

    setup_style()
    df = load_data(n_points)
    df = add_size_bin(df)

    if df.height == 0:
        rprint("[yellow]No SR data found[/yellow]")
        return

    suffix = "_sasa" if time_col == "sasa_time_ms" else ""
    plot_dir = PLOTS_DIR.joinpath(f"efficiency{suffix}")
    plot_dir.mkdir(parents=True, exist_ok=True)

    # Baseline: threads=1
    df_t1 = df.filter(pl.col("threads") == 1).select(
        [
            "tool_label",
            "structure",
            "n_atoms",
            "size_bin",
            pl.col(time_col).alias("t1_ms"),
        ]
    )

    df_eff = df.join(df_t1, on=["tool_label", "structure", "n_atoms", "size_bin"])
    df_eff = df_eff.with_columns(
        (pl.col("t1_ms") / (pl.col(time_col) * pl.col("threads"))).alias("efficiency")
    )

    tool_labels = ["zsasa_f64", "freesasa", "rustsasa"]
    max_threads = df["threads"].max()

    df_tmax = df_eff.filter(pl.col("threads") == max_threads)

    summary_table = Table(
        title=f"Parallel Efficiency by Size (SR, threads={max_threads})"
    )
    summary_table.add_column("Size Bin", style="cyan")
    summary_table.add_column("Count", justify="right")
    summary_table.add_column("zsasa (f64)", justify="right")
    summary_table.add_column("FreeSASA", justify="right")
    summary_table.add_column("RustSASA", justify="right")

    bin_order = [b[2] for b in BINS]

    for bin_name in bin_order:
        df_bin = df_tmax.filter(pl.col("size_bin") == bin_name)
        if df_bin.height == 0:
            continue

        count = df_bin.select("structure").unique().height
        row_values = [bin_name, f"{count:,}"]

        for tool in tool_labels:
            df_tool = df_bin.filter(pl.col("tool_label") == tool)
            if df_tool.height > 0:
                eff = df_tool["efficiency"].median()
                row_values.append(f"{eff:.2f}")
            else:
                row_values.append("-")

        summary_table.add_row(*row_values)

    rprint(summary_table)
    rprint("[dim]Efficiency = T1 / (TN * N). Higher = better thread utilization[/dim]")

    # Plot efficiency by size bin
    fig, ax = plt.subplots(figsize=(14, 6))

    x_labels = [b for b in bin_order if b in df_tmax["size_bin"].unique().to_list()]
    x_pos = list(range(len(x_labels)))
    bar_width = 0.25

    for i, tool in enumerate(tool_labels):
        tool_eff = []
        for bin_name in x_labels:
            df_bin_tool = df_tmax.filter(
                (pl.col("size_bin") == bin_name) & (pl.col("tool_label") == tool)
            )
            if df_bin_tool.height > 0:
                tool_eff.append(df_bin_tool["efficiency"].median())
            else:
                tool_eff.append(0)

        offset = (i - 1) * bar_width
        ax.bar(
            [x + offset for x in x_pos],
            tool_eff,
            bar_width,
            label=display_name(tool),
            color=COLORS.get(tool, "#95a5a6"),
        )

    ax.set_xlabel("Structure Size (atoms)")
    ax.set_ylabel(f"Parallel Efficiency (threads={max_threads})")
    ax.set_title("Parallel Efficiency by Structure Size")
    ax.set_xticks(x_pos)
    ax.set_xticklabels(x_labels, rotation=45, ha="right")
    ax.legend()
    ax.grid(True, alpha=0.3, axis="y")
    ax.set_ylim(0, 0.5)

    fig.tight_layout()
    out_path = plot_dir.joinpath("summary.png")
    fig.savefig(out_path)
    plt.close(fig)
    rprint(f"[green]Saved:[/green] {out_path}")


def plot_speedup(
    min_atoms: int = 50000,
    top_n: int = 5,
    n_points: int = 100,
    time_col: str = "time_ms",
):
    """Find structures with best zsasa speedup at any thread count."""
    from rich.table import Table

    setup_style()
    df = load_data(n_points)
    df = add_size_bin(df)

    df_large = df.filter(pl.col("n_atoms") >= min_atoms)

    if df_large.height == 0:
        rprint("[yellow]No SR data with 50k+ atoms found[/yellow]")
        return

    pivot = (
        df_large.select(["structure", "n_atoms", "threads", "tool_label", time_col])
        .pivot(
            on="tool_label", index=["structure", "n_atoms", "threads"], values=time_col
        )
        .drop_nulls()
    )

    zsasa_col = "zsasa_f64" if "zsasa_f64" in pivot.columns else "zsasa"
    if zsasa_col not in pivot.columns:
        rprint("[red]No zsasa data found[/red]")
        return

    speedup_cols = []
    if "freesasa" in pivot.columns:
        pivot = pivot.with_columns(
            (pl.col("freesasa") / pl.col(zsasa_col)).alias("vs_freesasa")
        )
        speedup_cols.append("vs_freesasa")
    if "rustsasa" in pivot.columns:
        pivot = pivot.with_columns(
            (pl.col("rustsasa") / pl.col(zsasa_col)).alias("vs_rustsasa")
        )
        speedup_cols.append("vs_rustsasa")

    suffix = "_sasa" if time_col == "sasa_time_ms" else ""
    plot_dir = PLOTS_DIR.joinpath(f"speedup{suffix}")
    plot_dir.mkdir(parents=True, exist_ok=True)

    results = {}

    for col in speedup_cols:
        comparison = col.replace("vs_", "")
        top_entries = pivot.sort(col, descending=True).head(top_n)

        table = Table(
            title=f"Top {top_n} zsasa Speedup vs {display_name(comparison)} ({min_atoms:,}+ atoms)"
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

    # Generate comparison plot
    if len(results) < 2:
        return

    best_fs = results["freesasa"][0]
    best_rs = results["rustsasa"][0]

    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    for idx, (comparison, best) in enumerate(
        [("freesasa", best_fs), ("rustsasa", best_rs)]
    ):
        ax = axes[idx]
        struct = best["structure"]
        n_atoms = best["n_atoms"]

        df_struct = pivot.filter(pl.col("structure") == struct).sort("threads")

        threads_list = df_struct["threads"].to_list()
        vs_fs = df_struct["vs_freesasa"].to_list()
        vs_rs = df_struct["vs_rustsasa"].to_list()

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
            vs_rs,
            marker="s",
            label="vs RustSASA",
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

        ax.annotate(
            f"{vs_fs[-1]:.2f}x",
            xy=(threads_list[-1], vs_fs[-1]),
            xytext=(5, 0),
            textcoords="offset points",
            color="#3498db",
            fontsize=10,
        )
        ax.annotate(
            f"{vs_rs[-1]:.2f}x",
            xy=(threads_list[-1], vs_rs[-1]),
            xytext=(5, 0),
            textcoords="offset points",
            color="#e74c3c",
            fontsize=10,
        )

        best_fs_val = max(vs_fs)
        best_fs_t = threads_list[vs_fs.index(best_fs_val)]
        best_rs_val = max(vs_rs)
        best_rs_t = threads_list[vs_rs.index(best_rs_val)]

        if comparison == "freesasa":
            title = f"{struct} ({n_atoms:,} atoms)\nBest vs FreeSASA: {best_fs_val:.2f}x @ threads={best_fs_t}"
        else:
            title = f"{struct} ({n_atoms:,} atoms)\nBest vs RustSASA: {best_rs_val:.2f}x @ threads={best_rs_t}"
        ax.set_title(title)

    fig.suptitle("zsasa Speedup vs FreeSASA and RustSASA (SR Algorithm)", fontsize=14)
    fig.tight_layout()

    out_path = plot_dir.joinpath("comparison.png")
    fig.savefig(out_path)
    plt.close(fig)
    rprint(f"[green]Saved:[/green] {out_path}")
