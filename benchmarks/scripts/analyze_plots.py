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
    detect_outlier_rows,
    display_name,
    load_data,
    metric_label,
    metric_suffix,
    setup_style,
    split_outliers,
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


_COMPARISONS_F64 = [
    ("zsasa_f64_vs_freesasa", "f64 vs FreeSASA", "o", "#3498db"),
    ("zsasa_f64_vs_rustsasa", "f64 vs RustSASA", "s", "#e74c3c"),
]

_COMPARISONS_BITMASK = [
    ("zsasa_f32_bitmask_vs_freesasa", "bitmask vs FreeSASA", "^", "#2980b9"),
    ("zsasa_f32_bitmask_vs_rustsasa", "bitmask vs RustSASA", "v", "#c0392b"),
]


def _plot_speedup_single(
    df_sr: pl.DataFrame,
    threads: int,
    ax,
    show_legend: bool = True,
    time_col: str = "time_ms",
    comparisons: list | None = None,
):
    """Plot speedup for a single thread count on given axes."""
    bin_labels = [b[2] for b in BINS]
    speedup_data = compute_speedup_by_bin(df_sr, threads=threads, time_col=time_col)
    data_dict = {r["size_bin"]: r for r in speedup_data.iter_rows(named=True)}

    x_labels = [b for b in bin_labels if b in data_dict]
    x_pos = list(range(len(x_labels)))

    if comparisons is None:
        comparisons = _COMPARISONS_F64 + _COMPARISONS_BITMASK

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
    ax.set_ylim(bottom=0)
    ax.grid(True, alpha=0.3)

    if show_legend:
        ax.legend(fontsize=8)


# === Plot Command Implementations ===


def plot_scatter(n_points: int = 100, time_col: str = "time_ms"):
    """Generate atoms vs time scatter plot."""
    setup_style()
    df = load_data(n_points)

    suffix = metric_suffix(time_col)
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

    # Grid of all threads (2x2 layout for 4 thread counts)
    n_threads = len(thread_counts)
    n_cols = min(2, n_threads)
    n_rows = (n_threads + n_cols - 1) // n_cols

    fig, axes = plt.subplots(n_rows, n_cols, figsize=(7 * n_cols, 6 * n_rows))
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

    suffix = metric_suffix(time_col)
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

    df, df_outliers = split_outliers(df, time_col=time_col)
    if df_outliers.height > 0:
        n_out = df_outliers.select("structure").unique().height
        rprint(f"[dim]Excluded {n_out} outlier structure(s) from speedup grid[/dim]")

    suffix = metric_suffix(time_col)
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

    # Generate combined grid (all comparisons)
    n_threads = len(thread_counts)
    n_cols = min(2, n_threads)
    n_rows = (n_threads + n_cols - 1) // n_cols

    fig, axes = plt.subplots(
        n_rows, n_cols, figsize=(7 * n_cols, 5 * n_rows), sharey=True
    )
    if n_threads == 1:
        axes = [[axes]]
    elif n_rows == 1:
        axes = [axes]
    elif n_cols == 1:
        axes = [[ax] for ax in axes]

    for idx, threads in enumerate(thread_counts):
        row, col = idx // n_cols, idx % n_cols
        ax = axes[row][col]
        _plot_speedup_single(
            df,
            threads,
            ax,
            show_legend=(idx == 0),
            time_col=time_col,
        )
        ax.set_title(f"threads={threads}")
        if col > 0:
            ax.set_ylabel("")

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

    df, df_outliers = split_outliers(df, time_col=time_col)
    if df_outliers.height > 0:
        n_out = df_outliers.select("structure").unique().height
        rprint(f"[dim]Excluded {n_out} outlier structure(s) from samples[/dim]")

    df = add_size_bin(df)

    suffix = metric_suffix(time_col)
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

    df, df_outliers = split_outliers(df, time_col=time_col)
    if df_outliers.height > 0:
        n_out = df_outliers.select("structure").unique().height
        rprint(f"[dim]Excluded {n_out} outlier structure(s) from large analysis[/dim]")

    df = add_size_bin(df)

    suffix = metric_suffix(time_col)
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

    max_err = max(r["q75"] - r["speedup"] for r in results)
    for i, s in enumerate(speedups):
        ax.text(s + max_err + 0.08, i, f"{s:.2f}x", va="center", fontsize=10)

    n_structs = df_large.select("structure").unique().height
    ax.set_title(
        f"zsasa Speedup on Large Structures (50k+ atoms, n={n_structs}, threads={max_threads})"
    )
    ax.set_xlim(0, max(s + max_err for s in speedups) + 0.4)
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


def plot_memory(n_points: int = 100):
    """Generate peak memory (RSS) comparison plots."""
    setup_style()
    df = load_data(n_points)
    df = add_size_bin(df)

    # Filter to single-threaded and drop rows without memory data
    df_t1 = df.filter((pl.col("threads") == 1) & pl.col("memory_bytes").is_not_null())

    if df_t1.height == 0:
        rprint("[yellow]No memory data found[/yellow]")
        return

    df_t1 = df_t1.with_columns(
        (pl.col("memory_bytes") / (1024 * 1024)).alias("memory_mb")
    )

    plot_dir = PLOTS_DIR.joinpath("memory")
    plot_dir.mkdir(parents=True, exist_ok=True)

    # Exclude f32 variants for cleaner comparison
    tool_labels = [
        t for t in sorted(df_t1["tool_label"].unique().to_list()) if "f32" not in t
    ]

    # --- Scatter: atoms vs memory ---
    fig, ax = plt.subplots(figsize=(10, 6))

    for tool_label in tool_labels:
        df_tool = df_t1.filter(pl.col("tool_label") == tool_label)
        ax.scatter(
            df_tool["n_atoms"].to_list(),
            df_tool["memory_mb"].to_list(),
            label=display_name(tool_label),
            alpha=0.4,
            s=10,
            color=COLORS.get(tool_label, "#95a5a6"),
        )

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Number of Atoms")
    ax.set_ylabel("Peak RSS (MB)")
    ax.set_title("SR: Peak Memory vs Structure Size (threads=1)")
    ax.legend()
    ax.grid(True, alpha=0.3)

    fig.tight_layout()
    out_path = plot_dir.joinpath("scatter.png")
    fig.savefig(out_path)
    plt.close(fig)
    rprint(f"[green]Saved:[/green] {out_path}")

    # --- Bar chart: median memory by size bin ---
    bin_order = [b[2] for b in BINS]
    available_bins = [b for b in bin_order if b in df_t1["size_bin"].unique().to_list()]

    if not available_bins:
        return

    fig, ax = plt.subplots(figsize=(14, 6))

    x_pos = list(range(len(available_bins)))
    bar_width = 0.8 / max(len(tool_labels), 1)

    for i, tool in enumerate(tool_labels):
        medians = []
        for bin_name in available_bins:
            df_bt = df_t1.filter(
                (pl.col("size_bin") == bin_name) & (pl.col("tool_label") == tool)
            )
            medians.append(df_bt["memory_mb"].median() if df_bt.height > 0 else 0)

        offset = (i - len(tool_labels) / 2 + 0.5) * bar_width
        ax.bar(
            [x + offset for x in x_pos],
            medians,
            bar_width,
            label=display_name(tool),
            color=COLORS.get(tool, "#95a5a6"),
        )

    ax.set_xlabel("Structure Size (atoms)")
    ax.set_ylabel("Median Peak RSS (MB)")
    ax.set_title("SR: Peak Memory by Structure Size (threads=1)")
    ax.set_xticks(x_pos)
    ax.set_xticklabels(available_bins, rotation=45, ha="right")
    ax.legend()
    ax.grid(True, alpha=0.3, axis="y")

    fig.tight_layout()
    out_path = plot_dir.joinpath("by_size.png")
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

    df, df_outliers = split_outliers(df, time_col=time_col)
    if df_outliers.height > 0:
        n_out = df_outliers.select("structure").unique().height
        rprint(
            f"[dim]Excluded {n_out} outlier structure(s) from speedup analysis[/dim]"
        )

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

    suffix = metric_suffix(time_col)
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


def plot_outliers(n_points: int = 100, time_col: str = "sasa_time_ms"):
    """Generate plots highlighting outlier structures where competitors struggle.

    Shows structures where FreeSASA or RustSASA exhibit pathological behavior
    (>10x slower) while zsasa handles them without issues.
    """
    from rich.table import Table

    setup_style()
    df = load_data(n_points)
    outlier_rows = detect_outlier_rows(df, time_col=time_col)

    if not outlier_rows:
        rprint("[yellow]No outlier structures detected[/yellow]")
        return

    outlier_names = sorted({name for name, _ in outlier_rows})
    df_outliers = df.filter(pl.col("structure").is_in(outlier_names))
    rprint(
        f"[bold]Outlier structures ({len(outlier_names)}):[/bold] {', '.join(outlier_names)}"
    )

    suffix = metric_suffix(time_col)
    plot_dir = PLOTS_DIR.joinpath(f"outliers{suffix}")
    plot_dir.mkdir(parents=True, exist_ok=True)

    # --- Table: outlier summary (only flagged thread counts) ---
    pivot_all = (
        df_outliers.select(["structure", "tool_label", "n_atoms", "threads", time_col])
        .pivot(
            on="tool_label",
            index=["structure", "n_atoms", "threads"],
            values=time_col,
        )
        .sort(["n_atoms", "threads"])
    )
    # Filter to only the (structure, threads) pairs that were actually flagged
    # and have zsasa data (drop rows where all tool values are null)
    zsasa_col = "zsasa_f64" if "zsasa_f64" in pivot_all.columns else "zsasa"
    pivot = pivot_all.filter(
        pl.struct(["structure", "threads"]).map_elements(
            lambda row: (row["structure"], row["threads"]) in outlier_rows,
            return_dtype=pl.Boolean,
        )
    )
    table = Table(title="Outlier Structures: zsasa Advantage")
    table.add_column("Structure", style="cyan")
    table.add_column("Atoms", justify="right")
    table.add_column("T", justify="right")
    table.add_column("zsasa (ms)", justify="right")
    for competitor in ("freesasa", "rustsasa"):
        if competitor in pivot.columns:
            table.add_column(f"{display_name(competitor)} (ms)", justify="right")
            table.add_column("Ratio", justify="right")

    for row in pivot.iter_rows(named=True):
        zsasa_ms = row.get(zsasa_col)
        cells = [row["structure"], f"{row['n_atoms']:,}", str(row["threads"])]
        cells.append(f"{zsasa_ms:.1f}" if zsasa_ms is not None else "-")
        for competitor in ("freesasa", "rustsasa"):
            if competitor not in pivot.columns:
                continue
            comp_ms = row.get(competitor)
            if comp_ms is not None and zsasa_ms:
                ratio = comp_ms / zsasa_ms
                color = "green" if ratio > 1.0 else "red"
                cells.append(f"{comp_ms:.1f}")
                cells.append(f"[{color}]{ratio:.1f}x[/{color}]")
            else:
                cells.extend(["-", "-"])
        table.add_row(*cells)

    rprint(table)

    # --- Per-competitor per-thread bar charts ---
    if zsasa_col not in pivot.columns:
        return

    thread_counts = sorted(df_outliers["threads"].unique().to_list())
    for threads in thread_counts:
        pivot_t = pivot_all.filter(pl.col("threads") == threads)

        ratio_data: dict[str, list[tuple[str, int, float]]] = {}
        for competitor in ("freesasa", "rustsasa"):
            if competitor not in pivot_t.columns:
                continue
            entries = []
            for row in pivot_t.iter_rows(named=True):
                zsasa_ms = row.get(zsasa_col)
                comp_ms = row.get(competitor)
                if zsasa_ms and comp_ms:
                    ratio = comp_ms / zsasa_ms
                    if ratio > 1.5:
                        entries.append((row["structure"], row["n_atoms"], ratio))
            if entries:
                ratio_data[competitor] = sorted(entries, key=lambda x: x[1])

        for competitor, entries in ratio_data.items():
            if len(entries) < 2:
                continue

            comp_name = display_name(competitor)
            structures = [e[0] for e in entries]
            n_atoms_list = [e[1] for e in entries]
            ratios = [e[2] for e in entries]
            labels = [f"{s}\n({n:,})" for s, n in zip(structures, n_atoms_list)]

            fig, ax = plt.subplots(figsize=(max(8, len(entries) * 0.8), 6))

            x_pos = list(range(len(entries)))
            bars = ax.bar(
                x_pos,
                ratios,
                0.6,
                color=COLORS.get(competitor, "#95a5a6"),
                alpha=0.8,
            )

            for bar, ratio in zip(bars, ratios):
                ax.text(
                    bar.get_x() + bar.get_width() / 2,
                    bar.get_height() + 0.5,
                    f"{ratio:.1f}x",
                    ha="center",
                    va="bottom",
                    fontsize=8,
                    fontweight="bold",
                )

            ax.axhline(y=1.0, color="gray", linestyle="--", linewidth=1.5)
            ax.set_xticks(x_pos)
            ax.set_xticklabels(labels, fontsize=8)
            ax.set_ylabel(f"{comp_name} / zsasa Ratio ({metric_label(time_col)})")
            ax.set_title(
                f"Outlier Structures: {comp_name} Pathological Cases"
                f" (t={threads}, n={len(entries)})"
            )
            ax.grid(True, alpha=0.3, axis="y")

            fig.tight_layout()
            out_path = plot_dir.joinpath(f"{competitor}_bar_t{threads}.png")
            fig.savefig(out_path)
            plt.close(fig)
            rprint(f"[green]Saved:[/green] {out_path}")

    # --- Thread scaling for each outlier structure ---
    thread_counts = sorted(df_outliers["threads"].unique().to_list())
    for struct in outlier_names:
        df_struct = df_outliers.filter(pl.col("structure") == struct)
        if df_struct.height == 0:
            continue

        n_atoms = df_struct["n_atoms"][0]
        fig, ax = plt.subplots(figsize=(8, 5))

        for tool_label in sorted(df_struct["tool_label"].unique().to_list()):
            if "f32" in tool_label:
                continue
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
                linewidth=2,
                markersize=6,
            )

        ax.set_xlabel("Threads")
        ax.set_ylabel(f"{metric_label(time_col)} (ms)")
        ax.set_xticks(thread_counts)
        ax.set_title(f"{struct} ({n_atoms:,} atoms) — Outlier")
        ax.set_yscale("log")
        ax.grid(True, alpha=0.3)
        ax.legend()

        fig.tight_layout()
        out_path = plot_dir.joinpath(f"{struct}.png")
        fig.savefig(out_path)
        plt.close(fig)
        rprint(f"[green]Saved:[/green] {out_path}")
