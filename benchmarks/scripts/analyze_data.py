"""Data loading and utility functions for benchmark analysis."""

from pathlib import Path

import matplotlib.pyplot as plt
import polars as pl

# === Path Constants ===

RESULTS_DIR = Path(__file__).parent.parent.joinpath("results")
PLOTS_DIR = RESULTS_DIR.joinpath("plots")

# === Visual Constants ===

COLORS = {
    "zsasa": "#f39c12",  # Orange (LR uses plain "zsasa")
    "zsasa_f32": "#e67e22",  # Dark orange
    "zsasa_f64": "#f39c12",  # Light orange
    "freesasa": "#3498db",  # Blue
    "rustsasa": "#e74c3c",  # Red
}

LINESTYLES = {
    "zsasa": "-",
    "zsasa_f32": "--",
    "zsasa_f64": "-",
    "freesasa": "-",
    "rustsasa": "-",
}

DISPLAY_NAMES = {
    "zsasa_f64": "zsasa (f64)",
    "zsasa_f32": "zsasa (f32)",
    "zsasa": "zsasa",
    "freesasa": "FreeSASA",
    "rustsasa": "RustSASA",
}


def display_name(tool_label: str) -> str:
    """Get display name for a tool label."""
    return DISPLAY_NAMES.get(tool_label, tool_label)


# === Size Bin Definitions (v2: 14 bins) ===

BINS = [
    (0, 500, "0-500"),
    (500, 1000, "500-1k"),
    (1000, 2000, "1k-2k"),
    (2000, 3000, "2k-3k"),
    (3000, 5000, "3k-5k"),
    (5000, 10000, "5k-10k"),
    (10000, 20000, "10k-20k"),
    (20000, 50000, "20k-50k"),
    (50000, 75000, "50k-75k"),
    (75000, 100000, "75k-100k"),
    (100000, 150000, "100k-150k"),
    (150000, 200000, "150k-200k"),
    (200000, 500000, "200k-500k"),
    (500000, float("inf"), "500k+"),
]


# === Data Loading Functions ===


def load_data() -> pl.DataFrame:
    """Load benchmark results from all result directories.

    Auto-discovers directories under RESULTS_DIR that contain results.csv.
    Derives tool_label from the CSV data (tool + precision for zig).
    Aggregates multiple runs into mean values per (structure, threads).
    """
    csv_files = sorted(RESULTS_DIR.glob("*/results.csv"))

    if not csv_files:
        raise FileNotFoundError(f"No results.csv files found in {RESULTS_DIR}")

    dfs = []
    for f in csv_files:
        df = pl.read_csv(f)
        dfs.append(df)

    df = pl.concat(dfs, how="diagonal")

    # Ensure precision column exists (older CSVs might lack it)
    if "precision" not in df.columns:
        df = df.with_columns(pl.lit(None).cast(pl.Utf8).alias("precision"))

    # Remap tool names: zig → zsasa, rust → rustsasa
    df = df.with_columns(
        pl.col("tool").replace({"zig": "zsasa", "rust": "rustsasa"}).alias("tool"),
    )

    # Create tool_label: zsasa_f64, zsasa_f32, freesasa, rustsasa
    df = df.with_columns(
        pl.when((pl.col("tool") == "zsasa") & pl.col("precision").is_not_null())
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
            pl.col("sasa_time_ms").std().alias("time_std"),
            pl.col("sasa_time_ms").len().alias("n_runs"),
            pl.col("total_sasa").first(),
        )
        .sort(["algorithm", "tool_label", "n_atoms"])
    )


# === Utility Functions ===


def add_size_bin(df: pl.DataFrame) -> pl.DataFrame:
    """Add size bin column to DataFrame."""
    expr = pl.lit("500k+")
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

    cols = pivot.columns

    # zsasa_f64 vs FreeSASA
    if "zsasa_f64" in cols and "freesasa" in cols:
        pivot = pivot.with_columns(
            (pl.col("freesasa") / pl.col("zsasa_f64")).alias("zsasa_f64_vs_freesasa")
        )
    # Legacy: zsasa vs FreeSASA (LR without precision suffix)
    elif "zsasa" in cols and "freesasa" in cols:
        pivot = pivot.with_columns(
            (pl.col("freesasa") / pl.col("zsasa")).alias("zsasa_vs_freesasa")
        )

    # zsasa_f64 vs RustSASA
    if "zsasa_f64" in cols and "rustsasa" in cols:
        pivot = pivot.with_columns(
            (pl.col("rustsasa") / pl.col("zsasa_f64")).alias("zsasa_f64_vs_rustsasa")
        )

    # zsasa_f32 vs zsasa_f64
    if "zsasa_f32" in cols and "zsasa_f64" in cols:
        pivot = pivot.with_columns(
            (pl.col("zsasa_f64") / pl.col("zsasa_f32")).alias("zsasa_f32_vs_zsasa_f64")
        )

    # FreeSASA vs RustSASA
    if "freesasa" in cols and "rustsasa" in cols:
        pivot = pivot.with_columns(
            (pl.col("rustsasa") / pl.col("freesasa")).alias("freesasa_vs_rustsasa")
        )

    pivot = add_size_bin(pivot)

    # Aggregate by bin - median with IQR for error bars
    speedup_cols = [
        c
        for c in pivot.columns
        if c.endswith(("_vs_freesasa", "_vs_rustsasa", "_vs_zsasa_f64"))
    ]

    agg_exprs = [pl.len().alias("count")]
    for col in speedup_cols:
        agg_exprs.extend(
            [
                pl.col(col).median().alias(col),
                pl.col(col).quantile(0.25).alias(f"{col}_q25"),
                pl.col(col).quantile(0.75).alias(f"{col}_q75"),
            ]
        )

    return pivot.group_by("size_bin").agg(agg_exprs).sort("size_bin")


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
