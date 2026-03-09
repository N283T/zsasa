"""Data loading and utility functions for benchmark analysis."""

from pathlib import Path

import matplotlib.pyplot as plt
import polars as pl

# === Path Constants ===

_BENCHMARKS_DIR = Path(__file__).parent.parent
RESULTS_BASE = _BENCHMARKS_DIR.joinpath("results", "single")
PLOTS_DIR = _BENCHMARKS_DIR.joinpath("results", "plots", "single")

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


def _remap_tool_names(df: pl.DataFrame) -> pl.DataFrame:
    """Remap raw tool names to canonical display names."""
    return df.with_columns(
        pl.col("tool")
        .str.replace(r"^zig", "zsasa")
        .str.replace(r"^rust$", "rustsasa")
        .alias("tool"),
    )


def _merge_timing_data(df: pl.DataFrame, sr_dir: Path) -> pl.DataFrame:
    """Merge timing.csv data into the main DataFrame.

    Looks for timing.csv files alongside results.csv in each tool directory.
    Aggregates duplicate rows (from multiple sasa runs) by mean, then
    left-joins on (tool, structure, threads).
    """
    timing_files = sorted(sr_dir.glob("*/timing.csv"))
    if not timing_files:
        return df

    timing_dfs = [pl.read_csv(f) for f in timing_files]
    timing = pl.concat(timing_dfs, how="diagonal")
    timing = _remap_tool_names(timing)

    # Aggregate to one row per (tool, structure, threads) to prevent row duplication
    timing = (
        timing.with_columns(pl.col("threads").cast(pl.Int64))
        .group_by(["tool", "structure", "threads"])
        .agg(
            pl.col("parse_time_ms").mean().alias("parse_time_ms"),
            pl.col("sasa_time_ms").mean().alias("sasa_time_ms"),
        )
    )

    # Drop existing empty timing columns from df before join
    drop_cols = [c for c in ("parse_time_ms", "sasa_time_ms") if c in df.columns]
    if drop_cols:
        df = df.drop(drop_cols)

    # Left join: keep all rows from df, add timing where available
    return df.join(timing, on=["tool", "structure", "threads"], how="left")


def load_data(n_points: int = 100) -> pl.DataFrame:
    """Load SR benchmark results from RESULTS_BASE/<n_points>/**/results.csv.

    Derives tool_label from the CSV data (tool + precision for zig).
    Aggregates multiple runs into mean values per (structure, threads).
    """
    sr_dir = RESULTS_BASE.joinpath(str(n_points))
    csv_files = sorted(sr_dir.glob("*/results.csv"))

    if not csv_files:
        raise FileNotFoundError(f"No results.csv files found in {sr_dir}")

    dfs = []
    for f in csv_files:
        df = pl.read_csv(f)
        dfs.append(df)

    df = pl.concat(dfs, how="diagonal")

    # Ensure precision column exists (older CSVs might lack it)
    if "precision" not in df.columns:
        df = df.with_columns(pl.lit(None).cast(pl.Utf8).alias("precision"))

    # Remap tool names: zig* → zsasa, rust → rustsasa
    df = _remap_tool_names(df)

    # Create tool_label from tool column (already contains precision/variant info)
    df = df.with_columns(pl.col("tool").alias("tool_label"))

    # Convert current CSV format (mean_s etc.) to time_ms used downstream.
    # Old format had sasa_time_ms per-run rows; new format has pre-aggregated
    # mean_s/stddev_s/median_s in seconds.
    if "mean_s" in df.columns:
        df = df.with_columns(
            (pl.col("mean_s") * 1000).alias("time_ms"),
            (pl.col("stddev_s") * 1000).alias("time_std"),
            (pl.col("median_s") * 1000).alias("median_ms"),
        )
        # total_sasa not available in new format
        if "total_sasa" not in df.columns:
            df = df.with_columns(pl.lit(None).cast(pl.Float64).alias("total_sasa"))

        # Internal timing: merge from timing.csv if available
        df = _merge_timing_data(df, sr_dir)

        # Ensure timing columns exist (may still be absent after merge)
        for col in ("parse_time_ms", "sasa_time_ms"):
            if col not in df.columns:
                df = df.with_columns(pl.lit(None).cast(pl.Float64).alias(col))

        # Memory column (from hyperfine, may be absent in old CSVs)
        if "memory_bytes" not in df.columns:
            df = df.with_columns(pl.lit(None).cast(pl.Int64).alias("memory_bytes"))

        return (
            df.select(
                [
                    "tool",
                    "tool_label",
                    "structure",
                    "n_atoms",
                    "algorithm",
                    "precision",
                    "threads",
                    "time_ms",
                    "time_std",
                    "total_sasa",
                    "parse_time_ms",
                    "sasa_time_ms",
                    "memory_bytes",
                ]
            )
            .with_columns(pl.lit(1).cast(pl.UInt32).alias("n_runs"))
            .sort(["algorithm", "tool_label", "n_atoms"])
        )

    # Legacy format: per-run rows with sasa_time_ms
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


def compute_speedup_by_bin(
    df: pl.DataFrame, threads: int = 1, time_col: str = "time_ms"
) -> pl.DataFrame:
    """Compute speedup ratios by size bin for given thread count."""
    df_t = df.filter(pl.col("threads") == threads)

    pivot = (
        df_t.select(["structure", "tool_label", "n_atoms", time_col])
        .pivot(on="tool_label", index=["structure", "n_atoms"], values=time_col)
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

    # zsasa_f32_bitmask vs FreeSASA
    if "zsasa_f32_bitmask" in cols and "freesasa" in cols:
        pivot = pivot.with_columns(
            (pl.col("freesasa") / pl.col("zsasa_f32_bitmask")).alias(
                "zsasa_f32_bitmask_vs_freesasa"
            )
        )

    # zsasa_f32_bitmask vs RustSASA
    if "zsasa_f32_bitmask" in cols and "rustsasa" in cols:
        pivot = pivot.with_columns(
            (pl.col("rustsasa") / pl.col("zsasa_f32_bitmask")).alias(
                "zsasa_f32_bitmask_vs_rustsasa"
            )
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


# === Outlier Detection ===

# Per-bin IQR outlier detection with a minimum ratio floor.
# k=1.5 is the standard Tukey "outlier" fence (Q3 + k*IQR).
# The floor prevents flagging structures in tight-IQR bins where the actual
# ratio is modest (e.g. 1.6x in a bin where IQR is 0.07).
OUTLIER_IQR_K = 1.5
OUTLIER_RATIO_FLOOR = 5.0


def detect_outliers(
    df: pl.DataFrame,
    time_col: str = "sasa_time_ms",
    iqr_k: float = OUTLIER_IQR_K,
    ratio_floor: float = OUTLIER_RATIO_FLOOR,
) -> set[str]:
    """Return set of structure names where any competitor is a per-bin outlier.

    For each size bin × competitor, computes the Tukey fence Q3 + k*IQR on the
    ratio (competitor_time / zsasa_time).  A structure is flagged only if it
    exceeds BOTH the per-bin fence AND the absolute ratio floor, at ANY thread
    count.
    """
    return {name for name, _ in detect_outlier_rows(df, time_col, iqr_k, ratio_floor)}


def detect_outlier_rows(
    df: pl.DataFrame,
    time_col: str = "sasa_time_ms",
    iqr_k: float = OUTLIER_IQR_K,
    ratio_floor: float = OUTLIER_RATIO_FLOOR,
) -> set[tuple[str, int]]:
    """Return set of (structure, threads) pairs that are per-bin outliers.

    For each size bin × competitor, computes the Tukey fence Q3 + k*IQR on the
    ratio (competitor_time / zsasa_time).  A row is flagged only if it exceeds
    BOTH the per-bin fence AND the absolute ratio floor.
    """
    pivot = (
        df.select(["structure", "threads", "tool_label", "n_atoms", time_col])
        .pivot(
            on="tool_label", index=["structure", "n_atoms", "threads"], values=time_col
        )
        .drop_nulls()
    )
    pivot = add_size_bin(pivot)

    zsasa_col = "zsasa_f64" if "zsasa_f64" in pivot.columns else "zsasa"
    if zsasa_col not in pivot.columns:
        return set()

    outlier_rows: set[tuple[str, int]] = set()

    for competitor in ("freesasa", "rustsasa"):
        if competitor not in pivot.columns:
            continue
        ratio_col = f"{competitor}_ratio"
        with_ratio = pivot.with_columns(
            (pl.col(competitor) / pl.col(zsasa_col)).alias(ratio_col)
        )

        # Compute per-bin Tukey fence
        fences = (
            with_ratio.group_by("size_bin")
            .agg(
                pl.col(ratio_col).quantile(0.75).alias("q75"),
                (
                    pl.col(ratio_col).quantile(0.75) - pl.col(ratio_col).quantile(0.25)
                ).alias("iqr"),
            )
            .with_columns((pl.col("q75") + iqr_k * pl.col("iqr")).alias("fence"))
            .select(["size_bin", "fence"])
        )

        joined = with_ratio.join(fences, on="size_bin")
        outliers = joined.filter(
            (pl.col(ratio_col) > pl.col("fence")) & (pl.col(ratio_col) > ratio_floor)
        )
        for row in outliers.iter_rows(named=True):
            outlier_rows.add((row["structure"], row["threads"]))

    return outlier_rows


def split_outliers(
    df: pl.DataFrame,
    time_col: str = "sasa_time_ms",
) -> tuple[pl.DataFrame, pl.DataFrame]:
    """Split DataFrame into (clean, outliers) based on outlier detection."""
    outlier_names = detect_outliers(df, time_col=time_col)
    if not outlier_names:
        return df, df.head(0)
    clean = df.filter(~pl.col("structure").is_in(outlier_names))
    outliers = df.filter(pl.col("structure").is_in(outlier_names))
    return clean, outliers


METRIC_LABELS = {
    "time_ms": "Wall-clock Time",
    "sasa_time_ms": "SASA-only Time",
}


def metric_label(time_col: str) -> str:
    """Human-readable label for a time column."""
    return METRIC_LABELS.get(time_col, time_col)


def metric_suffix(time_col: str) -> str:
    """File-name suffix for a time column (empty for wall-clock)."""
    return "_sasa" if time_col == "sasa_time_ms" else ""


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
