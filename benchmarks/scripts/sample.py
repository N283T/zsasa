#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.11"
# dependencies = ["rich>=13.0", "typer>=0.9.0", "matplotlib>=3.0"]
# ///
"""Stratified sampling of benchmark structures by atom count.

Reads an index file and performs stratified sampling based on atom count bins.
Uses logarithmic scale bins to ensure representation across all size ranges.

Usage:
    # Analyze distribution
    ./benchmarks/scripts/sample.py benchmarks/inputs/index.json --analyze

    # Generate sample
    ./benchmarks/scripts/sample.py benchmarks/inputs/index.json \\
        --target 75000 --seed 42 --output benchmarks/samples/stratified_75k.json
"""

from __future__ import annotations

import json
import random
from datetime import datetime, timezone
from pathlib import Path
from typing import Annotated

import typer
from rich.console import Console
from rich.table import Table

app = typer.Typer(help="Stratified sampling of benchmark structures")
console = Console()

# Logarithmic scale bins (atom count ranges)
BINS = [
    ("0-500", 0, 500),
    ("500-2k", 500, 2_000),
    ("2k-10k", 2_000, 10_000),
    ("10k-50k", 10_000, 50_000),
    ("50k-200k", 50_000, 200_000),
    ("200k+", 200_000, float("inf")),
]

# Bins considered "rare" - take all samples from these
RARE_BINS = {"50k-200k", "200k+"}


def load_index(index_path: Path) -> dict[str, int]:
    """Load index file and return entries dict."""
    data = json.loads(index_path.read_text())
    return data.get("entries", {})


def categorize_by_bin(entries: dict[str, int]) -> dict[str, list[str]]:
    """Categorize entries into bins based on atom count."""
    bins: dict[str, list[str]] = {name: [] for name, _, _ in BINS}

    for pdb_id, n_atoms in entries.items():
        for name, low, high in BINS:
            if low <= n_atoms < high:
                bins[name].append(pdb_id)
                break

    return bins


def analyze_distribution(entries: dict[str, int]) -> None:
    """Display distribution analysis."""
    bins = categorize_by_bin(entries)

    table = Table(title="Atom Count Distribution")
    table.add_column("Bin", style="cyan")
    table.add_column("Range", style="dim")
    table.add_column("Count", justify="right")
    table.add_column("Percent", justify="right")
    table.add_column("Status", style="dim")

    total = len(entries)

    for name, low, high in BINS:
        count = len(bins[name])
        pct = count / total * 100 if total > 0 else 0
        high_str = f"{high:,}" if high != float("inf") else "+"
        range_str = f"{low:,} - {high_str}"
        status = "rare" if name in RARE_BINS else ""
        table.add_row(name, range_str, f"{count:,}", f"{pct:.1f}%", status)

    table.add_section()
    table.add_row("Total", "", f"{total:,}", "100.0%", "")

    console.print(table)


def stratified_sample(
    entries: dict[str, int],
    target: int,
    seed: int,
) -> tuple[list[str], dict[str, dict[str, int]]]:
    """Perform stratified sampling.

    Strategy:
    - Rare bins (50k+): Take all samples
    - Other bins: Proportional sampling from remaining target

    Returns:
        (samples, distribution) where distribution maps bin name to
        {"total": n, "sampled": m}
    """
    rng = random.Random(seed)
    bins = categorize_by_bin(entries)

    # Calculate rare bin contribution
    rare_samples: list[str] = []
    for name in RARE_BINS:
        rare_samples.extend(bins[name])

    remaining_target = max(0, target - len(rare_samples))

    # Calculate proportional sampling for non-rare bins
    non_rare_bins = {name: ids for name, ids in bins.items() if name not in RARE_BINS}
    non_rare_total = sum(len(ids) for ids in non_rare_bins.values())

    samples: list[str] = list(rare_samples)
    distribution: dict[str, dict[str, int]] = {}

    for name, _, _ in BINS:
        bin_ids = bins[name]
        total = len(bin_ids)

        if name in RARE_BINS:
            # Take all from rare bins
            sampled = total
        else:
            # Proportional sampling
            if non_rare_total > 0:
                proportion = len(bin_ids) / non_rare_total
                sample_count = int(remaining_target * proportion)
                # Ensure we don't sample more than available
                sample_count = min(sample_count, len(bin_ids))
                sampled_ids = rng.sample(bin_ids, sample_count)
                samples.extend(sampled_ids)
                sampled = sample_count
            else:
                sampled = 0

        distribution[name] = {"total": total, "sampled": sampled}

    # Shuffle final list for randomized order
    rng.shuffle(samples)

    return samples, distribution


@app.command()
def main(
    index_path: Annotated[
        Path,
        typer.Argument(help="Path to index.json file"),
    ],
    analyze: Annotated[
        bool,
        typer.Option("--analyze", "-a", help="Show distribution analysis only"),
    ] = False,
    target: Annotated[
        int,
        typer.Option("--target", "-n", help="Target sample size"),
    ] = 75_000,
    seed: Annotated[
        int,
        typer.Option("--seed", "-s", help="Random seed for reproducibility"),
    ] = 42,
    output: Annotated[
        Path | None,
        typer.Option("--output", "-o", help="Output sample file path"),
    ] = None,
) -> None:
    """Perform stratified sampling based on atom count bins."""
    if not index_path.exists():
        console.print(f"[red]Error:[/red] Index file not found: {index_path}")
        raise typer.Exit(1)

    # Load index
    console.print(f"Loading [cyan]{index_path}[/cyan]...")
    entries = load_index(index_path)
    console.print(f"Loaded [cyan]{len(entries):,}[/cyan] entries\n")

    # Analyze mode
    if analyze:
        analyze_distribution(entries)
        return

    # Validate target
    if target <= 0:
        console.print("[red]Error:[/red] Target must be positive")
        raise typer.Exit(1)

    if target > len(entries):
        console.print(
            f"[yellow]Warning:[/yellow] Target {target:,} exceeds total {len(entries):,}"
        )
        target = len(entries)

    # Perform sampling
    samples, distribution = stratified_sample(entries, target, seed)

    # Display distribution
    table = Table(title=f"Stratified Sample (target={target:,}, seed={seed})")
    table.add_column("Bin", style="cyan")
    table.add_column("Total", justify="right")
    table.add_column("Sampled", justify="right")
    table.add_column("Rate", justify="right")

    for name, _, _ in BINS:
        dist = distribution[name]
        rate = dist["sampled"] / dist["total"] * 100 if dist["total"] > 0 else 0
        table.add_row(
            name,
            f"{dist['total']:,}",
            f"{dist['sampled']:,}",
            f"{rate:.1f}%",
        )

    table.add_section()
    total_entries = sum(d["total"] for d in distribution.values())
    table.add_row("Total", f"{total_entries:,}", f"{len(samples):,}", "")

    console.print(table)

    # Save output
    if output is None:
        console.print("\n[dim]Use --output to save samples to file[/dim]")
        return

    # Ensure output directory exists
    output.parent.mkdir(parents=True, exist_ok=True)

    sample_data = {
        "version": 1,
        "created": datetime.now(timezone.utc).isoformat(),
        "parameters": {
            "target": target,
            "seed": seed,
            "source": str(index_path),
        },
        "distribution": distribution,
        "total_sampled": len(samples),
        "samples": samples,
    }

    output.write_text(json.dumps(sample_data, indent=2))
    console.print(f"\n[green]Saved:[/green] {output}")
    console.print(f"  Samples: [cyan]{len(samples):,}[/cyan]")


@app.command()
def plot(
    sample_path: Annotated[
        Path,
        typer.Argument(help="Path to sample JSON file"),
    ],
    output: Annotated[
        Path | None,
        typer.Option(
            "--output", "-o", help="Output PNG path (default: plots/dataset/<name>.png)"
        ),
    ] = None,
) -> None:
    """Generate distribution bar chart for a sample file."""
    import matplotlib.pyplot as plt

    if not sample_path.exists():
        console.print(f"[red]Error:[/red] Sample file not found: {sample_path}")
        raise typer.Exit(1)

    # Load sample data
    data = json.loads(sample_path.read_text())
    dist = data.get("distribution", {})
    total = data.get("total_sampled", 0)

    if not dist:
        console.print("[red]Error:[/red] No distribution data in sample file")
        raise typer.Exit(1)

    # Extract data for plotting
    bins = list(dist.keys())
    counts = [dist[b]["sampled"] for b in bins]

    # Create plot
    fig, ax = plt.subplots(figsize=(10, 6))
    x = range(len(bins))
    bars = ax.bar(x, counts, color="#f39c12", alpha=0.8, edgecolor="#e67e22")

    ax.set_xlabel("Structure Size (atoms)", fontsize=11)
    ax.set_ylabel("Number of Structures", fontsize=11)
    ax.set_title(f"Dataset Distribution (n={total:,})", fontsize=13, fontweight="bold")
    ax.set_xticks(x)
    ax.set_xticklabels(bins, rotation=45, ha="right")

    # Add value labels
    for bar in bars:
        height = bar.get_height()
        if height > 0:
            ax.annotate(
                f"{int(height):,}",
                xy=(bar.get_x() + bar.get_width() / 2, height),
                xytext=(0, 3),
                textcoords="offset points",
                ha="center",
                va="bottom",
                fontsize=9,
            )

    plt.tight_layout()

    # Determine output path
    if output is None:
        plots_dir = Path(__file__).parent.parent / "results" / "plots" / "dataset"
        plots_dir.mkdir(parents=True, exist_ok=True)
        output = plots_dir / f"{sample_path.stem}.png"

    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=150, bbox_inches="tight")
    plt.close(fig)

    console.print(f"[green]Saved:[/green] {output}")


if __name__ == "__main__":
    app()
