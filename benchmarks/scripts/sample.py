#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.12"
# dependencies = ["rich>=13.0", "typer>=0.9.0", "matplotlib>=3.0"]
# ///
"""Stratified sampling of benchmark structures by atom count.

Reads a v2 index file (with source info) and performs stratified sampling
using 14 atom count bins. Prefers AFDB (human) entries, fills with PDB.

Usage:
    # Analyze distribution
    ./benchmarks/scripts/sample.py benchmarks/inputs/index.json --analyze

    # Generate sample (150 per bin, seed 42)
    ./benchmarks/scripts/sample.py benchmarks/inputs/index.json \\
        --output benchmarks/dataset/sample.json

    # Plot distribution
    ./benchmarks/scripts/sample.py plot benchmarks/dataset/sample.json
"""

from __future__ import annotations

import json
import random
from pathlib import Path
from typing import Annotated

import typer
from rich.console import Console
from rich.table import Table

app = typer.Typer(help="Stratified sampling of benchmark structures")
console = Console()

# 14-bin scheme matching current sample.json
BINS = [
    ("0-500", 0, 500),
    ("500-1k", 500, 1_000),
    ("1k-2k", 1_000, 2_000),
    ("2k-3k", 2_000, 3_000),
    ("3k-5k", 3_000, 5_000),
    ("5k-10k", 5_000, 10_000),
    ("10k-20k", 10_000, 20_000),
    ("20k-50k", 20_000, 50_000),
    ("50k-75k", 50_000, 75_000),
    ("75k-100k", 75_000, 100_000),
    ("100k-150k", 100_000, 150_000),
    ("150k-200k", 150_000, 200_000),
    ("200k-500k", 200_000, 500_000),
    ("500k+", 500_000, float("inf")),
]


def load_index(index_path: Path) -> dict[str, dict]:
    """Load v2 index file and return entries dict.

    Supports both v1 (flat int values) and v2 (dict with n_atoms/source).
    """
    data = json.loads(index_path.read_text())
    raw_entries = data.get("entries", {})

    entries: dict[str, dict] = {}
    for entry_id, value in raw_entries.items():
        if isinstance(value, int):
            # v1 format: just atom count
            entries[entry_id] = {"n_atoms": value, "source": "pdb"}
        else:
            entries[entry_id] = value

    return entries


def categorize_by_bin(
    entries: dict[str, dict],
) -> dict[str, list[dict]]:
    """Categorize entries into bins. Each bin entry is {id, n_atoms, source}."""
    bins: dict[str, list[dict]] = {name: [] for name, _, _ in BINS}

    for entry_id, info in entries.items():
        n_atoms = info["n_atoms"]
        for name, low, high in BINS:
            if low <= n_atoms < high:
                bins[name].append(
                    {
                        "id": entry_id,
                        "n_atoms": n_atoms,
                        "source": info.get("source", "pdb"),
                    }
                )
                break

    return bins


def analyze_distribution(entries: dict[str, dict]) -> None:
    """Display distribution analysis with source breakdown."""
    bins = categorize_by_bin(entries)

    table = Table(title="Atom Count Distribution")
    table.add_column("Bin", style="cyan")
    table.add_column("Range", style="dim")
    table.add_column("Count", justify="right")
    table.add_column("Human", justify="right", style="green")
    table.add_column("PDB", justify="right", style="blue")
    table.add_column("Percent", justify="right")

    total = len(entries)

    for name, low, high in BINS:
        bin_entries = bins[name]
        count = len(bin_entries)
        human = sum(1 for e in bin_entries if e["source"] == "human")
        pdb = sum(1 for e in bin_entries if e["source"] == "pdb")
        pct = count / total * 100 if total > 0 else 0
        high_str = f"{high:,}" if high != float("inf") else "+"
        range_str = f"{low:,} - {high_str}"
        table.add_row(
            name, range_str, f"{count:,}", f"{human:,}", f"{pdb:,}", f"{pct:.1f}%"
        )

    table.add_section()
    total_human = sum(1 for e in entries.values() if e.get("source") == "human")
    total_pdb = sum(1 for e in entries.values() if e.get("source") == "pdb")
    table.add_row(
        "Total", "", f"{total:,}", f"{total_human:,}", f"{total_pdb:,}", "100.0%"
    )

    console.print(table)


def stratified_sample(
    entries: dict[str, dict],
    n_per_bin: int,
    seed: int,
) -> dict[str, list[dict]]:
    """Perform stratified sampling with AFDB-first preference.

    For each bin:
    1. Prefer source=human entries
    2. Fill remaining slots with source=pdb
    3. If bin has fewer than n_per_bin, take all

    Returns dict mapping bin name to list of {id, n_atoms, source}.
    """
    rng = random.Random(seed)
    bins = categorize_by_bin(entries)
    samples: dict[str, list[dict]] = {}

    for name, _, _ in BINS:
        bin_entries = bins[name]

        if len(bin_entries) <= n_per_bin:
            # Take all, sorted by n_atoms for determinism
            selected = sorted(bin_entries, key=lambda e: (e["n_atoms"], e["id"]))
        else:
            # Split by source
            human = [e for e in bin_entries if e["source"] == "human"]
            pdb = [e for e in bin_entries if e["source"] == "pdb"]

            rng.shuffle(human)
            rng.shuffle(pdb)

            if len(human) >= n_per_bin:
                # Enough human entries
                selected = human[:n_per_bin]
            else:
                # Take all human, fill with PDB
                remaining = n_per_bin - len(human)
                selected = human + pdb[:remaining]

            # Sort for deterministic output
            selected = sorted(selected, key=lambda e: (e["n_atoms"], e["id"]))

        samples[name] = selected

    return samples


@app.command()
def sample(
    index_path: Annotated[
        Path,
        typer.Argument(help="Path to v2 index.json file"),
    ],
    analyze: Annotated[
        bool,
        typer.Option("--analyze", "-a", help="Show distribution analysis only"),
    ] = False,
    n_per_bin: Annotated[
        int,
        typer.Option("--n-per-bin", "-n", help="Maximum samples per bin"),
    ] = 150,
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

    # Perform sampling
    samples = stratified_sample(entries, n_per_bin, seed)
    total_sampled = sum(len(v) for v in samples.values())

    # Display distribution
    table = Table(title=f"Stratified Sample (n_per_bin={n_per_bin}, seed={seed})")
    table.add_column("Bin", style="cyan")
    table.add_column("Available", justify="right")
    table.add_column("Sampled", justify="right")
    table.add_column("Human", justify="right", style="green")
    table.add_column("PDB", justify="right", style="blue")

    bins_categorized = categorize_by_bin(entries)
    for name, _, _ in BINS:
        available = len(bins_categorized[name])
        sampled = samples[name]
        human = sum(1 for e in sampled if e["source"] == "human")
        pdb = sum(1 for e in sampled if e["source"] == "pdb")
        table.add_row(
            name, f"{available:,}", f"{len(sampled):,}", f"{human:,}", f"{pdb:,}"
        )

    table.add_section()
    total_available = len(entries)
    table.add_row("Total", f"{total_available:,}", f"{total_sampled:,}", "", "")
    console.print(table)

    # Save output
    if output is None:
        console.print("\n[dim]Use --output to save samples to file[/dim]")
        return

    output.parent.mkdir(parents=True, exist_ok=True)

    sample_data = {
        "seed": seed,
        "n_per_bin": n_per_bin,
        "bins": {name: low for name, low, _ in BINS},
        "total": total_sampled,
        "samples": samples,
    }

    output.write_text(json.dumps(sample_data, indent=2))
    console.print(f"\n[green]Saved:[/green] {output}")
    console.print(f"  Samples: [cyan]{total_sampled:,}[/cyan]")


# ---------------------------------------------------------------------------
# Plot subcommand
# ---------------------------------------------------------------------------

# Preset colors for plots
COLOR_PRESETS = {
    "orange": ("#f39c12", "#e67e22"),
    "blue": ("#3498db", "#2980b9"),
    "green": ("#2ecc71", "#27ae60"),
    "red": ("#e74c3c", "#c0392b"),
    "purple": ("#9b59b6", "#8e44ad"),
}


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
    color: Annotated[
        str,
        typer.Option(
            "--color", "-c", help="Bar color (orange/blue/green/red/purple or hex code)"
        ),
    ] = "orange",
) -> None:
    """Generate distribution bar chart for a sample file."""
    import matplotlib.pyplot as plt

    if not sample_path.exists():
        console.print(f"[red]Error:[/red] Sample file not found: {sample_path}")
        raise typer.Exit(1)

    # Load sample data
    data = json.loads(sample_path.read_text())
    samples = data.get("samples", {})
    total = data.get("total", 0)

    if not samples:
        console.print("[red]Error:[/red] No samples data in file")
        raise typer.Exit(1)

    # Resolve color
    if color in COLOR_PRESETS:
        fill_color, edge_color = COLOR_PRESETS[color]
    else:
        fill_color = color
        edge_color = color

    # Extract data for plotting
    bin_names = list(samples.keys())
    counts = [len(samples[b]) for b in bin_names]

    # Create plot
    fig, ax = plt.subplots(figsize=(12, 6))
    x = range(len(bin_names))
    bars = ax.bar(x, counts, color=fill_color, alpha=0.8, edgecolor=edge_color)

    ax.set_xlabel("Structure Size (atoms)", fontsize=11)
    ax.set_ylabel("Number of Structures", fontsize=11)
    ax.set_title(f"Dataset Distribution (n={total:,})", fontsize=13, fontweight="bold")
    ax.set_xticks(x)
    ax.set_xticklabels(bin_names, rotation=45, ha="right")

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
        plots_dir = Path(__file__).parent.parent.joinpath("results", "plots", "dataset")
        plots_dir.mkdir(parents=True, exist_ok=True)
        output = plots_dir.joinpath(f"{sample_path.stem}.png")

    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=150, bbox_inches="tight")
    plt.close(fig)

    console.print(f"[green]Saved:[/green] {output}")


if __name__ == "__main__":
    app()
