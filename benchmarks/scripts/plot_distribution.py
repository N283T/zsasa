#!/usr/bin/env python3
# /// script
# requires-python = ">=3.11"
# dependencies = ["matplotlib", "rich"]
# ///
"""Generate dataset distribution plots for benchmark documentation."""

import json
from pathlib import Path

import matplotlib.pyplot as plt


def load_sample(path: Path) -> dict:
    """Load sample JSON file."""
    with open(path) as f:
        return json.load(f)


def plot_distribution(samples_dir: Path, output_dir: Path):
    """Create distribution bar chart comparing SR (100k) and LR (30k) datasets."""
    sr_data = load_sample(samples_dir / "stratified_100k.json")
    lr_data = load_sample(samples_dir / "stratified_30k.json")

    # Extract distribution - use more granular bins for display
    # The sample files use 6 bins, but we want 10 bins for the chart
    # We'll show what we have from the sample files
    sr_dist = sr_data["distribution"]
    lr_dist = lr_data["distribution"]

    # Bins from sample files
    bins = list(sr_dist.keys())
    sr_counts = [sr_dist[b]["sampled"] for b in bins]
    lr_counts = [lr_dist[b]["sampled"] for b in bins]

    # Create figure
    fig, ax = plt.subplots(figsize=(10, 6))

    x = range(len(bins))
    width = 0.35

    bars1 = ax.bar([i - width/2 for i in x], sr_counts, width, label=f'SR (~100k)', color='#2ecc71', alpha=0.8)
    bars2 = ax.bar([i + width/2 for i in x], lr_counts, width, label=f'LR (~30k)', color='#3498db', alpha=0.8)

    ax.set_xlabel('Structure Size (atoms)', fontsize=12)
    ax.set_ylabel('Number of Structures', fontsize=12)
    ax.set_title('Benchmark Dataset Distribution', fontsize=14, fontweight='bold')
    ax.set_xticks(x)
    ax.set_xticklabels(bins, rotation=45, ha='right')
    ax.legend()

    # Add value labels on bars
    def add_labels(bars):
        for bar in bars:
            height = bar.get_height()
            if height > 0:
                ax.annotate(f'{int(height):,}',
                           xy=(bar.get_x() + bar.get_width() / 2, height),
                           xytext=(0, 3),
                           textcoords="offset points",
                           ha='center', va='bottom', fontsize=8)

    add_labels(bars1)
    add_labels(bars2)

    # Add totals
    sr_total = sum(sr_counts)
    lr_total = sum(lr_counts)
    ax.text(0.98, 0.95, f'SR Total: {sr_total:,}\nLR Total: {lr_total:,}',
            transform=ax.transAxes, ha='right', va='top',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5),
            fontsize=10)

    plt.tight_layout()

    # Save
    output_dir.mkdir(parents=True, exist_ok=True)
    output_path = output_dir / "dataset_distribution.png"
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close()

    print(f"✓ Saved: {output_path}")
    return output_path


def main():
    base = Path(__file__).parent.parent
    samples_dir = base / "samples"
    output_dir = base / "results" / "plots" / "dataset"

    plot_distribution(samples_dir, output_dir)


if __name__ == "__main__":
    main()
