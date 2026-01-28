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

    sr_dist = sr_data["distribution"]
    lr_dist = lr_data["distribution"]

    bins = list(sr_dist.keys())
    sr_counts = [sr_dist[b]["sampled"] for b in bins]
    lr_counts = [lr_dist[b]["sampled"] for b in bins]

    # Create 2-grid figure
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))

    x = range(len(bins))

    # SR plot (left)
    bars1 = ax1.bar(x, sr_counts, color="#2ecc71", alpha=0.8, edgecolor="#27ae60")
    ax1.set_xlabel("Structure Size (atoms)", fontsize=11)
    ax1.set_ylabel("Number of Structures", fontsize=11)
    ax1.set_title(
        f"Shrake-Rupley (n={sum(sr_counts):,})", fontsize=13, fontweight="bold"
    )
    ax1.set_xticks(x)
    ax1.set_xticklabels(bins, rotation=45, ha="right")

    # LR plot (right)
    bars2 = ax2.bar(x, lr_counts, color="#3498db", alpha=0.8, edgecolor="#2980b9")
    ax2.set_xlabel("Structure Size (atoms)", fontsize=11)
    ax2.set_ylabel("Number of Structures", fontsize=11)
    ax2.set_title(
        f"Lee-Richards (n={sum(lr_counts):,})", fontsize=13, fontweight="bold"
    )
    ax2.set_xticks(x)
    ax2.set_xticklabels(bins, rotation=45, ha="right")

    # Add value labels on bars
    def add_labels(ax, bars):
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

    add_labels(ax1, bars1)
    add_labels(ax2, bars2)

    plt.tight_layout()

    # Save
    output_dir.mkdir(parents=True, exist_ok=True)
    output_path = output_dir / "dataset_distribution.png"
    plt.savefig(output_path, dpi=150, bbox_inches="tight")
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
