#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.11"
# dependencies = ["rich>=13.0", "typer>=0.9.0", "polars", "matplotlib>=3.8"]
# ///
"""IPC (Instructions Per Cycle) benchmark.

Measures CPU efficiency metrics for representative structures per size bin.

Usage:
    ./benchmarks/scripts/ipc.py run
    ./benchmarks/scripts/ipc.py run --threads 1,4,10
"""

from __future__ import annotations

import gzip
import json
import re
import shutil
import subprocess
import tempfile
from pathlib import Path
from typing import Annotated

import polars as pl
import typer
from rich import print as rprint
from rich.table import Table

app = typer.Typer(help="IPC benchmark")

ROOT = Path(__file__).parent.parent.parent
DATASET_DIR = ROOT / "benchmarks" / "dataset"
INPUTS_DIR = ROOT / "benchmarks" / "inputs"
RESULTS_DIR = ROOT / "benchmarks" / "results" / "ipc"

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


def get_size_bin(n_atoms: int) -> str:
    for low, high, label in BINS:
        if low <= n_atoms < high:
            return label
    return "200k+"


def get_binary_path(tool: str) -> Path:
    if tool == "zig":
        return ROOT / "zig-out" / "bin" / "zsasa"
    elif tool == "freesasa":
        return ROOT / "benchmarks" / "external" / "freesasa-bench" / "src" / "freesasa"
    elif tool == "rust":
        return (
            ROOT
            / "benchmarks"
            / "external"
            / "rustsasa-bench"
            / "target"
            / "release"
            / "rust-sasa"
        )
    raise ValueError(f"Unknown tool: {tool}")


def select_representatives() -> list[dict]:
    """Select one representative structure per bin from inputs."""
    index_path = INPUTS_DIR / "index.json"
    if not index_path.exists():
        # Fallback to dataset
        rprint("[yellow]No index.json, using dataset structures[/yellow]")
        structures = []
        for f in DATASET_DIR.glob("*.json.gz"):
            with gzip.open(f, "rt") as fp:
                data = json.load(fp)
            n_atoms = len(data.get("x", []))
            structures.append(
                {
                    "pdb_id": f.stem.replace(".json", ""),
                    "n_atoms": n_atoms,
                    "path": f,
                }
            )
        return structures

    with open(index_path) as f:
        index = json.load(f)

    # Group by bin - index format is {"entries": {"pdb_id": n_atoms, ...}}
    by_bin: dict[str, list] = {b[2]: [] for b in BINS}
    for pdb_id, n_atoms in index["entries"].items():
        bin_name = get_size_bin(n_atoms)
        by_bin[bin_name].append({"pdb_id": pdb_id, "n_atoms": n_atoms})

    # Select median-sized structure from each bin
    representatives = []
    for bin_name in [b[2] for b in BINS]:
        entries = by_bin[bin_name]
        if not entries:
            continue
        entries.sort(key=lambda x: x["n_atoms"])
        mid = entries[len(entries) // 2]
        # Find the JSON file - could be .json or .json.gz
        json_path = INPUTS_DIR / f"{mid['pdb_id']}.json.gz"
        if not json_path.exists():
            json_path = INPUTS_DIR / f"{mid['pdb_id']}.json"
        if not json_path.exists():
            continue
        representatives.append(
            {
                "pdb_id": mid["pdb_id"],
                "n_atoms": mid["n_atoms"],
                "bin": bin_name,
                "path": json_path,
            }
        )

    return representatives


def measure_ipc(tool: str, json_path: Path, threads: int) -> dict | None:
    """Run benchmark with /usr/bin/time -l and extract IPC metrics."""
    binary = get_binary_path(tool)
    if not binary.exists():
        return None

    # Prepare input (decompress if needed)
    input_path = json_path
    cleanup = False
    if str(json_path).endswith(".gz"):
        with tempfile.NamedTemporaryFile(suffix=".json", delete=False) as f:
            input_path = Path(f.name)
        cleanup = True
        with gzip.open(json_path, "rb") as f_in:
            with open(input_path, "wb") as f_out:
                shutil.copyfileobj(f_in, f_out)

    try:
        with tempfile.NamedTemporaryFile(suffix=".json", delete=False) as f:
            output_path = Path(f.name)

        if tool == "zig":
            cmd = [
                str(binary),
                "--algorithm=sr",
                f"--threads={threads}",
                str(input_path),
                str(output_path),
            ]
        elif tool == "freesasa":
            cmd = [
                str(binary),
                "--json-input",
                str(input_path),
                "--shrake-rupley",
                "--resolution=100",
                f"--n-threads={threads}",
            ]
        elif tool == "rust":
            cmd = [
                str(binary),
                "-J",
                str(input_path),
                "-n",
                "100",
                "-t",
                str(threads),
            ]
        else:
            return None

        result = subprocess.run(
            ["/usr/bin/time", "-l"] + cmd,
            capture_output=True,
            text=True,
            timeout=600,
        )

        # Parse metrics from stderr (time -l output)
        stderr = result.stderr
        instructions = 0
        cycles = 0

        for line in stderr.split("\n"):
            if "instructions retired" in line:
                match = re.search(r"(\d+)", line)
                if match:
                    instructions = int(match.group(1))
            elif "cycles elapsed" in line:
                match = re.search(r"(\d+)", line)
                if match:
                    cycles = int(match.group(1))

        if cycles == 0:
            return None

        return {
            "instructions": instructions,
            "cycles": cycles,
            "ipc": instructions / cycles,
        }

    finally:
        if cleanup:
            input_path.unlink(missing_ok=True)
        output_path.unlink(missing_ok=True)


@app.callback(invoke_without_command=True)
def run(
    threads: Annotated[str, typer.Option("--threads", "-t")] = "1,10",
    tools: Annotated[str, typer.Option("--tools")] = "zig,freesasa,rust",
):
    """Run IPC benchmark for representative structures."""
    thread_list = [int(t) for t in threads.split(",")]
    tool_list = [t.strip() for t in tools.split(",")]

    rprint("[bold]Selecting representative structures per bin...[/bold]")
    representatives = select_representatives()
    rprint(f"Selected {len(representatives)} structures\n")

    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    all_results = []

    for struct in representatives:
        pdb_id = struct["pdb_id"]
        n_atoms = struct["n_atoms"]
        bin_name = struct.get("bin", get_size_bin(n_atoms))
        path = struct["path"]

        rprint(f"[cyan]{bin_name}[/cyan]: {pdb_id} ({n_atoms:,} atoms)")

        for t in thread_list:
            for tool in tool_list:
                if tool == "rust" and t > 1:
                    # Rust's threading is different, skip multi-thread for now
                    pass
                metrics = measure_ipc(tool, path, t)
                if metrics:
                    all_results.append(
                        {
                            "bin": bin_name,
                            "pdb_id": pdb_id,
                            "n_atoms": n_atoms,
                            "tool": tool,
                            "threads": t,
                            **metrics,
                        }
                    )
                    rprint(f"  {tool} t={t}: IPC={metrics['ipc']:.2f}")

    # Save results
    if all_results:
        df = pl.DataFrame(all_results)
        csv_path = RESULTS_DIR / "results.csv"
        df.write_csv(csv_path)
        rprint(f"\n[green]Saved:[/green] {csv_path}")

        # Summary table
        _print_summary(all_results, thread_list, tool_list)


def _print_summary(results: list[dict], thread_list: list[int], tool_list: list[str]):
    """Print summary table."""
    for t in thread_list:
        table = Table(title=f"IPC Summary (threads={t})")
        table.add_column("Bin", style="cyan")
        table.add_column("Atoms", justify="right")
        for tool in tool_list:
            table.add_column(f"{tool.capitalize()} IPC", justify="right")
            table.add_column("Inst (B)", justify="right")

        bin_order = [b[2] for b in BINS]
        by_bin = {}
        for r in results:
            if r["threads"] != t:
                continue
            key = r["bin"]
            if key not in by_bin:
                by_bin[key] = {"n_atoms": r["n_atoms"]}
            by_bin[key][r["tool"]] = r

        for bin_name in bin_order:
            if bin_name not in by_bin:
                continue
            data = by_bin[bin_name]
            row = [bin_name, f"{data['n_atoms']:,}"]
            for tool in tool_list:
                if tool in data:
                    ipc = data[tool]["ipc"]
                    inst_b = data[tool]["instructions"] / 1e9
                    # Color IPC
                    if ipc >= 3.5:
                        row.append(f"[green]{ipc:.2f}[/green]")
                    elif ipc >= 3.0:
                        row.append(f"[yellow]{ipc:.2f}[/yellow]")
                    else:
                        row.append(f"[red]{ipc:.2f}[/red]")
                    row.append(f"{inst_b:.2f}")
                else:
                    row.extend(["-", "-"])
            table.add_row(*row)

        rprint("\n")
        rprint(table)

    rprint("\n[dim]IPC = Instructions Per Cycle. Higher = better CPU utilization[/dim]")
    rprint(
        "[dim]Inst = Total instructions (billions). Lower = more efficient code[/dim]"
    )


@app.command()
def plot():
    """Generate IPC and instruction count plots."""
    import matplotlib.pyplot as plt

    csv_path = RESULTS_DIR / "results.csv"
    if not csv_path.exists():
        rprint("[red]No results.csv found. Run the benchmark first.[/red]")
        raise typer.Exit(1)

    df = pl.read_csv(csv_path)

    # Use t=1 for cleaner comparison
    df_t1 = df.filter(pl.col("threads") == 1)

    bin_order = [b[2] for b in BINS]
    tools = ["zig", "freesasa", "rust"]
    colors = {"zig": "#E85D04", "freesasa": "#3D5A80", "rust": "#8B4513"}
    labels = {"zig": "Zig", "freesasa": "FreeSASA (C)", "rust": "RustSASA"}

    # Pivot data
    inst_data = {}
    ipc_data = {}
    atoms_data = {}
    for tool in tools:
        tool_df = df_t1.filter(pl.col("tool") == tool).sort(
            pl.col("bin").cast(pl.Categorical)
        )
        # Sort by bin order
        rows = {r["bin"]: r for r in tool_df.to_dicts()}
        inst_data[tool] = [
            rows.get(b, {}).get("instructions", 0) / 1e9 for b in bin_order
        ]
        ipc_data[tool] = [rows.get(b, {}).get("ipc", 0) for b in bin_order]
        if tool == "zig":
            atoms_data = {b: rows.get(b, {}).get("n_atoms", 0) for b in bin_order}

    # Filter to bins with data
    valid_bins = [b for b in bin_order if atoms_data.get(b, 0) > 0]

    # --- Plot 1: Instruction Count Comparison (Main) ---
    fig, ax = plt.subplots(figsize=(12, 7))

    x = range(len(valid_bins))
    width = 0.25

    for i, tool in enumerate(tools):
        values = [inst_data[tool][bin_order.index(b)] for b in valid_bins]
        bars = ax.bar(
            [xi + i * width for xi in x],
            values,
            width,
            label=labels[tool],
            color=colors[tool],
        )
        # Add ratio labels for FreeSASA
        if tool == "freesasa":
            zig_values = [inst_data["zig"][bin_order.index(b)] for b in valid_bins]
            for bar, zig_v, fs_v in zip(bars, zig_values, values):
                if zig_v > 0:
                    ratio = fs_v / zig_v
                    ax.annotate(
                        f"{ratio:.1f}x",
                        xy=(bar.get_x() + bar.get_width() / 2, bar.get_height()),
                        ha="center",
                        va="bottom",
                        fontsize=8,
                        color="#3D5A80",
                    )

    ax.set_xlabel("Size Bin (atoms)")
    ax.set_ylabel("Instructions (Billions)")
    ax.set_title(
        "Instruction Count Comparison (Single Thread)\nLower = More Efficient Code",
        fontsize=14,
    )
    ax.set_xticks([xi + width for xi in x])
    ax.set_xticklabels(valid_bins, rotation=45, ha="right")
    ax.legend()
    ax.grid(axis="y", alpha=0.3)

    plt.tight_layout()
    out_path = RESULTS_DIR / "instructions.png"
    plt.savefig(out_path, dpi=150)
    rprint(f"[green]Saved:[/green] {out_path}")
    plt.close()

    # --- Plot 2: Instruction Ratio (Zig = 1.0 baseline) ---
    fig, ax = plt.subplots(figsize=(10, 6))

    x = range(len(valid_bins))
    width = 0.35

    zig_values = [inst_data["zig"][bin_order.index(b)] for b in valid_bins]
    freesasa_ratios = []
    rust_ratios = []

    for b in valid_bins:
        idx = bin_order.index(b)
        zig_v = inst_data["zig"][idx]
        if zig_v > 0:
            freesasa_ratios.append(inst_data["freesasa"][idx] / zig_v)
            rust_ratios.append(inst_data["rust"][idx] / zig_v)
        else:
            freesasa_ratios.append(0)
            rust_ratios.append(0)

    ax.bar(
        [xi - width / 2 for xi in x],
        freesasa_ratios,
        width,
        label="FreeSASA (C)",
        color=colors["freesasa"],
    )
    ax.bar(
        [xi + width / 2 for xi in x],
        rust_ratios,
        width,
        label="RustSASA",
        color=colors["rust"],
    )
    ax.axhline(
        y=1.0, color=colors["zig"], linestyle="--", linewidth=2, label="Zig (baseline)"
    )

    ax.set_xlabel("Size Bin (atoms)")
    ax.set_ylabel("Instruction Ratio (vs Zig)")
    ax.set_title("Zig Executes Fewer Instructions\n(Zig = 1.0 baseline)", fontsize=14)
    ax.set_xticks(x)
    ax.set_xticklabels(valid_bins, rotation=45, ha="right")
    ax.legend()
    ax.grid(axis="y", alpha=0.3)

    # Add percentage labels
    for i, (fs, rs) in enumerate(zip(freesasa_ratios, rust_ratios)):
        if fs > 0:
            ax.annotate(
                f"{fs:.1f}x",
                xy=(i - width / 2, fs),
                ha="center",
                va="bottom",
                fontsize=9,
            )
        if rs > 0:
            ax.annotate(
                f"{rs:.1f}x",
                xy=(i + width / 2, rs),
                ha="center",
                va="bottom",
                fontsize=9,
            )

    plt.tight_layout()
    out_path = RESULTS_DIR / "instruction_ratio.png"
    plt.savefig(out_path, dpi=150)
    rprint(f"[green]Saved:[/green] {out_path}")
    plt.close()

    # --- Plot 3: IPC Comparison ---
    fig, ax = plt.subplots(figsize=(10, 6))

    x = range(len(valid_bins))
    width = 0.25

    for i, tool in enumerate(tools):
        values = [ipc_data[tool][bin_order.index(b)] for b in valid_bins]
        ax.bar(
            [xi + i * width for xi in x],
            values,
            width,
            label=labels[tool],
            color=colors[tool],
        )

    ax.set_xlabel("Size Bin (atoms)")
    ax.set_ylabel("IPC (Instructions Per Cycle)")
    ax.set_title(
        "IPC Comparison (Single Thread)\nHigher = Better CPU Utilization", fontsize=14
    )
    ax.set_xticks([xi + width for xi in x])
    ax.set_xticklabels(valid_bins, rotation=45, ha="right")
    ax.legend()
    ax.grid(axis="y", alpha=0.3)
    ax.set_ylim(0, 4.5)

    plt.tight_layout()
    out_path = RESULTS_DIR / "ipc.png"
    plt.savefig(out_path, dpi=150)
    rprint(f"[green]Saved:[/green] {out_path}")
    plt.close()

    rprint("\n[bold]Generated plots:[/bold]")
    rprint("  • instructions.png - Raw instruction counts (main)")
    rprint("  • instruction_ratio.png - Zig baseline comparison")
    rprint("  • ipc.png - IPC comparison")


if __name__ == "__main__":
    app()
