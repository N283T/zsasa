#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.11"
# dependencies = ["rich>=13.0", "typer>=0.9.0"]
# ///
"""Build atom count index for benchmark structures.

Scans all .json.gz files in the input directory and creates an index
mapping PDB ID to atom count. Uses parallel processing for speed.

Usage:
    ./benchmarks/scripts/build_index.py benchmarks/inputs
    ./benchmarks/scripts/build_index.py benchmarks/inputs --output index.json
"""

from __future__ import annotations

import gzip
import json
import os
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone
from pathlib import Path
from typing import Annotated

import typer
from rich.console import Console
from rich.progress import (
    BarColumn,
    MofNCompleteColumn,
    Progress,
    SpinnerColumn,
    TextColumn,
    TimeElapsedColumn,
    TimeRemainingColumn,
)

app = typer.Typer(help="Build atom count index for benchmark structures")
console = Console()


def get_atom_count(file_path: Path) -> tuple[str, int]:
    """Read a JSON file and return (pdb_id, atom_count).

    Returns atom_count of 0 on error.
    """
    try:
        # Extract PDB ID from filename
        name = file_path.name
        if name.endswith(".json.gz"):
            pdb_id = name[:-8]
        elif name.endswith(".json"):
            pdb_id = name[:-5]
        else:
            return ("", 0)

        # Read and parse JSON
        if name.endswith(".gz"):
            with gzip.open(file_path, "rt") as f:
                data = json.load(f)
        else:
            with open(file_path) as f:
                data = json.load(f)

        # Count atoms (length of x coordinate array)
        n_atoms = len(data.get("x", []))
        return (pdb_id, n_atoms)

    except Exception:
        return ("", 0)


def scan_files(input_dir: Path) -> list[Path]:
    """Scan directory for JSON files using os.scandir for performance."""
    files = []
    with os.scandir(input_dir) as it:
        for entry in it:
            if entry.is_file() and (
                entry.name.endswith(".json.gz") or entry.name.endswith(".json")
            ):
                files.append(Path(entry.path))
    return sorted(files, key=lambda p: p.name)


@app.command()
def main(
    input_dir: Annotated[
        Path,
        typer.Argument(help="Directory containing .json.gz files"),
    ],
    output: Annotated[
        Path | None,
        typer.Option("--output", "-o", help="Output index file path"),
    ] = None,
    workers: Annotated[
        int,
        typer.Option("--workers", "-w", help="Number of worker threads"),
    ] = 8,
) -> None:
    """Build atom count index from JSON structure files."""
    if not input_dir.exists():
        console.print(f"[red]Error:[/red] Directory not found: {input_dir}")
        raise typer.Exit(1)

    # Default output path
    if output is None:
        output = input_dir / "index.json"

    # Scan for files
    console.print(f"Scanning [cyan]{input_dir}[/cyan]...")
    files = scan_files(input_dir)

    if not files:
        console.print(f"[red]Error:[/red] No .json.gz files found in {input_dir}")
        raise typer.Exit(1)

    console.print(f"Found [cyan]{len(files):,}[/cyan] files")
    console.print(f"Using [cyan]{workers}[/cyan] worker threads\n")

    # Process files in parallel
    entries: dict[str, int] = {}
    errors = 0

    with Progress(
        SpinnerColumn(),
        TextColumn("[progress.description]{task.description}"),
        BarColumn(),
        MofNCompleteColumn(),
        TimeElapsedColumn(),
        TimeRemainingColumn(),
        console=console,
    ) as progress:
        task = progress.add_task("Building index", total=len(files))

        with ThreadPoolExecutor(max_workers=workers) as executor:
            futures = {executor.submit(get_atom_count, f): f for f in files}

            for future in as_completed(futures):
                pdb_id, n_atoms = future.result()
                if pdb_id and n_atoms > 0:
                    entries[pdb_id] = n_atoms
                else:
                    errors += 1
                progress.advance(task)

    # Build index
    index = {
        "version": 1,
        "created": datetime.now(timezone.utc).isoformat(),
        "total_files": len(entries),
        "entries": dict(sorted(entries.items())),
    }

    # Write output
    output.write_text(json.dumps(index, indent=2))

    console.print("\n[green]Done![/green]")
    console.print(f"  Indexed: [cyan]{len(entries):,}[/cyan] structures")
    if errors > 0:
        console.print(f"  Errors: [yellow]{errors:,}[/yellow]")
    console.print(f"  Output: [cyan]{output}[/cyan]")


if __name__ == "__main__":
    app()
