#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.12"
# dependencies = ["gemmi>=0.7.4", "rich>=13.0", "typer>=0.9.0"]
# ///
"""Build atom count index from AFDB PDB files and/or PDB mmCIF mirror.

Scans actual data sources to build an index mapping structure ID to atom count
and source. Supports two scan modes that can run separately or together:

  --afdb-dir: Scan AlphaFold DB PDB files (fast ATOM record counting)
  --pdb-dir:  Scan PDB mmCIF mirror via gemmi (clean + count)

Usage:
    # Both sources (typical)
    ./benchmarks/scripts/build_index.py --afdb-dir $AFDB_DIR --pdb-dir $PDB_DIR

    # AFDB only
    ./benchmarks/scripts/build_index.py --afdb-dir $AFDB_DIR -o afdb_index.json

    # PDB only (large structures)
    ./benchmarks/scripts/build_index.py --pdb-dir $PDB_DIR --min-atoms 20000
"""

from __future__ import annotations

import json
import os
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor, as_completed
from datetime import datetime, timezone
from pathlib import Path
from typing import Annotated

import gemmi
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

app = typer.Typer(help="Build atom count index from AFDB and/or PDB sources")
console = Console()


# ---------------------------------------------------------------------------
# AFDB scanning
# ---------------------------------------------------------------------------


def _count_afdb_atoms(pdb_path: str) -> tuple[str, int]:
    """Count ATOM records in an AFDB PDB file.

    Returns (entry_id, n_atoms). entry_id is the lowercased stem:
    AF-A0A0A0MT70-F1-model_v6.pdb -> af-a0a0a0mt70-f1-model_v6
    """
    name = os.path.basename(pdb_path)
    entry_id = name.removesuffix(".pdb").lower()
    count = 0
    try:
        with open(pdb_path) as f:
            for line in f:
                if line.startswith("ATOM  "):
                    count += 1
    except Exception:
        return (entry_id, 0)
    return (entry_id, count)


def _scan_afdb(afdb_dir: Path, workers: int) -> dict[str, int]:
    """Scan AFDB directory for PDB files, return {id: n_atoms}.

    Only scans the top-level directory (assumes flat layout).
    """
    files: list[str] = []
    with os.scandir(afdb_dir) as it:
        for entry in it:
            if entry.is_file() and entry.name.endswith(".pdb"):
                files.append(entry.path)

    if not files:
        return {}

    console.print(f"  Found [cyan]{len(files):,}[/cyan] AFDB PDB files")

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
        task = progress.add_task("  Scanning AFDB", total=len(files))

        with ThreadPoolExecutor(max_workers=workers) as executor:
            futures = {executor.submit(_count_afdb_atoms, f): f for f in files}
            for future in as_completed(futures):
                entry_id, n_atoms = future.result()
                if n_atoms > 0:
                    entries[entry_id] = n_atoms
                else:
                    errors += 1
                progress.advance(task)

    if errors > 0:
        console.print(f"  [yellow]AFDB errors: {errors:,}[/yellow]")
    console.print(f"  [green]AFDB indexed: {len(entries):,}[/green]")
    return entries


# ---------------------------------------------------------------------------
# PDB mmCIF scanning
# ---------------------------------------------------------------------------


def _process_cif(args: tuple[str, int]) -> tuple[str, int]:
    """Process a single mmCIF file: read, clean, count atoms.

    Args is (cif_path, min_atoms) for pickling compatibility.
    Returns (entry_id, n_atoms). n_atoms=0 means skip.
    """
    cif_path_str, min_atoms = args
    try:
        st = gemmi.read_structure(cif_path_str)

        # Pre-filter: raw atom count before cleaning
        if len(st) == 0:
            return ("", 0)
        raw_count = st[0].count_atom_sites()
        if raw_count < min_atoms:
            return ("", 0)

        # Clean
        st.setup_entities()
        st.remove_hydrogens()
        st.remove_alternative_conformations()
        st.remove_ligands_and_waters()
        st.remove_empty_chains()

        if len(st) == 0 or len(st[0]) == 0:
            return ("", 0)

        # Filter to L-peptide chains only
        model = st[0]
        chains_to_remove = []
        for chain in model:
            polymer = chain.get_polymer()
            if (
                not polymer
                or polymer.check_polymer_type() != gemmi.PolymerType.PeptideL
            ):
                chains_to_remove.append(chain.name)
        for name in chains_to_remove:
            model.remove_chain(name)

        n_atoms = sum(1 for chain in model for res in chain for _ in res)
        if n_atoms == 0:
            return ("", 0)

        # Extract ID from path: .../3lml.cif.gz -> 3lml
        basename = os.path.basename(cif_path_str)
        entry_id = basename.split(".")[0].lower()
        return (entry_id, n_atoms)

    except Exception:
        return ("", 0)


def _scan_pdb(pdb_dir: Path, min_atoms: int, workers: int) -> dict[str, int]:
    """Scan PDB mmCIF mirror using gemmi.CifWalk, return {id: n_atoms}."""
    # gemmi.CifWalk needs $PDB_DIR set
    os.environ["PDB_DIR"] = str(pdb_dir)

    console.print(f"  Walking mmCIF mirror at [cyan]{pdb_dir}[/cyan]...")
    cif_files: list[str] = []
    for path in gemmi.CifWalk(str(pdb_dir)):
        cif_files.append(path)

    if not cif_files:
        return {}

    console.print(f"  Found [cyan]{len(cif_files):,}[/cyan] CIF files")
    console.print(f"  Pre-filter: raw atoms >= [cyan]{min_atoms:,}[/cyan]")

    entries: dict[str, int] = {}
    errors = 0

    work_items = [(path, min_atoms) for path in cif_files]

    with Progress(
        SpinnerColumn(),
        TextColumn("[progress.description]{task.description}"),
        BarColumn(),
        MofNCompleteColumn(),
        TimeElapsedColumn(),
        TimeRemainingColumn(),
        console=console,
    ) as progress:
        task = progress.add_task("  Scanning PDB", total=len(work_items))

        with ProcessPoolExecutor(max_workers=workers) as executor:
            futures = {executor.submit(_process_cif, item): item for item in work_items}
            for future in as_completed(futures):
                entry_id, n_atoms = future.result()
                if entry_id and n_atoms > 0:
                    entries[entry_id] = n_atoms
                else:
                    errors += 1
                progress.advance(task)

    skipped = errors  # includes pre-filtered + actual errors
    console.print(f"  [dim]PDB skipped/filtered: {skipped:,}[/dim]")
    console.print(f"  [green]PDB indexed: {len(entries):,}[/green]")
    return entries


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


@app.command()
def main(
    afdb_dir: Annotated[
        Path | None,
        typer.Option("--afdb-dir", help="AlphaFold DB PDB directory"),
    ] = None,
    pdb_dir: Annotated[
        Path | None,
        typer.Option("--pdb-dir", help="PDB mmCIF mirror root directory"),
    ] = None,
    output: Annotated[
        Path,
        typer.Option("--output", "-o", help="Output index file path"),
    ] = Path("benchmarks/inputs/index.json"),
    min_atoms: Annotated[
        int,
        typer.Option("--min-atoms", help="Minimum raw atom count for PDB pre-filter"),
    ] = 20_000,
    workers: Annotated[
        int,
        typer.Option("--workers", "-w", help="Number of parallel workers"),
    ] = 8,
) -> None:
    """Build atom count index from AFDB and/or PDB mmCIF sources."""
    if afdb_dir is None and pdb_dir is None:
        console.print(
            "[red]Error:[/red] Specify at least one of --afdb-dir or --pdb-dir"
        )
        raise typer.Exit(1)

    all_entries: dict[str, dict] = {}
    source_counts: dict[str, int] = {}

    # Scan AFDB
    if afdb_dir is not None:
        if not afdb_dir.exists():
            console.print(f"[red]Error:[/red] AFDB directory not found: {afdb_dir}")
            raise typer.Exit(1)
        console.print("[bold]Scanning AFDB...[/bold]")
        afdb_entries = _scan_afdb(afdb_dir, workers)
        for entry_id, n_atoms in afdb_entries.items():
            all_entries[entry_id] = {"n_atoms": n_atoms, "source": "human"}
        source_counts["human"] = len(afdb_entries)
        console.print()

    # Scan PDB
    if pdb_dir is not None:
        if not pdb_dir.exists():
            console.print(f"[red]Error:[/red] PDB directory not found: {pdb_dir}")
            raise typer.Exit(1)
        console.print("[bold]Scanning PDB mmCIF mirror...[/bold]")
        pdb_entries = _scan_pdb(pdb_dir, min_atoms, workers)
        pdb_added = 0
        for entry_id, n_atoms in pdb_entries.items():
            # Don't overwrite AFDB entries if there's a collision
            if entry_id not in all_entries:
                all_entries[entry_id] = {"n_atoms": n_atoms, "source": "pdb"}
                pdb_added += 1
        source_counts["pdb"] = pdb_added
        if pdb_added < len(pdb_entries):
            collisions = len(pdb_entries) - pdb_added
            console.print(
                f"  [dim]PDB collisions with AFDB (skipped): {collisions:,}[/dim]"
            )
        console.print()

    if not all_entries:
        console.print("[red]Error:[/red] No entries found")
        raise typer.Exit(1)

    # Build v2 index
    index = {
        "version": 2,
        "created": datetime.now(timezone.utc).isoformat(),
        "sources": source_counts,
        "total": len(all_entries),
        "entries": dict(sorted(all_entries.items())),
    }

    # Write output
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(index, indent=2))

    console.print("[green]Done![/green]")
    console.print(f"  Total: [cyan]{len(all_entries):,}[/cyan] entries")
    for src, count in sorted(source_counts.items()):
        console.print(f"  {src}: [cyan]{count:,}[/cyan]")
    console.print(f"  Output: [cyan]{output}[/cyan]")


if __name__ == "__main__":
    app()
