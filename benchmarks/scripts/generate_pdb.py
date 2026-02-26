#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.12"
# dependencies = [
#     "gemmi>=0.7.4",
#     "rich>=13.0",
#     "typer>=0.9",
# ]
# ///
"""Convert mmCIF files to cleaned PDB files for benchmarking.

Accepts either a sample.json file or a directory of mmCIF files.
When using sample.json, resolves each entry by source:
- source=human (AFDB): Copy PDB from AlphaFold DB directory (already clean)
- source=pdb: Read mmCIF, clean up (remove H, altconf, ligands/waters, L-peptide only), write PDB

Usage:
    # From sample file (resolves AFDB + PDB sources)
    ./benchmarks/scripts/generate_pdb.py benchmarks/dataset/sample.json benchmarks/dataset/pdb/

    # From a directory of mmCIF files (all treated as PDB source)
    ./benchmarks/scripts/generate_pdb.py /path/to/cif/dir benchmarks/dataset/pdb/
"""

from __future__ import annotations

import multiprocessing as mp
import os
import shutil
from concurrent.futures import ProcessPoolExecutor, as_completed
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
)
from rich.table import Table

console = Console()

# Default AFDB directory (AlphaFold Human Proteome v6)
# Overridden at runtime by --afdb-dir CLI option
_DEFAULT_AFDB_DIR = Path(
    os.environ.get(
        "AFDB_DIR",
        "/data/alphafold/UP000005640_9606_HUMAN_v6/pdb",
    )
)


def _count_pdb_atoms(pdb_path: Path) -> int:
    """Count ATOM records in a PDB file."""
    count = 0
    with open(pdb_path) as f:
        for line in f:
            if line.startswith("ATOM  "):
                count += 1
    return count


def afdb_id_to_filename(entry_id: str) -> str:
    """Convert AFDB ID to PDB filename.

    af-a0a0a0mt70-f1-model_v6 -> AF-A0A0A0MT70-F1-model_v6.pdb
    """
    parts = entry_id.split("-")
    converted = "-".join(p.upper() if i < 3 else p for i, p in enumerate(parts))
    return f"{converted}.pdb"


def resolve_cif_path(entry_id: str) -> Path:
    """Resolve PDB code to mmCIF path using gemmi and $PDB_DIR.

    Uses gemmi.expand_if_pdb_code which follows the BioJava convention:
    3lml -> $PDB_DIR/structures/divided/mmCIF/lm/3lml.cif.gz

    Requires $PDB_DIR environment variable to be set.
    """
    return Path(gemmi.expand_if_pdb_code(entry_id))


def clean_cif_to_pdb(cif_path: Path, output_path: Path) -> tuple[str, int]:
    """Read mmCIF, clean up, write PDB. Returns (id, n_atoms).

    If the structure has no atoms after cleaning (H removal, L-peptide filter),
    no output file is written and n_atoms is 0.
    """
    st = gemmi.read_structure(str(cif_path))
    st.setup_entities()
    st.remove_hydrogens()
    st.remove_alternative_conformations()
    st.remove_ligands_and_waters()
    st.remove_empty_chains()

    if len(st) == 0:
        return (cif_path.stem.replace(".cif", ""), 0)

    # Filter to L-peptide chains only
    model = st[0]
    chains_to_remove = []
    for chain in model:
        polymer = chain.get_polymer()
        if not polymer or polymer.check_polymer_type() != gemmi.PolymerType.PeptideL:
            chains_to_remove.append(chain.name)

    for name in chains_to_remove:
        model.remove_chain(name)

    # Count atoms in cleaned model (all HETATM/ligands already removed above)
    n_atoms = sum(1 for chain in model for res in chain for _ in res)
    if n_atoms == 0:
        return (cif_path.stem.replace(".cif", ""), 0)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    gemmi.shorten_chain_names(st)
    st.write_pdb(str(output_path))

    return (cif_path.stem.replace(".cif", ""), n_atoms)


def process_entry(
    args: tuple[str, str, Path, Path],
) -> tuple[str, int, str]:
    """Process a single sample entry. Returns (id, n_atoms, status).

    Status is one of: "skipped" (already exists), "copied" (AFDB source),
    "processed" (CIF converted), "empty" (no atoms after cleaning),
    or "error: <detail>".

    PDB source entries use gemmi.expand_if_pdb_code (reads $PDB_DIR env var).
    """
    entry_id, source, output_path, afdb_dir = args

    try:
        if output_path.exists():
            n_atoms = _count_pdb_atoms(output_path)
            return (entry_id, n_atoms, "skipped")

        if source == "human":
            # AFDB: copy PDB directly
            filename = afdb_id_to_filename(entry_id)
            src_path = afdb_dir.joinpath(filename)
            if not src_path.exists():
                return (entry_id, 0, f"error: AFDB file not found: {src_path}")
            output_path.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(src_path, output_path)
            n_atoms = _count_pdb_atoms(output_path)
            return (entry_id, n_atoms, "copied")

        elif source == "pdb":
            # PDB: resolve via gemmi ($PDB_DIR), preprocess mmCIF -> PDB
            cif_path = resolve_cif_path(entry_id)
            if not cif_path.exists():
                return (entry_id, 0, f"error: CIF not found: {cif_path}")
            _, n_atoms = clean_cif_to_pdb(cif_path, output_path)
            if n_atoms == 0:
                return (entry_id, 0, "empty")
            return (entry_id, n_atoms, "processed")

        else:
            return (entry_id, 0, f"error: unknown source: {source}")

    except Exception as e:
        return (entry_id, 0, f"error: {type(e).__name__}: {e}")


def _process_dir_entry(args: tuple[Path, Path]) -> tuple[str, int, str]:
    """Process a single CIF file for directory mode. Returns (id, n_atoms, status)."""
    cif_path, output_path = args
    try:
        if output_path.exists():
            n_atoms = _count_pdb_atoms(output_path)
            return (cif_path.stem, n_atoms, "skipped")
        entry_id, n_atoms = clean_cif_to_pdb(cif_path, output_path)
        if n_atoms == 0:
            return (entry_id, 0, "empty")
        return (entry_id, n_atoms, "processed")
    except Exception as e:
        return (cif_path.stem, 0, f"error: {type(e).__name__}: {e}")


app = typer.Typer(help=__doc__)


@app.command()
def generate(
    input_path: Annotated[
        Path,
        typer.Argument(help="sample.json file or directory of mmCIF files"),
    ],
    output_dir: Annotated[
        Path,
        typer.Argument(help="Output directory for PDB files"),
    ],
    workers: Annotated[
        int,
        typer.Option("--workers", "-w", help="Number of parallel workers (0 = auto)"),
    ] = 0,
    afdb_dir: Annotated[
        Path | None,
        typer.Option("--afdb-dir", help="AFDB PDB directory (overrides $AFDB_DIR)"),
    ] = None,
    pdb_dir: Annotated[
        Path | None,
        typer.Option(
            "--pdb-dir",
            help="PDB mirror root directory (overrides $PDB_DIR). "
            "Uses gemmi convention: $PDB_DIR/structures/divided/mmCIF/{mid2}/{id}.cif.gz",
        ),
    ] = None,
) -> None:
    """Generate cleaned PDB files from sample.json or mmCIF directory."""

    resolved_afdb_dir = afdb_dir if afdb_dir is not None else _DEFAULT_AFDB_DIR

    # Set $PDB_DIR for gemmi.expand_if_pdb_code (inherited by worker processes)
    if pdb_dir is not None:
        os.environ["PDB_DIR"] = str(pdb_dir)

    if not input_path.exists():
        console.print(f"[red]Error: Path not found: {input_path}[/red]")
        raise typer.Exit(1)

    if workers <= 0:
        workers = max(1, mp.cpu_count() - 1)

    output_dir.mkdir(parents=True, exist_ok=True)

    if input_path.is_file() and input_path.suffix == ".json":
        _process_sample_file(input_path, output_dir, workers, resolved_afdb_dir)
    elif input_path.is_dir():
        _process_directory(input_path, output_dir, workers)
    else:
        console.print(f"[red]Error: Unsupported input: {input_path}[/red]")
        raise typer.Exit(1)


def _process_sample_file(
    sample_path: Path,
    output_dir: Path,
    workers: int,
    afdb_dir: Path,
) -> None:
    """Process entries from sample.json."""
    import json

    with open(sample_path) as f:
        data = json.load(f)

    if "samples" not in data:
        console.print("[red]Error: Invalid sample file (missing 'samples' key)[/red]")
        raise typer.Exit(1)

    # Collect all entries with source info
    work_items: list[tuple[str, str, Path, Path]] = []
    samples = data["samples"]

    if isinstance(samples, list):
        # v1 format: flat list of IDs (no source info, defaults to mmCIF->PDB conversion)
        for entry_id in samples:
            output_path = output_dir.joinpath(f"{entry_id}.pdb")
            work_items.append((entry_id, "pdb", output_path, afdb_dir))
    else:
        # v2 format: dict of bins with entries
        for bin_name, entries in samples.items():
            for i, entry in enumerate(entries):
                if "id" not in entry:
                    console.print(
                        f"[red]Error: Entry {i} in bin '{bin_name}' missing 'id' key[/red]"
                    )
                    raise typer.Exit(1)
                entry_id = entry["id"]
                source = entry.get("source", "pdb")
                output_path = output_dir.joinpath(f"{entry_id}.pdb")
                work_items.append((entry_id, source, output_path, afdb_dir))

    if not work_items:
        console.print("[yellow]No entries to process[/yellow]")
        return

    source_counts: dict[str, int] = {}
    for _, source, _, _ in work_items:
        source_counts[source] = source_counts.get(source, 0) + 1

    console.print(
        f"[bold]Processing {len(work_items):,} entries from {sample_path}[/bold]"
    )
    for source, count in sorted(source_counts.items()):
        console.print(f"  {source}: {count:,}")
    console.print(f"\nWorkers: [cyan]{workers}[/cyan]")
    console.print(f"Output: [cyan]{output_dir}[/cyan]\n")

    _run_parallel(work_items, workers, output_dir)


def _process_directory(input_dir: Path, output_dir: Path, workers: int) -> None:
    """Process all mmCIF files in directory."""
    console.print(f"[bold]Scanning {input_dir}...[/bold]")

    cif_files = []
    for path in gemmi.CifWalk(str(input_dir)):
        cif_files.append(Path(path))

    if not cif_files:
        console.print("[yellow]No CIF files found[/yellow]")
        return

    console.print(f"Found [cyan]{len(cif_files):,}[/cyan] CIF files")
    console.print(f"Workers: [cyan]{workers}[/cyan]\n")

    # Directory mode uses clean_cif_to_pdb directly (no source routing via process_entry)
    _run_parallel_dir(cif_files, output_dir, workers)


def _run_parallel_dir(cif_files: list[Path], output_dir: Path, workers: int) -> None:
    """Process CIF files from directory in parallel."""
    success_count = 0
    skip_count = 0
    error_count = 0
    total_atoms = 0

    work_items: list[tuple[Path, Path]] = []
    for cif_path in cif_files:
        stem = cif_path.stem.replace(".cif", "").lower()
        output_path = output_dir.joinpath(f"{stem}.pdb")
        work_items.append((cif_path, output_path))

    with Progress(
        SpinnerColumn(),
        TextColumn("[progress.description]{task.description}"),
        BarColumn(),
        MofNCompleteColumn(),
        TimeElapsedColumn(),
        console=console,
    ) as progress:
        task = progress.add_task("Processing", total=len(work_items))

        with ProcessPoolExecutor(max_workers=workers) as executor:
            futures = {
                executor.submit(_process_dir_entry, item): item for item in work_items
            }

            for future in as_completed(futures):
                try:
                    entry_id, n_atoms, status = future.result()
                except Exception as exc:
                    item = futures[future]
                    console.print(f"  [red]{item[0].stem}: worker crashed: {exc}[/red]")
                    error_count += 1
                    progress.advance(task)
                    continue

                progress.advance(task)

                if status in ("processed", "copied"):
                    success_count += 1
                    total_atoms += n_atoms
                elif status == "skipped":
                    skip_count += 1
                    total_atoms += n_atoms
                elif status == "empty":
                    skip_count += 1
                else:
                    error_count += 1
                    console.print(f"  [red]{entry_id}: {status}[/red]")

    _print_summary(success_count, skip_count, error_count, total_atoms, output_dir)


def _run_parallel(
    work_items: list[tuple[str, str, Path, Path]],
    workers: int,
    output_dir: Path,
) -> None:
    """Run parallel processing for sample-based entries."""
    success_count = 0
    skip_count = 0
    error_count = 0
    total_atoms = 0

    with Progress(
        SpinnerColumn(),
        TextColumn("[progress.description]{task.description}"),
        BarColumn(),
        MofNCompleteColumn(),
        TimeElapsedColumn(),
        console=console,
    ) as progress:
        task = progress.add_task("Processing", total=len(work_items))

        with ProcessPoolExecutor(max_workers=workers) as executor:
            futures = {
                executor.submit(process_entry, item): item for item in work_items
            }

            for future in as_completed(futures):
                try:
                    entry_id, n_atoms, status = future.result()
                except Exception as exc:
                    item = futures[future]
                    console.print(f"  [red]{item[0]}: worker crashed: {exc}[/red]")
                    error_count += 1
                    progress.advance(task)
                    continue

                progress.advance(task)

                if status in ("processed", "copied"):
                    success_count += 1
                    total_atoms += n_atoms
                elif status == "skipped":
                    skip_count += 1
                    total_atoms += n_atoms
                elif status == "empty":
                    skip_count += 1
                else:
                    error_count += 1
                    console.print(f"  [red]{entry_id}: {status}[/red]")

    _print_summary(success_count, skip_count, error_count, total_atoms, output_dir)


def _print_summary(
    success: int, skipped: int, errors: int, total_atoms: int, output_dir: Path
) -> None:
    """Print processing summary."""
    console.print()
    table = Table(title="Summary")
    table.add_column("Metric", style="cyan")
    table.add_column("Value", justify="right")

    table.add_row("Processed", f"{success:,}")
    table.add_row("Skipped (existing/empty)", f"{skipped:,}")
    table.add_row("Errors", f"{errors:,}")
    table.add_row("Total atoms", f"{total_atoms:,}")

    console.print(table)
    console.print(f"\n[green]Done![/green] Output saved to: {output_dir}")


if __name__ == "__main__":
    app()
