#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.11"
# dependencies = [
#     "freesasa>=2.0",
#     "gemmi>=0.7.4",
#     "rich>=13.0",
#     "typer>=0.9",
# ]
# ///
"""Convert CIF files to minimal JSON with ProtOr radii for SASA calculation.

Creates minimal JSON files (x, y, z, r only) suitable for benchmarking.
- Protein atoms only (HETATM excluded)
- Hydrogens excluded
- ProtOr radii pre-applied (via freesasa classifier)
- Output is gzip compressed

Usage:
    # Process directory (recursive, parallel)
    ./benchmarks/scripts/generate_json.py /path/to/cif/files /path/to/output

    # Create tar archive for distribution (no double compression)
    ./benchmarks/scripts/generate_json.py /path/to/cif/files output.tar --archive

    # Split into multiple archives (10,000 files each)
    ./benchmarks/scripts/generate_json.py /path/to/cif/files output --archive --split 10000

    # Process single file
    ./benchmarks/scripts/generate_json.py --file input.cif output.json.gz
"""

from __future__ import annotations

import gzip
import json
import multiprocessing as mp
import tarfile
import tempfile
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import Annotated

import freesasa
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


def get_protor_radius(residue: str, atom_name: str) -> float:
    """Get ProtOr radius for an atom."""
    # Create classifier per call (for multiprocessing safety)
    classifier = freesasa.Classifier.getStandardClassifier("protor")
    radius = classifier.radius(residue, atom_name)

    # Fallback to element-based if not found
    if radius <= 0:
        element = atom_name.strip()[0] if atom_name.strip() else "C"
        fallback = {
            "C": 1.70,
            "N": 1.55,
            "O": 1.52,
            "S": 1.80,
            "P": 1.80,
        }
        radius = fallback.get(element, 1.70)

    return radius


def cif_to_benchmark_json(cif_path: Path, output_path: Path) -> tuple[str, int, int]:
    """Convert CIF to minimal benchmark JSON with ProtOr radii.

    Returns (pdb_id, n_atoms, file_size_bytes).
    """
    try:
        st = gemmi.read_structure(str(cif_path))

        # Setup entities to populate entity_type
        st.setup_entities()

        # Clean up
        st.remove_hydrogens()
        st.remove_alternative_conformations()
        st.remove_ligands_and_waters()
        st.remove_empty_chains()

        # Skip if no models
        if len(st) == 0:
            return (cif_path.stem, 0, 0)

        xs: list[float] = []
        ys: list[float] = []
        zs: list[float] = []
        rs: list[float] = []

        model = st[0]

        for chain in model:
            polymer = chain.get_polymer()
            if not polymer:
                continue

            # Only process L-peptides (standard proteins)
            if polymer.check_polymer_type() != gemmi.PolymerType.PeptideL:
                continue

            for residue in polymer:
                for atom in residue:
                    radius = get_protor_radius(residue.name, atom.name)
                    xs.append(atom.pos.x)
                    ys.append(atom.pos.y)
                    zs.append(atom.pos.z)
                    rs.append(radius)

        # Skip empty structures (no protein atoms)
        if len(xs) == 0:
            return (cif_path.stem, 0, 0)

        # Remove duplicate coordinates (some PDB entries have exact duplicates)
        seen: set[tuple[float, float, float]] = set()
        unique_xs: list[float] = []
        unique_ys: list[float] = []
        unique_zs: list[float] = []
        unique_rs: list[float] = []
        for i in range(len(xs)):
            coord = (xs[i], ys[i], zs[i])
            if coord not in seen:
                seen.add(coord)
                unique_xs.append(xs[i])
                unique_ys.append(ys[i])
                unique_zs.append(zs[i])
                unique_rs.append(rs[i])
        xs, ys, zs, rs = unique_xs, unique_ys, unique_zs, unique_rs

        # Minimal JSON (no residue/atom_name - radii are pre-computed)
        data = {"x": xs, "y": ys, "z": zs, "r": rs}
        json_bytes = json.dumps(data, separators=(",", ":")).encode("utf-8")

        # Write gzipped
        output_path.parent.mkdir(parents=True, exist_ok=True)
        with gzip.open(output_path, "wb", compresslevel=6) as f:
            f.write(json_bytes)

        pdb_id = cif_path.stem.replace(".cif", "")
        return (pdb_id, len(xs), output_path.stat().st_size)

    except Exception as e:
        # Return error info
        return (f"{cif_path.stem} (error: {e})", -1, 0)


def process_single_file(args: tuple[Path, Path]) -> tuple[str, int, int]:
    """Worker function for parallel processing."""
    cif_path, output_path = args
    return cif_to_benchmark_json(cif_path, output_path)


def find_cif_files(input_dir: Path) -> list[Path]:
    """Find all CIF files using gemmi.CifWalk."""
    cif_files = []
    for path in gemmi.CifWalk(str(input_dir)):
        cif_files.append(Path(path))
    return cif_files


app = typer.Typer(help=__doc__)


@app.command()
def generate(
    input_path: Annotated[
        Path | None,
        typer.Argument(
            help="Input directory (recursive) or use --file for single file"
        ),
    ] = None,
    output_path: Annotated[
        Path | None,
        typer.Argument(help="Output directory (or .tar.gz file with --archive)"),
    ] = None,
    file_mode: Annotated[
        bool,
        typer.Option("--file", "-f", help="Single file mode"),
    ] = False,
    archive: Annotated[
        bool,
        typer.Option("--archive", "-a", help="Create tar archive for distribution"),
    ] = False,
    split: Annotated[
        int,
        typer.Option("--split", "-s", help="Split archive into chunks of N files"),
    ] = 0,
    workers: Annotated[
        int,
        typer.Option("--workers", "-w", help="Number of parallel workers"),
    ] = 0,
) -> None:
    """Convert CIF file(s) to JSON with ProtOr radii."""

    if input_path is None:
        console.print("[red]Error: Input path required[/red]")
        raise typer.Exit(1)

    if not input_path.exists():
        console.print(f"[red]Error: Path not found: {input_path}[/red]")
        raise typer.Exit(1)

    # Single file mode
    if file_mode:
        if archive:
            console.print("[red]Error: --archive not supported with --file[/red]")
            raise typer.Exit(1)

        if output_path is None:
            stem = input_path.stem.replace(".cif", "")
            output_path = input_path.parent / f"{stem}.json.gz"

        pdb_id, n_atoms, size = cif_to_benchmark_json(input_path, output_path)
        if n_atoms > 0:
            size_kb = size / 1024
            console.print(
                f"[green]Generated {output_path.name}[/green] "
                f"({n_atoms} atoms, {size_kb:.1f} KB)"
            )
        else:
            console.print(
                f"[yellow]Skipped {input_path.name}[/yellow] (no protein atoms)"
            )
        return

    # Directory mode
    if output_path is None:
        console.print("[red]Error: Output path required for directory mode[/red]")
        raise typer.Exit(1)

    # Archive mode: process to temp dir, then create tar
    if archive:
        process_directory_archive(input_path, output_path, workers, split)
    else:
        process_directory(input_path, output_path, workers)


def process_directory(input_dir: Path, output_dir: Path, workers: int) -> None:
    """Process all CIF files in directory recursively."""

    console.print(f"[bold]Scanning {input_dir}...[/bold]")
    cif_files = find_cif_files(input_dir)
    console.print(f"Found [cyan]{len(cif_files):,}[/cyan] CIF files\n")

    if len(cif_files) == 0:
        console.print("[yellow]No CIF files found[/yellow]")
        return

    output_dir.mkdir(parents=True, exist_ok=True)

    # Prepare work items
    work_items: list[tuple[Path, Path]] = []
    for cif_path in cif_files:
        # Extract PDB ID from path (handle various naming conventions)
        stem = cif_path.stem.replace(".cif", "").lower()
        output_file = output_dir / f"{stem}.json.gz"
        work_items.append((cif_path, output_file))

    # Determine worker count
    if workers <= 0:
        workers = max(1, mp.cpu_count() - 1)

    console.print(f"Processing with [cyan]{workers}[/cyan] workers...\n")

    # Process in parallel
    success_count = 0
    skip_count = 0
    error_count = 0
    total_atoms = 0
    total_size = 0

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
                executor.submit(process_single_file, item): item for item in work_items
            }

            for future in as_completed(futures):
                pdb_id, n_atoms, size = future.result()
                progress.advance(task)

                if n_atoms > 0:
                    success_count += 1
                    total_atoms += n_atoms
                    total_size += size
                elif n_atoms == 0:
                    skip_count += 1
                else:
                    error_count += 1

    # Summary
    console.print()
    table = Table(title="Summary")
    table.add_column("Metric", style="cyan")
    table.add_column("Value", justify="right")

    table.add_row("Processed", f"{success_count:,}")
    table.add_row("Skipped (no protein)", f"{skip_count:,}")
    table.add_row("Errors", f"{error_count:,}")
    table.add_row("Total atoms", f"{total_atoms:,}")
    table.add_row("Total size", f"{total_size / 1024 / 1024:.1f} MB")

    console.print(table)
    console.print(f"\n[green]Done![/green] Output saved to: {output_dir}")


def process_directory_archive(
    input_dir: Path, archive_path: Path, workers: int, split: int = 0
) -> None:
    """Process all CIF files and create a tar archive."""

    # Ensure archive path ends with .tar
    base_path = archive_path
    if str(archive_path).endswith(".tar"):
        base_path = Path(str(archive_path)[:-4])
    archive_path = Path(str(base_path) + ".tar")

    # Create temp directory for processing
    with tempfile.TemporaryDirectory() as temp_dir:
        temp_output = Path(temp_dir) / "dataset"
        temp_output.mkdir()

        # Process files to temp directory
        console.print(f"[bold]Scanning {input_dir}...[/bold]")
        cif_files = find_cif_files(input_dir)
        console.print(f"Found [cyan]{len(cif_files):,}[/cyan] CIF files\n")

        if len(cif_files) == 0:
            console.print("[yellow]No CIF files found[/yellow]")
            return

        # Prepare work items
        work_items: list[tuple[Path, Path]] = []
        for cif_path in cif_files:
            stem = cif_path.stem.replace(".cif", "").lower()
            output_file = temp_output / f"{stem}.json.gz"
            work_items.append((cif_path, output_file))

        # Determine worker count
        if workers <= 0:
            workers = max(1, mp.cpu_count() - 1)

        console.print(f"Processing with [cyan]{workers}[/cyan] workers...\n")

        # Process in parallel
        success_count = 0
        skip_count = 0
        error_count = 0
        total_atoms = 0
        total_size = 0

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
                    executor.submit(process_single_file, item): item
                    for item in work_items
                }

                for future in as_completed(futures):
                    pdb_id, n_atoms, size = future.result()
                    progress.advance(task)

                    if n_atoms > 0:
                        success_count += 1
                        total_atoms += n_atoms
                        total_size += size
                    elif n_atoms == 0:
                        skip_count += 1
                    else:
                        error_count += 1

        # Create tar archive(s)
        console.print("\n[bold]Creating archive...[/bold]")
        base_path.parent.mkdir(parents=True, exist_ok=True)

        json_files = sorted(temp_output.glob("*.json.gz"))
        archive_paths: list[Path] = []
        total_archive_size = 0

        if split > 0 and len(json_files) > split:
            # Split into multiple archives
            n_chunks = (len(json_files) + split - 1) // split
            with Progress(
                SpinnerColumn(),
                TextColumn("[progress.description]{task.description}"),
                BarColumn(),
                MofNCompleteColumn(),
                TimeElapsedColumn(),
                console=console,
            ) as progress:
                task = progress.add_task("Creating archives", total=n_chunks)

                for i in range(0, len(json_files), split):
                    chunk = json_files[i : i + split]
                    chunk_num = i // split + 1
                    chunk_path = Path(f"{base_path}_{chunk_num:03d}.tar")
                    archive_paths.append(chunk_path)

                    progress.update(task, description=f"Creating {chunk_path.name}")

                    with tarfile.open(chunk_path, "w") as tar:
                        for json_file in chunk:
                            tar.add(json_file, arcname=f"dataset/{json_file.name}")

                    total_archive_size += chunk_path.stat().st_size
                    progress.advance(task)
        else:
            # Single archive
            with tarfile.open(archive_path, "w") as tar:
                for json_file in json_files:
                    tar.add(json_file, arcname=f"dataset/{json_file.name}")
            archive_paths.append(archive_path)
            total_archive_size = archive_path.stat().st_size

        # Summary
        console.print()
        table = Table(title="Summary")
        table.add_column("Metric", style="cyan")
        table.add_column("Value", justify="right")

        table.add_row("Processed", f"{success_count:,}")
        table.add_row("Skipped (no protein)", f"{skip_count:,}")
        table.add_row("Errors", f"{error_count:,}")
        table.add_row("Total atoms", f"{total_atoms:,}")
        table.add_row("Uncompressed size", f"{total_size / 1024 / 1024:.1f} MB")
        table.add_row("Archive size", f"{total_archive_size / 1024 / 1024:.1f} MB")
        if len(archive_paths) > 1:
            table.add_row("Archives", f"{len(archive_paths)}")

        console.print(table)
        if len(archive_paths) == 1:
            console.print(
                f"\n[green]Done![/green] Archive saved to: {archive_paths[0]}"
            )
        else:
            console.print(f"\n[green]Done![/green] Archives saved to:")
            for p in archive_paths:
                console.print(f"  - {p}")


if __name__ == "__main__":
    app()
