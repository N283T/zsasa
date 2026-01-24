#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.11"
# dependencies = ["rich>=13.0", "typer>=0.9.0"]
# ///
"""Validate Zig SASA implementation against FreeSASA reference values.

Compares the Zig implementation output with FreeSASA reference data
and reports accuracy metrics.
"""

from __future__ import annotations

import json
import re
import subprocess
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Annotated

import typer
from rich.console import Console
from rich.table import Table

app = typer.Typer(help="Validate Zig SASA against FreeSASA reference")
console = Console()


@dataclass
class ValidationResult:
    """Result of validating one structure."""

    pdb_id: str
    n_atoms: int
    reference_sasa: float
    zig_sasa: float
    difference_percent: float
    passed: bool
    time_ms: float | None = None


def run_zig_sasa(
    input_path: Path,
    algorithm: str = "sr",
    n_points: int = 100,
    n_slices: int = 20,
    classifier: str | None = None,
) -> tuple[float, float | None]:
    """Run Zig SASA implementation and return (total_area, time_ms)."""
    zig_binary = Path(__file__).parent.parent / "zig-out" / "bin" / "freesasa_zig"

    if not zig_binary.exists():
        raise FileNotFoundError(f"Zig binary not found: {zig_binary}")

    with tempfile.NamedTemporaryFile(suffix=".json", delete=False) as f:
        output_path = Path(f.name)

    try:
        cmd = [
            str(zig_binary),
            f"--algorithm={algorithm}",
        ]
        if algorithm == "sr":
            cmd.append(f"--n-points={n_points}")
        else:
            cmd.append(f"--n-slices={n_slices}")

        if classifier:
            cmd.append(f"--classifier={classifier}")

        cmd.extend([str(input_path), str(output_path)])

        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=300,
        )

        if result.returncode != 0:
            raise RuntimeError(f"Zig failed: {result.stderr}")

        # Parse time from stderr (e.g., "Total area: 18908.90 Å²")
        time_ms = None
        for line in result.stderr.split("\n"):
            if "ms" in line.lower():
                # Try to extract time
                match = re.search(r"(\d+(?:\.\d+)?)\s*ms", line, re.IGNORECASE)
                if match:
                    time_ms = float(match.group(1))

        with open(output_path) as f:
            data = json.load(f)

        return data["total_area"], time_ms

    finally:
        output_path.unlink(missing_ok=True)


def validate_structure(
    pdb_id: str,
    inputs_dir: Path,
    references_dir: Path,
    algorithm: str,
    tolerance: float,
    classifier: str | None = None,
) -> ValidationResult | None:
    """Validate one structure against reference."""
    input_path = inputs_dir / f"{pdb_id}.json"
    reference_path = references_dir / f"{pdb_id}_n100_p1.4.json"

    if not input_path.exists() or not reference_path.exists():
        return None

    with open(reference_path) as f:
        ref = json.load(f)

    reference_sasa = ref["total_area"]
    n_atoms = ref["n_atoms"]

    try:
        zig_sasa, time_ms = run_zig_sasa(input_path, algorithm, classifier=classifier)
    except Exception as e:
        console.print(f"  [red]{pdb_id}: ERROR[/red] - {e}")
        return None

    if reference_sasa == 0:
        diff_percent = 0.0 if zig_sasa == 0 else float("inf")
    else:
        diff_percent = abs(zig_sasa - reference_sasa) / reference_sasa * 100
    passed = diff_percent <= tolerance

    return ValidationResult(
        pdb_id=pdb_id,
        n_atoms=n_atoms,
        reference_sasa=reference_sasa,
        zig_sasa=zig_sasa,
        difference_percent=diff_percent,
        passed=passed,
        time_ms=time_ms,
    )


@app.command()
def main(
    algorithm: Annotated[
        str, typer.Option("--algorithm", "-a", help="Algorithm: sr or lr")
    ] = "sr",
    tolerance: Annotated[
        float, typer.Option("--tolerance", "-t", help="Max allowed difference (%)")
    ] = 2.0,
    classifier: Annotated[
        str | None,
        typer.Option("--classifier", "-c", help="Classifier (protor, naccess)"),
    ] = "protor",
    no_classifier: Annotated[
        bool, typer.Option("--no-classifier", help="Disable classifier (element-based)")
    ] = False,
) -> None:
    """Validate Zig SASA against FreeSASA reference values."""
    if no_classifier:
        classifier = None

    base_dir = Path(__file__).parent.parent / "benchmarks"
    inputs_dir = base_dir / "inputs"
    references_dir = base_dir / "references"

    structures = ["1crn", "1ubq", "1a0q", "3hhb", "1aon", "4v6x"]

    clf_str = classifier if classifier else "none (element-based)"
    console.rule("[bold]SASA Validation[/bold]")
    console.print(
        f"algorithm={algorithm}, classifier={clf_str}, tolerance={tolerance}%"
    )

    table = Table(show_header=True, header_style="bold cyan")
    table.add_column("PDB", style="bold")
    table.add_column("Atoms", justify="right")
    table.add_column("FreeSASA", justify="right", style="yellow")
    table.add_column("Zig", justify="right", style="green")
    table.add_column("Diff%", justify="right")
    table.add_column("Status", justify="center")

    results: list[ValidationResult] = []
    for pdb_id in structures:
        result = validate_structure(
            pdb_id, inputs_dir, references_dir, algorithm, tolerance, classifier
        )
        if result:
            results.append(result)
            status = "[green]PASS[/green]" if result.passed else "[red]FAIL[/red]"
            table.add_row(
                result.pdb_id.upper(),
                f"{result.n_atoms:,}",
                f"{result.reference_sasa:.2f}",
                f"{result.zig_sasa:.2f}",
                f"{result.difference_percent:.3f}%",
                status,
            )
        else:
            table.add_row(pdb_id.upper(), "-", "-", "-", "-", "[dim]SKIPPED[/dim]")

    console.print()
    console.print(table)

    if results:
        passed = sum(1 for r in results if r.passed)
        total = len(results)
        avg_diff = sum(r.difference_percent for r in results) / len(results)

        console.print()
        console.print(f"Summary: {passed}/{total} passed (avg diff: {avg_diff:.3f}%)")

        if passed == total:
            console.print("[bold green]All validations PASSED![/bold green]")
        else:
            console.print("[bold red]Some validations FAILED![/bold red]")
            raise typer.Exit(1)
    else:
        raise typer.Exit(1)


if __name__ == "__main__":
    app()
