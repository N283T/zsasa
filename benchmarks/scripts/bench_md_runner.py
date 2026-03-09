#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.12"
# dependencies = [
#     "cffi",
#     "MDAnalysis",
#     "mdsasa-bolt",
#     "mdtraj",
#     "numpy",
#     "typer>=0.9.0",
# ]
# ///
"""Runner script for MD trajectory SASA benchmarks.

Called by bench_md.py via hyperfine to time individual tool executions.
Each invocation loads trajectory data, computes SASA, then exits.

Usage:
    ./bench_md_runner.py --tool mdtraj --xtc traj.xtc --pdb top.pdb
    ./bench_md_runner.py --tool zsasa_mdtraj --xtc traj.xtc --pdb top.pdb --threads 4
"""

from __future__ import annotations

import sys
from enum import Enum
from pathlib import Path
from typing import Annotated

import typer

# Add zsasa Python package to path
sys.path.insert(0, str(Path(__file__).parent.parent.parent.joinpath("python")))

app = typer.Typer()


class Tool(str, Enum):
    mdtraj = "mdtraj"
    zsasa_mdtraj = "zsasa_mdtraj"
    zsasa_mdanalysis = "zsasa_mdanalysis"
    mdsasa_bolt = "mdsasa_bolt"


def run_mdtraj(xtc: Path, pdb: Path, n_points: int, stride: int) -> None:
    """Run native MDTraj shrake_rupley (single-threaded)."""
    import mdtraj as md

    traj = md.load(str(xtc), top=str(pdb), stride=stride)
    md.shrake_rupley(traj, n_sphere_points=n_points)


def run_zsasa_mdtraj(
    xtc: Path,
    pdb: Path,
    n_threads: int,
    n_points: int,
    stride: int,
    use_bitmask: bool = False,
) -> None:
    """Run zsasa via MDTraj wrapper."""
    import mdtraj as md

    from zsasa.mdtraj import compute_sasa

    traj = md.load(str(xtc), top=str(pdb), stride=stride)
    compute_sasa(
        traj,
        n_points=n_points,
        n_threads=n_threads,
        mode="total",
        use_bitmask=use_bitmask,
    )


def run_zsasa_mdanalysis(
    xtc: Path,
    pdb: Path,
    n_threads: int,
    n_points: int,
    stride: int,
    use_bitmask: bool = False,
) -> None:
    """Run zsasa via MDAnalysis wrapper."""
    import MDAnalysis as mda

    from zsasa.mdanalysis import SASAAnalysis

    u = mda.Universe(str(pdb), str(xtc))
    analysis = SASAAnalysis(u)
    analysis.run(
        step=stride, n_points=n_points, n_threads=n_threads, use_bitmask=use_bitmask
    )


def run_mdsasa_bolt(xtc: Path, pdb: Path, n_points: int, stride: int) -> None:
    """Run mdsasa-bolt (RustSASA, all cores)."""
    try:
        import MDAnalysis as mda
        from mdsasa_bolt import SASAAnalysis
    except ImportError as e:
        print(f"Error: mdsasa-bolt not installed: {e}", file=sys.stderr)
        sys.exit(1)

    u = mda.Universe(str(pdb), str(xtc))
    sasa = SASAAnalysis(u.atoms, n_points=n_points)
    sasa.run(step=stride)


@app.command()
def main(
    tool: Annotated[Tool, typer.Option("--tool", help="Tool to benchmark")],
    xtc: Annotated[Path, typer.Option("--xtc", help="XTC trajectory file")],
    pdb: Annotated[Path, typer.Option("--pdb", help="Topology PDB file")],
    threads: Annotated[int, typer.Option("--threads", help="Number of threads")] = 0,
    n_points: Annotated[
        int, typer.Option("--n-points", help="Test points per atom")
    ] = 100,
    stride: Annotated[int, typer.Option("--stride", help="Frame stride")] = 1,
    use_bitmask: Annotated[
        bool, typer.Option("--use-bitmask", help="Use bitmask LUT optimization")
    ] = False,
) -> None:
    """Execute a single SASA benchmark tool."""
    if tool == Tool.mdtraj:
        run_mdtraj(xtc, pdb, n_points, stride)
    elif tool == Tool.zsasa_mdtraj:
        run_zsasa_mdtraj(xtc, pdb, threads, n_points, stride, use_bitmask=use_bitmask)
    elif tool == Tool.zsasa_mdanalysis:
        run_zsasa_mdanalysis(
            xtc, pdb, threads, n_points, stride, use_bitmask=use_bitmask
        )
    elif tool == Tool.mdsasa_bolt:
        run_mdsasa_bolt(xtc, pdb, n_points, stride)


if __name__ == "__main__":
    app()
