#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.12"
# dependencies = [
#     "MDAnalysis>=2.7.0",
# ]
# ///
"""Convert an XTC trajectory to DCD format using MDAnalysis.

Usage:
    ./scripts/convert_xtc_to_dcd.py test_data/1l2y.xtc test_data/1l2y.pdb test_data/1l2y.dcd
"""

import sys

import MDAnalysis as mda


def main() -> None:
    if len(sys.argv) != 4:
        print(f"Usage: {sys.argv[0]} <xtc> <topology> <output.dcd>")
        sys.exit(1)

    xtc_path = sys.argv[1]
    top_path = sys.argv[2]
    dcd_path = sys.argv[3]

    print(f"Loading: {xtc_path} with topology {top_path}")
    u = mda.Universe(top_path, xtc_path)

    print(f"Atoms: {u.atoms.n_atoms}")
    print(f"Frames: {u.trajectory.n_frames}")

    # Write DCD
    with mda.Writer(dcd_path, u.atoms.n_atoms) as w:
        for ts in u.trajectory:
            w.write(u.atoms)

    print(f"Written: {dcd_path}")

    # Verify by reading back
    u2 = mda.Universe(top_path, dcd_path)
    print(f"Verification - Frames: {u2.trajectory.n_frames}, Atoms: {u2.atoms.n_atoms}")


if __name__ == "__main__":
    main()
