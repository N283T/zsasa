# CLI Progress Bars Design

## Goal

Replace the CLI's hand-written carriage-return progress messages for `batch` and `traj` with Zig's standard `std.Progress` reporting while preserving current quiet-mode and output behavior.

## Scope

- Applies to `zsasa batch` and `zsasa traj` only.
- Keeps `zsasa calc` unchanged.
- Keeps `-q, --quiet` as the single user-facing way to suppress progress.
- Does not add new CLI flags.

## Architecture

Use `std.Progress.start(...)` to create a progress root in each subcommand path that currently prints periodic progress. Create child nodes for file/item processing in `batch` and frame processing in `traj`. Existing status summaries, timing output, result files, CSV, JSON, and JSONL output stay unchanged.

For parallel batch processing, workers continue to update the existing atomic `processed_count`; the monitor loop owns the progress node and advances it by the observed delta. This avoids cross-thread progress node mutation. For sequential batch and trajectory paths, the processing loop advances the node after each completed item or frame.

## Behavior

- Non-quiet `batch` shows a standard progress bar/count for processed work items.
- Non-quiet `traj` shows standard progress for processed frames when the total is known from `--end`; otherwise it shows an indeterminate/processed-count style progress node.
- Quiet mode suppresses progress bars and keeps script-friendly output behavior.
- Existing final summaries remain visible in non-quiet mode.

## Testing

- Unit tests should cover argument parsing remains unchanged.
- Build/test with Zig 0.16.0 through the repository's Nix/flake tooling.
- Smoke test `batch` on `examples/` and `traj` on `test_data/1l2y.xtc` + `test_data/1l2y.pdb` with a short frame range.
