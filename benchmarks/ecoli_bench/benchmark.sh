#!/usr/bin/env bash
# E. coli Proteome Benchmark (RustSASA-style)
#
# All tools use PDB input for fair comparison.
# Methodology matches RustSASA paper: hyperfine, warmup 3, runs 3
#
# Usage:
#   ./benchmark.sh              # Run all
#   ./benchmark.sh --tool zig   # Run specific tool
#   ./benchmark.sh --dry-run    # Show commands only

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "$SCRIPT_DIR/../.." && pwd)"

# Input: PDB files
ECOLI_PDB="$ROOT_DIR/benchmarks/UP000000625_83333_ECOLI_v6/pdb"

# Output
RESULTS_DIR="$SCRIPT_DIR/results"
TEMP_OUT="$SCRIPT_DIR/temp_out"

# Binaries
ZSASA="$ROOT_DIR/zig-out/bin/zsasa"
FREESASA_BATCH="$SCRIPT_DIR/sasa_batch"
RUSTSASA="$ROOT_DIR/benchmarks/external/rustsasa-bench/target/release/rust-sasa"

# Parameters (RustSASA paper: warmup 3, runs 3)
WARMUP=3
RUNS=3
THREADS=8

# Parse args
TOOL=""
DRY_RUN=""

while [[ $# -gt 0 ]]; do
    case $1 in
        --tool|-t) TOOL="$2"; shift 2 ;;
        --threads|-T) THREADS="$2"; shift 2 ;;
        --warmup) WARMUP="$2"; shift 2 ;;
        --runs) RUNS="$2"; shift 2 ;;
        --dry-run) DRY_RUN=1; shift ;;
        *) echo "Unknown: $1"; exit 1 ;;
    esac
done

echo "=== E. coli Proteome Benchmark ==="
echo "Input: $ECOLI_PDB (PDB)"
echo "Warmup: $WARMUP, Runs: $RUNS, Threads: $THREADS"
echo ""

mkdir -p "$RESULTS_DIR" "$TEMP_OUT"

run_bench() {
    local name="$1"
    local cmd="$2"
    local json_out="$RESULTS_DIR/bench_${name}.json"

    echo ">>> $name"
    if [[ -n "$DRY_RUN" ]]; then
        echo "    $cmd"
        return
    fi

    # Note: Each tool clears its own output directory before running
    hyperfine --warmup "$WARMUP" --runs "$RUNS" --export-json "$json_out" "$cmd"
    echo ""
}

# === zsasa ===
run_zig() {
    [[ -x "$ZSASA" ]] || { echo "[SKIP] zsasa not found"; return; }

    local out_dir="$TEMP_OUT/zig"
    rm -rf "$out_dir" && mkdir -p "$out_dir"

    # f64 (default, double precision)
    run_bench "zsasa_f64_${THREADS}t" \
        "$ZSASA $ECOLI_PDB $out_dir --threads=$THREADS --parallelism=file --precision=f64"

    run_bench "zsasa_f64_1t" \
        "$ZSASA $ECOLI_PDB $out_dir --threads=1 --precision=f64"

    # f32 (single precision)
    run_bench "zsasa_f32_${THREADS}t" \
        "$ZSASA $ECOLI_PDB $out_dir --threads=$THREADS --parallelism=file --precision=f32"

    run_bench "zsasa_f32_1t" \
        "$ZSASA $ECOLI_PDB $out_dir --threads=1 --precision=f32"
}

# === FreeSASA (sasa_batch, single-threaded) ===
run_freesasa() {
    [[ -x "$FREESASA_BATCH" ]] || { echo "[SKIP] sasa_batch not found"; return; }

    local out_dir="$TEMP_OUT/freesasa"
    rm -rf "$out_dir" && mkdir -p "$out_dir"

    run_bench "freesasa_1t" \
        "$FREESASA_BATCH $ECOLI_PDB $out_dir"
}

# === RustSASA ===
run_rustsasa() {
    [[ -x "$RUSTSASA" ]] || { echo "[SKIP] RustSASA not found"; return; }

    local out_dir="$TEMP_OUT/rustsasa"
    rm -rf "$out_dir" && mkdir -p "$out_dir"

    run_bench "rustsasa_${THREADS}t" \
        "$RUSTSASA $ECOLI_PDB $out_dir --format json -t $THREADS"

    run_bench "rustsasa_1t" \
        "$RUSTSASA $ECOLI_PDB $out_dir --format json -t 1"
}

# === Main ===
case "${TOOL:-all}" in
    zig|zsasa) run_zig ;;
    freesasa) run_freesasa ;;
    rust|rustsasa) run_rustsasa ;;
    all) run_zig; run_freesasa; run_rustsasa ;;
    *) echo "Unknown tool: $TOOL"; exit 1 ;;
esac

echo "=== Done! Results: $RESULTS_DIR ==="

# Summary
if [[ -z "$DRY_RUN" ]]; then
    echo ""
    for f in "$RESULTS_DIR"/bench_*.json; do
        [[ -f "$f" ]] || continue
        name=$(basename "$f" .json | sed 's/bench_//')
        mean=$(jq -r '.results[0].mean' "$f" 2>/dev/null || echo "N/A")
        stddev=$(jq -r '.results[0].stddev' "$f" 2>/dev/null || echo "N/A")
        [[ "$mean" != "N/A" ]] && printf "  %-20s %8.3f s (± %.3f)\n" "$name" "$mean" "$stddev"
    done
fi
