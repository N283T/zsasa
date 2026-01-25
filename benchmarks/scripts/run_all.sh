#!/bin/bash
# Full benchmark suite for all tools and algorithms
# Runs with 10-minute intervals between each benchmark
#
# Usage:
#   ./run_all.sh                                    # All 238K structures
#   ./run_all.sh benchmarks/samples/stratified_75k.json  # With sample file

set -e

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
INPUT_DIR="benchmarks/inputs"
SAMPLE_FILE="${1:-}"  # Optional: first argument is sample file
THREADS="1,2,4,8,10"
RUNS=1
INTERVAL=600  # 10 minutes in seconds

log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1"
}

run_benchmark() {
    local tool=$1
    local algorithm=$2

    log "Starting: $tool $algorithm (threads=$THREADS, runs=$RUNS)"

    local cmd=("$SCRIPT_DIR/run.py"
        --tool "$tool"
        --algorithm "$algorithm"
        --input-dir "$INPUT_DIR"
        --threads "$THREADS"
        --runs "$RUNS"
    )

    # Add sample file if specified
    if [[ -n "$SAMPLE_FILE" ]]; then
        cmd+=(--sample-file "$SAMPLE_FILE")
    fi

    "${cmd[@]}"

    log "Completed: $tool $algorithm"
}

wait_interval() {
    log "Waiting ${INTERVAL}s (10 min) before next benchmark..."
    sleep "$INTERVAL"
}

# Main
log "=== Full Benchmark Suite ==="
log "Input: $INPUT_DIR"
if [[ -n "$SAMPLE_FILE" ]]; then
    log "Sample: $SAMPLE_FILE"
fi
log "Threads: $THREADS, Runs: $RUNS"
log ""

# SR algorithm
run_benchmark freesasa sr
wait_interval

run_benchmark rust sr
wait_interval

run_benchmark zig sr
wait_interval

# LR algorithm
run_benchmark freesasa lr
wait_interval

run_benchmark zig lr

log ""
log "=== All benchmarks completed ==="
