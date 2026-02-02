#!/usr/bin/env bash
# E. coli proteome benchmark runner
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "$SCRIPT_DIR/../.." && pwd)"

INPUT_DIR="$ROOT_DIR/benchmarks/UP000000625_83333_ECOLI_v6/json"
OUTPUT_BASE="$ROOT_DIR/benchmarks/results/ecoli"
RUN_SCRIPT="$ROOT_DIR/benchmarks/scripts/run.py"

THREADS="1,2,4,8,10"
RUNS=10

# Parse arguments
TOOLS="zig,freesasa,rust"
ALGORITHM="sr"

while [[ $# -gt 0 ]]; do
    case $1 in
        --tools)
            TOOLS="$2"
            shift 2
            ;;
        --algorithm)
            ALGORITHM="$2"
            shift 2
            ;;
        --threads)
            THREADS="$2"
            shift 2
            ;;
        --runs)
            RUNS="$2"
            shift 2
            ;;
        --dry-run)
            DRY_RUN=1
            shift
            ;;
        *)
            echo "Unknown option: $1"
            exit 1
            ;;
    esac
done

echo "E. coli Proteome Benchmark"
echo "=========================="
echo "Input: $INPUT_DIR"
echo "Output: $OUTPUT_BASE"
echo "Tools: $TOOLS"
echo "Algorithm: $ALGORITHM"
echo "Threads: $THREADS"
echo "Runs: $RUNS"
echo ""

mkdir -p "$OUTPUT_BASE"

IFS=',' read -ra TOOL_ARRAY <<< "$TOOLS"
for tool in "${TOOL_ARRAY[@]}"; do
    # Skip rust for LR (not supported)
    if [[ "$tool" == "rust" && "$ALGORITHM" == "lr" ]]; then
        echo "Skipping rust (LR not supported)"
        continue
    fi

    output_dir="$OUTPUT_BASE/${tool}_${ALGORITHM}"

    # Add precision suffix for zig
    if [[ "$tool" == "zig" ]]; then
        output_dir="${output_dir}_f64"
    fi

    echo "Running: $tool $ALGORITHM -> $output_dir"

    cmd="$RUN_SCRIPT \
        --tool $tool \
        --algorithm $ALGORITHM \
        --input-dir $INPUT_DIR \
        --output-dir $output_dir \
        --threads $THREADS \
        --runs $RUNS"

    if [[ -n "${DRY_RUN:-}" ]]; then
        echo "  [dry-run] $cmd"
    else
        $cmd
    fi

    echo ""
done

echo "Done! Results saved to: $OUTPUT_BASE"
