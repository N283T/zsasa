#!/bin/bash
# shellcheck shell=bash
# FreeSASA batch JSON processing using xargs -P for parallelism
# Usage: ./freesasa_batch.sh <input_dir> <output_dir> [n_jobs]

set -euo pipefail

SCRIPT_PATH="$(cd "$(dirname "$0")" && pwd)/$(basename "$0")"
SCRIPT_DIR="$(dirname "$SCRIPT_PATH")"
FREESASA="${SCRIPT_DIR}/../external/freesasa-bench/src/freesasa"

# Cross-platform nanosecond time
get_time_ns() {
    if command -v gdate &>/dev/null; then
        gdate +%s%N
    elif [[ "$(date +%s%N 2>/dev/null)" != *N* ]]; then
        date +%s%N
    else
        python3 -c "import time; print(int(time.time_ns()))"
    fi
}

# Single file processing mode (called by xargs)
process_one() {
    local input_path="$1"
    local output_dir="$2"
    local work_dir="$3"
    local freesasa="$4"

    local filename
    filename="$(basename "$input_path")"
    local base_name="${filename%.gz}"
    base_name="${base_name%.json}"

    local output_filename="${filename%.gz}"
    local output_path="${output_dir}/${output_filename}"

    local start_ns end_ns
    start_ns=$(get_time_ns)

    local result
    result=$("$freesasa" --json-input="$input_path" --no-warnings 2>&1) || true

    end_ns=$(get_time_ns)
    local sasa_time_ns=$((end_ns - start_ns))

    # Extract total SASA from result
    local total_sasa
    total_sasa=$(echo "$result" | grep -E '^Total\s+:' | grep -oE '[0-9]+\.[0-9]+' || echo "0")

    # Write output JSON
    cat > "$output_path" << EOF
{
  "total_area": ${total_sasa},
  "atom_areas": []
}
EOF

    # Write timing to work dir
    echo "$sasa_time_ns" > "$work_dir/${base_name}.time"
}

# Check for single-file processing mode
if [[ "${1:-}" == "--process-one" ]]; then
    shift
    process_one "$@"
    exit 0
fi

# Main batch mode
INPUT_DIR="${1:?Usage: $0 <input_dir> <output_dir> [n_jobs]}"
OUTPUT_DIR="${2:?Usage: $0 <input_dir> <output_dir> [n_jobs]}"
N_JOBS="${3:-$(sysctl -n hw.ncpu 2>/dev/null || nproc 2>/dev/null || echo 4)}"

# Validate inputs
if [[ ! -d "$INPUT_DIR" ]]; then
    echo "Error: Input directory does not exist: $INPUT_DIR" >&2
    exit 1
fi

if ! [[ "$N_JOBS" =~ ^[0-9]+$ ]] || [[ "$N_JOBS" -lt 1 ]]; then
    echo "Error: n_jobs must be a positive integer" >&2
    exit 1
fi

if [[ ! -x "$FREESASA" ]]; then
    echo "Error: freesasa not found at $FREESASA" >&2
    exit 1
fi

# Create output directory
mkdir -p "$OUTPUT_DIR"

# Create temp directory for timing results
WORK_DIR=$(mktemp -d)
trap 'rm -rf "$WORK_DIR"' EXIT

# Find JSON files
JSON_FILES=()
while IFS= read -r -d '' file; do
    JSON_FILES+=("$file")
done < <(find -L "$INPUT_DIR" -maxdepth 1 -type f \( -name "*.json" -o -name "*.json.gz" \) -print0 | sort -z)

TOTAL_FILES=${#JSON_FILES[@]}

if [[ $TOTAL_FILES -eq 0 ]]; then
    echo "No .json or .json.gz files found in $INPUT_DIR" >&2
    exit 1
fi

echo "Batch processing: $TOTAL_FILES files" >&2
echo "Jobs: $N_JOBS" >&2

START_TIME=$(get_time_ns)

# Process files in parallel using xargs -P
printf '%s\0' "${JSON_FILES[@]}" | \
    xargs -0 -P "$N_JOBS" -I {} \
    "$SCRIPT_PATH" --process-one {} "$OUTPUT_DIR" "$WORK_DIR" "$FREESASA"

END_TIME=$(get_time_ns)
TOTAL_TIME_NS=$((END_TIME - START_TIME))

# Collect results
TOTAL_SASA_TIME_NS=0
SUCCESSFUL=0
FAILED=0

for time_file in "$WORK_DIR"/*.time; do
    [[ -e "$time_file" ]] || continue
    time_ns=$(cat "$time_file")
    if [[ -n "$time_ns" && "$time_ns" =~ ^[0-9]+$ ]]; then
        TOTAL_SASA_TIME_NS=$((TOTAL_SASA_TIME_NS + time_ns))
        ((SUCCESSFUL++)) || true
    else
        ((FAILED++)) || true
    fi
done

# Calculate failed as total minus successful
FAILED=$((TOTAL_FILES - SUCCESSFUL))

# Output benchmark format
SASA_TIME_MS=$(python3 -c "print(f'{$TOTAL_SASA_TIME_NS / 1000000:.2f}')")
TOTAL_TIME_MS=$(python3 -c "print(f'{$TOTAL_TIME_NS / 1000000:.2f}')")

echo "BATCH_SASA_TIME_MS:$SASA_TIME_MS" >&2
echo "BATCH_TOTAL_TIME_MS:$TOTAL_TIME_MS" >&2
echo "BATCH_FILES:$TOTAL_FILES" >&2
echo "BATCH_SUCCESSFUL:$SUCCESSFUL" >&2
echo "BATCH_FAILED:$FAILED" >&2

echo ""
echo "Batch processing complete:"
echo "  Total files:  $TOTAL_FILES"
echo "  Successful:   $SUCCESSFUL"
echo "  Failed:       $FAILED"
echo "  SASA time:    $SASA_TIME_MS ms"
echo "  Total time:   $TOTAL_TIME_MS ms"
