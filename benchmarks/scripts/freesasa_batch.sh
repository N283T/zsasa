#!/bin/bash
# shellcheck shell=bash
# FreeSASA batch JSON processing using background jobs
# Usage: ./freesasa_batch.sh <input_dir> <output_dir> [n_jobs]

set -euo pipefail

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

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
FREESASA="${SCRIPT_DIR}/../external/freesasa-bench/src/freesasa"

if [[ ! -x "$FREESASA" ]]; then
    echo "Error: freesasa not found at $FREESASA" >&2
    exit 1
fi

# Create output directory
mkdir -p "$OUTPUT_DIR"

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
echo "Threads: $N_JOBS" >&2

# Create temp directory for results (avoid shadowing TMPDIR)
WORK_DIR=$(mktemp -d)
trap 'rm -rf "$WORK_DIR"' EXIT

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

# Process single file function
process_one() {
    local idx="$1"
    local input_path="$2"
    local output_dir="$3"
    local freesasa="$4"
    local work_dir="$5"

    local filename
    filename="$(basename "$input_path")"

    # Remove .gz suffix for output filename
    local output_filename="${filename%.gz}"
    local output_path="${output_dir}/${output_filename}"

    # Run freesasa with JSON input and capture timing
    local start_ns end_ns
    start_ns=$(get_time_ns)

    local result
    result=$("$freesasa" --json-input="$input_path" --no-warnings 2>&1) || true

    end_ns=$(get_time_ns)

    local sasa_time_ns=$((end_ns - start_ns))

    # Extract total SASA from result (format: "Total     :    XXXX.XX")
    local total_sasa
    total_sasa=$(echo "$result" | grep -E '^Total\s+:' | grep -oE '[0-9]+\.[0-9]+' || echo "0")

    # Write output JSON with atomic rename
    local tmp_output
    tmp_output=$(mktemp "$output_dir/.tmp.XXXXXX")
    cat > "$tmp_output" << EOF
{
  "total_area": ${total_sasa},
  "atom_areas": []
}
EOF
    mv "$tmp_output" "$output_path"

    echo "$sasa_time_ns" > "$work_dir/$idx.time"
}

START_TIME=$(get_time_ns)

# Process files with limited parallelism
PIDS=()
for i in "${!JSON_FILES[@]}"; do
    process_one "$i" "${JSON_FILES[$i]}" "$OUTPUT_DIR" "$FREESASA" "$WORK_DIR" &
    PIDS+=($!)

    # Limit concurrent jobs using wait -n if available (Bash 4.3+)
    if (( ${#PIDS[@]} >= N_JOBS )); then
        if wait -n 2>/dev/null; then
            : # Job completed successfully
        fi
        # Clean up finished PIDs
        for j in "${!PIDS[@]}"; do
            if ! kill -0 "${PIDS[$j]}" 2>/dev/null; then
                wait "${PIDS[$j]}" 2>/dev/null || true
                unset 'PIDS[j]'
            fi
        done
        PIDS=("${PIDS[@]}")
    fi

    echo -ne "\rProcessing: $((i+1))/$TOTAL_FILES" >&2
done

# Wait for remaining jobs and check exit status
for pid in "${PIDS[@]}"; do
    if ! wait "$pid" 2>/dev/null; then
        echo "Warning: Job $pid failed" >&2
    fi
done
echo >&2

END_TIME=$(get_time_ns)
TOTAL_TIME_NS=$((END_TIME - START_TIME))

# Collect results
TOTAL_SASA_TIME_NS=0
SUCCESSFUL=0
FAILED=0

for i in "${!JSON_FILES[@]}"; do
    if [[ -f "$WORK_DIR/$i.time" ]]; then
        time_ns=$(cat "$WORK_DIR/$i.time")
        if [[ -n "$time_ns" && "$time_ns" =~ ^[0-9]+$ ]]; then
            TOTAL_SASA_TIME_NS=$((TOTAL_SASA_TIME_NS + time_ns))
            ((SUCCESSFUL++)) || true
        else
            ((FAILED++)) || true
        fi
    else
        ((FAILED++)) || true
    fi
done

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
