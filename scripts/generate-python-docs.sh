#!/usr/bin/env bash
# Generate Python API documentation using pdoc.
#
# Usage:
#   ./scripts/generate-python-docs.sh
#
# Output: docs/python-api/auto/
#
# Prerequisites:
#   cd python && uv pip install -e ".[dev]"

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
OUTPUT_DIR="$PROJECT_ROOT/docs/python-api/auto"

cd "$PROJECT_ROOT/python"

echo "Generating Python API docs with pdoc..."
uv run pdoc zsasa/ -o "$OUTPUT_DIR"

echo "Done. Output: $OUTPUT_DIR"
