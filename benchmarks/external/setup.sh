#!/usr/bin/env bash
# Setup external benchmark dependencies.
#
# Clones, builds, and symlinks all tool binaries into bin/.
#
# Usage:
#   cd benchmarks/external
#   ./setup.sh          # build all
#   ./setup.sh freesasa # build one tool
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
BIN_DIR="$SCRIPT_DIR/bin"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"

mkdir -p "$BIN_DIR"

# --- helpers ---

info()  { printf '\033[1;34m==>\033[0m %s\n' "$*"; }
ok()    { printf '\033[1;32m  ✓\033[0m %s\n' "$*"; }
skip()  { printf '\033[1;33m  ⊘\033[0m %s (already built)\n' "$*"; }
err()   { printf '\033[1;31m  ✗\033[0m %s\n' "$*" >&2; }

symlink() {
    local src="$1" name="$2"
    ln -sf "$src" "$BIN_DIR/$name"
    ok "$name -> $src"
}

# --- tools ---

build_freesasa() {
    info "FreeSASA (vanilla)"
    cd "$SCRIPT_DIR"
    if [ ! -d freesasa ]; then
        git clone https://github.com/mittinatten/freesasa.git
    fi
    cd freesasa
    if [ ! -f src/freesasa ]; then
        ./configure --enable-threads --disable-json --disable-xml
        make -j"$(nproc 2>/dev/null || sysctl -n hw.ncpu)"
    else
        skip "freesasa"
    fi
    symlink "$SCRIPT_DIR/freesasa/src/freesasa" freesasa
}

build_freesasa_batch() {
    info "freesasa_batch"
    if [ ! -f "$SCRIPT_DIR/freesasa/src/.libs/libfreesasa.a" ]; then
        err "Build freesasa first (need libfreesasa.a)"
        return 1
    fi
    cd "$SCRIPT_DIR/freesasa_batch"
    if [ ! -f freesasa_batch ]; then
        make
    else
        skip "freesasa_batch"
    fi
    symlink "$SCRIPT_DIR/freesasa_batch/freesasa_batch" freesasa_batch
}

build_rustsasa() {
    info "RustSASA (vanilla)"
    cd "$SCRIPT_DIR"
    if [ ! -d rustsasa ]; then
        git clone --recursive https://github.com/mcisb/rustsasa.git
    fi
    cd rustsasa
    if [ ! -f target/release/rust-sasa ]; then
        cargo build --release --features cli
    else
        skip "rust-sasa"
    fi
    symlink "$SCRIPT_DIR/rustsasa/target/release/rust-sasa" rust-sasa
}

build_lahuta() {
    info "Lahuta"
    cd "$SCRIPT_DIR"
    if [ ! -d lahuta ]; then
        git clone https://github.com/DominikSko/lahuta.git
    fi
    cd lahuta
    if [ ! -f build/cli/lahuta ]; then
        cmake -B build -DCMAKE_BUILD_TYPE=Release
        cmake --build build --config Release -j"$(nproc 2>/dev/null || sysctl -n hw.ncpu)"
    else
        skip "lahuta"
    fi
    symlink "$SCRIPT_DIR/lahuta/build/cli/lahuta" lahuta
}

link_zsasa() {
    info "zsasa"
    local zsasa="$PROJECT_ROOT/zig-out/bin/zsasa"
    if [ -f "$zsasa" ]; then
        symlink "$zsasa" zsasa
    else
        err "zsasa not found at $zsasa (run 'zig build' first)"
    fi
}

# --- main ---

if [ $# -eq 0 ]; then
    # Build all
    build_freesasa
    build_freesasa_batch
    build_rustsasa
    build_lahuta
    link_zsasa
    echo ""
    info "Done! Binaries in $BIN_DIR:"
    ls -la "$BIN_DIR"/
else
    for tool in "$@"; do
        case "$tool" in
            freesasa)       build_freesasa ;;
            freesasa_batch) build_freesasa_batch ;;
            rustsasa)       build_rustsasa ;;
            lahuta)         build_lahuta ;;
            zsasa)          link_zsasa ;;
            *)              err "Unknown tool: $tool"; exit 1 ;;
        esac
    done
fi
