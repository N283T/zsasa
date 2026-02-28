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

# Auto-enter nix-shell if not already inside one
if [ -z "${IN_NIX_SHELL:-}" ] && command -v nix-shell &>/dev/null; then
    exec nix-shell "$SCRIPT_DIR/shell.nix" --run "$0 $*"
fi

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
        autoreconf -i
        ./configure --enable-threads --disable-json --disable-xml
        make -j"$(nproc 2>/dev/null || sysctl -n hw.ncpu)"
    else
        skip "freesasa"
    fi
    symlink "$SCRIPT_DIR/freesasa/src/freesasa" freesasa
}

build_freesasa_batch() {
    info "freesasa_batch"
    if [ ! -f "$SCRIPT_DIR/freesasa/src/libfreesasa.a" ]; then
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
        git clone --recursive https://github.com/maxall41/RustSASA.git rustsasa
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
        git clone --recursive https://github.com/bisejdiu/lahuta.git
    fi
    cd lahuta
    if [ ! -f build/cli/lahuta ]; then
        cmake -B build -DCMAKE_BUILD_TYPE=Release -DLAHUTA_BUILD_PYTHON=OFF -DLAHUTA_BUILD_SHARED_CORE=OFF
        cmake --build build --config Release -j"$(nproc 2>/dev/null || sysctl -n hw.ncpu)"
    else
        skip "lahuta"
    fi
    symlink "$SCRIPT_DIR/lahuta/build/cli/lahuta" lahuta
}

build_zsasa() {
    info "zsasa (ReleaseFast)"
    (cd "$PROJECT_ROOT" && zig build --release=fast)
    local zsasa="$PROJECT_ROOT/zig-out/bin/zsasa"
    if [ -f "$zsasa" ]; then
        symlink "$zsasa" zsasa
    else
        err "zsasa build failed (binary not found at $zsasa)"
        return 1
    fi
}

# --- verify ---

verify_tools() {
    info "Verifying tools..."
    local testdata="$SCRIPT_DIR/testdata"
    local tmpdir
    tmpdir="$(mktemp -d)"
    trap 'rm -rf "$tmpdir"' RETURN
    local test_pdb
    test_pdb="$(ls "$testdata"/*.pdb | head -1)"
    local failed=0

    # zsasa
    if "$BIN_DIR/zsasa" calc "$test_pdb" > "$tmpdir/zsasa.out" 2>&1; then
        ok "zsasa"
    else
        err "zsasa"; failed=1
    fi

    # zsasa batch
    if "$BIN_DIR/zsasa" batch "$testdata" --format=jsonl -o "$tmpdir/zsasa_batch.jsonl" --threads=2 > /dev/null 2>&1; then
        ok "zsasa batch"
    else
        err "zsasa batch"; failed=1
    fi

    # freesasa
    if "$BIN_DIR/freesasa" "$test_pdb" > "$tmpdir/freesasa.out" 2>&1; then
        ok "freesasa"
    else
        err "freesasa"; failed=1
    fi

    # freesasa_batch
    if "$BIN_DIR/freesasa_batch" "$testdata" "$tmpdir/freesasa_batch" > /dev/null 2>&1; then
        ok "freesasa_batch"
    else
        err "freesasa_batch"; failed=1
    fi

    # rustsasa
    if "$BIN_DIR/rust-sasa" "$test_pdb" "$tmpdir/rustsasa_out" > /dev/null 2>&1; then
        ok "rust-sasa"
    else
        err "rust-sasa"; failed=1
    fi

    # lahuta
    if "$BIN_DIR/lahuta" sasa-sr -f "$test_pdb" -o "$tmpdir/lahuta.jsonl" --progress 0 --is_af2_model > /dev/null 2>&1; then
        ok "lahuta"
    else
        err "lahuta"; failed=1
    fi

    if [ "$failed" -eq 0 ]; then
        echo ""
        ok "All tools verified successfully"
    else
        echo ""
        err "Some tools failed verification"
        return 1
    fi
}

# --- main ---

if [ $# -eq 0 ]; then
    # Build all
    build_freesasa
    build_freesasa_batch
    build_rustsasa
    build_lahuta
    build_zsasa
    echo ""
    info "Done! Binaries in $BIN_DIR:"
    ls -la "$BIN_DIR"/
    echo ""
    verify_tools
else
    for tool in "$@"; do
        case "$tool" in
            freesasa)       build_freesasa ;;
            freesasa_batch) build_freesasa_batch ;;
            rustsasa)       build_rustsasa ;;
            lahuta)         build_lahuta ;;
            zsasa)          build_zsasa ;;
            verify)         verify_tools ;;
            *)              err "Unknown tool: $tool"; exit 1 ;;
        esac
    done
fi
