#!/bin/sh
# install.sh - zsasa CLI installer
#
# Usage:
#   curl -fsSL https://raw.githubusercontent.com/N283T/zsasa/main/install.sh | sh
#   ./install.sh
#
# Environment variables:
#   VERSION     - Pin a specific version (default: latest release)
#   INSTALL_DIR - Override install directory (default: ~/.local/bin)

set -eu

REPO="N283T/zsasa"
REQUIRED_ZIG_VERSION="0.15.2"
DEFAULT_INSTALL_DIR="${HOME}/.local/bin"

# ---------------------------------------------------------------------------
# Color helpers (no-op when stdout is not a terminal)
# ---------------------------------------------------------------------------
_supports_color() {
    [ -t 1 ] && command -v tput > /dev/null 2>&1 && tput colors > /dev/null 2>&1
}

if _supports_color; then
    COLOR_RESET="$(tput sgr0)"
    COLOR_BOLD="$(tput bold)"
    COLOR_GREEN="$(tput setaf 2)"
    COLOR_RED="$(tput setaf 1)"
    COLOR_CYAN="$(tput setaf 6)"
else
    COLOR_RESET=""
    COLOR_BOLD=""
    COLOR_GREEN=""
    COLOR_RED=""
    COLOR_CYAN=""
fi

info()    { printf "%s==>%s %s%s\n"  "${COLOR_CYAN}"  "${COLOR_RESET}" "$*" "${COLOR_RESET}"; }
success() { printf "%s✓%s %s%s\n"   "${COLOR_GREEN}" "${COLOR_RESET}" "$*" "${COLOR_RESET}"; }
warn()    { printf "%sWarning:%s %s%s\n" "${COLOR_BOLD}" "${COLOR_RESET}" "$*" "${COLOR_RESET}" >&2; }
error()   { printf "%sError:%s %s%s\n" "${COLOR_RED}" "${COLOR_RESET}" "$*" "${COLOR_RESET}" >&2; }
bold()    { printf "%s%s%s\n" "${COLOR_BOLD}" "$*" "${COLOR_RESET}"; }

die() {
    error "$*"
    exit 1
}

# ---------------------------------------------------------------------------
# Temp dir + cleanup trap
# ---------------------------------------------------------------------------
TMPDIR_WORK=""

cleanup() {
    if [ -n "${TMPDIR_WORK}" ] && [ -d "${TMPDIR_WORK}" ]; then
        rm -rf "${TMPDIR_WORK}"
    fi
}

trap cleanup EXIT INT TERM

make_tmpdir() {
    TMPDIR_WORK="$(mktemp -d)"
}

# ---------------------------------------------------------------------------
# Download helper: curl preferred, wget fallback
# ---------------------------------------------------------------------------
download() {
    _url="$1"
    _dest="$2"
    if command -v curl > /dev/null 2>&1; then
        curl -fsSL --retry 3 --retry-delay 2 -o "${_dest}" "${_url}"
    elif command -v wget > /dev/null 2>&1; then
        wget -q --tries=3 -O "${_dest}" "${_url}"
    else
        die "Neither curl nor wget found. Please install one and retry."
    fi
}

download_stdout() {
    _url="$1"
    if command -v curl > /dev/null 2>&1; then
        curl -fsSL --retry 3 --retry-delay 2 "${_url}"
    elif command -v wget > /dev/null 2>&1; then
        wget -q --tries=3 -O - "${_url}"
    else
        die "Neither curl nor wget found. Please install one and retry."
    fi
}

# ---------------------------------------------------------------------------
# OS / Arch detection
# ---------------------------------------------------------------------------
detect_target() {
    _os="$(uname -s)"
    _arch="$(uname -m)"

    case "${_os}" in
        Linux)  _os_name="linux" ;;
        Darwin) _os_name="macos" ;;
        *)      die "Unsupported operating system: ${_os}. Only Linux and macOS are supported." ;;
    esac

    case "${_arch}" in
        x86_64 | amd64)  _arch_name="x86_64" ;;
        aarch64 | arm64) _arch_name="aarch64" ;;
        *)                die "Unsupported architecture: ${_arch}." ;;
    esac

    echo "${_os_name}-${_arch_name}"
}

# ---------------------------------------------------------------------------
# Zig version check
# ---------------------------------------------------------------------------
zig_version() {
    # Returns the version string from `zig version`, or empty if zig not found
    if command -v zig > /dev/null 2>&1; then
        zig version 2>/dev/null || true
    fi
}

has_required_zig() {
    _ver="$(zig_version)"
    [ "${_ver}" = "${REQUIRED_ZIG_VERSION}" ]
}

# ---------------------------------------------------------------------------
# SHA256 helper: sha256sum (Linux) or shasum -a 256 (macOS)
# ---------------------------------------------------------------------------
_sha256() {
    if command -v sha256sum > /dev/null 2>&1; then
        sha256sum "$1" | awk '{print $1}'
    elif command -v shasum > /dev/null 2>&1; then
        shasum -a 256 "$1" | awk '{print $1}'
    else
        return 1
    fi
}

# ---------------------------------------------------------------------------
# Fetch latest version from GitHub API
# ---------------------------------------------------------------------------
fetch_latest_version() {
    _api_url="https://api.github.com/repos/${REPO}/releases/latest"
    _json="$(download_stdout "${_api_url}")" \
        || die "Could not fetch release information from GitHub."

    if command -v jq > /dev/null 2>&1; then
        _tag="$(printf '%s' "${_json}" | jq -r '.tag_name // empty' | sed 's/^v//')"
    else
        # Extract tag_name without jq -- use sed with POSIX BRE
        _tag="$(printf '%s' "${_json}" | sed -n 's/.*"tag_name"[[:space:]]*:[[:space:]]*"v\([^"]*\)".*/\1/p' | head -n 1)"
    fi

    if [ -z "${_tag}" ]; then
        die "Could not determine latest version. Try: VERSION=x.y.z $0"
    fi
    printf '%s' "${_tag}"
}

# ---------------------------------------------------------------------------
# Build from source
# ---------------------------------------------------------------------------
build_from_source() {
    _version="$1"
    _install_dir="$2"

    info "Cloning zsasa v${_version}..."
    _src_dir="${TMPDIR_WORK}/zsasa-src"
    git clone --depth 1 --branch "v${_version}" \
        "https://github.com/${REPO}.git" "${_src_dir}" \
        || die "git clone failed. Ensure git is installed and you have network access."

    info "Building with zig build -Doptimize=ReleaseFast..."
    ( cd "${_src_dir}" && zig build -Doptimize=ReleaseFast ) \
        || die "Build failed. Check that Zig ${REQUIRED_ZIG_VERSION} is correctly installed."

    _binary="${_src_dir}/zig-out/bin/zsasa"
    [ -f "${_binary}" ] || die "Built binary not found at ${_binary}."

    info "Installing zsasa to ${_install_dir}/zsasa..."
    mkdir -p "${_install_dir}"
    install -m 755 "${_binary}" "${_install_dir}/zsasa"
}

# ---------------------------------------------------------------------------
# Download pre-built binary
# ---------------------------------------------------------------------------
download_binary() {
    _version="$1"
    _target="$2"
    _install_dir="$3"

    _asset_name="zsasa-${_version}-${_target}"
    _url="https://github.com/${REPO}/releases/download/v${_version}/${_asset_name}"
    _tmp_bin="${TMPDIR_WORK}/zsasa"
    _checksum_url="https://github.com/${REPO}/releases/download/v${_version}/SHA256SUMS"
    _tmp_checksums="${TMPDIR_WORK}/SHA256SUMS"

    info "Downloading ${_asset_name}..."
    download "${_url}" "${_tmp_bin}" \
        || die "Download failed. Check your network or visit https://github.com/${REPO}/releases."

    info "Verifying checksum..."
    if download "${_checksum_url}" "${_tmp_checksums}" 2>/dev/null; then
        _expected="$(grep "${_asset_name}" "${_tmp_checksums}" | awk '{print $1}')"
        if [ -n "${_expected}" ]; then
            _actual="$(_sha256 "${_tmp_bin}")"
            if [ "${_expected}" != "${_actual}" ]; then
                die "Checksum mismatch! Expected ${_expected}, got ${_actual}"
            fi
            info "Checksum verified."
        fi
    else
        warn "Checksum file not available. Skipping verification."
    fi

    info "Installing zsasa to ${_install_dir}/zsasa..."
    mkdir -p "${_install_dir}"
    install -m 755 "${_tmp_bin}" "${_install_dir}/zsasa"
}

# ---------------------------------------------------------------------------
# PATH guidance
# ---------------------------------------------------------------------------
print_path_guidance() {
    _install_dir="$1"

    # Check if _install_dir is already in PATH
    _in_path=0
    _saved_IFS="${IFS}"
    IFS=":"
    for _dir in ${PATH}; do
        if [ "${_dir}" = "${_install_dir}" ]; then
            _in_path=1
            break
        fi
    done
    IFS="${_saved_IFS}"

    if [ "${_in_path}" -eq 0 ]; then
        printf "\n"
        info "Add to your PATH (if not already):"
        printf "  export PATH=\"%s:\$PATH\"\n" "${_install_dir}"
        printf "Add this to your %s~/.zshrc%s or %s~/.bashrc%s to make it permanent.\n" \
            "${COLOR_BOLD}" "${COLOR_RESET}" "${COLOR_BOLD}" "${COLOR_RESET}"
    fi
}

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
main() {
    bold ""
    bold "zsasa installer"
    bold ""

    make_tmpdir

    # Resolve install directory
    INSTALL_DIR="${INSTALL_DIR:-${DEFAULT_INSTALL_DIR}}"

    # Validate INSTALL_DIR is an absolute path
    case "${INSTALL_DIR}" in
        /*) ;;
        *) die "INSTALL_DIR must be an absolute path, got: '${INSTALL_DIR}'" ;;
    esac

    # Detect platform
    TARGET="$(detect_target)"
    info "Detected: ${TARGET}"
    info "Install directory: ${INSTALL_DIR}"
    printf "\n"

    # Resolve version
    if [ -z "${VERSION:-}" ]; then
        info "Fetching latest version from GitHub..."
        VERSION="$(fetch_latest_version)"
    fi

    # Validate version format
    case "${VERSION}" in
        [0-9]*.[0-9]*.[0-9]*) ;;
        *) die "Invalid version format: '${VERSION}'. Expected: X.Y.Z" ;;
    esac

    info "Version: ${VERSION}"
    printf "\n"

    # Decide build strategy
    if has_required_zig; then
        command -v git > /dev/null 2>&1 \
            || die "git is required for building from source. Install git or remove Zig from PATH to use binary download."
        info "Zig ${REQUIRED_ZIG_VERSION} found — building from source..."
        build_from_source "${VERSION}" "${INSTALL_DIR}"
    else
        _zig_ver="$(zig_version)"
        if [ -n "${_zig_ver}" ]; then
            info "Zig ${_zig_ver} found (need ${REQUIRED_ZIG_VERSION}) — downloading pre-built binary..."
        else
            info "Zig not found — downloading pre-built binary..."
        fi
        download_binary "${VERSION}" "${TARGET}" "${INSTALL_DIR}"
    fi

    printf "\n"

    # Verify installation
    _installed_bin="${INSTALL_DIR}/zsasa"
    if ! "${_installed_bin}" --version > /dev/null 2>&1; then
        die "Installation verification failed. Binary at ${_installed_bin} did not run successfully."
    fi

    success "zsasa ${VERSION} installed successfully!"

    print_path_guidance "${INSTALL_DIR}"

    printf "\n"
}

main "$@"
