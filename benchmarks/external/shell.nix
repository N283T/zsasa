{ pkgs ? import <nixpkgs> {} }:
pkgs.mkShell {
  packages = with pkgs; [
    git
    # FreeSASA (C, autotools)
    autoconf
    automake
    libtool
    pkg-config
    # RustSASA (Rust)
    cargo
    rustc
    # Lahuta (C++, cmake)
    cmake
    zlib
  ];
}
