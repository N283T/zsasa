{
  description = "zsasa - fast SASA calculator written in Zig";

  inputs = {
    nixpkgs.url = "github:NixOS/nixpkgs/nixpkgs-unstable";
    flake-utils.url = "github:numtide/flake-utils";
    zig-overlay = {
      url = "github:mitchellh/zig-overlay";
      inputs.nixpkgs.follows = "nixpkgs";
    };
  };

  outputs =
    {
      self,
      nixpkgs,
      flake-utils,
      zig-overlay,
    }:
    flake-utils.lib.eachSystem
      [
        "x86_64-linux"
        "aarch64-linux"
        "x86_64-darwin"
        "aarch64-darwin"
      ]
      (
        system:
        let
          pkgs = nixpkgs.legacyPackages.${system};
          zig = zig-overlay.packages.${system}."0.15.2";

          # Pre-fetch Zig dependencies as a fixed-output derivation.
          # This runs `zig build --fetch` with network access and captures
          # the resulting package cache directory ($ZIG_GLOBAL_CACHE_DIR/p).
          zigDeps = pkgs.runCommand "zsasa-zig-deps"
            {
              src = ./.;
              nativeBuildInputs = [ zig ];
              outputHashAlgo = "sha256";
              outputHashMode = "recursive";
              outputHash = "sha256-wxgdxQiNj7hOTagippwrDkeiZ5RZsagvJI67T28jf04=";
            }
            ''
              export ZIG_GLOBAL_CACHE_DIR=$(mktemp -d)
              cp -r $src/. .
              zig build --fetch
              mv $ZIG_GLOBAL_CACHE_DIR/p $out
            '';

          zsasa = pkgs.stdenv.mkDerivation {
            pname = "zsasa";
            version = "0.2.3";

            src = ./.;

            nativeBuildInputs = [ zig ];

            dontConfigure = true;
            dontFixup = true;

            buildPhase = ''
              export ZIG_GLOBAL_CACHE_DIR=$(mktemp -d)
              zig build \
                -Doptimize=ReleaseFast \
                --system ${zigDeps} \
                --prefix $out \
                -j$NIX_BUILD_CORES
            '';

            installPhase = ''
              # zig build --prefix already installed the binary
              true
            '';

            meta = with pkgs.lib; {
              description = "Fast solvent accessible surface area (SASA) calculator written in Zig";
              homepage = "https://github.com/N283T/zsasa";
              license = licenses.mit;
              mainProgram = "zsasa";
              platforms = platforms.unix;
            };
          };
        in
        {
          packages = {
            default = zsasa;
            inherit zsasa;
          };

          apps.default = flake-utils.lib.mkApp {
            drv = zsasa;
          };
        }
      );
}
