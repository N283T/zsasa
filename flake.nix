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
          zig = zig-overlay.packages.${system}."0.16.0";

          # Pre-fetch Zig dependencies as a fixed-output derivation.
          # This runs `zig build --fetch` with network access and captures the
          # unpacked packages. Since Zig 0.16 they are unpacked into ./zig-pkg
          # (the global cache only keeps the tarballs under $ZIG_GLOBAL_CACHE_DIR/p),
          # and that directory is what `zig build --system` expects.
          #
          # outputHash must be refreshed whenever the dependencies in
          # build.zig.zon change (or the Zig version does): run
          # `scripts/check_nix_deps_hash.py --refresh`, see AGENTS.md.
          zigDeps = pkgs.runCommand "zsasa-zig-deps"
            {
              src = ./.;
              nativeBuildInputs = [ zig ];
              outputHashAlgo = "sha256";
              outputHashMode = "recursive";
              # zig-deps-fingerprint: 10e9e106d4e4a19a0e54b3161037408cdc0edfa42214966585d7cc65da2d075e
              outputHash = "sha256-JxhOkVdcM9wtvh5CGOxz/PU5RhnHoq7RwwlLXXeCHdY=";
            }
            ''
              export ZIG_GLOBAL_CACHE_DIR=$(mktemp -d)
              cp -r $src/. .
              zig build --fetch
              mv zig-pkg $out
            '';

          zsasa = pkgs.stdenv.mkDerivation {
            pname = "zsasa";
            version = "0.12.0";

            src = ./.;

            nativeBuildInputs = [ zig ];

            dontConfigure = true;
            dontFixup = true;

            # Without -Dcpu, Zig targets the CPU of the machine that builds the
            # package, so a binary from a cache or another machine could fault
            # on an older CPU. Match the release binaries: x86_64_v3 (AVX2,
            # FMA) on x86-64 and the baseline elsewhere.
            zigCpu = if pkgs.stdenv.hostPlatform.isx86_64 then "x86_64_v3" else "baseline";

            buildPhase = ''
              export ZIG_GLOBAL_CACHE_DIR=$(mktemp -d)
              zig build \
                -Doptimize=ReleaseFast \
                -Dcpu=$zigCpu \
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
