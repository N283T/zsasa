# Build stage
FROM alpine:3.21 AS builder

# Install Zig. The tarball is verified against the sha256 that ziglang.org
# publishes in https://ziglang.org/download/index.json for this exact version
# (the build fails if it does not match). When bumping ZIG_VERSION, update
# both checksums, python/pyproject.toml (cibuildwheel before-all) and the Zig
# version in .github/workflows/*.yml.
ARG ZIG_VERSION=0.16.0
ARG ZIG_SHA256_X86_64=70e49664a74374b48b51e6f3fdfbf437f6395d42509050588bd49abe52ba3d00
ARG ZIG_SHA256_AARCH64=ea4b09bfb22ec6f6c6ceac57ab63efb6b46e17ab08d21f69f3a48b38e1534f17
RUN apk add --no-cache curl xz git && \
    ARCH=$(uname -m) && \
    case "${ARCH}" in \
      x86_64)  ZIG_SHA256="${ZIG_SHA256_X86_64}" ;; \
      aarch64) ZIG_SHA256="${ZIG_SHA256_AARCH64}" ;; \
      *) echo "unsupported architecture: ${ARCH}" >&2; exit 1 ;; \
    esac && \
    curl -fsSL -o /tmp/zig.tar.xz \
      "https://ziglang.org/download/${ZIG_VERSION}/zig-${ARCH}-linux-${ZIG_VERSION}.tar.xz" && \
    { [ "$(sha256sum /tmp/zig.tar.xz | cut -d' ' -f1)" = "${ZIG_SHA256}" ] || \
      { echo "sha256 mismatch for the Zig tarball" >&2; exit 1; }; } && \
    tar -xJf /tmp/zig.tar.xz -C /usr/local && \
    rm /tmp/zig.tar.xz && \
    ln -s /usr/local/zig-${ARCH}-linux-${ZIG_VERSION}/zig /usr/local/bin/zig

# Copy source and build
WORKDIR /src
COPY build.zig build.zig.zon ./
COPY src/ src/
# Build for a fixed CPU instead of the build machine's: x86_64_v3 (AVX2, FMA)
# on x86-64 and the baseline on arm64, like the release binaries.
RUN case "$(uname -m)" in \
      x86_64) ZIG_CPU=x86_64_v3 ;; \
      *)      ZIG_CPU=baseline ;; \
    esac && \
    zig build -Doptimize=ReleaseFast -Dcpu="${ZIG_CPU}" && \
    cp zig-out/bin/zsasa /zsasa

# Runtime stage — statically linked, no OS needed
FROM scratch
COPY --from=builder /zsasa /usr/local/bin/zsasa
ENTRYPOINT ["/usr/local/bin/zsasa"]
