# Build stage
FROM alpine:3.21 AS builder

# Install Zig
RUN apk add --no-cache curl xz git && \
    ARCH=$(uname -m) && \
    curl -fsSL "https://ziglang.org/download/0.16.0/zig-${ARCH}-linux-0.16.0.tar.xz" | \
    tar -xJ -C /usr/local && \
    ln -s /usr/local/zig-${ARCH}-linux-0.16.0/zig /usr/local/bin/zig

# Copy source and build
WORKDIR /src
COPY build.zig build.zig.zon ./
COPY src/ src/
RUN zig build -Doptimize=ReleaseFast && \
    cp zig-out/bin/zsasa /zsasa

# Runtime stage — statically linked, no OS needed
FROM scratch
COPY --from=builder /zsasa /usr/local/bin/zsasa
ENTRYPOINT ["/usr/local/bin/zsasa"]
