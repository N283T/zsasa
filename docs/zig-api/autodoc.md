# Zig Autodoc (Interactive)

Auto-generated interactive API documentation from Zig source code.

The autodoc provides searchable, browsable documentation for all public types,
functions, and modules in zsasa.

**[Open Zig Autodoc](../zig-autodoc/){ .md-button .md-button--primary }**

## Features

- Search across public declarations and doc comments
- Interactive type browser with source links
- Doc comment rendering from `///` comments
- Built by the Zig compiler from `src/root.zig`

## Generate Locally

```bash
zig build docs
# Serve with HTTP server (WASM requires HTTP, not file://)
python3 -m http.server 8080 --directory zig-out/docs
# Open http://localhost:8080
```
