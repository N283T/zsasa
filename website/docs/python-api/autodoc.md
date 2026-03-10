---
sidebar_position: 6
---

# Python Autodoc (pdoc)

Auto-generated API documentation from Python docstrings using [pdoc](https://pdoc.dev/).

The autodoc provides searchable, browsable documentation for all public classes,
functions, and modules in the zsasa Python package.

**[Open Python Autodoc](pathname:///zsasa/python-autodoc/zsasa.html)**

## Features

- Complete API reference generated from source docstrings
- Type annotations and signatures
- Searchable across all modules
- Covers core API, integrations, and analysis utilities

## Generate Locally

```bash
cd python
uv run pdoc zsasa --output-dir /tmp/python-autodoc
# Open /tmp/python-autodoc/zsasa.html in browser
```

Or serve with live reload:

```bash
cd python
uv run pdoc zsasa
# Opens http://localhost:8080 automatically
```
