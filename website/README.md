# zsasa website

The site at <https://n283t.github.io/zsasa/> is built from this directory by a small Python script. There is no Node toolchain.

| Path | Contents |
| --- | --- |
| `docs/` | Documentation pages, in Markdown |
| `src/index.html`, `src/404.html` | Hand-written landing and error pages |
| `src/assets/` | Stylesheet and scripts (`site.css`, `site.js`, `charts.js`, hero figures) |
| `templates/` | Shared page chrome and the docs page template |
| `static/` | Files copied to the site root as they are |
| `data/benchmarks/` | Chart and table data for the benchmark pages (generated) |
| `build.py` | The builder |

## Build and preview

```bash
uv run website/build.py
python3 -m http.server 4321 --directory website/dist
```

The build fails on broken internal links, missing heading anchors, and unknown chart or table references. The sidebar order is the `NAV` list at the top of `build.py`; a new page must be added there.

## Markdown extensions

The docs use a few constructs beyond CommonMark, all handled by `build.py`:

- Admonitions: `:::note`, `:::tip`, `:::info`, `:::warning`, `:::danger`, with an optional `[Title]`.
- Tabs: `<Tabs>` with `<TabItem label="...">` children.
- Explicit heading ids: `## Heading {#custom-id}`.
- Benchmark figures: `<div data-chart="batch/map"></div>` and `<div data-table="batch/t10-ecoli-afdb"></div>`, where the first path segment is a file in `data/benchmarks/` and the second a key inside it.

## Updating the benchmark data

The JSON files under `data/benchmarks/` are exported from a [`zsasa-benchmarks`](https://github.com/N283T/zsasa-benchmarks) checkout with populated `results/`. Run the exporter with that repository's environment, then rebuild:

```bash
uv run --project ../zsasa-benchmarks python website/scripts/export_benchmarks.py ../zsasa-benchmarks
uv run website/build.py
```

Do not edit the JSON by hand. The headline numbers on the landing page come from the same files.
