"""Export benchmark charts and tables for the website from a zsasa-benchmarks checkout.

Reads the summary tables under ``results/tables`` and the validation database of
the benchmark repository and writes compact chart specifications to
``website/data/benchmarks/*.json``. The site builder inlines those specs into the
benchmark pages; ``website/src/assets/charts.js`` draws them.

Run it with the benchmark repository's environment, because the validation data
is loaded through that repository's own loader so the website uses exactly the
runs adopted for the published figures:

    uv run --project ../zsasa-benchmarks python website/scripts/export_benchmarks.py ../zsasa-benchmarks
"""

from __future__ import annotations

import argparse
import csv
import json
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

OUT = Path(__file__).resolve().parents[1] / "data" / "benchmarks"

ZSASA_VERSION = "0.9.0"

BATCH_DATASETS = {
    "E. coli AFDB": "E. coli AFDB",
    "Human AFDB": "Human AFDB",
    "UP000005640_9606_HUMAN_v6_cif": "Human AFDB, mmCIF",
    "SwissProt AFDB": "SwissProt AFDB",
}
# zsasa mmCIF parser-path variants are bitmask f32 runs; their names do not say so.
BITMASK_VARIANTS = {
    "zsasa_generic_read",
    "zsasa_generic_mmap",
    "zsasa_af_fast_read",
    "zsasa_af_fast_mmap",
}
TRAJECTORIES = ["5wvo_C", "6sup_A", "5vz0_A"]
POINTS = [64, 128, 256, 512, 1024]


def tool_of(variant: str) -> str:
    """Colour slot: zsasa, the reference implementation, the Rust implementation, or other."""
    if variant.startswith("zsasa"):
        return "zsasa"
    if variant.startswith(("freesasa", "mdtraj")):
        return "ref"
    if variant.startswith(("rustsasa", "mdsasa")):
        return "rust"
    return "other"


def mode_of(variant: str) -> str:
    return "bitmask" if "bitmask" in variant or variant in BITMASK_VARIANTS else "exact"


def sig(value: float, digits: int = 4) -> float:
    """Round to significant digits to keep the exported JSON small."""
    if value == 0 or not np.isfinite(value):
        return 0.0
    return float(f"{value:.{digits}g}")


def num(row: dict[str, str], key: str) -> float | None:
    text = row.get(key, "")
    return float(text) if text not in ("", None) else None


def fmt(value: float | None, digits: int = 1, suffix: str = "") -> str:
    if value is None:
        return "–"
    return f"{value:,.{digits}f}{suffix}"


class Tables:
    def __init__(self, bench: Path) -> None:
        self.dir = bench / "results" / "tables"

    def rows(self, name: str) -> list[dict[str, str]]:
        with (self.dir / f"{name}.csv").open() as handle:
            return list(csv.DictReader(handle))


def marker(variant: str, name: str) -> dict[str, str | bool]:
    """Marker shape, fill and line dash that tell apart series sharing a tool colour."""
    bitmask = mode_of(variant) == "bitmask"
    if variant in BITMASK_VARIANTS:  # mmCIF parser paths
        shape, hollow, dash = (
            ("triangle" if "af_fast" in variant else "square"),
            "mmap" in variant,
            "",
        )
    elif variant.startswith("zsasa_md"):  # Python integrations
        shape, hollow, dash = (
            ("square" if "mdtraj" in variant else "triangle"),
            bitmask,
            "dash" if bitmask else "",
        )
    else:
        shape, hollow, dash = (
            ("diamond" if bitmask else "circle"),
            "f32" in name,
            "dash" if bitmask else "",
        )
    if "corrected" in name:
        dash = "dot"
    return {"shape": shape, "hollow": hollow, "dash": dash}


def series_meta(
    row: dict[str, str], keep_version: bool = False
) -> dict[str, str | bool]:
    name = row["display_name"]
    if not keep_version:
        name = name.replace(f"zsasa {ZSASA_VERSION} ", "zsasa ")
    return {
        "name": name,
        "tool": tool_of(row["variant"]),
        "mode": mode_of(row["variant"]),
    } | marker(row["variant"], name)


def throughput_err(row: dict[str, str], key: str) -> float:
    """Propagate the runtime standard deviation onto a rate (items / runtime)."""
    return (
        num(row, key)
        * (num(row, "runtime_stddev_s") or 0.0)
        / num(row, "runtime_mean_s")
    )


# ---------------------------------------------------------------- batch


def export_batch(t: Tables) -> dict:
    t10 = defaultdict(list)
    for row in t.rows("batch_t10_summary"):
        t10[row["dataset_label"]].append(row)
    scaling = defaultdict(lambda: defaultdict(list))
    for row in t.rows("batch_thread_scaling"):
        scaling[row["dataset_label"]][row["variant"]].append(row)

    thr = "throughput_structures_per_sec"
    charts, tables = {}, {}

    charts["map"] = {
        "type": "scatter",
        "title": "Throughput against peak memory",
        "note": "10 threads, 128 sphere points. Points are three-run means; whiskers show one standard deviation. "
        "Higher and further left is better.",
        "x": {"label": "Peak RSS (MiB)", "digits": 1, "unit": "MiB"},
        "y": {"label": "Throughput (structures/s)", "digits": 0, "unit": "str/s"},
        "controls": [{"label": "Dataset", "options": ["E. coli AFDB", "Human AFDB"]}],
        "views": {
            str(i): {
                "points": [
                    series_meta(r)
                    | {
                        "x": sig(num(r, "peak_rss_mean_mib")),
                        "y": sig(num(r, thr)),
                        "ex": sig((num(r, "peak_rss_stddev_mib") or 0.0), 3),
                        "ey": sig(throughput_err(r, thr), 3),
                    }
                    for r in t10[label]
                ]
            }
            for i, label in enumerate(["E. coli AFDB", "Human AFDB"])
        },
    }

    metrics = [
        ("Throughput", thr, "Throughput (structures/s)", "str/s", 0),
        ("Peak RSS", "peak_rss_mean_mib", "Peak RSS (MiB)", "MiB", 1),
    ]
    views = {}
    for i, key in enumerate(BATCH_DATASETS):
        for j, (_, col, label, unit, digits) in enumerate(metrics):
            rows = []
            for r in t10[key]:
                err = (
                    throughput_err(r, thr)
                    if col == thr
                    else (num(r, "peak_rss_stddev_mib") or 0.0)
                )
                entry = series_meta(r, keep_version=key == "SwissProt AFDB") | {
                    "v": sig(num(r, col)),
                    "e": sig(err, 3),
                }
                speedup = num(r, "runtime_speedup_vs_freesasa")
                if speedup is not None and r["variant"] != "freesasa_batch":
                    entry["note"] = f"{speedup:.2f}× the FreeSASA batch throughput"
                rows.append(entry)
            views[f"{i}|{j}"] = {
                "label": label,
                "unit": unit,
                "digits": digits,
                "rows": rows,
            }
    charts["bars"] = {
        "type": "bar",
        "title": "Batch throughput and memory at 10 threads",
        "note": "128 sphere points. Bars are three-run means; whiskers show one standard deviation.",
        "controls": [
            {"label": "Dataset", "options": list(BATCH_DATASETS.values())},
            {"label": "Metric", "options": [m[0] for m in metrics]},
        ],
        "views": views,
    }

    story = [
        "zsasa_0_9_0_f64",
        "zsasa_0_9_0_bitmask_f32",
        "freesasa_batch",
        "rustsasa",
        "lahuta_bitmask",
    ]
    charts["scaling"] = {
        "type": "line",
        "title": "E. coli AFDB throughput by thread count",
        "note": "4,370 structures, 128 sphere points. Points are three-run means with one standard deviation.",
        "x": {"label": "Threads", "values": [1, 4, 8, 10]},
        "views": {
            "0": {
                "y": {
                    "label": "Throughput (structures/s)",
                    "digits": 0,
                    "unit": "str/s",
                },
                "series": [
                    series_meta(scaling["E. coli AFDB"][v][0])
                    | {
                        "pts": [
                            [
                                int(r["threads"]),
                                sig(num(r, thr)),
                                sig(throughput_err(r, thr), 3),
                            ]
                            for r in sorted(
                                scaling["E. coli AFDB"][v],
                                key=lambda r: int(r["threads"]),
                            )
                            if int(r["threads"]) <= 10
                        ]
                    }
                    for v in story
                ],
            }
        },
    }

    # Busy cores is (user + system CPU time) / wall time: how many cores the run kept working.
    overcommit_metrics = metrics + [
        (
            "Busy cores",
            "cpu_utilization_proxy",
            "Busy cores (CPU time / wall time)",
            "cores",
            1,
        )
    ]
    views = {}
    for i, key in enumerate(BATCH_DATASETS):
        for j, (_, col, label, unit, digits) in enumerate(overcommit_metrics):
            series = []
            for variant, rows in scaling[key].items():
                rows = sorted(
                    (r for r in rows if int(r["threads"]) >= 10),
                    key=lambda r: int(r["threads"]),
                )
                if len(rows) < 3 or tool_of(variant) != "zsasa":
                    continue
                errors = {
                    thr: lambda r: throughput_err(r, thr),
                    "peak_rss_mean_mib": lambda r: num(r, "peak_rss_stddev_mib") or 0.0,
                }
                err = errors.get(col, lambda r: 0.0)
                series.append(
                    series_meta(rows[0], keep_version=key == "SwissProt AFDB")
                    | {
                        "pts": [
                            [int(r["threads"]), sig(num(r, col)), sig(err(r), 3)]
                            for r in rows
                        ]
                    }
                )
            views[f"{i}|{j}"] = {
                "y": {"label": label, "digits": digits, "unit": unit},
                "series": series,
            }
    charts["overcommit"] = {
        "type": "line",
        "title": "zsasa beyond the core count",
        "note": "Worker counts above the 10 logical CPUs of the benchmark machine, at 128 sphere points. Three-run means, "
        "except SwissProt, which was measured once. Busy cores is CPU time divided by wall time.",
        "x": {"label": "Worker threads", "values": [10, 20, 40]},
        "controls": [
            {"label": "Dataset", "options": list(BATCH_DATASETS.values())},
            {"label": "Metric", "options": [m[0] for m in overcommit_metrics]},
        ],
        "views": views,
    }

    for key, label in BATCH_DATASETS.items():
        rows = t10[key]
        comparators = [
            (c, n)
            for c, n in (
                ("freesasa", "FreeSASA batch"),
                ("rustsasa", "RustSASA"),
                ("lahuta_bitmask", "Lahuta bitmask"),
            )
            if any(num(r, f"runtime_speedup_vs_{c}") is not None for r in rows)
        ]
        slug = label.lower().replace(". ", "").replace(", ", "-").replace(" ", "-")
        tables[f"t10-{slug}"] = {
            "caption": f"{label}: {int(rows[0]['expected_count']):,} structures, 10 threads, 128 sphere points.",
            "columns": ["Tool", "Runtime (s)", "Structures/s", "Peak RSS (MiB)"]
            + [f"vs {n}" for _, n in comparators],
            "rows": [
                [
                    series_meta(r, keep_version=key == "SwissProt AFDB")["name"],
                    fmt(num(r, "runtime_mean_s"), 3),
                    fmt(num(r, thr), 0),
                    fmt(num(r, "peak_rss_mean_mib")),
                ]
                + [
                    fmt(num(r, f"runtime_speedup_vs_{c}"), 2, "×")
                    for c, _ in comparators
                ]
                for r in rows
            ],
        }

    def headline(key: str) -> dict[str, str]:
        r = next(r for r in t10[key] if r["variant"] == "zsasa_0_9_0_bitmask_f32")
        return {
            "structures": f"{int(r['expected_count']):,}",
            "runtime": f"{num(r, 'runtime_mean_s'):.1f}",
            "throughput": f"{num(r, thr):,.0f}",
            "rss": f"{num(r, 'peak_rss_mean_mib'):.1f}",
            "rss_ceil": str(int(-(-num(r, "peak_rss_mean_mib") // 5) * 5)),
            "speedup": f"{num(r, 'runtime_speedup_vs_freesasa'):.2f}",
        }

    # Headline numbers quoted on the landing page: zsasa bitmask f32 at 10 threads.
    facts = {
        f"{p}_{k}": v
        for p, key in (("ecoli", "E. coli AFDB"), ("human", "Human AFDB"))
        for k, v in headline(key).items()
    }
    return {"facts": facts, "charts": charts, "tables": tables}


# ---------------------------------------------------------------- trajectories


def export_md(t: Tables) -> dict:
    summary = defaultdict(list)
    for row in t.rows("md_summary"):
        summary[row["dataset_label"]].append(row)
    scaling = defaultdict(lambda: defaultdict(list))
    for row in t.rows("md_thread_scaling"):
        scaling[row["dataset_label"]][(row["variant"], row["bitmask_variant"])].append(
            row
        )

    def view_label(name: str) -> str:
        r = summary[name][0]
        return (
            f"{name} ({int(r['frame_count']):,} frames, {int(r['atom_count']):,} atoms)"
        )

    labels = [view_label(n) for n in TRAJECTORIES]
    fps = "frames_per_sec"
    charts, tables = {}, {}

    charts["map"] = {
        "type": "scatter",
        "title": "Trajectory throughput against peak memory",
        "note": "10 threads, 128 sphere points. Points are three-run means; whiskers show one standard deviation. "
        "Both axes are logarithmic.",
        "x": {"label": "Peak RSS (MiB)", "digits": 0, "unit": "MiB", "scale": "log"},
        "y": {
            "label": "Throughput (frames/s)",
            "digits": 1,
            "unit": "frames/s",
            "scale": "log",
        },
        "controls": [{"label": "Trajectory", "options": labels}],
        "views": {
            str(i): {
                "points": [
                    series_meta(r)
                    | {
                        "x": sig(num(r, "peak_rss_mean_mib")),
                        "y": sig(num(r, fps)),
                        "ex": sig((num(r, "peak_rss_stddev_mib") or 0.0), 3),
                        "ey": sig(throughput_err(r, fps), 3),
                    }
                    for r in summary[name]
                ]
            }
            for i, name in enumerate(TRAJECTORIES)
        },
    }

    metrics = [
        ("Throughput", fps, "Throughput (frames/s)", "frames/s", 1),
        ("Peak RSS", "peak_rss_mean_mib", "Peak RSS (MiB)", "MiB", 0),
    ]
    views = {}
    for i, name in enumerate(TRAJECTORIES):
        for j, (_, col, label, unit, digits) in enumerate(metrics):
            rows = []
            for r in summary[name]:
                err = (
                    throughput_err(r, fps)
                    if col == fps
                    else (num(r, "peak_rss_stddev_mib") or 0.0)
                )
                entry = series_meta(r) | {"v": sig(num(r, col)), "e": sig(err, 3)}
                notes = [
                    f"{num(r, f'runtime_speedup_vs_{c}'):.1f}× {n}"
                    for c, n in (("mdtraj", "MDTraj"), ("mdsasa_bolt", "mdsasa-bolt"))
                    if num(r, f"runtime_speedup_vs_{c}") is not None
                    and not r["variant"].startswith(c)
                ]
                if notes and tool_of(r["variant"]) == "zsasa":
                    entry["note"] = "Throughput: " + ", ".join(notes)
                rows.append(entry)
            views[f"{i}|{j}"] = {
                "label": label,
                "unit": unit,
                "digits": digits,
                "rows": rows,
            }
    charts["bars"] = {
        "type": "bar",
        "title": "Trajectory throughput and memory at 10 threads",
        "note": "128 sphere points, stride 1. Bars are three-run means; whiskers show one standard deviation.",
        "controls": [
            {"label": "Trajectory", "options": labels},
            {"label": "Metric", "options": [m[0] for m in metrics]},
        ],
        "views": views,
    }

    native = [
        ("zsasa_cli_f64", ""),
        ("zsasa_cli_f32", ""),
        ("zsasa_cli_bitmask_f64", "single_corrected"),
        ("zsasa_cli_bitmask_f32", "single_corrected"),
    ]
    views = {}
    for i, name in enumerate(TRAJECTORIES):
        series = []
        for key in native:
            rows = sorted(scaling[name][key], key=lambda r: int(r["threads"]))
            if rows:
                series.append(
                    series_meta(rows[0])
                    | {
                        "pts": [
                            [
                                int(r["threads"]),
                                sig(num(r, fps)),
                                sig(throughput_err(r, fps), 3),
                            ]
                            for r in rows
                        ]
                    }
                )
        views[str(i)] = {
            "y": {"label": "Throughput (frames/s)", "digits": 1, "unit": "frames/s"},
            "series": series,
        }
    charts["overcommit"] = {
        "type": "line",
        "title": "Native zsasa beyond the core count",
        "note": "Worker counts above the 10 logical CPUs of the benchmark machine. 128 sphere points, three-run means. "
        "Bitmask runs use the single-LUT mode with the experimental bias correction.",
        "x": {"label": "Worker threads", "values": [10, 20, 40]},
        "controls": [{"label": "Trajectory", "options": labels}],
        "views": views,
    }

    for name in TRAJECTORIES:
        rows = summary[name]
        comparators = [
            (c, n)
            for c, n in (("mdtraj", "MDTraj"), ("mdsasa_bolt", "mdsasa-bolt"))
            if any(num(r, f"runtime_speedup_vs_{c}") is not None for r in rows)
        ]
        tables[f"summary-{name.lower().replace('_', '-')}"] = {
            "caption": f"{view_label(name)}: 128 sphere points; zsasa at 10 threads.",
            "columns": ["Tool", "Runtime (s)", "Frames/s", "Peak RSS (MiB)"]
            + [f"vs {n}" for _, n in comparators],
            "rows": [
                [
                    r["display_name"],
                    fmt(num(r, "runtime_mean_s"), 2),
                    fmt(num(r, fps)),
                    fmt(num(r, "peak_rss_mean_mib"), 0),
                ]
                + [
                    fmt(num(r, f"runtime_speedup_vs_{c}"), 1, "×")
                    for c, _ in comparators
                ]
                for r in rows
            ],
        }
    return {"charts": charts, "tables": tables}


# ---------------------------------------------------------------- single file


def export_single(t: Tables) -> dict:
    t10 = defaultdict(lambda: defaultdict(list))
    for row in t.rows("single_file_t10_summary"):
        t10[row["input_format"]][row["structure_id"]].append(row)
    scaling = defaultdict(lambda: defaultdict(list))
    for row in t.rows("single_file_thread_scaling"):
        scaling[(row["input_format"], row["structure_id"])][row["variant"]].append(row)

    formats = ["PDB", "mmCIF"]
    atoms = {s: int(rows[0]["n_atoms"]) for s, rows in t10["PDB"].items()}
    structures = sorted(atoms, key=atoms.get)
    structure_labels = [f"{s} ({atoms[s]:,} atoms)" for s in structures]
    story = ["zsasa_f64", "zsasa_bitmask_f32", "freesasa", "rustsasa", "pdbtools_jl"]
    metrics = [
        ("Runtime", "runtime_mean_s", "runtime_stddev_s", "Runtime (s)", "s", 3),
        (
            "Peak RSS",
            "peak_rss_mean_mib",
            "peak_rss_stddev_mib",
            "Peak RSS (MiB)",
            "MiB",
            1,
        ),
    ]
    charts, tables = {}, {}

    views = {}
    for i, fmt_ in enumerate(formats):
        for j, (_, col, err, label, unit, digits) in enumerate(metrics):
            series = []
            for variant in story:
                pts = []
                for s in structures:
                    row = next(
                        (r for r in t10[fmt_][s] if r["variant"] == variant), None
                    )
                    if row and num(row, col) is not None:
                        pts.append(
                            [
                                atoms[s],
                                sig(num(row, col)),
                                sig(num(row, err) or 0, 3),
                                s,
                            ]
                        )
                if pts:
                    series.append(series_meta(row) | {"pts": pts})
            views[f"{i}|{j}"] = {
                "y": {"label": label, "digits": digits, "unit": unit, "scale": "log"},
                "series": series,
            }
    charts["atoms"] = {
        "type": "line",
        "title": "Single-structure cost against size",
        "note": "Eight structures from 10,919 to 4,506,416 atoms. 10 threads, 100 sphere points, three-run means. "
        "Both axes are logarithmic.",
        "x": {"label": "Atoms", "scale": "log", "digits": 0},
        "controls": [
            {"label": "Format", "options": formats},
            {"label": "Metric", "options": [m[0] for m in metrics]},
        ],
        "views": views,
    }

    views = {}
    for k, s in enumerate(structures):
        for i, fmt_ in enumerate(formats):
            for j, (_, col, err, label, unit, digits) in enumerate(metrics):
                rows = [
                    series_meta(r)
                    | {"v": sig(num(r, col)), "e": sig(num(r, err) or 0, 3)}
                    for r in t10[fmt_][s]
                    if num(r, col) is not None
                ]
                views[f"{k}|{i}|{j}"] = {
                    "label": label,
                    "unit": unit,
                    "digits": digits,
                    "rows": rows,
                }
    charts["bars"] = {
        "type": "bar",
        "title": "Per-structure runtime and memory at 10 threads",
        "note": "100 sphere points. Bars are three-run means; whiskers show one standard deviation.",
        "controls": [
            {"label": "Structure", "options": structure_labels},
            {"label": "Format", "options": formats},
            {"label": "Metric", "options": [m[0] for m in metrics]},
        ],
        "views": views,
    }

    views = {}
    for k, s in enumerate(structures):
        for i, fmt_ in enumerate(formats):
            series = []
            for variant in story:
                rows = sorted(
                    scaling[(fmt_, s)].get(variant, []), key=lambda r: int(r["threads"])
                )
                if rows:
                    series.append(
                        series_meta(rows[0])
                        | {
                            "pts": [
                                [
                                    int(r["threads"]),
                                    sig(num(r, "runtime_mean_s")),
                                    sig(num(r, "runtime_stddev_s") or 0, 3),
                                ]
                                for r in rows
                            ]
                        }
                    )
            views[f"{k}|{i}"] = {
                "y": {"label": "Runtime (s)", "digits": 3, "unit": "s", "scale": "log"},
                "series": series,
            }
    charts["threads"] = {
        "type": "line",
        "title": "Single-structure runtime by thread count",
        "note": "100 sphere points, three-run means. The runtime axis is logarithmic.",
        "x": {"label": "Threads", "values": [1, 4, 8, 10]},
        "controls": [
            {"label": "Structure", "options": structure_labels},
            {"label": "Format", "options": formats},
        ],
        "views": views,
    }

    for fmt_ in formats:
        rows = []
        for s in structures:
            by = {r["variant"]: r for r in t10[fmt_][s]}
            z = by["zsasa_f64"]
            rows.append(
                [
                    s,
                    f"{atoms[s]:,}",
                    fmt(num(z, "runtime_mean_s"), 3),
                    fmt(num(z, "peak_rss_mean_mib")),
                ]
                + [
                    fmt(num(z, f"runtime_speedup_vs_{c}"), 1, "×")
                    for c in ("freesasa", "rustsasa", "pdbtools_jl")
                ]
            )
        tables[f"t10-{fmt_.lower()}"] = {
            "caption": f"{fmt_} input, zsasa f64 at 10 threads and 100 sphere points.",
            "columns": [
                "Structure",
                "Atoms",
                "Runtime (s)",
                "Peak RSS (MiB)",
                "vs FreeSASA",
                "vs RustSASA",
                "vs PDBTools.jl",
            ],
            "rows": rows,
        }
    return {"charts": charts, "tables": tables}


# ---------------------------------------------------------------- validation


def band(pv, tables, reference: str, candidates: list[tuple[str, str]]) -> list[dict]:
    series = []
    for column, name in candidates:
        pts = []
        for table in tables:
            errors = pv.signed_relative_errors(table.rows, reference, column)
            if errors:
                pts.append(
                    [table.points]
                    + [
                        sig(float(v), 3)
                        for v in (
                            np.median(errors),
                            np.percentile(errors, 5),
                            np.percentile(errors, 95),
                        )
                    ]
                )
        series.append(
            {"name": name, "tool": "zsasa", "mode": mode_of(column), "pts": pts}
            | marker(column, name)
        )
    return series


def export_validation(bench: Path, t: Tables) -> dict:
    sys.path.insert(0, str(bench / "scripts"))
    import plot_validation_figures as pv  # noqa: PLC0415 (benchmark repository module)

    db = bench / "results" / "benchmark.duckdb"
    static = sorted(
        pv.load_validation_tables_from_db(
            db,
            benchmark_kind="validation",
            reference_tool_id="freesasa_batch",
            reference_column="freesasa",
            algorithm="sr",
        ),
        key=lambda x: x.points,
    )
    md = sorted(
        pv.load_validation_tables_from_db(
            db,
            benchmark_kind="trajectory_validation",
            reference_tool_id="mdtraj",
            reference_column="mdtraj",
            algorithm="sr",
        ),
        key=lambda x: x.points,
    )
    by_points = {x.points: x for x in static}
    charts, tables = {}, {}

    modes = [("zsasa_f64", "zsasa f64"), ("zsasa_bitmask_f32", "zsasa bitmask f32")]
    dense_points = [128, 1024]
    rows0 = by_points[128].rows
    ids = [r["structure"].removeprefix("AF-").split("-F")[0] for r in rows0]
    views = {}
    for i, (column, name) in enumerate(modes):
        for j, points in enumerate(dense_points):
            rows = by_points[points].rows
            assert [r["structure"] for r in rows] == [r["structure"] for r in rows0]
            ref = [float(r["freesasa"]) for r in rows]
            obs = [float(r[column]) for r in rows]
            s = pv.summarize_pair(rows, "freesasa", column)
            views[f"{i}|{j}"] = {
                "x": [round(v) for v in ref],
                "y": [
                    round(100 * (o - v) / v, 3) for v, o in zip(ref, obs, strict=True)
                ],
                "tool": "zsasa",
                "mode": mode_of(column),
                "stats": f"n = {s.n:,} · R² = {s.r2:.6f} · mean |difference| {s.mean_error_percent:.3f}% · "
                f"max {s.max_error_percent:.2f}%",
            }
    charts["static-diff"] = {
        "type": "dense",
        "title": "zsasa against FreeSASA, structure by structure",
        "note": "4,370 E. coli AlphaFold structures. Each dot is one structure: its total SASA from FreeSASA, and how far "
        "zsasa is from that value. Negative values mean zsasa reports less area.",
        "x": {
            "label": "FreeSASA total SASA (Å²)",
            "scale": "log",
            "digits": 0,
            "unit": "Å²",
        },
        "y": {"label": "Difference from FreeSASA (%)", "digits": 3, "unit": "%"},
        "ids": ids,
        "controls": [
            {"label": "Mode", "options": [n for _, n in modes]},
            {"label": "Sphere points", "options": [str(p) for p in dense_points]},
        ],
        "views": views,
    }

    charts["static-points"] = {
        "type": "band",
        "title": "Difference from FreeSASA by sphere-point count",
        "note": "4,370 E. coli AlphaFold structures. Lines are medians of the signed relative difference; bands span the "
        "5th to 95th percentile.",
        "x": {"label": "Sphere points", "values": POINTS},
        "y": {"label": "Difference from FreeSASA (%)", "digits": 3, "unit": "%"},
        "views": {
            "0": {
                "series": band(
                    pv,
                    static,
                    "freesasa",
                    [
                        ("zsasa_f64", "zsasa f64"),
                        ("zsasa_bitmask_f32", "zsasa bitmask f32"),
                    ],
                )
            }
        },
    }
    charts["md-points"] = {
        "type": "band",
        "title": "Difference from MDTraj by sphere-point count",
        "note": "1,001 frames of the 5wvo_C trajectory, MDTraj computed frame by frame. Lines are medians of the signed "
        "relative difference; bands span the 5th to 95th percentile.",
        "x": {"label": "Sphere points", "values": POINTS},
        "y": {"label": "Difference from MDTraj (%)", "digits": 3, "unit": "%"},
        "views": {
            "0": {
                "series": band(
                    pv,
                    md,
                    "mdtraj",
                    [
                        ("zsasa_cli_f64", "zsasa CLI f64"),
                        ("zsasa_mdtraj", "zsasa + MDTraj"),
                    ],
                )
            }
        },
    }
    charts["md-bitmask"] = {
        "type": "band",
        "title": "Trajectory bitmask mode against exact zsasa",
        "note": "1,001 frames of the 5wvo_C trajectory, relative to zsasa CLI f32 at the same sphere-point count. "
        "Lines are medians; bands span the 5th to 95th percentile.",
        "x": {"label": "Sphere points", "values": POINTS},
        "y": {"label": "Difference from zsasa f32 (%)", "digits": 3, "unit": "%"},
        "views": {
            "0": {
                "series": band(
                    pv,
                    md,
                    "zsasa_cli_f32",
                    [
                        ("zsasa_cli_bitmask_f32_single", "zsasa CLI bitmask f32"),
                        (
                            "zsasa_cli_bitmask_f32_single_corrected",
                            "zsasa CLI bitmask f32, experimental correction",
                        ),
                    ],
                )
            }
        },
    }

    def summary_table(
        tables_, reference: str, columns: list[tuple[str, str]], caption: str
    ) -> dict:
        rows = []
        for column, name in columns:
            for table in tables_:
                if table.points not in (128, 1024):
                    continue
                s = pv.summarize_pair(table.rows, reference, column)
                rows.append(
                    [
                        name,
                        str(table.points),
                        f"{s.r2:.6f}",
                        f"{s.mean_error_percent:.3f}%",
                        f"{s.max_error_percent:.2f}%",
                    ]
                )
        return {
            "caption": caption,
            "columns": [
                "Tool",
                "Sphere points",
                "R²",
                "Mean |difference|",
                "Max |difference|",
            ],
            "rows": rows,
        }

    tables["static-summary"] = summary_table(
        static,
        "freesasa",
        [
            ("zsasa_f64", "zsasa f64"),
            ("zsasa_f32", "zsasa f32"),
            ("zsasa_bitmask_f64", "zsasa bitmask f64"),
            ("zsasa_bitmask_f32", "zsasa bitmask f32"),
            ("rustsasa", "RustSASA"),
            ("lahuta", "Lahuta"),
            ("lahuta_bitmask", "Lahuta bitmask"),
        ],
        "Agreement with FreeSASA over 4,370 E. coli AlphaFold structures.",
    )
    tables["md-summary"] = summary_table(
        md,
        "mdtraj",
        [
            ("zsasa_cli_f64", "zsasa CLI f64"),
            ("zsasa_cli_f32", "zsasa CLI f32"),
            ("zsasa_cli_bitmask_f32_single", "zsasa CLI bitmask f32"),
            (
                "zsasa_cli_bitmask_f32_single_corrected",
                "zsasa CLI bitmask f32, experimental correction",
            ),
            ("zsasa_mdtraj", "zsasa + MDTraj"),
            ("zsasa_mdanalysis", "zsasa + MDAnalysis"),
        ],
        "Agreement with MDTraj over 1,001 frames of the 5wvo_C trajectory.",
    )
    return {"charts": charts, "tables": tables}


def main() -> None:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "benchmarks",
        type=Path,
        help="path to a zsasa-benchmarks checkout with populated results/",
    )
    args = parser.parse_args()
    bench = args.benchmarks.resolve()
    t = Tables(bench)

    commit = subprocess.run(
        ["git", "-C", str(bench), "rev-parse", "--short", "HEAD"],
        capture_output=True,
        text=True,
        check=True,
    ).stdout.strip()
    source = {
        "repository": "N283T/zsasa-benchmarks",
        "commit": commit,
        "zsasa": ZSASA_VERSION,
    }

    OUT.mkdir(parents=True, exist_ok=True)
    pages = {
        "batch": export_batch(t),
        "md": export_md(t),
        "single": export_single(t),
        "validation": export_validation(bench, t),
    }
    for name, page in pages.items():
        path = OUT / f"{name}.json"
        path.write_text(
            json.dumps(
                {"source": source} | page, ensure_ascii=False, separators=(",", ":")
            )
            + "\n"
        )
        print(
            f"{path.relative_to(OUT.parents[2])}: {len(page['charts'])} charts, {len(page['tables'])} tables, {path.stat().st_size / 1024:.0f} KiB"
        )


if __name__ == "__main__":
    main()
