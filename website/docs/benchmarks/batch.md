# Batch Processing Benchmarks

Batch benchmarks measure complete directory processing: parsing, SASA calculation, output writing, and worker scheduling. The current current pinned results use `zsasa` v0.6.0, pinned comparator builds, 128 sphere points, and 10 threads unless noted.

## TL;DR

| Dataset | Structures | Best `zsasa` mode | Runtime | Throughput | RSS | Speedup vs FreeSASA batch |
| --- | ---: | --- | ---: | ---: | ---: | ---: |
| *E. coli* AFDB | 4,370 | bitmask f32 | 1.481 s | 2,951 str/s | 45.1 MiB | 8.77× |
| Human AFDB | 23,586 | bitmask f32 | 13.814 s | 1,707 str/s | 79.5 MiB | 9.70× |

FreeSASA has no native directory-batch command, so the FreeSASA batch rows use a thin `freesasa_batch` wrapper around the pinned FreeSASA C API.

## *E. coli* AFDB batch

The *E. coli* AFDB collection contains 4,370 structures and is used for both comparator benchmarking and thread scaling.

| Tool | Runtime | Structures/s | RSS | Speedup vs FreeSASA | Speedup vs RustSASA | Speedup vs Lahuta bitmask |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| `zsasa` f64 | 4.411 s | 991 | 45.5 MiB | 2.94× | 1.31× | 0.46× |
| `zsasa` f32 | 4.246 s | 1,029 | 43.5 MiB | 3.06× | 1.36× | 0.48× |
| `zsasa` bitmask f64 | 1.504 s | 2,905 | 48.5 MiB | 8.63× | 3.84× | 1.35× |
| `zsasa` bitmask f32 | **1.481 s** | **2,951** | **45.1 MiB** | **8.77×** | **3.90×** | **1.37×** |
| FreeSASA batch | 12.982 s | 337 | 212.6 MiB | baseline | 0.44× | 0.16× |
| RustSASA | 5.769 s | 757 | 170.5 MiB | 2.25× | baseline | 0.35× |
| Lahuta bitmask | 2.034 s | 2,148 | 180.7 MiB | 6.38× | 2.84× | baseline |

![E. coli throughput vs peak RSS](pathname:///zsasa/assets/benchmarks/paper/batch/ecoli_t10_throughput_vs_peak_rss.png)

**Figure 1. E. coli throughput versus peak RSS.** This map shows the practical trade-off directly: higher throughput and lower RSS are better.

![E. coli throughput and RSS bars](pathname:///zsasa/assets/benchmarks/paper/batch/ecoli_throughput_rss_bar_2grid.png)

**Figure 2. E. coli absolute throughput and peak RSS.** The bar view complements the map by showing the absolute values for each tool and mode.

### Thread scaling

| Mode | 1 thread | 4 threads | 8 threads | 10 threads | 10-thread speedup |
| --- | ---: | ---: | ---: | ---: | ---: |
| `zsasa` f64 | 147 str/s | 563 str/s | 871 str/s | 991 str/s | 6.72× |
| `zsasa` bitmask f32 | 435 str/s | 1,688 str/s | 2,623 str/s | 2,951 str/s | 6.78× |
| FreeSASA batch | 105 str/s | 366 str/s | 352 str/s | 337 str/s | 3.21× |
| RustSASA | 142 str/s | 538 str/s | 697 str/s | 757 str/s | 5.34× |
| Lahuta bitmask | 325 str/s | 1,300 str/s | 1,917 str/s | 2,148 str/s | 6.61× |

![E. coli throughput by thread count](pathname:///zsasa/assets/benchmarks/paper/details/ecoli_throughput_vs_threads.png)

**Figure 3. E. coli throughput by thread count.** `zsasa` exact and bitmask modes scale strongly up to 10 threads on the M4 benchmark machine.

## Human AFDB batch

The Human AFDB collection contains 23,586 structures, 5.4× more than the *E. coli* AFDB collection. `zsasa` retains low peak memory because batch mode streams structures instead of holding the full collection in memory.

| Tool | Runtime | Structures/s | RSS | Speedup vs FreeSASA | Speedup vs RustSASA | Speedup vs Lahuta bitmask |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| `zsasa` f64 | 45.508 s | 518 | 82.3 MiB | 2.94× | 1.48× | 0.43× |
| `zsasa` f32 | 44.352 s | 532 | 76.4 MiB | 3.02× | 1.52× | 0.44× |
| `zsasa` bitmask f64 | 14.150 s | 1,667 | 83.7 MiB | 9.47× | 4.77× | 1.38× |
| `zsasa` bitmask f32 | **13.814 s** | **1,707** | **79.5 MiB** | **9.70×** | **4.89×** | **1.42×** |
| FreeSASA batch | 133.960 s | 176 | 627.5 MiB | baseline | 0.50× | 0.15× |
| RustSASA | 67.531 s | 349 | 330.0 MiB | 1.98× | baseline | 0.29× |
| Lahuta bitmask | 19.566 s | 1,205 | 326.9 MiB | 6.85× | 3.45× | baseline |

![Human throughput vs peak RSS](pathname:///zsasa/assets/benchmarks/paper/batch/human_t10_throughput_vs_peak_rss.png)

**Figure 4. Human throughput versus peak RSS.** `zsasa` keeps the same low-memory/high-throughput position on the larger Human AFDB collection.

![Human throughput and RSS bars](pathname:///zsasa/assets/benchmarks/paper/batch/human_throughput_rss_bar_2grid.png)

**Figure 5. Human absolute throughput and peak RSS.** The absolute bars make the memory difference visible without relying on speedup ratios.

## Legacy SwissProt benchmark (pre-pinned)

The SwissProt result below is retained only as historical context from the earlier website benchmark set. It was collected before the current `zsasa-benchmarks` pinned v0.6.0 harness and should not be mixed with the current pinned headline claims above.

Dataset: SwissProt PDB v6, 550,122 structures, PDB format. Benchmark settings: warmup=3, runs=3, threads=10.

### M2 Max, 96 GB

| Tool | Time | files/s | RSS |
| --- | ---: | ---: | ---: |
| `zsasa` bitmask f32 | **4m 02s** | **2,269** | **157 MB** |
| `zsasa` bitmask f64 | 4m 07s | 2,229 | 162 MB |
| Lahuta bitmask | 5m 12s | 1,761 | 2,187 MB |
| RustSASA | 10m 58s | 835 | 1,131 MB |
| FreeSASA | 32m 21s | 283 | 2,875 MB |

### M4, 32 GB

| Tool | Time | files/s | RSS |
| --- | ---: | ---: | ---: |
| `zsasa` bitmask f32 | 11m 05s | 828 | 157 MB |
| `zsasa` bitmask f64 | 11m 07s | 824 | 161 MB |
| Lahuta bitmask | 11m 08s | 823 | 2,152 MB |
| RustSASA | 26m 16s | 349 | 1,131 MB |
| FreeSASA | 31m 42s | 289 | 2,440 MB |

On the 32 GB M4 system, the dataset exceeded available RAM and the run became I/O-bound; `zsasa` and Lahuta bitmask converged in wall-clock time, while `zsasa` retained much lower peak RSS.

## Reproducing the current batch results

The current pinned benchmark data is exported from `zsasa-benchmarks/results/tables/batch_t10_summary.csv` and `batch_thread_scaling.csv`. The corresponding full harness is in [`N283T/zsasa-benchmarks`](https://github.com/N283T/zsasa-benchmarks).
