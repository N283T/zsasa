# Batch Processing Benchmarks

Batch benchmarks time a complete directory run: parsing, SASA calculation, output writing and worker scheduling. The runs use `zsasa` 0.12.1, pinned comparator builds and 128 sphere points, and the comparators were measured in the same session. Rows measured with another `zsasa` release name it.

FreeSASA has no native directory mode, so its rows use a thin `freesasa_batch` wrapper around the pinned FreeSASA C API.

## Throughput and memory

<div data-chart="batch/map"></div>

On both proteomes the `zsasa` bitmask modes sit in the top-left corner: the highest throughput at the lowest peak memory. The exact f64 and f32 modes are faster than FreeSASA batch and RustSASA. The comparators use 3.1 to 4.2 times the peak memory of `zsasa` on *E. coli* and 3.4 to 7.2 times on Human.

<div data-chart="batch/bars"></div>

### *E. coli* AFDB

<div data-table="batch/t10-ecoli-afdb"></div>

### Human AFDB

<div data-table="batch/t10-human-afdb"></div>

The "vs" columns are runtime ratios: the comparator's runtime divided by the row's runtime.

## Thread scaling

<div data-chart="batch/scaling"></div>

From 1 to 10 threads, `zsasa` f64 gains a factor of 6.6 and bitmask f32 a factor of 6.7 on a machine with 4 performance and 6 efficiency cores. Parallel efficiency is 0.95 at 4 threads and falls once the threads exceed the performance cores. FreeSASA batch does not scale beyond 4 threads.

## Beyond the core count

`zsasa` accepts more worker threads than the machine has cores. The benchmark machine has 10 logical CPUs, so 20 and 40 workers are overcommitted.

<div data-chart="batch/overcommit"></div>

On the *E. coli* and Human runs, extra workers change throughput by less than 3% and raise peak memory with the worker count. On SwissProt, 20 workers raise bitmask throughput by 4% and 40 workers lower it by 2%.

The "Busy cores" metric shows why. The *E. coli* and Human collections occupy 0.8 and 9.1 GB on disk and fit in the file cache of the 32 GB machine, and `zsasa` keeps 8.6 to 9.6 of the 10 cores busy on them. SwissProt occupies 119 GB, and with 10 workers the bitmask run keeps only 4.8 cores busy on average. The idle workers wait inside `open()`: the benchmark machine runs endpoint-security software that scans every file when it is first opened (see [SwissProt AFDB](#swissprot)). Using more workers than cores is therefore an option for runs that wait on file access, not a general setting.

## Human AFDB, mmCIF input {#mmcif}

The mmCIF runs compare two parser paths, the generic mmCIF parser and a fast path for AlphaFold models ("AF fast"), with two input strategies, reading the whole file or memory-mapping it. All four are bitmask f32 runs of `zsasa` 0.12.1, and Lahuta was measured in the same session.

<div data-table="batch/t10-human-afdb-mmcif"></div>

The AF fast path is about 4% ahead of Lahuta bitmask at about a quarter of its peak memory, and the generic parser is 6 to 7% behind. Reading and memory-mapping the input take the same time within 2% on this warm-cache run; memory-mapping has the lower peak memory, 84 against 95 MiB for the AF fast path at 10 threads and 240 against 319 MiB at 40 workers.

## SwissProt AFDB {#swissprot}

SwissProt is the largest dataset, 550,122 structures. Each configuration was measured once, so this section reports single observations without an uncertainty estimate.

The benchmark machine scans every file when it is first opened, and this collection is too large for its files to stay cached and scanned between runs. macOS also slows that scan down while a tool reads from the same disk with the default I/O policy. The SwissProt commands were therefore run under the passive disk I/O policy (`taskpolicy -d passive`), for `zsasa` and Lahuta alike.

<div data-table="batch/t10-swissprot-afdb"></div>

At 10 threads `zsasa` 0.12.1 bitmask f32 processes the collection in 5.7 minutes at 1,599 structures/s, against 1,269 structures/s for Lahuta bitmask, at 44% of its peak memory. Exact f32 runs at 735 structures/s. Lahuta's bitmask SASA mode reads only AlphaFold models and was measured at 10 threads; it has no equivalent of worker overcommit. With 20 and 40 workers `zsasa` 0.12.1 bitmask f32 reaches 1,666 and 1,566 structures/s at 373 and 531 MiB.

Since 0.11.0, [`--passive-io`](../guide/batch.md#file-scanners-on-macos) requests the passive policy from within `zsasa`.

The same commands under the default I/O policy were measured in July 2026 with `zsasa` 0.9.0 and not repeated with 0.12.1. Under that policy both tools were slower on this machine and the order at 10 threads was reversed: `zsasa` 0.9.0 bitmask f32 ran at 821 structures/s and Lahuta bitmask at 935. The [investigation notes](https://github.com/N283T/zsasa-benchmarks/tree/main/investigations/swissprot-open-latency) in the benchmark repository have the measurements behind this.

## Reproducing these results

The charts and tables on this page are generated from `results/tables/batch_t10_summary.csv` and `batch_thread_scaling.csv` in [`N283T/zsasa-benchmarks`](https://github.com/N283T/zsasa-benchmarks), which also holds the full harness and manifests.
