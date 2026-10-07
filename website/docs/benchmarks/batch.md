# Batch Processing Benchmarks

Batch benchmarks time a complete directory run: parsing, SASA calculation, output writing and worker scheduling. All runs use `zsasa` 0.9.0, pinned comparator builds and 128 sphere points. The comparator builds did not change between `zsasa` releases, so their timings are retained from the 0.6.0 session on the same machine.

FreeSASA has no native directory mode, so its rows use a thin `freesasa_batch` wrapper around the pinned FreeSASA C API.

## Throughput and memory

<div data-chart="batch/map"></div>

On both proteomes the `zsasa` bitmask modes sit in the top-left corner: the highest throughput at the lowest peak memory. The exact f64 and f32 modes are faster than FreeSASA batch and RustSASA while using a quarter of their memory or less.

<div data-chart="batch/bars"></div>

### *E. coli* AFDB

<div data-table="batch/t10-ecoli-afdb"></div>

### Human AFDB

<div data-table="batch/t10-human-afdb"></div>

The "vs" columns are runtime ratios: the comparator's runtime divided by the row's runtime.

## Thread scaling

<div data-chart="batch/scaling"></div>

From 1 to 10 threads, `zsasa` f64 and bitmask f32 both gain a factor of about 6.9 on a machine with 4 performance and 6 efficiency cores. FreeSASA batch peaks at 4 threads, and Lahuta bitmask stops improving between 8 and 10.

## Beyond the core count

`zsasa` accepts more worker threads than the machine has cores. The benchmark machine has 10 logical CPUs, so 20 and 40 workers are overcommitted.

<div data-chart="batch/overcommit"></div>

On the *E. coli* and Human runs, extra workers leave throughput unchanged within about 2% and raise peak memory roughly in proportion to the worker count. SwissProt is the exception: throughput keeps rising at 20 and 40 workers.

The "Busy cores" metric shows why. The *E. coli* and Human collections occupy 0.8 and 9.1 GB on disk and fit in the file cache of the 32 GB machine, and `zsasa` keeps about 9 of the 10 cores busy on them. SwissProt occupies 119 GB. With 10 workers the bitmask run keeps only 2.3 cores busy on average while the workers wait for input, and 40 workers raise that to 4.1. File input limits the largest run at least as much as the SASA calculation does. Using more workers than cores is therefore an option for runs like this one, not a general setting.

## Human AFDB, mmCIF input {#mmcif}

The mmCIF runs compare two parser paths, the generic mmCIF parser and a fast path for AlphaFold models ("AF fast"), with two input strategies, reading the whole file or memory-mapping it. All four are bitmask f32 runs.

<div data-table="batch/t10-human-afdb-mmcif"></div>

The AF fast path with read-all input is the fastest configuration, about 10% ahead of Lahuta bitmask at under a third of its peak memory. The generic parser matches Lahuta bitmask on throughput.

## SwissProt AFDB {#swissprot}

SwissProt is the largest dataset, 550,122 structures. Each configuration was measured once, so this section reports single observations without an uncertainty estimate.

<div data-table="batch/t10-swissprot-afdb"></div>

At 10 threads Lahuta bitmask is the fastest tool on this dataset, about 14% ahead of `zsasa` 0.9.0 bitmask f32, while using about 2.5 times the memory. This is Lahuta's home ground: its bitmask SASA mode reads only AlphaFold models, and SwissProt AFDB is half a million of them.

Because this run is limited by file input (see [Beyond the core count](#beyond-the-core-count)), worker overcommit pays off here. With 20 workers `zsasa` 0.9.0 bitmask f32 is about 23% above Lahuta bitmask's 10-thread result, and with 40 workers about 31% above it, at half to 70% of Lahuta's peak memory. Lahuta has no equivalent of worker overcommit and was measured at 10 threads.

`zsasa` 0.9.0 and 0.6.0 are within 2% of each other at 10 threads. The "generic/read" row is `zsasa` 0.9.0 bitmask f32 with read-all input instead of memory mapping.

## Reproducing these results

The charts and tables on this page are generated from `results/tables/batch_t10_summary.csv` and `batch_thread_scaling.csv` in [`N283T/zsasa-benchmarks`](https://github.com/N283T/zsasa-benchmarks), which also holds the full harness and manifests.
