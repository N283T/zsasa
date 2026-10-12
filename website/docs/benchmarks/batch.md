# Batch Processing Benchmarks

Batch benchmarks time a complete directory run: parsing, SASA calculation, output writing and worker scheduling. The runs use `zsasa` 0.9.0, pinned comparator builds and 128 sphere points; rows measured with a later release name it. The comparator builds did not change between `zsasa` releases, so their timings are retained from the 0.6.0 session on the same machine.

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

On the *E. coli* and Human runs, extra workers leave throughput unchanged within about 2% and raise peak memory roughly in proportion to the worker count. SwissProt is the exception: bitmask throughput is about 10% higher at 20 workers.

The "Busy cores" metric shows why. The *E. coli* and Human collections occupy 0.8 and 9.1 GB on disk and fit in the file cache of the 32 GB machine, and `zsasa` keeps about 9 of the 10 cores busy on them. SwissProt occupies 119 GB, and with 10 workers the bitmask run keeps only 3.4 cores busy on average, 4.1 with 40 workers. The idle workers wait inside `open()`: the benchmark machine runs endpoint-security software that scans every file when it is first opened (see [SwissProt AFDB](#swissprot)). Using more workers than cores is therefore an option for runs that wait on file access, not a general setting.

## Human AFDB, mmCIF input {#mmcif}

The mmCIF runs compare two parser paths, the generic mmCIF parser and a fast path for AlphaFold models ("AF fast"), with two input strategies, reading the whole file or memory-mapping it. All four are bitmask f32 runs of `zsasa` 0.11.0: batch runs of 0.9.0 to 0.10.1 ignored `--input-io`, so the comparison needs the later release. Lahuta was measured in the same session.

<div data-table="batch/t10-human-afdb-mmcif"></div>

The AF fast path is about 5% ahead of Lahuta bitmask at under a third of its peak memory, and the generic parser is about 6% behind. Reading and memory-mapping the input take the same time within 2% on this warm-cache run; memory-mapping has the lower peak memory, 84 against 95 MiB for the AF fast path at 10 threads and 244 against 323 MiB at 40 workers.

## SwissProt AFDB {#swissprot}

SwissProt is the largest dataset, 550,122 structures. Each configuration was measured once, so this section reports single observations without an uncertainty estimate.

The benchmark machine scans every file when it is first opened, and this collection is too large for its files to stay cached and scanned between runs. macOS also slows that scan down while a tool reads from the same disk with the default I/O policy. The SwissProt commands were therefore run under the passive disk I/O policy (`taskpolicy -d passive`), for `zsasa` and Lahuta alike.

<div data-table="batch/t10-swissprot-afdb"></div>

At 10 threads `zsasa` 0.9.0 bitmask f32 processes the collection in 7.7 minutes, about 6% faster than Lahuta bitmask at 40% of its peak memory. Lahuta's bitmask SASA mode reads only AlphaFold models and was measured at 10 threads; it has no equivalent of worker overcommit. With 20 and 40 workers `zsasa` 0.9.0 bitmask f32 reaches 1,312 and 1,283 structures/s at 343 and 472 MiB.

The `zsasa` 0.11.0 row is one run of the later release under the same policy. From 0.11.0 on a batch reads each input with one read instead of memory-mapping it, which shortens the wait for files on collections larger than memory. The release also adds [`--passive-io`](../guide/batch.md#file-scanners-on-macos), which requests the passive policy from within `zsasa`.

Under the default I/O policy both tools are slower on this machine and the order at 10 threads is reversed: `zsasa` 0.9.0 bitmask f32 runs at 821 structures/s and Lahuta bitmask at 935. The [investigation notes](https://github.com/N283T/zsasa-benchmarks/tree/main/investigations/swissprot-open-latency) in the benchmark repository have the measurements behind this.

## Reproducing these results

The charts and tables on this page are generated from `results/tables/batch_t10_summary.csv` and `batch_thread_scaling.csv` in [`N283T/zsasa-benchmarks`](https://github.com/N283T/zsasa-benchmarks), which also holds the full harness and manifests.
