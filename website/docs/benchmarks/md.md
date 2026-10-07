# MD Trajectory Benchmarks

Trajectory benchmarks time frame-by-frame SASA over a whole trajectory, including reading the trajectory. All runs use `zsasa` 0.9.0, 128 sphere points, stride 1, the NACCESS classifier and explicit hydrogens.

## Throughput summary

<div data-chart="md/map"></div>

Native `zsasa` occupies the top-left corner on every trajectory: the highest frame rate at the lowest peak memory. The Python integrations, zsasa + MDTraj and zsasa + MDAnalysis, reach a similar frame rate, and their memory is set by the Python trajectory loader.

<div data-chart="md/bars"></div>

## Per-dataset comparison

### 5wvo_C

<div data-table="md/summary-5wvo-c"></div>

### 6sup_A

<div data-table="md/summary-6sup-a"></div>

### 5vz0_A

MDTraj was not measured on this trajectory.

<div data-table="md/summary-5vz0-a"></div>

The "vs" columns are runtime ratios: the comparator's runtime divided by the row's runtime.

:::warning[Bitmask rows on these trajectories]
The bitmask rows use the single-LUT mode with an experimental bias-correction option, which may change or be removed. Without it, bitmask mode is about 2% below MDTraj on an explicit-hydrogen trajectory, so it is not recommended for such trajectories. See [SASA Validation](validation.md#trajectory-validation-against-mdtraj) before relying on the bitmask frame rates.
:::

## Beyond the core count

The native trajectory path accepts more worker threads than the machine has cores. The benchmark machine has 10 logical CPUs.

<div data-chart="md/overcommit"></div>

Extra workers help a little here: 40 workers are 2% to 11% faster than 10, with the largest gain on the smallest system. Per-frame SASA values were identical at 10, 20 and 40 workers.

## Workloads

| Dataset | Frames | Atoms | Source | Use |
| --- | ---: | ---: | --- | --- |
| 5wvo_C | 1,001 | 3,858 | ATLAS | validation and throughput |
| 6sup_A | 1,001 | 33,377 | ATLAS | large-system throughput |
| 5vz0_A | 10,001 | 17,910 | ATLAS | long-trajectory throughput |

## Memory interpretation

`zsasa` streams trajectory frames and keeps memory close to the current-frame working set. In contrast, the mdsasa-bolt path uses an MDAnalysis front-end that materializes atom data for every frame before the Rust SASA core runs, so its peak memory grows with trajectory length.

## Validation pointer

Agreement with MDTraj on 5wvo_C is covered in [SASA Validation](validation.md#trajectory-validation-against-mdtraj).

## Evidence source

The charts and tables on this page are generated from `results/tables/md_summary.csv` and `md_thread_scaling.csv` in [`N283T/zsasa-benchmarks`](https://github.com/N283T/zsasa-benchmarks).
