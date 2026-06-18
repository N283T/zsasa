# SASA Validation

Validation checks whether `zsasa` agrees with established tools under matched settings. It is not an implementation-independent ground truth: SASA depends on the algorithm, radii, probe radius, sampling convention, hydrogens, and parser policy.

## Summary

| Comparison | Workload | Mode | Points | R² | Mean relative difference | Max relative difference |
| --- | --- | --- | ---: | ---: | ---: | ---: |
| FreeSASA | 4,370 *E. coli* AFDB structures | f64 | 100 | 1.000000 | 0.0000206% | 0.000205% |
| FreeSASA | 4,370 *E. coli* AFDB structures | f32 | 100 | 1.000000 | 0.000140% | 0.0150% |
| FreeSASA | 4,370 *E. coli* AFDB structures | bitmask f64 | 128 | 0.999811 | 0.662% | 2.02% |
| MDTraj | 5wvo_C, 1,001 frames | `zsasa` + MDTraj f64 | 500 | 0.9938 | 0.198% | 0.531% |
| MDTraj | 5wvo_C, 1,001 frames | `zsasa` + MDTraj f64 | 1000 | 0.9983 | 0.0998% | 0.288% |

## Static structure validation against FreeSASA

The exact Shrake--Rupley path closely reproduces FreeSASA total SASA on the *E. coli* AlphaFold Database validation set. At 100 sphere points, f64 has a mean relative difference of `2.06e-5%`; f32 remains very close, with a mean relative difference of `0.000140%`.

The bitmask path is intentionally different: it trades exact numerical identity for throughput and bounded approximation error. At 128 points, bitmask f64 and f32 both show mean relative differences of about `0.662%` versus FreeSASA.

[![Static validation mean relative error](pathname:///zsasa/assets/benchmarks/paper/details/static_sr_mean_relative_error.png)](/assets/benchmarks/paper/details/static_sr_mean_relative_error.png)

**Figure 1. Static validation error across point counts.** This summary view is easier to scan than the full scatter-grid output. Exact `zsasa` modes remain nearly identical to FreeSASA; bitmask mode has a visible but quantified approximation envelope.

## Trajectory validation against MDTraj

Trajectory validation uses the 5wvo_C ATLAS trajectory with 1,001 frames. Agreement improves as sphere-point count increases because MDTraj and `zsasa` use different sphere-point conventions at low point counts but converge toward the same surface area.

| Path | Points | R² | Mean relative difference | Max relative difference |
| --- | ---: | ---: | ---: | ---: |
| `zsasa` + MDTraj f64 | 100 | 0.872 | 0.946% | 1.90% |
| `zsasa` + MDTraj f64 | 500 | 0.9938 | 0.198% | 0.531% |
| `zsasa` + MDTraj f64 | 1000 | 0.9983 | 0.0998% | 0.288% |
| CLI f64 | 100 | 0.7725 | 1.27% | 2.66% |
| CLI f64 | 1000 | 0.9722 | 0.416% | 1.03% |

[![MD validation R2](pathname:///zsasa/assets/benchmarks/paper/details/md_r2.png)](/assets/benchmarks/paper/details/md_r2.png)

**Figure 2. MD validation R² across point counts.** Higher point counts reduce the implementation-specific sampling difference against MDTraj. The full scatter grids are available from the benchmark repository outputs, but the website keeps the validation page focused on summary readouts.

## Reproducibility notes

- Static validation used 4,370 *E. coli* K-12 AFDB structures.
- Static validation points: 100, 128, 200, 500, and 1,000.
- MD validation points: 100, 200, 500, and 1,000.
- Validation records were imported into DuckDB and exported through `zsasa-benchmarks/results/tables/validation_pairwise_summary.csv`.
