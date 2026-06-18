# SASA Validation

Validation is a consistency check rather than a comparison to an external ground truth. As in the paper, SASA is defined operationally by the chosen algorithm, radii, probe radius, sampling convention, hydrogen policy, and parser decisions, so there is no single implementation-independent reference value for these inputs. Close agreement under matched settings shows that `zsasa` stays in line with established tools, not that any one tool is uniquely correct.

## Summary

| Comparison | Workload | Mode | Points | R² | Mean relative difference | Max relative difference |
| --- | --- | --- | ---: | ---: | ---: | ---: |
| FreeSASA | 4,370 *E. coli* AFDB structures | f64 | 100 | 1.000000 | 0.0000206% | 0.000205% |
| FreeSASA | 4,370 *E. coli* AFDB structures | f32 | 100 | 1.000000 | 0.000140% | 0.0150% |
| FreeSASA | 4,370 *E. coli* AFDB structures | bitmask f64 | 128 | 0.999811 | 0.662% | 2.02% |
| MDTraj | 5wvo_C, 1,001 frames | `zsasa` + MDTraj f64 | 500 | 0.9938 | 0.198% | 0.531% |
| MDTraj | 5wvo_C, 1,001 frames | `zsasa` + MDTraj f64 | 1000 | 0.9983 | 0.0998% | 0.288% |

## Static structure validation against FreeSASA

The exact Shrake--Rupley path closely reproduces FreeSASA total SASA on the *E. coli* AlphaFold Database validation set. At 100 sphere points, f64 has a mean relative difference of `2.06e-5%`; f32 remains very close, with a mean relative difference of `0.000140%`. This near-identity is expected because exact `zsasa` follows FreeSASA's golden-spiral sphere-point convention.

The bitmask path is intentionally different: it trades exact numerical identity for throughput and bounded approximation error. At 128 points, bitmask f64 and f32 both show mean relative differences of about `0.662%` versus FreeSASA. These rows should be read as measurements of a lookup-table approximation, not as a simple sampling-density convergence curve: fixed direction and angle quantization introduces quantization error and a systematic offset, while finite sphere-point sampling error can either cancel or reinforce that offset. Increasing the requested point count can reduce worst-case deviations, but it does not guarantee a monotonic decrease in mean relative difference.

The two visible scatter plots below match the two validation panels used in paper Figure 2: an exact f64 row that demonstrates near numerical identity with FreeSASA, and the bitmask f32 row used for the throughput-oriented approximation claim. The full scatter grid is still available below as supporting detail.

<div className="benchmarkFigureGrid benchmarkFigureGrid--2">
  <a href="/zsasa/assets/benchmarks/paper/validation/100p_zsasa_f64_vs_freesasa.png" target="_blank" rel="noopener noreferrer">
    <img src="/zsasa/assets/benchmarks/paper/validation/100p_zsasa_f64_vs_freesasa.png" alt="100-point zsasa f64 versus FreeSASA scatter plot" />
  </a>
  <a href="/zsasa/assets/benchmarks/paper/validation/128p_zsasa_bitmask_f32_vs_freesasa.png" target="_blank" rel="noopener noreferrer">
    <img src="/zsasa/assets/benchmarks/paper/validation/128p_zsasa_bitmask_f32_vs_freesasa.png" alt="128-point zsasa bitmask f32 versus FreeSASA scatter plot" />
  </a>
</div>

<p className="benchmarkFigureCaption"><strong>Figures 1-2. Representative static-validation scatter plots.</strong> Left: exact <code>zsasa</code> f64 at 100 sphere points (R² = 1.000000, mean relative difference 2.1×10⁻⁵%). Right: <code>zsasa</code> bitmask f32 at 128 sphere points (R² = 0.999811, mean relative difference 0.66%). Click either panel to open the full-size PNG.</p>

[![Static validation mean relative error](pathname:///zsasa/assets/benchmarks/paper/details/static_sr_mean_relative_error.png)](/assets/benchmarks/paper/details/static_sr_mean_relative_error.png)

**Figure 3. Static validation error across point counts.** This summary view shows how agreement changes across the point-count sweep. Exact `zsasa` modes remain nearly identical to FreeSASA because they share the same golden-spiral sphere-point convention. `zsasa` bitmask mode has a visible but quantified approximation envelope from lookup-table quantization error. RustSASA and Lahuta use the same point-placement convention in this sweep, so their curves are nearly overlapping. In the paper-aligned comparison, Lahuta bitmask is applied only at the 128-point, Lahuta-compatible setting used for the paper/batch comparison rather than as the main point-count trend.

<details>
<summary>Show static-validation scatter grid</summary>

[![Static validation scatter grid](pathname:///zsasa/assets/benchmarks/paper/details/static_sr_scatter_grid.png)](/assets/benchmarks/paper/details/static_sr_scatter_grid.png)

The scatter grid is retained as supporting evidence. Click the image to inspect the full-size grid.

</details>

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

**Figure 4. MD validation R² across point counts.** Higher point counts reduce the implementation-specific sampling difference against MDTraj. This summary is easier to scan than the full scatter-grid output.

<details>
<summary>Show MD-validation scatter grid</summary>

[![MD validation scatter grid](pathname:///zsasa/assets/benchmarks/paper/details/md_scatter_grid.png)](/assets/benchmarks/paper/details/md_scatter_grid.png)

The scatter grid is retained as supporting evidence. Click the image to inspect the full-size grid.

</details>

## Reproducibility notes

- Static validation used 4,370 *E. coli* K-12 AFDB structures.
- Static validation points: 100, 128, 200, 500, and 1,000.
- MD validation points: 100, 200, 500, and 1,000.
- Validation records were imported into DuckDB and exported through `zsasa-benchmarks/results/tables/validation_pairwise_summary.csv`.
