# SASA Validation

Validation compares total SASA from `zsasa` with established implementations on the same inputs: FreeSASA for static structures and MDTraj for trajectory frames. All values are from `zsasa` 0.9.0.

## Static structure validation against FreeSASA

The static set is the *E. coli* K-12 AlphaFold proteome, 4,370 structures, evaluated at 64, 128, 256, 512 and 1,024 sphere points.

<div data-chart="validation/static-diff"></div>

In exact f64 mode `zsasa` and FreeSASA agree to within floating-point noise on every structure, so the points form a flat line at zero. Bitmask f32 reports slightly less area than FreeSASA, with a spread that narrows as the structure grows.

<div data-chart="validation/static-points"></div>

The bitmask difference is a systematic offset rather than sampling noise. The lookup table quantizes directions and angles in fixed steps, so the median stays near −0.7% at every point count, and more sphere points narrow the spread without moving the median.

For scale, two exact implementations differ by a similar amount: in the table below RustSASA differs from FreeSASA by 0.28% on average at 128 points, with a maximum of 2.15%. The bitmask offset is small enough for screening, ranking and feature generation, where SASA is used as a relative descriptor. Use an exact mode when values must match FreeSASA.

<div data-table="validation/static-summary"></div>

"Mean" and "max" are over the absolute relative difference from FreeSASA.

## Trajectory validation against MDTraj

Trajectory validation uses 1,001 frames of the 5wvo_C ATLAS trajectory, with explicit hydrogens and the NACCESS classifier. The MDTraj reference is computed one frame at a time.

<div data-chart="validation/md-points"></div>

Called through its MDTraj integration, which uses MDTraj's element radii, `zsasa` agrees with MDTraj without bias: the median difference is within 0.02% at every point count. The band narrows as the point count rises, because the two tools place their sphere points along different axes and that difference averages out.

The CLI differs from MDTraj by a constant −0.3%. That offset comes from the radius set, not the algorithm: the CLI runs use the NACCESS classifier, and MDTraj uses its own element radii.

<div data-chart="validation/md-bitmask"></div>

On this explicit-hydrogen trajectory the raw bitmask result is about 1.7% below exact `zsasa`, a larger offset than for static structures. Bitmask mode is therefore not recommended for explicit-hydrogen trajectories. An experimental bias-correction option removes most of the offset; it may change or be removed in a later release.

<div data-table="validation/md-summary"></div>

## Reproducibility notes

- Static validation: 4,370 *E. coli* K-12 AFDB structures, 52 conditions.
- Trajectory validation: 5wvo_C, 1,001 frames, 75 conditions.
- Sphere-point counts: 64, 128, 256, 512 and 1,024.
- The charts are computed from the per-structure and per-frame records in the `zsasa-benchmarks` result database, using the same run selection as the published figures.
