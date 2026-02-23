# SASA Validation

Accuracy comparison of zsasa against reference implementations, independent of timing benchmarks.

## Test Environment

| Item | Value |
|------|-------|
| Machine | MacBook Pro |
| Chip | Apple M4 (10 cores: 4P + 6E) |
| Memory | 32 GB |
| OS | macOS |

## Static Structure (PDB)

Compares total SASA per structure across an entire proteome dataset.

### Shrake-Rupley (E. coli Proteome, 4,370 structures)

Dataset: AlphaFold E. coli K-12 proteome, n_points=100.

| Tool | N | R² | Mean Error % | Max Error % |
|------|--:|---:|--------------:|------------:|
| zsasa f64 | 4,370 | 1.000000 | 0.0000 | 0.0000 |
| zsasa f32 | 4,370 | 1.000000 | 0.0001 | 0.0150 |

- **f64**: Bit-identical to FreeSASA (both use the same algorithm parameters)
- **f32**: Max 0.015% error from floating-point rounding — negligible for practical use

![SR validation scatter](pathname://../../benchmarks/results/validation/ecoli/sr/validation_sr.png)

### Lee-Richards (E. coli Proteome, 4,370 structures)

Dataset: AlphaFold E. coli K-12 proteome, n_slices=20.

| Tool | N | R² | Mean Error % | Max Error % |
|------|--:|---:|--------------:|------------:|
| zsasa f32 | 4,370 | 1.000000 | 0.2214 | 0.3097 |
| zsasa f64 | 4,370 | 1.000000 | 0.2214 | 0.3102 |

- **f32 and f64 are identical** — LR error comes from slice discretization differences, not floating-point precision
- Mean error ~0.22% is higher than SR (<0.001%) due to different slicing implementations between zsasa and FreeSASA
- R² = 1.000000 confirms strong linear agreement

![LR validation scatter](pathname://../../benchmarks/results/validation/ecoli/lr/validation_lr.png)

## MD Trajectory

Compares per-frame total SASA across trajectory tools to validate:

1. **SASA accuracy**: zsasa vs MDTraj at varying n_points
2. **XTC reader consistency**: Python bindings (mdtraj I/O) vs CLI (native Zig XTC reader)

### Tools

| Tool | SASA Engine | XTC Reader | Notes |
|------|-------------|------------|-------|
| **mdtraj** | MDTraj (C) | MDTraj | Reference implementation |
| **zsasa_mdtraj** | zsasa (Zig) | MDTraj | Python bindings, mdtraj loads trajectory |
| **zsasa_cli** | zsasa (Zig) | Zig (native) | CLI with built-in XTC reader |

### 5wvo_C (DNMT1, 250 residues, 1,001 frames)

#### SASA Accuracy vs MDTraj (n_points convergence)

Reference: mdtraj at n_points=960.

| Tool | n_points | R² | Mean Error % | Max Error % |
|------|--------:|----|-------------:|------------:|
| mdtraj | 100 | 0.992886 | 0.8707 | 1.6661 |
| zsasa_mdtraj | 100 | 0.992725 | 0.2020 | 0.9218 |
| zsasa_cli | 100 | 0.992727 | 0.2020 | 0.9221 |
| mdtraj | 500 | 0.999051 | 0.1073 | 0.3291 |
| zsasa_mdtraj | 500 | 0.999054 | 0.1111 | 0.3235 |
| zsasa_cli | 500 | 0.999054 | 0.1110 | 0.3233 |
| zsasa_mdtraj | 960 | 0.999444 | 0.1050 | 0.3222 |
| zsasa_cli | 960 | 0.999444 | 0.1050 | 0.3229 |

**Observations:**
- At n_points=100, MDTraj native shows ~0.87% mean error vs n_points=960, while zsasa shows ~0.20% — both use Golden Section Spiral ([MDTraj `sasa.cpp:146-176`](https://github.com/mdtraj/mdtraj/blob/2f4b592f/mdtraj/geometry/src/sasa.cpp#L146-L176)) but with slightly different implementations (starting point, longitude accumulation), causing divergence at low n_points
- At n_points=500+, all tools converge to <0.12% mean error
- zsasa_mdtraj and zsasa_cli produce nearly identical results at each n_points (same SASA engine)

#### XTC Reader Consistency (Python Bindings vs CLI)

| Reference | Tool | R² | Mean Error % | Max Error % |
|-----------|------|----|-------------:|------------:|
| zsasa_mdtraj (100) | zsasa_cli (100) | 1.000000 | 0.0003 | 0.0132 |
| zsasa_mdtraj (500) | zsasa_cli (500) | 1.000000 | 0.0003 | 0.0050 |
| zsasa_mdtraj (960) | zsasa_cli (960) | 1.000000 | 0.0002 | 0.0017 |

- **R² = 1.000000** across all n_points — Python bindings and CLI read XTC identically
- Max error <0.014% comes from floating-point coordinate precision differences between MDTraj's C reader and zsasa's Zig XTC reader
- Error decreases with higher n_points (more points smooth out coordinate-level noise)

![MD validation n_points=100](pathname://../../benchmarks/results/validation_md/5wvo_C_R1/validation_md_100.png)
![MD validation n_points=500](pathname://../../benchmarks/results/validation_md/5wvo_C_R1/validation_md_500.png)
![MD validation n_points=960](pathname://../../benchmarks/results/validation_md/5wvo_C_R1/validation_md_960.png)

## Running Validation

### Static Structure (PDB)

```bash
# Shrake-Rupley (both f32 and f64 vs FreeSASA)
./benchmarks/scripts/validation.py run \
    -i benchmarks/UP000000625_83333_ECOLI_v6/pdb \
    -n ecoli --algorithm sr
# -> benchmarks/results/validation/ecoli/sr/

# Lee-Richards
./benchmarks/scripts/validation.py run \
    -i benchmarks/UP000000625_83333_ECOLI_v6/pdb \
    -n ecoli --algorithm lr
# -> benchmarks/results/validation/ecoli/lr/

# Re-analyze existing results
./benchmarks/scripts/validation.py compare \
    -d benchmarks/results/validation/ecoli/sr
```

### MD Trajectory

```bash
# Single n_points
./benchmarks/scripts/validation_md.py run \
    --xtc trajectory.xtc --pdb topology.pdb \
    -n my_test --n-points 100

# n_points convergence comparison
./benchmarks/scripts/validation_md.py run \
    --xtc trajectory.xtc --pdb topology.pdb \
    -n my_test --n-points 100,500,960

# Re-analyze existing results
./benchmarks/scripts/validation_md.py compare \
    -d benchmarks/results/validation_md/my_test
```

## Related Documents

- [single-file.md](single-file.md) - Single-file performance benchmarks
- [batch.md](batch.md) - Batch processing benchmarks
- [md.md](md.md) - MD trajectory performance benchmarks
