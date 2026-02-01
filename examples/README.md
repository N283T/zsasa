# Example Files

Sample structure files for testing zsasa.

| File | Format | Structure | Atoms |
|------|--------|-----------|-------|
| `1crn.pdb` | PDB | Crambin | 327 |
| `1crn.cif.gz` | mmCIF (gzip) | Crambin | 327 |
| `1ubq.pdb` | PDB | Ubiquitin | 660 |
| `1ubq.cif` | mmCIF | Ubiquitin | 660 |
| `3hhb.cif.gz` | mmCIF (gzip) | Hemoglobin | 4,779 |

## Usage

```bash
# Build
zig build -Doptimize=ReleaseFast

# Run on example files
./zig-out/bin/zsasa examples/1crn.pdb
./zig-out/bin/zsasa examples/1ubq.cif
./zig-out/bin/zsasa examples/3hhb.cif.gz

# With options
./zig-out/bin/zsasa examples/1crn.pdb --algorithm=lr --n-threads=4
./zig-out/bin/zsasa examples/1ubq.cif --probe-radius=1.2
```
