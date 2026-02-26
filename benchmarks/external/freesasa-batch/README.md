# FreeSASA Batch Wrapper

C++ wrappers for FreeSASA batch SASA computation.

Originally from [RustSASA benchmark](https://github.com/OWissett/rustsasa).

## Batch Binary (recommended)

The multi-threaded batch binary lives in [freesasa-bench](https://github.com/N283T/freesasa-bench) at `src/freesasa_batch.cc`.

### Build

```bash
cd ../freesasa-bench/src
c++ -O3 -std=c++17 -I. \
  -o freesasa_batch freesasa_batch.cc json_input.o \
  libfreesasa.a -lz -lpthread
```

### Usage

```bash
./freesasa_batch <input_dir> <output_dir> --n-threads=N --n-points=N
```

- Processes `.json`, `.json.gz`, `.pdb`, and `.cif` files
- `--n-threads`: Number of parallel file-processing threads (default: 1)
- `--n-points`: Shrake-Rupley sphere points per atom (default: 100)

## Legacy PDB-only Wrapper (deprecated)

The original `sasa_batch.cpp` is a single-threaded PDB-only wrapper.

### Build

```bash
c++ -O3 -std=c++17 \
  -I../freesasa-bench/src \
  -o sasa_batch sasa_batch.cpp \
  ../freesasa-bench/src/libfreesasa.a
```

### Usage

```bash
./sasa_batch <input_dir> <output_dir> [n_points]
```
