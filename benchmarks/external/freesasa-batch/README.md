# FreeSASA Batch Wrapper

C++ wrapper for FreeSASA that processes all PDB files in a directory sequentially.

Originally from [RustSASA benchmark](https://github.com/OWissett/rustsasa).

## Build

```bash
c++ -O3 -std=c++17 \
  -I../freesasa-bench/src \
  -o sasa_batch sasa_batch.cpp \
  ../freesasa-bench/src/libfreesasa.a
```

## Usage

```bash
./sasa_batch <input_dir> <output_dir> [n_points]
```

- Processes all `.pdb`/`.cif` files in `input_dir`
- Writes total SASA per file to `output_dir`
- `n_points`: Shrake-Rupley sphere points (default: 100)
