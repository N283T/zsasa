# FreeSASA Batch Wrapper

C++ wrapper for FreeSASA that processes all PDB files in a directory sequentially.

Originally from [RustSASA benchmark](https://github.com/OWissett/rustsasa).

## Build

```bash
c++ -O3 -o sasa_batch sasa_batch.cpp -lfreesasa
```

Requires FreeSASA library to be installed.

## Usage

```bash
./sasa_batch <input_dir> <output_dir>
```

Processes all `.pdb` files in `input_dir` and writes results to `output_dir`.
