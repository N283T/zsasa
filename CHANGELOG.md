# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

- **Parse-once SASA selection maps**: allow ordinary workflow chain maps to request multiple globally identified chain selections per structure, reuse one parsed/classified input and duplicate chain-set calculations, emit per-selection JSONL success/error rows, and optionally include stable source-indexed atom identity metadata.
- **Multi-interface BSA workflows**: allow CSV and JSON interface maps to request multiple stable-ID interfaces per structure while reusing one parsed/classified input, emitting per-interface success/error JSONL rows, auditable residue SASA components, and opt-in atom-level detail. (#408)

### Changed

- **HETATM records are excluded by default at every entry point (breaking)**: the CLI no longer includes HETATM automatically when the CCD classifier is active, and the gemmi, BioPython and Biotite integrations now default to `include_hetatm=False`. Previously the default result depended on the entry point: `zsasa calc`/`batch` with the default classifier counted ligands, ions and waters, while `process_directory()` did not (1ubq: 5656.65 vs 4834.72 Å²). Modified residues recorded as HETATM inside a chain (for example `MSE`) are excluded too. Pass `--include-hetatm`, `include_hetatm = true` (workflow) or `include_hetatm=True` (Python) to get the earlier CLI behavior. (#431)
- **Selection-map tail scheduling**: claim ordinary multi-selection chain-map files in deterministic longest-processing-time-first order using selection count, while preserving generic batch, BSA, and legacy one-row map ordering.
- **Faster test suite**: `zig build test` runs each test once instead of two or three times (the module and executable test artifacts are filtered to the tests no other artifact reaches, and `scripts/check_test_partition.py`, run in CI, fails when a test would run in no artifact or in more than one), and the batch directory tests in `src/c_api.zig` and `python/tests/test_batch_dir.py` process small temporary inputs with exact reference values instead of running Lee-Richards over the 38-model `test_data/1l2y.pdb`. The `json_writer`, `scanDirectory` and `readAtomInputFromFile` file tests now use temporary directories. (#436)

### Fixed

- **BSA workflow parallelism and progress**: apply workflow thread counts to concurrent structure-file workers, keep partner and complex SASA calculations single-threaded in the multi-file path, serialize concurrent JSONL rows safely, and report completed structure files through the standard terminal progress display.
- **Batch output name collisions**: reject directories whose inputs share an output stem (for example `1crn.pdb` and `1crn.cif.gz`) before processing, instead of reporting every file as successful while later results silently overwrite earlier ones. Applies to per-file output from `zsasa batch`, workflow jobs, and `process_directory(output_dir=...)`; JSONL output is unaffected. The C API reports the condition as `ZSASA_ERROR_OUTPUT_NAME_COLLISION`. (#415)
- **SDF/MOL batch output names**: give every molecule of an SDF or MOL input its own name and per-file output in `zsasa batch`, workflow jobs and `process_directory()`. Molecules of one file that share a title get their 1-based position appended (`dup_ethanol_1`, `dup_ethanol_2`) instead of overwriting each other's output and sharing one JSONL `filename`. The output file name is the molecule name with the extension appended, so dots in the file stem or the title are no longer cut off (`lig.v2.sdf` wrote every molecule to `lig.json`; now `lig.v2_methane.json`). `/`, `\` and control characters in a title are replaced by `_` in the output file name (the title `a/b` in `lig.sdf` failed with `FileNotFound`; now `lig_a_b.json`), so a title cannot write outside the output directory; the JSONL `filename` keeps the title as written. Names change for molecules whose name is shared within their file (the name and the output file) and for stems or titles that contain dots, path separators or control characters (the output file only); all other names are unchanged. The collision check compares the real molecule output names: `lig.sdf` and `lig.mol` are accepted unless their molecules clash, an SDF molecule output that equals another input's output (`lig_1.json` from an unnamed molecule in `lig.sdf` and from `lig_1.pdb`) is rejected, and output names that differ only in ASCII case (`PROT.pdb` and `prot.cif.gz`) are rejected on every platform because they are one file on case-insensitive filesystems. SDF inputs are parsed once per batch, and an SDF file that cannot be read or parsed is reported the same way whatever the thread count: with `--threads=1` it now gets a `read/parse failed` JSONL error row too, where it was only listed in the summary as `SDF parse failed`. (#417)
- **Batch error-path crashes**: report errors raised late in a batch run instead of crashing. A workflow that failed after job setup (for example `use_bitmask = true` with an unsupported `n_points`) freed its job states twice, a single-input batch whose JSONL output could not be created freed its result buffer twice, and a JSONL write failure after parallel processing joined the worker threads a second time. BinaryCIF files whose inline chemical component data failed to load (for example more than 65,535 atoms in one component) also freed decoded columns twice. (#418)
- **CCD-derived radii for components that list hydrogens**: count the hydrogens listed in a component definition together with those implied by unfilled valence when deriving ProtOr radii. Carbons were typed as hydrogen-free whenever the definition listed its hydrogens, as wwPDB CCD entries and typical SDF files do, so sp3 carbons got 1.61 Å instead of 1.88 Å and sp2/aromatic CH carbons got 1.61 Å instead of 1.76 Å. Hydrogen-free aromatic CH carbons written with aromatic bond orders (`AROM`, SDF bond type 4) are corrected the same way. This changes SASA values for ligands and modified residues whose radii come from `--ccd`, inline `_chem_comp_*` data, `--sdf`, or direct SDF input; residues in the built-in table (standard amino acids, nucleotides) are unaffected. (#421)
- **Multi-member gzip input**: decompress every member of a `.gz` input and concatenate the output, as RFC 1952 requires, instead of silently stopping at the end of the first member. Concatenated gzip files (`cat a.gz b.gz`) and `bgzip`/BGZF files previously yielded a truncated structure with exit status 0. Each member's CRC32 and size are verified, the decompressed-size limit applies to the total, and trailing bytes that are not a gzip member (including zero padding) are now rejected as `GzipReadFailed` instead of being ignored. (#426)
- **Trajectory `--algorithm=lr` at the default precision**: `zsasa traj --algorithm=lr` now runs Lee-Richards at `--precision=f32` (the `traj` default), in both the single-threaded and the multi-threaded path, and `--n-slices` takes effect there. Earlier versions silently computed Shrake-Rupley unless `--precision=f64` was also given, so results of `--algorithm=lr` at the default precision change. (#423)
- **Trajectory `--no-hydrogens`**: remove hydrogens from the topology and from every frame, so the option works on trajectories that contain hydrogens instead of failing with an atom count mismatch. The trajectory is still checked against the whole topology, so a hydrogen-free trajectory needs a hydrogen-free topology. (#423)
- **Trajectory topologies with HETATM records or several models**: read every ATOM and HETATM record of the first model of the topology, in file order, instead of dropping HETATM atoms and concatenating all models. (#423)
- **Trajectory output safety**: validate options before creating the output file, so an invalid option no longer truncates an existing results file, and write the frames completed before a read or calculation failure instead of leaving an empty file (the command still exits non-zero). (#423)
- **Trajectory argument handling**: accept `-o FILE`, `-o=FILE`, `--output=FILE` and `--output FILE`, and reject other dash arguments such as `-oFILE` as unknown options (`-o=out.csv` used to write a file named `=out.csv`); reject `--stride=0` instead of treating it as 1; range-check `--probe-radius`, `--n-points` and `--n-slices` like `calc`; drop the positional `[output]` from the documented `traj` synopsis, which the parser never accepted. (#423)
- **NACCESS/OONS atom names read as metals**: atoms outside the NACCESS and OONS tables now take their element from the input's element column instead of from a guess on the atom name, which read any name starting with a two-letter element symbol as that element: gamma hydrogens (`HG`, `HG2`, `HG21`, ...) as mercury (1.55 Å instead of 1.10 Å), HEM `NA` as sodium, HEM `CAA`-`CAD` as calcium, ATP `PB` as lead, and `CD`/`CD1`/`CD2` in unlisted residues as cadmium. **This changes SASA values for NACCESS and OONS runs that include hydrogens or ligands**, including the default `zsasa traj` configuration (for example, model 1 of 1L2Y with hydrogens goes from 1866.44 Å² to 1840.88 Å² with NACCESS); runs over standard amino acid and nucleotide heavy atoms are unchanged. The CLI summary now counts these atoms as "fallback" instead of "classified". Where no element is available (JSON input without `element` under any classifier, Python `classify_atoms` and `get_radius` with NACCESS or OONS), the guess now reads names starting with H, C, N, O, P or S as that element and recognizes an ion by a residue name equal to its atom name; the Python structure-library integrations use the element they already have. (#422)
- **PDB files without an element column**: the element inferred from the atom name now follows the PDB column rule instead of taking the first letter of the name. A name that starts in column 13 begins with a two-letter element (`FE`, `ZN`, `MG`, `MN`, `CL`, `BR`, `SE`, and ions such as `NA`, `CA`, `CU`, `CD`, `HG`), a name that starts in column 14 and a four-character hydrogen name (`HG21`) are one-letter elements, and an ion is also recognized by a residue name equal to its atom name. Previously `ZN`, `MG` and `MN` were unknown elements (1.00 Å), heme `FE` was fluorine, `NA`, `CL`, `CU`, `CD` and calcium ions were nitrogen or carbon, and a mercury atom named `HG` was removed as a hydrogen. Files that left-justify or center their atom names are detected and read by name alone. Text in columns 77-78 that is not an element symbol (the ID code and line number of files written before the element column existed) no longer suppresses the inference, which left every element unknown and kept hydrogens in the calculation. **This changes SASA values for such files with every classifier** when they contain ions, ligands or, in the last case, hydrogens; files with a valid element column are unaffected. (#422)
- **mmCIF/BinaryCIF non-polymer residue numbers**: number waters, ligands and glycans by `auth_seq_id` when `label_seq_id` is null, instead of giving all of them residue number 0. Per-residue rows for waters and ligands change: each residue now has its own row (with `--include-hetatm`, 58 `HOH` rows instead of one `HOH 0` row for `1ubq.cif`), while polymer rows keep their `label_seq_id` numbers by default. Atom counts and total SASA can change for mmCIF and BinaryCIF inputs with alternate-location waters or ligands, because altloc selection no longer treats all waters of a chain as one residue and drops those that have an altloc ID. **Behavior change for `--auth-chain`** (and workflow maps with `asym_id_type = auth`): residue numbers reported with auth chain IDs change from `label_seq_id` to `auth_seq_id` for every residue, so chain ID, residue number and insertion code match the PDB-format file. Polymer residue numbers in `--auth-chain` output therefore change wherever the two numberings differ. (#424)
- **Thread limit during a multi-threaded calculation**: a SASA calculation no longer crashes when the process cannot start all of its worker threads, for example at the per-process thread limit or in a container whose `pids.max` is below the CPU count. The thread pool returned the spawn error without waiting for the workers it had already started, and those kept reading and writing memory the caller had freed; through the Python API this surfaced as `RuntimeError` followed by a segmentation fault. The calculation now finishes on the threads that did start plus the calling thread, and returns the same result. Lee-Richards also reports a failed scratch-buffer allocation as out of memory, instead of leaking the per-atom array (single-threaded) or returning success with part of that array uncomputed (multi-threaded). (#429)
- **Malformed SDF and BinaryCIF input**: report a parse error instead of corrupting memory, exhausting the stack, or allocating memory for lengths that cannot be right. An SDF V3000 molecule whose atom or bond block has more or fewer entries than its `M  V30 COUNTS` line declares is rejected as `InvalidV3000`; extra entries used to be written past the end of the allocated lists, which crashed or returned a wrong atom count, and missing entries were accepted. BinaryCIF input nested deeper than 32 MessagePack levels (real files nest 12) is rejected as `InvalidMessagePack` instead of overflowing the stack. A BinaryCIF `RunLength` or `IntegerPacking` `srcSize` that does not agree with the category's `rowCount` is rejected before the column is decoded, and columns that zsasa does not read are no longer decoded, so a file of about 1 kB can no longer make the parser allocate hundreds of megabytes. As a result, a malformed column that zsasa does not read no longer fails the parse, and an encoding chain with `RunLength` below `IntegerPacking`, which no writer produces, is rejected as `UnsupportedEncoding`. (#427)

### Documentation

- **Documentation brought back in line with the CLI and Python API**: the CLI reference marks which of `calc`, `batch` and `traj` accepts each option and lists the previously undocumented ones (`--sdf`, `--mol`, `--input-io`, `--af-model-fast`, `--profile-stages`, `--exclude-hydrogens`), the error-message table matches the binary, and `--help` lists every supported input extension and describes the structure CSV columns; the workflow guide has a reference for every key and for legacy flat-root manifests; the Python docs match the signatures (`classifier_type` in the classifier utility table, `chunk_size` and `store_atom_areas` for MDAnalysis, `compute_sasa_trajectory_summary` and `zsasa.dcd`, `ClassifierType.PROTOR` described as its own value rather than an alias of CCD, default CCD radii in the examples); and the batch, trajectory and classifier guides, `CONTRIBUTING.md` and `README.md` were corrected. (#438)

## [0.9.1] - 2026-07-28

### Added

- **Per-file workflow chain maps**: add CSV and JSON maps for file-specific single- or multi-chain SASA selections and BSA/ΔSASA partner groups, with per-entry `label` or `auth` asym ID selection for mmCIF and BinaryCIF inputs. (#406)

## [0.9.0] - 2026-07-15

### Added

- **Batch AF model fast parser**: add experimental `zsasa batch --af-model-fast` for AlphaFold-like mmCIF batches. The fast path preserves multiple chains and atom metadata, and falls back to the generic parser for unsupported layouts.
- **Batch input I/O selection**: add `--input-io=auto|mmap|read` to compare or select file input strategies where supported.
- **Batch stage profiling**: add `--profile-stages` with `--timing` to report aggregate read/parse, classifier, and JSONL write times.

### Changed

- **Batch JSONL streaming**: buffer normal-sized JSONL records while streaming large records directly to avoid regressions for PDB batches with large atom-area payloads.

## [0.8.0] - 2026-07-01

### Added

- **AltLoc parsing modes**: add shared `--altloc` policy support for mmCIF and BinaryCIF inputs across `calc`, `batch`, and `traj`, including `auto`, `none`, `all`, `highest-occupancy`, and single-ID selection. (#394)
- **Workflow JSONL output controls**: add workflow `[output.jsonl]` controls for atom areas, total areas, decimal rounding, and metadata sidecars. (#396)

### Changed

- **Batch JSONL output**: add explicit success/error status rows, configurable decimal rounding, and stricter serialization/write error handling for JSONL output. (#395)
- **mmCIF parser performance**: reduce atom_site parsing overhead by batching tokenizer position updates, skipping unnecessary AltLoc selection, and lazily backfilling extended chain IDs. (#393)
- **ProtOr classifier**: make `--classifier=protor` use static ProtOr-compatible radii without loading inline, external, or SDF-derived CCD resources. This provides a faster protein-only path for mmCIF/PDB inputs such as AFDB models. (#392)

### Fixed

- **Inline CCD parsing**: gate inline CCD extraction so batch and trajectory topology parsing skip unnecessary CCD work, while `batch --classifier=ccd` still uses per-file inline CCD data when present. (#392)

## [0.7.1] - 2026-06-29

### Fixed

- **RSA text compatibility warning**: warn when `--format=rsa` output exceeds legacy NACCESS fixed-width columns, and document JSON/JSONL as the recommended machine-readable formats.

## [0.7.0] - 2026-06-28

### Added

- **FreeSASA-compatible outputs**: add FreeSASA-style `calc` output formats for interoperability with existing SASA workflows. (#382)
- **BSA workflow analysis**: add batch workflow support for two-partner buried surface area / ΔSASA analysis. (#388)
- **Trajectory bitmask LUT modes**: add configurable trajectory bitmask lookup-table reuse modes. (#387)

### Changed

- **Benchmark website**: refresh benchmark documentation and website assets. (#381)
- **Batch scheduling**: allow explicit batch thread overcommit for I/O-bound file sets. (#384)
- **Bitmask SR accuracy controls**: add experimental bitmask bias correction controls. (#386)

### Fixed

- **Workflow parsing performance**: optimize multi-chain workflow parsing. (#385)
- **Review hardening backlog**: harden C ABI test coverage, input validation, parser correctness, Python reliability, release workflows, source-install checks, and docs/package metadata. (#389)

## [0.6.0] - 2026-05-21

### Added

- **ztraj-backed trajectory formats**: support TRR and AMBER NetCDF in `zsasa traj` alongside XTC and DCD, with ztraj readers normalizing coordinates to Å before SASA calculation. (#378)

### Changed

- **Trajectory reader backend**: use `ztraj` for trajectory readers, centralizing XTC/DCD handling and dependency notices. (#377)
- **Python trajectory guidance**: recommend `pyztraj` for direct Python trajectory-file I/O while keeping `zsasa.xtc` and `zsasa.dcd` as compatibility APIs. (#379)

### Fixed

- **Batch JSONL and FFI concurrency coverage**: add regression coverage for batch JSONL output and concurrent FFI calls. (#376)

## [0.5.0] - 2026-05-20

### Added

- **Adaptive batch bitmask SR**: add experimental `zsasa batch --adaptive-sr` mode with coarse/fine point controls for two-stage bitmask Shrake-Rupley runs. (#371)
- **BinaryCIF input support**: accept `.bcif`, `.bcif.gz`, and `.bcif.zst` inputs in `calc` and `batch` by decoding `_atom_site` directly in Zig. (#372)

### Changed

- **Website documentation**: restructure the documentation site around task-oriented guides, CLI references, Python APIs, and integrations. (#370)
- **Workflow batch execution**: reuse parsed/classified structures across eligible workflow jobs so named chain analyses such as chain A, chain B, and complex AB avoid repeated parsing while preserving output schemas and compatibility fallbacks. (#374)

### Fixed

- **BinaryCIF CCD classification**: use inline CCD components from BinaryCIF inputs when available. (#373)

## [0.4.0] - 2026-05-19

### Added

- **Workflow files**: add TOML workflow support via `zsasa calc --workflow` and `zsasa batch --workflow`; `batch --manifest` remains a compatibility alias. (#366)

### Changed

- **Breaking**: custom classifier configs are TOML-only; legacy FreeSASA-style custom classifier files are no longer supported. (#366)
- **Breaking/API**: the public `batch_manifest` module was removed and replaced by `workflow_manifest`. (#366)
- **Repository automation**: add agent instructions and ignore local planning docs. (#367)

### Fixed

- **CI**: ignore gitignore-only changes in automation checks. (#368)

## [0.3.2] - 2026-05-17

### Added

- **Batch JSONL residue maps**: add opt-in `--residue-map` / `residue_map = true` output with compact columnar residue identifiers, atom ranges, and residue SASA for large-scale JSONL workflows. (#364)

## [0.3.1] - 2026-05-17

### Added

- **Batch TOML manifests**: add the original named multi-job workflow file support for chain A, chain B, and AB complex SASA; includes batch `--chain` / `--auth-chain` support, per-job JSONL and per-file output layout, and CLI/docs coverage. (#361)

## [0.3.0] - 2026-05-17

### Added

- **Zstandard decompression for structure inputs**: `.json.zst`, `.pdb.zst`, `.cif.zst`, `.mmcif.zst`, `.ent.zst`, `.sdf.zst`, and `.mol.zst` are now detected and transparently decompressed via native `std.compress.zstd`. (#357)

### Changed

- **CLI progress bars** now use Zig's standard progress API for `batch` and `traj`, reducing custom rendering logic while preserving progress reporting. (#358)
- **Project branding documentation**: add the project logo to README and refresh logo SVG assets. (#355, #356)

## [0.2.11] - 2026-04-26

### Fixed

- **Build & Publish workflow**: Bump remaining Zig **0.15.2 → 0.16.0** pins missed by PR1 (#345). The `Dockerfile` and the `cibuildwheel` `before-all` in `python/pyproject.toml` were both pulling the old toolchain and rejecting 0.16-only API; the result is that wheels and the Docker image have not been published since v0.2.8. Also updates `install.sh`, `python/hatch_build.py`, `python/zsasa/_ffi.py` error messages, `CONTRIBUTING.md`, the website docs, and the bug-report ISSUE template. (#353)
- **Windows wheel test**: `examples/1ubq.cif` was failing with `StreamTooLong` because the Windows `mmapFile` path used `.limited64(stat.size)` and Windows occasionally reports a byte count divergent from the buffered reader. Switch to `.unlimited` — `mmapFile`'s callers already trust the input path. (#353)
- **aarch64 wheel build**: `cibuildwheel`'s short alias `manylinux_2_28` resolved to the retired quay.io tag `2026.03.01-1`. Pin both `manylinux-x86_64-image` and `manylinux-aarch64-image` to an explicit current tag (`2026.04.25-0`). (#353)
- **`flake.nix` derivation version**: Catch up to `0.2.11` (was stale at `0.2.4` until v0.2.10).

## [0.2.10] - 2026-04-26

### Changed

- **Gzip decompression: C-zlib → native `std.compress.flate`**. Reverts the workaround introduced in #320 now that the upstream panic ([ziglang/zig#25035](https://github.com/ziglang/zig/issues/25035)) is fixed in Zig 0.16. Public API (`readGzip`, `readGzipLimited`, `DEFAULT_MAX_SIZE`, `GzipError`) and the 4 GB decompression-bomb cap unchanged; callers untouched. (#351)

### Added

- **Explicit gzip CRC32 + ISIZE trailer verification** in `gzip.zig`. `std.compress.flate` parses the gzip trailer fields but does not compare them against the decompressed bytes; without this check a corrupt `.cif.gz` could silently decompress to wrong bytes. Adds a regression test (`readGzip rejects gzip with corrupted CRC`). (#351)

### Removed

- **`zlib` dependency** — removed from `build.zig.zon`, `build.zig` (`b.dependency`, `linkLibrary`, `b.addTranslateC`), and `src/c/zlib_wrapper.h`. The shared library no longer links zlib. (#351)

## [0.2.9] - 2026-04-26

### Changed

- **Zig 0.16.0 migration**: minimum Zig version bumped to 0.16.0. Toolchain (build.zig.zon, flake.nix, CI workflows, README badge) updated. (#345)
- **`std.fs.*` → `std.Io.*`**: ~48 file I/O sites migrated to the new `std.Io` interface. "Juicy Main" pattern: a single `std.Io.Threaded` constructed in `main()` is threaded through subcommand handlers; `c_api.zig` constructs per-FFI-call `Io.Threaded` for batch entries (correct concurrency) and uses the global single-threaded variant for single-call entries. C ABI is preserved. (#345)
- **`std.heap.GeneralPurposeAllocator` → `std.heap.DebugAllocator`**, **`std.time.Timer` → `std.Io.Timestamp`**, **`std.Thread.Mutex` → `std.Io.Mutex`** (with `lockUncancelable`), **`std.Thread.sleep` → `std.Io.sleep`**. (#345)
- **`std.Io.Writer.Allocating` / `Writer.fixed` / `Reader.fixed`** replace `ArrayList(u8).writer()` / `fixedBufferStream`. (#345)
- **`mem.indexOf*` → `mem.find*`** rename across 18 sites; **Managed → Unmanaged HashMaps** across 11 files (`StringHashMapUnmanaged.empty` + op-time allocator); **`std.mem.trimStart/trimEnd`** replace deprecated `trimLeft/trimRight`; **`std.posix.PROT`** packed-struct syntax. (#345)
- **`@cImport(@cInclude("zlib.h"))`** moved to `b.addTranslateC` in `build.zig` (uses `src/c/zlib_wrapper.h` because dep lazy paths cannot serve as `root_source_file`). (#345)
- **`@Vector` annotation** for `@sqrt` coercion in `simd.zig` (8 sites — array → vector); `@bitCast(@Vector(N, bool))` → `@bitCast(@as(@Vector(N, u1), @intFromBool(...)))` for cross-platform vector-mask packing (Linux x86_64 fix). (#345)
- **`@floor`/`@ceil`/`@round`/`@trunc` audit**: simplify `@as(usize, @intFromFloat(@ceil(x)))` → `@as(usize, @ceil(x))` where 0.16 supports direct int return. (#345)
- **`std.testing.Smith`** for fuzz tests in `mmcif_parser`, `pdb_parser`, `cif_tokenizer`. (#345)
- **`zxdrfile` v0.1.1 → v0.4.0** (Zig 0.16 compatible release, 2026-04-26). `XtcReader.open` signature gained `io_handle` parameter; both call sites updated. (#345)

### Fixed

- **`Io.Writer.Allocating` errdefer leak** in `json_writer.zig`: `sasaResultToCsv` and `sasaResultToRichCsv` now properly clean up the buffer on write errors. (#345)
- **`mmap_reader.zig` Windows path**: restored `size` cap (`.limited64(size)`) lost during migration; added TOCTOU assertion. (#345)
- **`JsonlStreamWriter.writeResult`** now records write failures via `hasError()` for downstream propagation; `mutex.lock(io) catch unreachable` replaced with `lockUncancelable(io)`. (#345)

## [0.2.8] - 2026-04-13

### Added

- **SDF/MOL file support**: New SDF parser supporting V2000 and V3000 formats. Calculate SASA for small molecules directly from SDF files (`zsasa calc molecule.sdf`) (#338)
- **`--sdf` option**: Provide bond topology for CCD-unregistered compounds (e.g., Boltz-predicted ligands) via `--sdf=ligand.sdf`. Available in `calc`, `batch`, and `traj` subcommands (#338)
- **`--mol` option**: Select a specific molecule from multi-molecule SDF by name or 1-based index (`--mol=water` or `--mol=2`) (#339)
- **Batch SDF expansion**: Multi-molecule SDF files in batch mode are expanded into individual items, each molecule calculated independently (#339)

### Changed

- **ProtOr is now an alias for CCD**: `--classifier=protor` uses the same CCD bond-topology classifier. Default classifier for Python bindings changed from NACCESS to CCD (#335, #337)
- **Trajectory classifier default**: `traj` subcommand defaults to NACCESS for MD trajectories (#337)

### Fixed

- **SDF per-molecule SASA**: Multi-molecule SDF files now calculate each molecule independently instead of combining them into one structure (#339)
- **Python classifier tests**: Updated expected values for CCD default (ALA:O = 1.42 Å) (#338)
- **Memory safety**: Fixed StoredComponent leaks on ComponentDict insertion failure, toAtomInput over-allocation with >26 molecules (#338, #339)

### Documentation

- CCD classifier documentation rewritten with accurate architecture details (#334, #336)

## [0.2.7] - 2026-04-12

### Added

- **CCD classifier** (`--classifier=ccd`): New classifier that derives ProtOr-compatible radii from CCD (Chemical Component Dictionary) bond topology, enabling accurate radius assignment for any chemical component — not just standard amino acids (#326, #327)
- **External CCD dictionary** (`--ccd=<path>`): Load external CCD dictionary for non-standard residues. Supports both CIF text (`.cif`, `.cif.gz`) and binary ZSDC format (#328)
- **`compile-dict` subcommand**: Convert CCD dictionary from CIF text to compact binary ZSDC format for faster loading (`zsasa compile-dict components.cif.gz -o components.zsdc`) (#328)
- **Python**: `ClassifierType.CCD` added to Python bindings (#330)

### Changed

- **CCD classifier auto-includes HETATM**: When using `--classifier=ccd`, HETATM records are included automatically without needing `--include-hetatm` (#329)
- **Python**: Split monolithic `core.py` into focused modules (`_ffi.py`, `sasa.py`, `classifier.py`, `rsa.py`, `batch.py`). Public API unchanged (#332)

## [0.2.6] - 2026-03-22

### Fixed

- **Docker build**: use PIC-enabled zlib for shared library to fix `R_X86_64_32` relocation errors in Docker builds (#323)
- **CI**: vendor zlib from source (allyourcodebase/zlib) instead of linking system library, fixing builds on Windows and manylinux containers (#322)

### Changed

- **Publish workflow**: decouple job dependencies to prevent cascading failures — PyPI, GitHub Release, Docker, and package managers now run independently. Added `workflow_dispatch` with selective job re-runs (#324)

## [0.2.5] - 2026-03-22

### Fixed

- **Gzip decompression for mmCIF/PDB files**: `.cif.gz`, `.pdb.gz`, `.mmcif.gz`, `.ent.gz` files are now transparently decompressed. Previously only `.json.gz` was supported, and mmCIF/PDB parsers passed raw gzip data to the tokenizer, causing `NoAtomSiteLoop` errors (#319)

### Changed

- **Switch from Zig native flate to C zlib**: gzip decompression now uses C zlib (`gzopen`/`gzread`/`gzclose`) instead of Zig 0.15's `std.compress.flate`, which panics on certain valid gzip files (e.g. PDB entry 2OXD). This adds `libz` as a build dependency. See [ziglang/zig#25035](https://github.com/ziglang/zig/issues/25035). Will revert to native flate when the upstream bug is fixed (#320)
- **Decompression bomb protection**: `readGzip` enforces a 4 GB max decompressed size limit to prevent memory exhaustion from malicious `.gz` files

## [0.2.4] - 2026-03-11

### Added

- **CLI installer**: `install.sh` script for one-line CLI installation — builds from source with Zig or downloads pre-built binary from GitHub Releases (#307)
- **CLI release binaries**: pre-built CLI binaries for linux-x86_64, linux-aarch64, macos-x86_64, macos-aarch64, windows-x86_64 attached to GitHub Releases (#307)
- **Nix flake**: `nix run github:N283T/zsasa` or `nix profile install github:N283T/zsasa` for Nix-based installation (#311)

### Changed

- **XTC reader**: replaced local `src/xtc.zig` with external [zxdrfile](https://github.com/N283T/zxdrfile) package dependency (#309)
- **CI**: update cibuildwheel v2.22 to v3.4.0, remove obsolete skip selectors (#305, #306)

## [0.2.3] - 2026-03-10

### Fixed

- **Windows build**: fall back to heap allocation on Windows where POSIX mmap is unavailable (#301)
- **Windows build**: replace `std.mem.zeroes` with `undefined` for `JsonlStreamWriter` placeholder — zeroes fails on non-nullable pointers on Windows (#301)
- **Third-party license compliance**: add BSD-2-Clause notice for libxdrfile (via chemfiles/xdrfile) in `src/xtc.zig` and create `THIRD_PARTY_NOTICES.md` (#300)
- Correct Lahuta author name (Besian I. Sejdiu) in docs and README (#299)

### Changed

- Remove legacy `docs/` directory — all documentation migrated to website (#298)
- Add acknowledgments section to README and comparison page (#297)

## [0.2.2] - 2026-03-10

### Added

- **Bitmask variants for Python MD wrappers**: `use_bitmask` option for MDTraj, MDAnalysis, XTC, and DCD integrations (#275)
- **JSONL streaming batch output**: `--format=jsonl` for memory-efficient batch results (#227, #235)
- **Documentation site overhaul**:
  - Comparison page vs FreeSASA, RustSASA, and Lahuta with source code references (#292)
  - Landing page with hero section and feature cards (#291)
  - Python autodoc generation with pdoc (#290)
  - Split CLI reference into Commands, Input, and Output pages (#289)
  - Benchmarks overview page and changelog (#294)
  - Rewrote all benchmark pages with new results (#280, #282, #283, #284)

### Performance

- **mmap file reading**: replaced `readToEndAlloc` with memory-mapped I/O for structure files (#229)
- **Flat buffer NeighborList**: replaced dynamic `ArrayList` with pre-allocated flat buffers (#230)
- **Trajectory parallel workers**: aligned with batch allocator pattern for lower overhead (#231)
- **64KB write buffer** for JSONL output (#228)

### Fixed

- Corrected Lahuta metadata (URL, language) and RustSASA precision (f64 → f32) in benchmark docs (#293)
- Various PDB generation fixes: chain name shortening, serial number wrapping, CRYST1 Z value handling (#251, #252, #253)

## [0.2.1] - 2026-02-25

### Changed

- Relaxed bitmask `n_points` constraint from fixed 64/128/256 to any value 1..1024, with internal storage expanded from `[4]u64` to `[16]u64` (#210)
- Updated Python bindings, C API, CLI help text, and documentation to reflect new n_points range (#210)

### Added

- `--use-bitmask` and `--n-points` flags to benchmark scripts for bitmask LUT benchmarking (#208)

## [0.2.0] - 2026-02-25

### Added

- **Bitmask-optimized Shrake-Rupley algorithm** (`--use-bitmask`): precomputed occlusion bitmask LUT with O(1) octahedral encoding for direction lookup, replacing per-point neighbor testing (#197)
- **SIMD 4-neighbor batching** for bitmask SR: processes 4 neighbors simultaneously with `@Vector(4, T)`, branchless octahedral encoding via `@select`, and combined mask accumulation (#198)
- **Batch-mode `--use-bitmask` support**: shared LUT across all files in batch processing, avoiding redundant ~20ms LUT construction per file (#198)
- **Trajectory-mode `--use-bitmask` support**: build LUT once, reuse across all frames in both sequential and batch-parallel paths
- **Python bindings `use_bitmask` parameter**: added `use_bitmask=True` option to `calculate_sasa()` and `calculate_sasa_batch()`, with pass-through to MDTraj, MDAnalysis, XTC, and DCD integrations
- **C API bitmask exports**: `zsasa_calc_sr_bitmask`, `zsasa_calc_sr_batch_bitmask`, `zsasa_calc_sr_batch_bitmask_f32` for FFI access to bitmask LUT optimization

### Changed

- **BREAKING**: CLI now requires subcommands: `zsasa calc`, `zsasa batch`, `zsasa traj`
- **BREAKING**: Removed `--parallelism` option (`calc` uses atom-level, `batch` uses file-level parallelism)
- **Trajectory mode**: `--include-hydrogens` is now the default (hydrogen atoms included). Use `--no-hydrogens` to exclude. MD trajectories typically include all atoms.
- CI: removed Windows from PR checks (linux + macOS only); Windows builds remain in release workflow (#199)

### Removed

- Pipeline parallelism mode (`--parallelism=pipeline`)
- Atom-level batch parallelism (`--parallelism=atom`)

### Performance

- E.coli proteome batch (f32, 10 threads, 128 points): **4.92s → 2.57s** with bitmask SR (1.9x speedup, within 12% of lahuta reference)

## [0.1.3] - 2026-02-25

### Added

- **Directory batch processing C API**: `zsasa_batch_dir_*` functions for processing all structure files in a directory (#191)
- **Python bindings for directory batch processing**: `process_directory()` function and `BatchDirResult` dataclass wrapping the C API (#193)
  - Process all supported structure files in a directory from Python
  - Support for SR/LR algorithms, all classifiers, threading, output directory
  - Per-file results: filename, atom count, total SASA, status
  - Error mapping: `ValueError`, `FileNotFoundError`, `MemoryError`, `RuntimeError`
- Docusaurus documentation site with GitHub Pages deployment (#186)

### Changed

- Slimmed down README, added uv install instructions (#190)
- CI: skip workflow for website-only changes (#189)

### Fixed

- **json_writer performance regression**: Reverted from unbuffered streaming writes to in-memory string building + single writeAll, fixing ~8x slowdown in batch processing caused by millions of write syscalls (#157 regression, #194)
- Updated benchmark docs for new dataset (#188)
- Fixed absolute paths for benchmark images (#187)

### Removed

- **Streaming output** (`--stream`, `--stream-format`, `--stream-output`): Removed StreamWriter module and CLI options to reduce code complexity (#194)

## [0.1.2] - 2026-02-22

### Added

- JSON streaming output for batch processing (`--stream`): stream results as NDJSON or JSON array as each file completes (#157)
- Stream format selection (`--stream-format`): choose between `ndjson` (default) and `json` array format (#157)
- Stream output destination (`--stream-output`): write stream to file instead of stdout (#157)
- Writer-based streaming JSON output for per-file results, reducing memory usage (#157)
- Zig package manager (zon) distribution: zsasa can now be used as a library dependency via `zig fetch` (#160)
- Public library API in `root.zig`: `shrake_rupley`, `lee_richards`, `types`, `pdb_parser`, `mmcif_parser`, `json_parser`, `classifier`, `analysis`
- Fuzz tests for CIF tokenizer, PDB parser, and mmCIF parser using Zig's built-in `std.testing.fuzz()` (#161)
- TOML format support for custom classifier configs (`--config=file.toml`): human-friendly alternative to the legacy classifier text syntax with auto-detection by file extension (#158)
- **DCD trajectory reader** (native Zig): read NAMD/CHARMM DCD binary trajectories without external dependencies (#154)
  - Zig DCD reader (`src/dcd.zig`) with endianness auto-detection and CHARMM extension support
  - C API: `zsasa_dcd_open`, `zsasa_dcd_close`, `zsasa_dcd_read_frame`, `zsasa_dcd_get_natoms`
  - Python DCD reader (`zsasa.dcd`): `DcdReader` class and `compute_sasa_trajectory()` function
  - CLI: `zsasa traj` now supports `.dcd` files with auto-detection by file extension
- Zig library API reference documentation (`docs/zig-api/`): types, algorithms, parsers, classifier, analysis (#156)
- Python pdoc auto-generated API documentation (`scripts/generate-python-docs.sh`) (#156)
- Documentation site with MkDocs Material and GitHub Pages deployment (#184)
- `zig build docs` step for interactive Zig autodoc generation (#184)

### Changed

- `build.zig`: Removed boilerplate template comments (176 → 63 lines)
- `build.zig.zon`: Synced version with `build.zig`
- Homepage and Documentation URLs now point to GitHub Pages (#184)

### Fixed

- `simd.zig`: Fixed `std.math.atan2` comptime_float errors with explicit `@as(f64, ...)` casts
- `simd.zig`: Corrected `fastAtan2` test tolerance for negative quadrant inputs (0.005 → 0.07)
- `root.zig`: Re-enabled `shrake_rupley` and `lee_richards` in test block (previously excluded due to transitive simd test failure)

## [0.1.1] - 2026-02-22

### Added

- `CODE_OF_CONDUCT.md` (Contributor Covenant v2.1)
- `CITATION.cff` (CFF 1.2.0 format for academic citation)
- GitHub issue templates (bug report, feature request) in YAML form
- Pull request template
- `speedup_by_threads.png` plot generation in `analyze.py large` (thread scaling for 50k+ atoms)
- CI status, license, Zig, and Python badges to READMEs
- Pre-built wheel distribution via cibuildwheel (Linux x86_64/aarch64, macOS x86_64/arm64, Windows x86_64)
- Windows support: Zig build, tests, CLI, Python bindings
- PyPI publish workflow (`.github/workflows/publish.yml`) with OIDC trusted publishing
- `python -m ziglang` fallback in `hatch_build.py` for Zig discovery

### Changed

- `python/pyproject.toml`: Updated author name, added Documentation/Issues/Changelog URLs

### Fixed

- README: Corrected MD trajectory benchmark data (6sup_A_analysis: 4.3x speedup, verified against actual data)
- README: Fixed broken documentation link (`docs/python.md` → `docs/python-api/`)
- README: Fixed broken image reference for thread scaling plot

## [0.1.0] - 2026-01-31

### Added

- **PDB file format support**
  - Fixed-width PDB parser (`src/pdb_parser.zig`)
  - Auto-detection of `.pdb` and `.ent` files
  - MODEL/ENDMDL, alternate location, chain filtering
  - Element inference from atom names

- **mmCIF parser** (internal)
  - Native mmCIF parser (`src/mmcif_parser.zig`)
  - No external dependencies (replaces gemmi for CLI)

- **Gzip support**
  - Transparent decompression for `.json.gz` and `.cif.gz`
  - Streaming decompression with zlib

- **Batch processing** (directory input)
  - Process entire directories: `zsasa ./input_dir/ ./output_dir/`
  - File-level parallelism with work stealing
  - Per-thread arena allocators for memory efficiency
  - Progress bar with file count
  - `--parallelism` option for concurrent file processing
  - Duplicate coordinate detection and warning

- **f32 precision option** (`--precision=f32`)
  - Single-precision mode for reduced memory usage
  - Comptime generics for zero-cost abstraction
  - ~Same speed, slightly lower accuracy

- **AVX-512 auto-optimization**
  - 16-wide SIMD on supported CPUs
  - Automatic detection and fallback (16 → 8 → 4 → scalar)

- **Benchmark infrastructure**
  - `benchmarks/scripts/run.py` - Unified benchmark runner
  - `benchmarks/scripts/analyze.py` - Results analysis and plotting
  - `benchmarks/scripts/sample.py` - Stratified sampling for large datasets
  - `benchmarks/scripts/build_index.py` - Dataset indexing
  - Full PDB dataset benchmark (238,124 structures)
  - Batch benchmark: Zig +7% faster than RustSASA

- **Example files** (`examples/`)
  - Sample PDB and mmCIF files for quick testing
  - README with usage examples

- **8-wide SIMD optimization** for Shrake-Rupley algorithm
  - `@Vector(8, f64)` for processing 8 atoms in parallel
  - Tiered processing: 8-wide → 4-wide → scalar for remaining atoms
  - ~16% speedup on large structures (4V6X: 237k atoms)

- **Fast trigonometry** for Lee-Richards algorithm
  - Polynomial approximations for `acos` and `atan2`
  - ~37% speedup on large structures (4V6X: 1021ms → 743ms)
  - Accuracy within 0.3% of reference (well within 2% tolerance)

- **Area difference column** in benchmark output
  - Shows percentage difference between Zig and FreeSASA C results

- **Analysis options** (CLI)
  - `--per-residue` - Per-residue SASA aggregation
  - `--rsa` - Relative Solvent Accessibility calculation
  - `--polar` - Polar/nonpolar SASA classification

- **Python bindings** (`python/zsasa`)
  - C ABI shared library (`libzsasa.dylib/.so/.dll`)
  - NumPy-based Python API with ctypes bindings
  - Both SR and LR algorithms supported
  - `calculate_sasa(coords, radii, algorithm="sr"|"lr", ...)` function
  - RSA functions: `calculate_rsa()`, `calculate_rsa_batch()`, `get_max_sasa()`
  - Per-residue aggregation: `aggregate_by_residue()`, `ResidueResult` class
  - 161 unit tests with pytest

- **Python integrations** (optional dependencies)
  - Gemmi integration for mmCIF/PDB file loading
  - BioPython integration for structure file support
  - Biotite integration for structure analysis workflows
  - MDTraj integration (`zsasa.mdtraj`) - drop-in replacement for `mdtraj.shrake_rupley()`
  - MDAnalysis integration (`zsasa.mdanalysis`) - `SASAAnalysis` class compatible with `AnalysisBase`

- **XTC trajectory reader** (native Zig)
  - Zig port of GROMACS libxdrfile (BSD-2-Clause)
  - Python XTC reader (`zsasa.xtc`) - no MDTraj/MDAnalysis dependency required

- **Trajectory subcommand** (`zsasa traj`)
  - `zsasa traj trajectory.xtc topology.pdb` - CLI trajectory mode
  - Frame-level batch parallelism with work-stealing
  - `--stride=N`, `--start=N`, `--end=N` frame filtering
  - `--batch-size=N` for controlling parallel batch size
  - Default f32 precision for speed in trajectory mode

- **Hydrogen and HETATM filtering**
  - `--include-hydrogens` - Include H/D atoms (default: excluded)
  - `--include-hetatm` - Include HETATM records (default: excluded)
  - Applied to both PDB and mmCIF parsers

- **SASA validation infrastructure**
  - `benchmarks/scripts/validation.py` - Accuracy validation vs FreeSASA C
  - `benchmarks/scripts/validation_md.py` - MD trajectory validation across implementations
  - Lee-Richards validation with E. coli proteome (R²=1.0)

- **MD trajectory benchmarks**
  - `benchmarks/scripts/bench_md.py` - Hyperfine-based MD benchmark
  - `benchmarks/scripts/analyze_md.py` - MD analysis and plots
  - E. coli proteome batch benchmark

- **Timing breakdown** (`--timing` flag)
  - Reports detailed timing for each phase: parsing, classification, SASA calculation, output
  - Enables fair performance comparison by measuring SASA-only time

- **Benchmark dataset** (6 structures from tiny to xlarge)
  - 1CRN (327 atoms), 1UBQ (602), 1A0Q (3,183), 3HHB (4,384), 1AON (58,674), 4V6X (237,685)
  - `benchmarks/inputs_protor/` - Pre-generated inputs with ProtOr radii
  - `scripts/data/generate_protor.py` - Generate inputs with ProtOr radii
  - `scripts/benchmark.py` - Unified benchmark comparing Zig vs FreeSASA

- **Lee-Richards algorithm** (`--algorithm=lr`)
  - Slice-based method with exact arc integration
  - `--n-slices=N` option (default: 20)
  - Multi-threading and SIMD support
  - 1.1x-1.7x faster than FreeSASA C

- **Atom classifier module** with CLI integration
  - `classifier.zig` - Core data structures, element-based radius guessing, ClassifierType enum
  - `classifier_naccess.zig` - NACCESS-compatible built-in classifier
  - `classifier_protor.zig` - ProtOr classifier (hybridization-based, Tsai et al. 1999)
  - `classifier_oons.zig` - OONS classifier (older FreeSASA default)
  - O(1) compile-time hash lookup using `StaticStringMap`
  - Support for 20 standard amino acids + SEC/MSE/PYL/ASX/GLX
  - Support for RNA/DNA nucleotides (A, C, G, I, T, U, DA, DC, DG, DI, DT, DU)
  - ANY fallback for backbone atoms (NACCESS/OONS)
  - Element-based radius guessing from atom names
  - `--classifier=naccess|protor|oons` - Use built-in classifier
  - `--config=FILE` - Use custom classifier config file

- Extended input format with optional `residue` and `atom_name` fields

### Changed

- **Updated benchmark**: Now compares against FreeSASA C (native binary) instead of Python
  - SR: 1.2x-2.3x faster than FreeSASA C
  - LR: 1.1x-1.7x faster than FreeSASA C
- **Scripts reorganization**
  - Created `scripts/data/` subdirectory for data preparation scripts
  - Renamed: `benchmark_all.py` → `benchmark.py`, `validate_accuracy.py` → `validate.py`
  - Moved data scripts to `scripts/data/` with shorter names
- **Project renamed**: `freesasa-zig` → `zsasa` (repository, binary, Python package)
- **Default ProtOr classifier** for PDB/mmCIF input (no `--classifier` flag needed)
- **Internal**: Refactored `AtomInput.r` from `[]const f64` to `[]f64` to properly support classifier mutations
- **Internal**: `FixedString5` for mmCIF 5-character `comp_id` (modified residue support)
- Added LICENSE (MIT) and CONTRIBUTING.md
- Enabled GitHub Actions CI/CD (format, build, test, Python)

## [0.0.5] - 2025-01-23

### Added

- **Extended CLI options**
  - `--probe-radius=R` - Configure probe radius (default: 1.4 Å)
  - `--n-points=N` - Configure test points per atom (default: 100)
  - `--quiet` / `-q` - Suppress progress output
  - `--help` / `-h` - Show help message
  - `--version` / `-V` - Show version

- **Output format options**
  - `--format=json` - Pretty-printed JSON (default)
  - `--format=compact` - Single-line JSON
  - `--format=csv` - CSV format with header

- **Input validation**
  - `--validate` - Validate input without calculation
  - Array length consistency check
  - Coordinate finiteness check (NaN/Inf detection)
  - Radius range validation (positive, ≤ 100 Å)
  - Detailed error messages with atom index and value

## [0.0.4] - 2025-01-22

### Added

- Multi-threading support with configurable thread count
- `--threads=N` CLI option (auto-detect by default)
- Generic thread pool implementation with work-stealing

### Changed

- Default execution mode is now multi-threaded
- 6.4x faster than FreeSASA (Python) on 3,183 atoms

## [0.0.3] - 2025-01-21

### Added

- SIMD optimization using `@Vector(4, f64)` for batch distance calculations
- 4.5x faster than FreeSASA (Python) single-threaded

## [0.0.2] - 2025-01-20

### Added

- Neighbor list optimization for O(N) neighbor lookup
- Spatial hashing with configurable cell size
- 3.9x faster than FreeSASA (Python) single-threaded

### Changed

- Reduced algorithmic complexity from O(N²) to O(N)

## [0.0.1] - 2025-01-19

### Added

- Initial Shrake-Rupley algorithm implementation
- Golden Section Spiral test point generation
- JSON input/output format
- Basic CLI with input/output file arguments
- Python scripts for structure conversion and benchmarking
  - `cif_to_input_json.py` - Convert mmCIF/PDB to input JSON
  - `calc_reference_sasa.py` - Generate reference SASA
  - `benchmark.py` - Performance benchmarking

[Unreleased]: https://github.com/N283T/zsasa/compare/v0.9.1...HEAD
[0.9.1]: https://github.com/N283T/zsasa/compare/v0.9.0...v0.9.1
[0.9.0]: https://github.com/N283T/zsasa/compare/v0.8.0...v0.9.0
[0.8.0]: https://github.com/N283T/zsasa/compare/v0.7.1...v0.8.0
[0.7.1]: https://github.com/N283T/zsasa/compare/v0.7.0...v0.7.1
[0.7.0]: https://github.com/N283T/zsasa/compare/v0.6.0...v0.7.0
[0.6.0]: https://github.com/N283T/zsasa/compare/v0.5.0...v0.6.0
[0.5.0]: https://github.com/N283T/zsasa/compare/v0.4.0...v0.5.0
[0.4.0]: https://github.com/N283T/zsasa/compare/v0.3.2...v0.4.0
[0.3.2]: https://github.com/N283T/zsasa/compare/v0.3.1...v0.3.2
[0.3.1]: https://github.com/N283T/zsasa/compare/v0.3.0...v0.3.1
[0.3.0]: https://github.com/N283T/zsasa/compare/v0.2.11...v0.3.0
[0.2.2]: https://github.com/N283T/zsasa/compare/v0.2.1...v0.2.2
[0.2.1]: https://github.com/N283T/zsasa/compare/v0.2.0...v0.2.1
[0.2.0]: https://github.com/N283T/zsasa/compare/v0.1.3...v0.2.0
[0.1.3]: https://github.com/N283T/zsasa/compare/v0.1.2...v0.1.3
[0.1.2]: https://github.com/N283T/zsasa/compare/v0.1.1...v0.1.2
[0.1.1]: https://github.com/N283T/zsasa/compare/v0.1.0...v0.1.1
[0.1.0]: https://github.com/N283T/zsasa/compare/v0.0.5...v0.1.0
[0.0.5]: https://github.com/N283T/zsasa/compare/v0.0.4...v0.0.5
[0.0.4]: https://github.com/N283T/zsasa/compare/v0.0.3...v0.0.4
[0.0.3]: https://github.com/N283T/zsasa/compare/v0.0.2...v0.0.3
[0.0.2]: https://github.com/N283T/zsasa/compare/v0.0.1...v0.0.2
[0.0.1]: https://github.com/N283T/zsasa/releases/tag/v0.0.1
