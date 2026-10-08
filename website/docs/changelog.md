---
sidebar_position: 10
---

# Changelog

All notable changes to zsasa. See [GitHub Releases](https://github.com/N283T/zsasa/releases) for full details.

## Unreleased

### Added

- **Parse-once SASA selection maps**: allow ordinary workflow chain maps to request multiple globally identified chain selections per structure, reuse one parsed/classified input and duplicate chain-set calculations, emit per-selection JSONL success/error rows, and optionally include stable source-indexed atom identity metadata.
- **Multi-interface BSA workflows**: allow CSV and JSON interface maps to request multiple stable-ID interfaces per structure while reusing one parsed/classified input, emitting per-interface success/error JSONL rows, auditable residue SASA components, and opt-in atom-level detail. (#408)
- **C ABI version**: the library exports `zsasa_abi_version()`, an integer that changes whenever an existing exported signature or struct layout changes. The Python package declares the signatures by hand and now refuses a library with another ABI version, or without the function, with an `ImportError` that names the library file, instead of calling through a stale declaration. The C API gains the error codes `ZSASA_ERROR_INVALID_FORMAT` (-7: a file that exists but is not valid for the format it was opened as) and `ZSASA_ERROR_OUTPUT_DIR` (-8: the output directory of `zsasa_batch_dir_process` cannot be created). (#435)

### Changed

- **HETATM records are excluded by default at every entry point (breaking)**: the CLI no longer includes HETATM automatically when the CCD classifier is active, and the gemmi, BioPython and Biotite integrations now default to `include_hetatm=False`. Previously the default result depended on the entry point: `zsasa calc`/`batch` with the default classifier counted ligands, ions and waters, while `process_directory()` did not (1ubq: 5656.65 vs 4834.72 Å²). Modified residues recorded as HETATM inside a chain (for example `MSE`) are excluded too. Pass `--include-hetatm`, `include_hetatm = true` (workflow) or `include_hetatm=True` (Python) to get the earlier CLI behavior. (#431)
- **Selection-map tail scheduling**: claim ordinary multi-selection chain-map files in deterministic longest-processing-time-first order using selection count, while preserving generic batch, BSA, and legacy one-row map ordering.
- **Faster test suite**: `zig build test` runs each test once instead of two or three times (the module and executable test artifacts are filtered to the tests no other artifact reaches, and `scripts/check_test_partition.py`, run in CI, fails when a test would run in no artifact or in more than one), and the batch directory tests in `src/c_api.zig` and `python/tests/test_batch_dir.py` process small temporary inputs with exact reference values instead of running Lee-Richards over the 38-model `test_data/1l2y.pdb`. The `json_writer`, `scanDirectory` and `readAtomInputFromFile` file tests now use temporary directories. (#436)
- **Test suite cleanup**: a passing `zig build test` prints nothing from the tests (they mute stderr through `src/test_support.zig`; set `ZSASA_TEST_STDERR=1` to see it), every test binary also passes when run directly (the BSA workflow test no longer starts `std.Progress` a second time), the DCD tests fail instead of passing silently when `test_data/1l2y.dcd` is missing, tests without an assertion got one, and the `zsasa traj` frame selection and DCD input, `compile-dict` and the corrected bitmask C exports are tested. A new CI job, `python_integrations`, installs gemmi, BioPython, Biotite, MDAnalysis and mdtraj and fails when a Python test is skipped; the `test_1crn_pdb` integration tests (which asserted `2000 < total_area < 3000` against a result of 3001.13 Å²) and the other real-structure tests now compare with `zsasa calc` on the same file. (#436)
- **Lee-Richards computes its arc angles exactly by default (results change)**: every neighbor now gets `acos` and `atan2`. Until now the neighbors handled in the SIMD batches of 8 and 4 got polynomial approximations that are far less accurate than their comments said (up to 0.064 rad for `atan2` and 0.0042 rad for `acos`, against the documented 0.0015 and 0.0003), and only the last 0 to 3 neighbors of an atom got exact angles. This raised Lee-Richards totals by a few tenths of a percent whatever the slice count (1ubq: 4812.47 → 4804.06 Å², 3hhb: 25583.13 → 25515.99 Å² with default options; up to 1.3 Å² on a single atom) and made the area of an atom depend on the order of the atoms in the input. With exact angles, per-atom areas agree with an independent reference to rounding error, do not depend on the atom order, and Lee-Richards at a high `--n-slices` converges to Shrake-Rupley at a high `--n-points`. **Every Lee-Richards result changes by a few tenths of a percent**: `--algorithm=lr` in `zsasa calc`, `batch` and `traj`, workflows, and the Python and C APIs. The new option `--lr-trig=exact|fast` (workflow key `lr_trig` under `[calculation]`) selects the mode, and `--lr-trig=fast` restores the previous values. The Python and C APIs always use exact angles. A Lee-Richards calculation takes about 1.3 times as long as before. The error bounds in the `src/simd.zig` comments are corrected to measured values. (#430)
- **The rich CSV has an `insertion_code` column (schema change)**: `zsasa calc --format=csv` for structure input now writes the header `chain,residue,resnum,insertion_code,atom_name,x,y,z,radius,area`. The new column is inserted directly after `resnum`, so `atom_name`, `x`, `y`, `z`, `radius` and `area` each move one position to the right: **readers that take columns by position must be updated**, readers that take them by name are unaffected. The column is empty for a residue without an insertion code, and the trailing total row has one more empty field (`,,,,,,,,,<total>`). Until now residues `10`, `10A` and `10B` all had `resnum` 10 and could not be told apart. The basic CSV (`atom_index,area`), which `calc` writes for JSON input and `batch` writes for every input, is unchanged. (#433)
- **RSA text follows the NACCESS fixed columns (`--format=rsa` output changes)**: a `RES` row now has the residue name in columns 5-7, the chain ID in column 9, the residue number right-justified in columns 10-13, the insertion code in column 14 and the five pairs of absolute and relative values in columns 16-80, as NACCESS and FreeSASA write them and as fixed-column readers such as Biopython's `Bio.PDB.NACCESS` expect them. `CHAIN` rows have the chain ID in column 10, and `CHAIN` and `TOTAL` rows have their sums in columns 12-21, 25-34, 38-47, 51-60 and 64-73. Rows were written as `RES`, the residue name, the chain ID right-justified in three columns and the residue number left-justified in four, which put the chain ID in column 11, the residue number from column 13 and every value two columns to the right; the `REM RES _ NUM` column header did not line up with the rows below it. Before: `RES MET   A 1      50.77  22.7  17.13   N/A ...` and `CHAIN  1   A     4834.7 ...`; now: `RES MET A   1    50.77  22.7  17.13   N/A ...` and `CHAIN  1 A     4834.7 ...`. **Readers written for the earlier zsasa layout, and readers that split rows at blanks, must be updated**: chain ID, residue number and insertion code now follow each other without a blank when the number has four characters (`A1000B`, `A-999`), as in a PDB file. Every value keeps a blank in front of it, so an absolute value of 1000.00 or more no longer touches the value before it (`N/A2222.27`); such a row, like one with a chain ID of more than one character, a residue name of more than three or a residue number of more than four, is written in full and moves the rest of the row to the right, with the existing warning. The warning now follows these columns: it is printed for chain IDs of two and three characters and for absolute values from 1000.00 to 9999.99, which it missed, and no longer for a four-digit residue number with an insertion code, which fits. The header line `REM  Atomic radii and reference values for relative SASA: <classifier>` is replaced by `REM  Atomic radii and polar/non-polar classes: <classifier>` and `REM  Reference values for relative SASA: Tien et al. 2013`: relative values have always been computed from the Tien et al. (2013) maximum SASA table, whatever the classifier. Values are unchanged. (#433)

### Fixed

- **BSA workflow parallelism**: apply workflow thread counts to concurrent structure-file workers, keep partner and complex SASA calculations single-threaded in the multi-file path, and serialize concurrent JSONL rows safely.
- **Batch output name collisions**: reject directories whose inputs share an output stem (for example `1crn.pdb` and `1crn.cif.gz`) before processing, instead of reporting every file as successful while later results silently overwrite earlier ones. Applies to per-file output from `zsasa batch`, workflow jobs, and `process_directory(output_dir=...)`; JSONL output is unaffected. The C API reports the condition as `ZSASA_ERROR_OUTPUT_NAME_COLLISION`. (#415)
- **SDF/MOL batch output names**: give every molecule of an SDF or MOL input its own name and per-file output in `zsasa batch`, workflow jobs and `process_directory()`. Molecules of one file that share a title get their 1-based position appended (`dup_ethanol_1`, `dup_ethanol_2`) instead of overwriting each other's output and sharing one JSONL `filename`. The output file name is the molecule name with the extension appended, so dots in the file stem or the title are no longer cut off (`lig.v2.sdf` wrote every molecule to `lig.json`; now `lig.v2_methane.json`). `/`, `\` and control characters in a title are replaced by `_` in the output file name (the title `a/b` in `lig.sdf` failed with `FileNotFound`; now `lig_a_b.json`), so a title cannot write outside the output directory; the JSONL `filename` keeps the title as written. Names change for molecules whose name is shared within their file (the name and the output file) and for stems or titles that contain dots, path separators or control characters (the output file only); all other names are unchanged. The collision check compares the real molecule output names: `lig.sdf` and `lig.mol` are accepted unless their molecules clash, an SDF molecule output that equals another input's output (`lig_1.json` from an unnamed molecule in `lig.sdf` and from `lig_1.pdb`) is rejected, and output names that differ only in ASCII case (`PROT.pdb` and `prot.cif.gz`) are rejected on every platform because they are one file on case-insensitive filesystems. SDF inputs are parsed once per batch, and an SDF file that cannot be read or parsed is reported the same way whatever the thread count: with `--threads=1` it now gets a `read/parse failed` JSONL error row too, where it was only listed in the summary as `SDF parse failed`. (#417)
- **Batch error-path crashes**: report errors raised late in a batch run instead of crashing. A workflow that failed after job setup (for example `use_bitmask = true` with an unsupported `n_points`) freed its job states twice, a single-input batch whose JSONL output could not be created freed its result buffer twice, and a JSONL write failure after parallel processing joined the worker threads a second time. BinaryCIF files whose inline chemical component data failed to load (for example more than 65,535 atoms in one component) also freed decoded columns twice. (#418)
- **CCD-derived radii for components that list hydrogens**: count the hydrogens listed in a component definition together with those implied by unfilled valence when deriving ProtOr radii. Carbons were typed as hydrogen-free whenever the definition listed its hydrogens, as wwPDB CCD entries and typical SDF files do, so sp3 carbons got 1.61 Å instead of 1.88 Å and sp2/aromatic CH carbons got 1.61 Å instead of 1.76 Å. Hydrogen-free aromatic CH carbons written with aromatic bond orders (`AROM`, SDF bond type 4) are corrected the same way. This changes SASA values for ligands and modified residues whose radii come from `--ccd`, inline `_chem_comp_*` data, `--sdf`, or direct SDF input; residues in the built-in table (standard amino acids, nucleotides) are unaffected. (#421)
- **Multi-member gzip input**: decompress every member of a `.gz` input and concatenate the output, as RFC 1952 requires, instead of silently stopping at the end of the first member. Concatenated gzip files (`cat a.gz b.gz`) and `bgzip`/BGZF files previously yielded a truncated structure with exit status 0. Each member's CRC32 and size are verified, the decompressed-size limit applies to the total, and trailing bytes that are not a gzip member (including zero padding) are now rejected as `GzipReadFailed` instead of being ignored. (#426)
- **Truncated gzip input**: report a `.gz` input that ends inside its compressed data, such as an incomplete download, as `GzipReadFailed` instead of crashing. For roughly a quarter of the possible cut positions in a typical structure file the process died with a segmentation fault, which also ended a whole `zsasa batch` run or the calling Python process; the fault is in the gzip decompressor of the Zig 0.16.0 standard library, which is now kept from reading past the end of the file. Every cut inside a member is rejected; a file cut exactly between two members of a multi-member file is itself a complete gzip file and is still read. (#442)
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
- **Neighbor grid for far-apart atoms**: store the neighbor-search grid sparsely when the bounding box of the atoms is mostly empty, instead of always allocating one slot per cell of the box. Finite coordinates that are far apart used to need memory proportional to the volume of the box (about 300 MB for two atoms 2,000 Å apart on the diagonal, tens of GB for a protein plus one atom at a `9999.999` sentinel coordinate), and at larger separations the cell count overflowed and the calculation crashed with a segmentation fault in the CLI, the C API and Python. A grid with more than 65,536 cells and more than 64 cells per atom now keeps only its occupied cells, so memory and time follow the number of atoms. Both layouts have the same cells and build each neighbor list in the same order, so SASA values do not depend on the layout and results for inputs that ran before are unchanged bit for bit. Beyond 2^21 cells along one axis (about 1.4e7 Å) the cells are made larger, which finds the same neighbors. A coordinate range too wide for the calculation precision (for example `1e300` with `--precision=f32`) is reported as `CoordinateRangeTooLarge` (`ZSASA_ERROR_INVALID_INPUT` in the C API, `ValueError` in Python). The C API calculation functions also return `ZSASA_ERROR_OUT_OF_MEMORY` instead of `ZSASA_ERROR_CALCULATION` when an allocation fails during a calculation, so Python raises `MemoryError`. (#428)
- **SDF/MOL records with a blank title**: accept a blank first line as the molecule title instead of failing with `InvalidCountsLine`. The parser skipped blank lines while looking for the title, so a blank title, which is what RDKit writes for a molecule without a name, shifted the header by one line; the same happened to every record after a `$$$$` separator. The header is now read as exactly three lines followed by the counts line. Blank lines after a `$$$$` separator and at the end of the file are still skipped: a blank line is taken as a title only when the fourth line from it is a counts line. Such a molecule has an empty name; in batch mode it is named `stem_N`. It is classified from its own bond topology like a molecule with a title (see the next entry). (#434)
- **SDF/MOL molecules are classified from their own bond topology whatever their title is**: with the `ccd` classifier, the atoms of an SDF or MOL input are matched to the molecule's own bond table directly, instead of through a component dictionary that was keyed by the molecule title. A molecule without a title was never entered in that dictionary and got element-based fallback radii only (ethanol 169.51 Å² instead of 181.07 Å², reported as `0 atoms classified, 3 fallback`). A title equal to a built-in residue name was looked up in the built-in table instead: `ALA` and `HOH` gave fallback radii, and `A`, `G` or `DT` gave atoms named `C2` or `N1` the radii of the nucleotide's atoms of that name. In `zsasa calc`, `--sdf` (or the `sdf` list of a workflow) replaced the input's own topology, so the molecule fell back to element radii unless the dictionary held a molecule of the same title. All of these now get the radii of the molecule's own bond table, in `calc`, `batch`, workflow jobs and `process_directory()`, and several untitled molecules in one file each use their own. `--sdf` and `--ccd` describe residues of PDB, mmCIF and BinaryCIF structures; they are not consulted for SDF/MOL input, and `calc` says so when `--sdf` is given. The residue name in the output is still the title, blank for a molecule without one. **This changes SASA values for SDF/MOL molecules with a blank title (including a title of spaces) or a title equal to a built-in residue name, and for `calc` on SDF/MOL input together with `--sdf`**; other SDF/MOL results are unchanged. An untitled molecule in an `--sdf` dictionary still cannot be matched to a residue and is skipped; the warning now names the file and the molecule. (#434)
- **Deuterium and tritium in SDF/MOL input**: treat atoms written with the symbols `D` and `T` as hydrogen, like the PDB, mmCIF and BinaryCIF parsers do for deuterium. They were kept as atoms of an unknown element (1.00 Å) when hydrogens were excluded, so ethanol with one deuterium gave 4 atoms instead of 3, and they counted as heavy atoms when radii were derived from the bond table, so a CD3 carbon got the hydrogen-free radius (1.61 Å instead of 1.88 Å). They are now excluded by default, included with `--include-hydrogens` as hydrogens (1.10 Å, named `H1`, `H2`, ... together with the other hydrogens instead of `X1`), and counted as hydrogens of the atom they are bonded to. Deuterium is counted the same way in component definitions from `--ccd`, inline `_chem_comp_*` data and `--sdf`. **This changes SASA values for SDF/MOL inputs and component definitions that contain deuterium or tritium.** (#434)
- **SDF V3000 molecules without atoms**: read atoms only between `M  V30 BEGIN ATOM` and `M  V30 END ATOM`, bonds only between `BEGIN BOND` and `END BOND`, and stop at `M  END`. A molecule with `COUNTS 0 0` and no atom block made the parser read on into the following record, so the file failed with `InvalidV3000` or lost molecules. Such a molecule is now read as a molecule without atoms, like its V2000 counterpart: `zsasa calc` on it fails with `NoAtoms`, and in batch mode it is reported as failed while the other molecules of the file are calculated. An atom or bond block that is still open at `M  END`, `$$$$` or the end of the file is rejected as `InvalidV3000`. (#434)
- **SDF V3000 continuation lines**: join a V3000 line that ends in `-` with the `M  V30` line that continues it before reading it, in every block. A continuation line used to be read as a line of its own, so a continued atom or bond line failed with `InvalidAtomLine`, `InvalidBondLine` or `InvalidFloat`, and a continuation that happened to look like an atom line was added as an extra atom. A continued line that is not followed by an `M  V30` line is rejected as `InvalidV3000`. (#434)
- **SDF/MOL atom names past 999 atoms of an element**: keep the atom names that zsasa generates for a molecule (`C1`, `C2`, `O1`, ...) unique. Names are limited to four characters, and the counter was cut to fit: from the 1,000th atom of a one-letter element every atom was named just `C`, and from the 100th atom of a two-letter element the last digit was dropped (`Cl100` became `Cl10`). The CCD classifier looks radii up by atom name, so in a chain of 1,050 carbons every carbon from the 1,000th on got the radius of the last one. The counter now continues in base 36 with a leading letter, as in hybrid-36: `CA00` follows `C999` (up to `CZZZ`, the 34,695th atom) and `ClA0` follows `Cl99` (up to `ClZZ`, the 1,035th). Atoms beyond that are named with a four-character base-36 number without the element symbol (`0000`, `0001`, ...). Names of molecules with fewer atoms per element are unchanged. **This changes radii and SASA values for molecules with 1,000 or more atoms of carbon, nitrogen, oxygen, sulfur or phosphorus** under the default `ccd` classifier. (#434)
- **Lee-Richards tangent circles**: a slice circle that touches a neighbor's circle from inside is now buried, and a neighbor circle that touches from inside no longer buries the slice. The containment tests were strict, so at exact tangency the covered arc was computed as the whole circle or as an arc of zero width, and bringing its ends into [0, 2π] could then turn the whole circle into an empty arc and the empty arc into the whole circle, depending on the direction of the neighbor. With an atom of radius 0.6 Å at the origin, one of radius 1.6 Å at x = ±1 Å and the default probe (the expanded spheres touch from inside) and an odd `--n-slices`, the small atom got up to 50.27 Å² instead of 0 and the large atom as little as 0 instead of 113.10 Å². At f32 the same could happen on ordinary input through rounding. The fix applies to both `--lr-trig` modes and changes results only where such a tangent pair occurs (nowhere in the bundled examples). (#430)
- **f32 totals**: with `--precision=f32` (the `traj` default) the total area is now the sum of the per-atom areas accumulated in f64, for Shrake-Rupley, bitmask Shrake-Rupley and Lee-Richards. It was accumulated in f32, so it differed from the sum of the per-atom areas it was reported with and changed with the thread count (3hhb, Lee-Richards: 25516.0254 Å² with one thread, 25515.9902 Å² with four, sum of the per-atom areas 25515.9910 Å²). Per-atom areas are unchanged; f32 totals change in about the seventh significant digit. In the Zig API, `SasaResultGen(T).total_area` is now `f64` for every `T`. (#430)
- **`batch --af-model-fast` with `--include-hetatm`**: read the file with the generic mmCIF parser when HETATM records are requested. The fast parser reads only the leading `ATOM` rows, so `--include-hetatm` was ignored for a file whose `ATOM` rows are followed by `HETATM` rows, which the documentation described as a case for the generic fallback. Runs without `--include-hetatm` are unchanged. (#432)
- **BSA residue rows with chain IDs longer than four characters**: `residue_chain` in residue-level BSA rows holds the full chain ID, like `atom_chain` in the same row. It was cut to four characters, so the two arrays disagreed for mmCIF and BinaryCIF chains with longer IDs. (#432)
- **Stale JSONL file after a batch on an empty directory**: `zsasa batch empty/ -o out.jsonl --format=jsonl` leaves an empty `out.jsonl` whatever the thread count. With more than one thread the run returned before creating the file, so the rows of an earlier run survived; `--threads=1` already truncated it. A per-file output directory is likewise created in both cases. (#432)
- **Empty chain lists are rejected**: a workflow job with `chains = []` is an error when the workflow is read (`EmptyJobChains`), and `zsasa batch --chain=,` (or any `--chain` value without a chain ID) is an error before the input directory is scanned. An empty `chains` array selected every chain when the workflow ran file by file and no atom at all when it ran job by job, and an empty `--chain` list failed every input with `NoAtomsFound` or `NoAtomSiteLoop`. Leave `chains` out to select every chain. (#432)
- **Workflow JSONL without atom areas on standard output**: a single-job workflow with `[output.jsonl] atom_areas = false` and no output directory writes its rows to standard output in every case. When the workflow ran job by job (for example because the input directory contains an SDF file or the job sets `auth_chain`) it wrote nothing at all and still exited with status 0, because the batch runners took "keep atom areas" to mean "JSONL mode". They now go by the output format, also for the per-file output decision and the output-name collision check. (#432)
- **Workflow JSONL error rows for unreadable inputs**: a workflow that runs file by file (every file parsed once and shared by the jobs) writes a `{"status":"err",...}` row to the JSONL output of every job for an input that cannot be read, parsed or classified, with the wording of the other runners (`read/parse failed: NoAtomsFound`). It only counted such inputs and printed them to standard error, so the input was missing from the JSONL files, while `zsasa batch` and workflows that run job by job wrote the row. Every job's JSONL file now has one row per input whichever way the workflow runs. (#432)
- **A workflow job that fails as a whole gives a non-zero exit status**: when a workflow runs job by job (an SDF file in the input directory, a job-level `auth_chain`, a chain map) and a job cannot run at all, because the input directory is missing, its inputs share per-file output names, or its output directory or JSONL file cannot be created, `zsasa batch --workflow` now exits with status 1 (`Error: WorkflowJobFailed`). It used to count the job as one more "failed" input in `Workflow complete: N successful, M failed` and exit with status 0, while the same condition was already an error for `zsasa batch` and for workflows that run file by file. The remaining jobs still run. Each failed job is reported with its cause when it fails (`Error running workflow job 'all': cannot read input directory 'structures': FileNotFound`), the summary line counts inputs only, and a closing line names the failed jobs (`1 of 3 jobs failed: all`). A chain-map job that cannot run no longer stops the jobs after it. **Scripts that relied on status 0 from such a run now see the failure.** (#432)
- **`-q`/`--quiet` no longer hides inputs that fail**: `zsasa batch` and `zsasa batch --workflow` now print the inputs they could not process to standard error at the end of every run, with or without `--quiet`: `2 of 40 inputs failed:` followed by one `  <input>: <reason>` line per input, and for a workflow one such block per job (`Job 'chain_a': 2 of 40 inputs failed:`). **`-q` now prints these lines to standard error whenever an input fails**; it still prints nothing when every input succeeds, and quiet mode keeps suppressing progress and the summary. Until now the list of failed inputs was part of the summary, so `zsasa batch -q in/ out/` on a directory with a corrupt file wrote no output file for it, printed nothing and exited with status 0. The exit status is unchanged: inputs that fail are reported, not fatal. At most 20 inputs are listed per run or job; a last line counts the rest and, for JSONL output, names the file whose `"status":"err"` rows hold all of them. Without `--quiet` the report replaces the `Failed files:` block of the summary, so no failure is printed twice, and that list is now limited to 20 entries as well. Workflows that run file by file no longer print an `Error running workflow job '<job>' on '<file>': ...` line per failure as they go; the same information is in the per-job report, which workflows that run job by job, chain-map jobs (counted in selections) and BSA analysis (counted in interfaces) now print too. (#432)
- **One residue definition for every per-residue output**: the `--per-residue` / `--rsa` table, the RSA file (`--format=rsa`) and the JSONL residue map now group atoms by the same rule, so they report the same residues with the same atom counts and areas. A residue is a run of consecutive atoms with the same chain ID, residue number, insertion code and residue name. Until now the table ignored the residue name, so `GLY A 10` followed by `LYS A 10` was printed as one `GLY` row holding the atoms of both, while the RSA file and the residue map had two. The table and the RSA file also merged atoms of one residue from anywhere in the input, where the residue map has always had one entry per run. **Rows change in two cases**: for a residue whose atoms are not contiguous in the input, the table and the RSA file now have one row per run instead of one merged row; and for a multi-model file read with all models superimposed (the default), they now have one row per residue and model instead of one row per residue with the areas of all models added up (`test_data/1l2y.pdb`: 760 rows instead of 20; use `--model=N` for one model). The JSONL residue map is unchanged. (#433)
- **Long chain IDs in CSV and RSA output**: chain IDs longer than four characters (mmCIF and BinaryCIF input) are written in full to the rich CSV and to the RSA file (`--format=rsa`), and the RSA file tells chains apart by the full ID. Both printed only the first four characters and the RSA file also grouped by them, so chains `AAAAA`, `AAAAB` and `AAAA` gave one `RES` row with the atoms of all three and one `CHAIN` line. The `--per-residue` table and the JSONL residue map already used the full ID. (#433)
- **CSV quoting**: text fields of the rich CSV (chain ID, residue name, insertion code, atom name) are quoted as RFC 4180 specifies when they contain a comma, a double quote, a carriage return or a line feed: the field is enclosed in double quotes and every double quote in it is doubled. They were written as they are, so a PDB chain ID of `,` produced a row with one field too many, a chain ID of `"` an unterminated quoted field, and an SDF title such as `a,b` (the residue name of the molecule) a shifted row. Fields that need no quoting are written exactly as before. (#433)
- **RSA polar and non-polar columns follow the classifier**: the `Non-polar` and `All polar` columns of `--format=rsa` are split by the polarity class that the active classifier gives each atom. They were split by element whatever the classifier (N, O, P and S polar, everything else non-polar), although NACCESS classes sulfur as non-polar and OONS classes carbonyl carbon as polar, so the columns contradicted the classifier's own classes. **The two columns change for NACCESS and OONS runs** (1ubq `TOTAL` row, non-polar / polar: NACCESS 2462.8 / 2360.5 → 2469.4 / 2353.9, OONS 2583.0 / 2196.5 → 2542.6 / 2236.9); with CCD and ProtOr they change only for atoms that these classifiers class differently from the element rule, such as selenium. Atoms that the classifier does not class (hydrogens, ligands outside its tables) are still classed by element. `--polar` now prints this split too, as a second block `Polar/Nonpolar SASA by atom class (classifier: <name>)` below the existing summary by residue type, which is unchanged; its areas equal the last two values of the RSA `TOTAL` row. The all-atom, side-chain and main-chain columns are unchanged. (#433)
- **Python per-residue aggregation and insertion codes**: `aggregate_from_result()` no longer merges residues that differ only in their insertion code. `AtomData` had no insertion-code field and the gemmi, BioPython and Biotite extractors dropped the code, so `aggregate_by_residue()` grouped by chain ID and residue number alone: with residue 2 of 1ubq renumbered to `1A`, it returned 75 residues instead of 76, with `MET 1` holding the 17 atoms of both. `AtomData` has a new optional field `insertion_codes` (default `None`), which the three extractors fill; `aggregate_by_residue()` takes a new optional argument `insertion_codes` and groups by chain ID, residue number and insertion code; and `ResidueResult` has a new field `insertion_code` (default `""`). The new fields and the new argument come last and are optional, so existing calls and positional construction keep working, and results for structures without insertion codes are unchanged. (#433)
- **Alternate locations with equal occupancies**: `--altloc=highest-occupancy` now keeps exactly one alternate per atom when occupancies tie, the one that comes first in the file. Each atom of a tie used to count as the best one, so a 0.50/0.50 pair, or a file without an occupancy column, kept every alternate superimposed, the same as `--altloc=all`. Atoms without an altloc ID are always kept in this mode (they used to compete on occupancy with other atoms of the same name in the same residue), and alternates of an atom that is also listed without an altloc ID are dropped, as with `auto`. **Atom counts and SASA change for inputs with tied occupancies under `highest-occupancy`** (2VB1: 1008 to 1001 atoms, 1US0: 2549 to 2492); the default `auto` already kept one alternate on a tie. (#425)
- **Microheterogeneity**: where the alternates of one residue position are different residues (residue 22 of 1EJG is PRO as altloc A and SER as altlocs B and C), one whole residue is kept instead of both superimposed. The alternate is decided per residue position (chain, residue number, insertion code) before it is decided per atom: `auto` keeps the residue that has altloc A, or without an A the residue with the highest occupancy; `highest-occupancy` keeps the residue with the highest occupancy; a single ID keeps the residue that carries it. The occupancy of a residue is that of its first atom with an altloc ID in the file, and a tie goes to the residue that comes first. Atoms without an altloc ID are kept under their own residue name. **Atom counts and SASA change for entries with microheterogeneous residues** with `auto` and `highest-occupancy` (1EJG with `auto`: 340 to 327 atoms, 2939.42 to 2930.98 Å²), for PDB, mmCIF and BinaryCIF input. (#425)
- **`--altloc` for PDB input**: the option now applies to PDB files in `calc`, `batch` and `traj` topologies exactly as it does to mmCIF and BinaryCIF. It was ignored there, so every mode behaved like `auto`, including `--altloc=none`, which exists to fail fast. **Atom counts change for PDB input when `--altloc` is given with a mode other than `auto`**: `all` keeps every alternate, a single ID and `highest-occupancy` select as documented, and `none` fails with `UnexpectedAltLoc` when an atom has an altloc ID. The PDB and mmCIF files of an entry are now read by the same rules in every mode. (#425)
- **`--altloc` with `--workflow`**: `zsasa batch --workflow wf.toml --altloc=B` no longer behaves like `auto`. The command-line value was dropped in every batch workflow path (jobs, chain maps and BSA analysis). Workflows also accept an `altloc` key under `[calculation]` with the values of the option, for `calc --workflow` and `batch --workflow`; `--altloc` on the command line takes precedence, and an unknown value is rejected when the workflow is read. (#425)
- **Alternate-location parsing time**: selecting alternates is linear in the number of atoms. Every atom with an altloc ID used to start a scan over all atoms of the file, so a 120,000-atom file in which every atom has altloc A took 8.5 s (PDB) and 15.8 s (mmCIF) to parse; it now takes 0.04 s and 0.05 s. (#425)
- **Gemmi integration and alternate conformers**: `zsasa.integrations.gemmi` keeps one conformer per site by the rules of `--altloc=auto` instead of counting the atoms of every alternate conformer, so it uses the same atoms as `zsasa calc`. **Atom counts and SASA from the Gemmi integration change for structures with alternate locations** (1EJG: 424 to 327 atoms). The Gemmi and BioPython integrations also treat deuterium (element `D`) as hydrogen, as the CLI and the Biotite integration do, so `include_hydrogens=False` removes it. (#425)
- **Stale bundled library in a development checkout (Python)**: when both `python/zsasa/libzsasa.*` (a copy made by the build hook) and `zig-out/lib/libzsasa.*` exist and differ, the newer file is loaded and a warning names the one that was skipped, instead of always loading the bundled one. `ZSASA_LIB` still has the highest priority. The current working directory is no longer searched for the library, because loading a shared library from wherever the process runs would let a stray file execute code. (#435)
- **Misleading exceptions from the Python bindings**: `process_directory()` with an output directory that cannot be created (for example below a regular file) raises `NotADirectoryError`, `FileExistsError` or `PermissionError` naming the output path instead of `FileNotFoundError` blaming the input directory, and an input path that is a file raises `NotADirectoryError`; the C API reports the two cases as `ZSASA_ERROR_OUTPUT_DIR` and `ZSASA_ERROR_FILE_IO`. `XtcReader` and `DcdReader` raise `ValueError` for an existing file that is empty, truncated or in another format (an XTC reader on a DCD file used to raise `FileNotFoundError`, an empty file `RuntimeError: Error opening XTC file: -3`), `FileNotFoundError` only for a missing file, and `IsADirectoryError` or `PermissionError` where those apply; the C open and read functions return `ZSASA_ERROR_INVALID_FORMAT` for malformed files, and out-of-memory in them maps to `ZSASA_ERROR_OUT_OF_MEMORY`. `compute_sasa_trajectory()`, `compute_sasa_trajectory_summary()` (XTC and DCD) and `SASAAnalysis.run()` raise `ValueError` for `step=0`, a negative `step`, `start` or `stop` before opening anything, where `step=0` was a `ZeroDivisionError` and negative values were silently accepted. (#435)
- **Concurrent `process_directory()` calls and signal handlers**: the batch C API created and destroyed a thread runtime per call, which replaces the process's SIGPIPE and SIGIO handlers on creation and restores the previous ones on destruction; overlapping calls restored each other's replacement handlers. The runtime is now created once, on first use, and kept for the life of the process, so the handlers are replaced at most once and never restored. (#435)
- **MDAnalysis radii without an element attribute**: the integration took the first character of `atom.type` (or the name) when the topology has no element column, which read `CL`, `ZN` and `FE` as carbon, an unknown element and fluorine. It now reads two-letter elements from the type and name together with the residue name: `CL`, `BR`, `FE`, `ZN`, `MG` and a few more are that element, `CA`, `CD`, `NA`, `HG` and similar names only when the residue is the ion itself (the alpha carbon `CA` of a protein and the `NA` nitrogen of a heme stay organic), and the CHARMM ion names `SOD`, `POT`, `CLA`, `CAL` give their elements. The radii table is read from `MDAnalysis.guesser.tables` and, before MDAnalysis 2.8 where that module does not exist, from `MDAnalysis.topology.tables`; the six-element built-in table that was used silently on those versions is only the fallback when MDAnalysis is not installed. (#435)

- **`install.sh` no longer installs a binary it could not verify**: until now the installer warned and installed anyway when `SHA256SUMS` could not be downloaded, printed "Checksum verified." when neither `sha256sum` nor `shasum` was installed, and also accepted a `SHA256SUMS` without an entry for the binary. It now stops with an error in each of these cases, and matches the asset name exactly instead of as a substring. **Behavior change**: setting `SKIP_CHECKSUM=1` is the only way to install without verification (`curl ... | SKIP_CHECKSUM=1 sh`). (#437)
- **Stale conda-forge, AUR and Nix packaging metadata**: the checksums in `packaging/conda-forge/meta.yaml` were those of 0.2.4 although the recipe pointed at 0.9.1, and `packaging/aur/PKGBUILD` and `.SRCINFO` were still at 0.2.4; both now match the published 0.9.1 assets. The Nix flake could not be built at all: its dependency hash was out of date, and `zig build --fetch` of Zig 0.16 unpacks the packages into `./zig-pkg` instead of the global cache, which the flake did not capture (`package not found`). `nix build` works again. `scripts/update_packaging_checksums.py` fills the conda and AUR files from the checksums of a published release (`--check` fails while they are stale), `scripts/release_bump.py` now marks the recipe's checksums as pending instead of leaving the previous release's, and `scripts/check_nix_deps_hash.py` detects a Nix dependency hash that was not refreshed after `build.zig.zon` changed (`--refresh` recomputes it). (#437)

### Documentation

- **Documentation brought back in line with the CLI and Python API**: the CLI reference marks which of `calc`, `batch` and `traj` accepts each option and lists the previously undocumented ones (`--sdf`, `--mol`, `--input-io`, `--af-model-fast`, `--profile-stages`, `--exclude-hydrogens`), the error-message table matches the binary, and `--help` lists every supported input extension and describes the structure CSV columns; the workflow guide has a reference for every key and for legacy flat-root manifests; the Python docs match the signatures (`classifier_type` in the classifier utility table, `chunk_size` and `store_atom_areas` for MDAnalysis, `compute_sasa_trajectory_summary` and `zsasa.dcd`, `ClassifierType.PROTOR` described as its own value rather than an alias of CCD, default CCD radii in the examples); and the batch, trajectory and classifier guides, `CONTRIBUTING.md` and `README.md` were corrected. (#438)
- **`ClassifierType.CCD` in the Python classify functions and integrations**: the docstrings and the Python API, integration and classifier guides now state that `classify_atoms`, `get_radius`, `get_atom_class` and the gemmi, BioPython and Biotite integrations use the built-in ProtOr table for `CCD`, with no bond-topology analysis of unknown components, unlike the CLI and `process_directory()`. The documented default radii were checked against the code. (#435)

## [v0.9.1](https://github.com/N283T/zsasa/releases/tag/v0.9.1) — 2026-07-28

### Added

- **Per-file workflow chain maps**: add CSV and JSON maps for file-specific single- or multi-chain SASA selections and BSA/ΔSASA partner groups, with per-entry `label` or `auth` asym ID selection for mmCIF and BinaryCIF inputs. (#406)

## [v0.9.0](https://github.com/N283T/zsasa/releases/tag/v0.9.0) — 2026-07-15

### Added

- **Batch AF model fast parser**: add experimental `zsasa batch --af-model-fast` for AlphaFold-like mmCIF batches. The fast path preserves multiple chains and atom metadata, and falls back to the generic parser for unsupported layouts.
- **Batch input I/O selection**: add `--input-io=auto|mmap|read` to compare or select file input strategies where supported.
- **Batch stage profiling**: add `--profile-stages` with `--timing` to report aggregate read/parse, classifier, and JSONL write times.

### Changed

- **Batch JSONL streaming**: buffer normal-sized JSONL records while streaming large records directly to avoid regressions for PDB batches with large atom-area payloads.

## [v0.8.0](https://github.com/N283T/zsasa/releases/tag/v0.8.0) — 2026-07-01

### Added

- **AltLoc parsing modes**: add shared `--altloc` policy support for mmCIF and BinaryCIF inputs across `calc`, `batch`, and `traj`, including `auto`, `none`, `all`, `highest-occupancy`, and single-ID selection. (#394)
- **Workflow JSONL output controls**: add workflow `[output.jsonl]` controls for atom areas, total areas, decimal rounding, and metadata sidecars. (#396)

### Changed

- **Batch JSONL output**: add explicit success/error status rows, configurable decimal rounding, and stricter serialization/write error handling for JSONL output. (#395)
- **mmCIF parser performance**: reduce atom_site parsing overhead by batching tokenizer position updates, skipping unnecessary AltLoc selection, and lazily backfilling extended chain IDs. (#393)
- **ProtOr classifier**: make `--classifier=protor` use static ProtOr-compatible radii without loading inline, external, or SDF-derived CCD resources. This provides a faster protein-only path for mmCIF/PDB inputs such as AFDB models. (#392)

### Fixed

- **Inline CCD parsing**: gate inline CCD extraction so batch and trajectory topology parsing skip unnecessary CCD work, while `batch --classifier=ccd` still uses per-file inline CCD data when present. (#392)

## [v0.7.1](https://github.com/N283T/zsasa/releases/tag/v0.7.1) — 2026-06-29

### Fixed

- **RSA text compatibility warning**: warn when `--format=rsa` output exceeds legacy NACCESS fixed-width columns, and document JSON/JSONL as the recommended machine-readable formats.

## [v0.7.0](https://github.com/N283T/zsasa/releases/tag/v0.7.0) — 2026-06-28

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

## [v0.6.0](https://github.com/N283T/zsasa/releases/tag/v0.6.0) — 2026-05-21

### Added

- **ztraj-backed trajectory formats**: `zsasa traj` now supports TRR and AMBER NetCDF alongside XTC and DCD, with coordinates normalized to Å before SASA calculation. (#378)

### Changed

- **Trajectory reader backend**: use `ztraj` for trajectory readers and dependency notices. (#377)
- **Python trajectory guidance**: recommend `pyztraj` for direct Python trajectory-file I/O while keeping `zsasa.xtc` and `zsasa.dcd` as compatibility APIs. (#379)

### Fixed

- **Batch JSONL and FFI concurrency coverage**: add regression coverage for batch JSONL output and concurrent FFI calls. (#376)

## [v0.5.0](https://github.com/N283T/zsasa/releases/tag/v0.5.0) — 2026-05-20

### Added

- **BinaryCIF input support**: calc and batch now accept `.bcif`, `.bcif.gz`, and `.bcif.zst` files by decoding `_atom_site` directly in Zig. Inline CCD extraction from BinaryCIF remains out of scope for this release. (#372)
- **Adaptive batch bitmask SR**: add experimental `zsasa batch --adaptive-sr` mode with coarse/fine point controls for two-stage bitmask Shrake-Rupley runs. (#371)

### Changed

- **Workflow batch execution**: reuse parsed/classified structures across eligible workflow jobs so named chain analyses such as chain A, chain B, and complex AB avoid repeated parsing while preserving output schemas and compatibility fallbacks. (#374)
- **Website documentation**: restructure the documentation site around task-oriented guides, CLI references, Python APIs, and integrations. (#370)

### Fixed

- **BinaryCIF CCD classification**: use inline CCD components from BinaryCIF inputs when available. (#373)

## [v0.4.0](https://github.com/N283T/zsasa/releases/tag/v0.4.0) — 2026-05-19

### Added

- **Workflow files**: add TOML workflow support via `zsasa calc --workflow` and `zsasa batch --workflow`; `batch --manifest` remains a compatibility alias. (#366)

### Changed

- **Breaking**: custom classifier configs are TOML-only; legacy FreeSASA-style custom classifier files are no longer supported. (#366)
- **Breaking/API**: the public `batch_manifest` module was removed and replaced by `workflow_manifest`. (#366)
- **Repository automation**: add agent instructions and ignore local planning docs. (#367)

### Fixed

- **CI**: ignore gitignore-only changes in automation checks. (#368)

## [v0.3.2](https://github.com/N283T/zsasa/releases/tag/v0.3.2) — 2026-05-17

### Added

- **Batch JSONL residue maps**: add opt-in `--residue-map` / `residue_map = true` output with compact columnar residue identifiers, atom ranges, and residue SASA for large-scale JSONL workflows. (#364)

## [v0.3.1](https://github.com/N283T/zsasa/releases/tag/v0.3.1) — 2026-05-17

### Added

- **Batch TOML manifests**: add the original named multi-job workflow file support for chain A, chain B, and AB complex SASA; includes batch `--chain` / `--auth-chain` support, per-job JSONL and per-file output layout, and CLI/docs coverage. (#361)

## [v0.3.0](https://github.com/N283T/zsasa/releases/tag/v0.3.0) — 2026-05-17

### Added

- **Zstandard decompression for structure inputs**: `.json.zst`, `.pdb.zst`, `.cif.zst`, `.mmcif.zst`, `.ent.zst`, `.sdf.zst`, and `.mol.zst` are now detected and transparently decompressed via native `std.compress.zstd`. (#357)

### Changed

- **CLI progress bars** now use Zig's standard progress API for `batch` and `traj`, reducing custom rendering logic while preserving progress reporting. (#358)
- **Project branding documentation**: add the project logo to README and refresh logo SVG assets. (#355, #356)

## [v0.2.11](https://github.com/N283T/zsasa/releases/tag/v0.2.11) — 2026-04-26

### Fixed

- **Build & Publish workflow**: Bump remaining Zig **0.15.2 → 0.16.0** pins missed by PR1 (#345). The `Dockerfile` and the `cibuildwheel` `before-all` in `python/pyproject.toml` were both pulling the old toolchain and rejecting 0.16-only API; the result is that wheels and the Docker image have not been published since v0.2.8. Also updates `install.sh`, `python/hatch_build.py`, `python/zsasa/_ffi.py` error messages, `CONTRIBUTING.md`, the website docs, and the bug-report ISSUE template. (#353)
- **Windows wheel test**: `examples/1ubq.cif` was failing with `StreamTooLong` because the Windows `mmapFile` path used `.limited64(stat.size)` and Windows occasionally reports a byte count divergent from the buffered reader. Switch to `.unlimited` — `mmapFile`'s callers already trust the input path. (#353)
- **aarch64 wheel build**: `cibuildwheel`'s short alias `manylinux_2_28` resolved to the retired quay.io tag `2026.03.01-1`. Pin both `manylinux-x86_64-image` and `manylinux-aarch64-image` to an explicit current tag (`2026.04.25-0`). (#353)
- **`flake.nix` derivation version**: Catch up to `0.2.11` (was stale at `0.2.4` until v0.2.10).

## [v0.2.10](https://github.com/N283T/zsasa/releases/tag/v0.2.10) — 2026-04-26

### Changed

- **Gzip decompression: C-zlib → native `std.compress.flate`**. Reverts the workaround introduced in #320 now that the upstream panic ([ziglang/zig#25035](https://github.com/ziglang/zig/issues/25035)) is fixed in Zig 0.16. Public API (`readGzip`, `readGzipLimited`, `DEFAULT_MAX_SIZE`, `GzipError`) and the 4 GB decompression-bomb cap unchanged; callers untouched. (#351)

### Added

- **Explicit gzip CRC32 + ISIZE trailer verification** in `gzip.zig`. `std.compress.flate` parses the gzip trailer fields but does not compare them against the decompressed bytes; without this check a corrupt `.cif.gz` could silently decompress to wrong bytes. Adds a regression test (`readGzip rejects gzip with corrupted CRC`). (#351)

### Removed

- **`zlib` dependency** — removed from `build.zig.zon`, `build.zig` (`b.dependency`, `linkLibrary`, `b.addTranslateC`), and `src/c/zlib_wrapper.h`. The shared library no longer links zlib. (#351)

## [v0.2.9](https://github.com/N283T/zsasa/releases/tag/v0.2.9) — 2026-04-26

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

## [v0.2.8](https://github.com/N283T/zsasa/releases/tag/v0.2.8) — 2026-04-13

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

## [v0.2.7](https://github.com/N283T/zsasa/releases/tag/v0.2.7) — 2026-04-12

### Added

- **CCD classifier** (`--classifier=ccd`): New classifier that derives ProtOr-compatible radii from CCD (Chemical Component Dictionary) bond topology, enabling accurate radius assignment for any chemical component — not just standard amino acids (#326, #327)
- **External CCD dictionary** (`--ccd=<path>`): Load external CCD dictionary for non-standard residues. Supports both CIF text (`.cif`, `.cif.gz`) and binary ZSDC format (#328)
- **`compile-dict` subcommand**: Convert CCD dictionary from CIF text to compact binary ZSDC format for faster loading (`zsasa compile-dict components.cif.gz -o components.zsdc`) (#328)
- **Python**: `ClassifierType.CCD` added to Python bindings (#330)

### Changed

- **CCD classifier auto-includes HETATM**: When using `--classifier=ccd`, HETATM records are included automatically without needing `--include-hetatm` (#329)
- **Python**: Split monolithic `core.py` into focused modules (`_ffi.py`, `sasa.py`, `classifier.py`, `rsa.py`, `batch.py`). Public API unchanged (#332)

## [v0.2.6](https://github.com/N283T/zsasa/releases/tag/v0.2.6) — 2026-03-22

### Fixed

- **Docker build**: use PIC-enabled zlib for shared library to fix `R_X86_64_32` relocation errors in Docker builds (#323)
- **CI**: vendor zlib from source (allyourcodebase/zlib) instead of linking system library, fixing builds on Windows and manylinux containers (#322)

### Changed

- **Publish workflow**: decouple job dependencies to prevent cascading failures — PyPI, GitHub Release, Docker, and package managers now run independently. Added `workflow_dispatch` with selective job re-runs (#324)

## [v0.2.5](https://github.com/N283T/zsasa/releases/tag/v0.2.5) — 2026-03-22

### Fixed

- **Gzip decompression for mmCIF/PDB files**: `.cif.gz`, `.pdb.gz`, `.mmcif.gz`, `.ent.gz` files are now transparently decompressed. Previously only `.json.gz` was supported, and mmCIF/PDB parsers passed raw gzip data to the tokenizer, causing `NoAtomSiteLoop` errors (#319)

### Changed

- **Switch from Zig native flate to C zlib**: gzip decompression now uses C zlib (`gzopen`/`gzread`/`gzclose`) instead of Zig 0.15's `std.compress.flate`, which panics on certain valid gzip files (e.g. PDB entry 2OXD). This adds `libz` as a build dependency. See [ziglang/zig#25035](https://github.com/ziglang/zig/issues/25035). Will revert to native flate when the upstream bug is fixed (#320)
- **Decompression bomb protection**: `readGzip` enforces a 4 GB max decompressed size limit to prevent memory exhaustion from malicious `.gz` files

## [v0.2.4](https://github.com/N283T/zsasa/releases/tag/v0.2.4) — 2026-03-11

### Added

- **CLI installer**: `install.sh` script for one-line CLI installation — builds from source with Zig or downloads pre-built binary from GitHub Releases (#307)
- **CLI release binaries**: pre-built CLI binaries for linux-x86_64, linux-aarch64, macos-x86_64, macos-aarch64, windows-x86_64 attached to GitHub Releases (#307)
- **Nix flake**: `nix run github:N283T/zsasa` or `nix profile install github:N283T/zsasa` for Nix-based installation (#311)

### Changed

- **XTC reader**: replaced local `src/xtc.zig` with external [zxdrfile](https://github.com/N283T/zxdrfile) package dependency (#309)
- **CI**: update cibuildwheel v2.22 to v3.4.0, remove obsolete skip selectors (#305, #306)

## [v0.2.3](https://github.com/N283T/zsasa/releases/tag/v0.2.3) — 2026-03-10

### Fixed

- **Windows build**: fall back to heap allocation on Windows where POSIX mmap is unavailable (#301)
- **Windows build**: replace `std.mem.zeroes` with `undefined` for `JsonlStreamWriter` placeholder — zeroes fails on non-nullable pointers on Windows (#301)
- **Third-party license compliance**: add BSD-2-Clause notice for libxdrfile (via chemfiles/xdrfile) in `src/xtc.zig` and create `THIRD_PARTY_NOTICES.md` (#300)
- Correct Lahuta author name (Besian I. Sejdiu) in docs and README (#299)

### Changed

- Remove legacy `docs/` directory — all documentation migrated to website (#298)
- Add acknowledgments section to README and comparison page (#297)

## [v0.2.2](https://github.com/N283T/zsasa/releases/tag/v0.2.2) — 2026-03-10

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

## [v0.2.1](https://github.com/N283T/zsasa/releases/tag/v0.2.1) — 2026-02-25

### Changed

- Relaxed bitmask `n_points` constraint from fixed 64/128/256 to any value 1..1024, with internal storage expanded from `[4]u64` to `[16]u64` (#210)
- Updated Python bindings, C API, CLI help text, and documentation to reflect new n_points range (#210)

### Added

- `--use-bitmask` and `--n-points` flags to benchmark scripts for bitmask LUT benchmarking (#208)

## [v0.2.0](https://github.com/N283T/zsasa/releases/tag/v0.2.0) — 2026-02-25

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

## [v0.1.3](https://github.com/N283T/zsasa/releases/tag/v0.1.3) — 2026-02-25

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

## [v0.1.2](https://github.com/N283T/zsasa/releases/tag/v0.1.2) — 2026-02-22

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

## [v0.1.1](https://github.com/N283T/zsasa/releases/tag/v0.1.1) — 2026-02-22

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

## [v0.1.0](https://github.com/N283T/zsasa/releases/tag/v0.1.0) — 2026-01-31

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

## [v0.0.5](https://github.com/N283T/zsasa/releases/tag/v0.0.5) — 2025-01-23

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

## [v0.0.4](https://github.com/N283T/zsasa/releases/tag/v0.0.4) — 2025-01-22

### Added

- Multi-threading support with configurable thread count
- `--threads=N` CLI option (auto-detect by default)
- Generic thread pool implementation with work-stealing

### Changed

- Default execution mode is now multi-threaded
- 6.4x faster than FreeSASA (Python) on 3,183 atoms

## [v0.0.3](https://github.com/N283T/zsasa/releases/tag/v0.0.3) — 2025-01-21

### Added

- SIMD optimization using `@Vector(4, f64)` for batch distance calculations
- 4.5x faster than FreeSASA (Python) single-threaded

## [v0.0.2](https://github.com/N283T/zsasa/releases/tag/v0.0.2) — 2025-01-20

### Added

- Neighbor list optimization for O(N) neighbor lookup
- Spatial hashing with configurable cell size
- 3.9x faster than FreeSASA (Python) single-threaded

### Changed

- Reduced algorithmic complexity from O(N²) to O(N)

## [v0.0.1](https://github.com/N283T/zsasa/releases/tag/v0.0.1) — 2025-01-19

### Added

- Initial Shrake-Rupley algorithm implementation
- Golden Section Spiral test point generation
- JSON input/output format
- Basic CLI with input/output file arguments
- Python scripts for structure conversion and benchmarking
  - `cif_to_input_json.py` - Convert mmCIF/PDB to input JSON
  - `calc_reference_sasa.py` - Generate reference SASA
  - `benchmark.py` - Performance benchmarking
