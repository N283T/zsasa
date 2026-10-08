//! C ABI interface for zsasa library.
//!
//! This module provides C-compatible functions that can be called from
//! other languages (Python, C, etc.) via FFI/ctypes.

const std = @import("std");
const shrake_rupley = @import("shrake_rupley.zig");
const shrake_rupley_bitmask = @import("shrake_rupley_bitmask.zig");
const bitmask_lut = @import("bitmask_lut.zig");
const lee_richards = @import("lee_richards.zig");
const types = @import("types.zig");
const classifier = @import("classifier.zig");
const classifier_naccess = @import("classifier_naccess.zig");
const classifier_oons = @import("classifier_oons.zig");
const classifier_ccd = @import("classifier_ccd.zig");
const analysis = @import("analysis.zig");
const ztraj = @import("ztraj");
const xtc = ztraj.io.xtc;
const dcd = ztraj.io.dcd;
const batch = @import("batch.zig");

const AtomInput = types.AtomInput;
const Config = types.Config;

/// No error - calculation completed successfully
pub const ZSASA_OK: c_int = 0;
/// Invalid input parameters (n_atoms=0, n_points=0, n_slices=0, invalid probe_radius,
/// non-finite values, or a coordinate range too wide for the calculation precision)
/// Note: Passing NULL pointers results in undefined behavior
pub const ZSASA_ERROR_INVALID_INPUT: c_int = -1;
/// Memory allocation failed during calculation
pub const ZSASA_ERROR_OUT_OF_MEMORY: c_int = -2;
/// Internal calculation error
pub const ZSASA_ERROR_CALCULATION: c_int = -3;
/// Directory/file I/O error (directory open failure)
pub const ZSASA_ERROR_FILE_IO: c_int = -4;
/// Unsupported n_points value for bitmask algorithm (must be 1..1024)
pub const ZSASA_ERROR_UNSUPPORTED_N_POINTS: c_int = -5;
/// Several inputs map to the same per-file output name in the output directory
pub const ZSASA_ERROR_OUTPUT_NAME_COLLISION: c_int = -6;
/// An input file exists but is not valid for the format it was opened as: another
/// format, an empty or truncated file, or corrupt data (trajectory readers)
pub const ZSASA_ERROR_INVALID_FORMAT: c_int = -7;
/// The output directory cannot be created or is not a directory
/// (`zsasa_batch_dir_process`); the input directory is reported as ZSASA_ERROR_FILE_IO
pub const ZSASA_ERROR_OUTPUT_DIR: c_int = -8;

// =============================================================================
// Algorithm Constants
// =============================================================================

/// Shrake-Rupley algorithm (test point method)
pub const ZSASA_ALGORITHM_SR: c_int = 0;
/// Lee-Richards algorithm (slice method)
pub const ZSASA_ALGORITHM_LR: c_int = 1;

// =============================================================================
// Classifier Types
// =============================================================================

/// NACCESS-compatible radii
pub const ZSASA_CLASSIFIER_NACCESS: c_int = 0;
/// Static ProtOr-compatible radii
pub const ZSASA_CLASSIFIER_PROTOR: c_int = 1;
/// OONS radii (older FreeSASA default)
pub const ZSASA_CLASSIFIER_OONS: c_int = 2;
/// CCD classifier (default). The CLI and `zsasa_batch_dir_process` derive radii from
/// bond topology for components that are not in the built-in table. The classify
/// functions (`zsasa_classifier_get_radius`, `zsasa_classifier_get_class`,
/// `zsasa_classify_atoms`) have no way to load components, so there it returns the
/// same built-in ProtOr radii as ZSASA_CLASSIFIER_PROTOR.
pub const ZSASA_CLASSIFIER_CCD: c_int = 3;

// =============================================================================
// Atom Classes
// =============================================================================

/// Polar atom class
pub const ZSASA_ATOM_CLASS_POLAR: c_int = 0;
/// Apolar atom class
pub const ZSASA_ATOM_CLASS_APOLAR: c_int = 1;
/// Unknown atom class
pub const ZSASA_ATOM_CLASS_UNKNOWN: c_int = 2;

// Version string
const VERSION = "0.9.1";

/// Version of the C ABI, returned by `zsasa_abi_version()`.
///
/// The Python bindings declare the signatures of this library by hand and refuse
/// to load a library whose ABI version differs from the one they expect, because
/// calling through a stale declaration is undefined behavior.
///
/// Bump it whenever an existing exported function changes its signature (argument
/// or return types, order or count) or an exported struct changes its layout, and
/// change `_EXPECTED_ABI_VERSION` in `python/zsasa/_ffi.py` together with it.
/// Adding a new export or a new error code does not require a bump.
const ABI_VERSION: c_int = 1;

/// Thread-safe allocator for C API (uses C allocator for simplicity)
const c_allocator = std.heap.c_allocator;

/// Returns a single-threaded Io for FFI entries that do not spawn threads.
/// Multi-threaded FFI entries (e.g., zsasa_batch_dir_process) MUST use
/// `sharedThreadedIo()` instead of constructing their own std.Io.Threaded.
fn cIo() std.Io {
    return std.Io.Threaded.global_single_threaded.io();
}

// The process-wide multi-threaded Io. `std.Io.Threaded.init` replaces the
// process's SIGIO and SIGPIPE handlers with a no-op one and `deinit` puts back
// whatever `init` found. Two overlapping instances (concurrent calls from several
// Python threads) restore each other's no-op handler instead of the original one.
// So the instance is created once, on first use, and never deinitialized: the
// handlers are replaced at most once per process and nothing is restored while
// another call may still be running. The worker threads Threaded starts on demand
// stay parked until the process exits.
const shared_io_uninitialized: u8 = 0;
const shared_io_initializing: u8 = 1;
const shared_io_ready: u8 = 2;
var shared_io_state = std.atomic.Value(u8).init(shared_io_uninitialized);
var shared_threaded: std.Io.Threaded = undefined;

/// Returns the process-wide multi-threaded Io, creating it on the first call.
/// Thread-safe: concurrent first calls create exactly one instance; the others
/// wait for it.
fn sharedThreadedIo() std.Io {
    while (true) {
        switch (shared_io_state.load(.acquire)) {
            shared_io_ready => return shared_threaded.io(),
            shared_io_uninitialized => {
                if (shared_io_state.cmpxchgStrong(shared_io_uninitialized, shared_io_initializing, .acquire, .monotonic) == null) {
                    shared_threaded = std.Io.Threaded.init(c_allocator, .{});
                    shared_io_state.store(shared_io_ready, .release);
                    return shared_threaded.io();
                }
            },
            else => std.Thread.yield() catch {},
        }
    }
}

fn isPositiveFinite(comptime T: type, value: T) bool {
    return std.math.isFinite(value) and value > 0.0;
}

fn isNonNegativeFinite(comptime T: type, value: T) bool {
    return std.math.isFinite(value) and value >= 0.0;
}

fn validateAtomArraysF64(
    x: [*]const f64,
    y: [*]const f64,
    z: [*]const f64,
    radii: [*]const f64,
    n_atoms: usize,
) bool {
    for (0..n_atoms) |i| {
        if (!std.math.isFinite(x[i]) or
            !std.math.isFinite(y[i]) or
            !std.math.isFinite(z[i]) or
            !isNonNegativeFinite(f64, radii[i]))
        {
            return false;
        }
    }
    return true;
}

fn validateBatchArraysF32(
    coordinates: [*]const f32,
    n_frames: usize,
    n_atoms: usize,
    radii: [*]const f32,
) bool {
    const frame_atoms = std.math.mul(usize, n_frames, n_atoms) catch return false;
    const coordinate_count = std.math.mul(usize, frame_atoms, 3) catch return false;

    for (coordinates[0..coordinate_count]) |coordinate| {
        if (!std.math.isFinite(coordinate)) return false;
    }
    for (radii[0..n_atoms]) |radius| {
        if (!isNonNegativeFinite(f32, radius)) return false;
    }
    return true;
}

/// Map an error returned by a SASA calculation to a C API error code.
fn calcErrorCode(err: anyerror) c_int {
    return switch (err) {
        error.OutOfMemory => ZSASA_ERROR_OUT_OF_MEMORY,
        // Finite coordinates whose range the calculation precision cannot represent
        error.CoordinateRangeTooLarge => ZSASA_ERROR_INVALID_INPUT,
        else => ZSASA_ERROR_CALCULATION,
    };
}

/// Map an error from opening or reading a trajectory file to a C API error code.
///
/// A missing file stays ZSASA_ERROR_INVALID_INPUT (the code the open functions
/// have always returned for it). A file that exists but cannot be parsed as the
/// requested format (another format, empty, truncated, corrupt) is
/// ZSASA_ERROR_INVALID_FORMAT, so that callers can tell it from a missing file.
fn trajectoryErrorCode(err: anyerror) c_int {
    return switch (err) {
        error.FileNotFound => ZSASA_ERROR_INVALID_INPUT,
        error.InvalidMagic,
        error.BadFormat,
        error.EndOfFile, // the header of an empty file; `next()` turns a frame-boundary end into null
        error.ReadError,
        error.DecompressionError,
        error.FixedAtomsNotSupported,
        => ZSASA_ERROR_INVALID_FORMAT,
        error.OutOfMemory => ZSASA_ERROR_OUT_OF_MEMORY,
        else => ZSASA_ERROR_CALCULATION,
    };
}

/// Get library version string.
export fn zsasa_version() callconv(.c) [*:0]const u8 {
    return VERSION;
}

/// Get the C ABI version (see `ABI_VERSION` for when it changes).
export fn zsasa_abi_version() callconv(.c) c_int {
    return ABI_VERSION;
}

/// Calculate SASA using Shrake-Rupley algorithm.
///
/// Parameters:
///   x, y, z: Atom coordinates (arrays of n_atoms elements)
///   radii: Atom radii (array of n_atoms elements)
///   n_atoms: Number of atoms
///   n_points: Number of test points per atom (e.g., 100)
///   probe_radius: Water probe radius in Angstroms (e.g., 1.4)
///   n_threads: Number of threads (0 = auto-detect)
///   atom_areas: Output buffer for per-atom SASA (must be pre-allocated, n_atoms elements)
///   total_area: Output pointer for total SASA
///
/// Returns:
///   ZSASA_OK (0) on success, negative error code on failure.
export fn zsasa_calc_sr(
    x: [*]const f64,
    y: [*]const f64,
    z: [*]const f64,
    radii: [*]const f64,
    n_atoms: usize,
    n_points: u32,
    probe_radius: f64,
    n_threads: usize,
    atom_areas: [*]f64,
    total_area: *f64,
) callconv(.c) c_int {
    // Validate input
    if (n_atoms == 0 or n_points == 0 or !isPositiveFinite(f64, probe_radius) or
        !validateAtomArraysF64(x, y, z, radii, n_atoms))
    {
        return ZSASA_ERROR_INVALID_INPUT;
    }

    // Duplicate radii (AtomInput.r requires []f64 for classifier support)
    const r_copy = c_allocator.dupe(f64, radii[0..n_atoms]) catch {
        return ZSASA_ERROR_OUT_OF_MEMORY;
    };
    defer c_allocator.free(r_copy);

    // Create AtomInput from raw arrays
    const input = AtomInput{
        .x = x[0..n_atoms],
        .y = y[0..n_atoms],
        .z = z[0..n_atoms],
        .r = r_copy,
        .allocator = c_allocator,
    };

    const config = Config{
        .n_points = n_points,
        .probe_radius = probe_radius,
    };

    // Calculate SASA
    const result = if (n_threads == 1)
        shrake_rupley.calculateSasa(c_allocator, input, config) catch |err| {
            return calcErrorCode(err);
        }
    else
        shrake_rupley.calculateSasaParallel(c_allocator, input, config, n_threads) catch |err| {
            return calcErrorCode(err);
        };
    defer {
        // Free the result's internal allocation
        c_allocator.free(result.atom_areas);
    }

    // Copy results to output buffers
    @memcpy(atom_areas[0..n_atoms], result.atom_areas);
    total_area.* = result.total_area;

    return ZSASA_OK;
}

/// Calculate SASA using Shrake-Rupley bitmask algorithm.
///
/// Uses a bitmask lookup-table optimization for faster SASA calculation.
/// Supports n_points values of 1..1024.
///
/// Parameters:
///   x, y, z: Atom coordinates (arrays of n_atoms elements)
///   radii: Atom radii (array of n_atoms elements)
///   n_atoms: Number of atoms
///   n_points: Number of test points per atom (must be 1..1024)
///   probe_radius: Water probe radius in Angstroms (e.g., 1.4)
///   n_threads: Number of threads (0 = auto-detect)
///   atom_areas: Output buffer for per-atom SASA (must be pre-allocated, n_atoms elements)
///   total_area: Output pointer for total SASA
///
/// Returns:
///   ZSASA_OK (0) on success, negative error code on failure.
export fn zsasa_calc_sr_bitmask(
    x: [*]const f64,
    y: [*]const f64,
    z: [*]const f64,
    radii: [*]const f64,
    n_atoms: usize,
    n_points: u32,
    probe_radius: f64,
    n_threads: usize,
    atom_areas: [*]f64,
    total_area: *f64,
) callconv(.c) c_int {
    if (n_atoms == 0 or n_points == 0 or !isPositiveFinite(f64, probe_radius) or
        !validateAtomArraysF64(x, y, z, radii, n_atoms))
    {
        return ZSASA_ERROR_INVALID_INPUT;
    }

    if (!bitmask_lut.isSupportedNPoints(n_points)) {
        return ZSASA_ERROR_UNSUPPORTED_N_POINTS;
    }

    const r_copy = c_allocator.dupe(f64, radii[0..n_atoms]) catch {
        return ZSASA_ERROR_OUT_OF_MEMORY;
    };
    defer c_allocator.free(r_copy);

    const input = AtomInput{
        .x = x[0..n_atoms],
        .y = y[0..n_atoms],
        .z = z[0..n_atoms],
        .r = r_copy,
        .allocator = c_allocator,
    };

    const config = Config{
        .n_points = n_points,
        .probe_radius = probe_radius,
    };

    const result = if (n_threads == 1)
        shrake_rupley_bitmask.calculateSasa(c_allocator, input, config) catch |err| {
            return calcErrorCode(err);
        }
    else
        shrake_rupley_bitmask.calculateSasaParallel(c_allocator, input, config, n_threads) catch |err| {
            return calcErrorCode(err);
        };
    defer c_allocator.free(result.atom_areas);

    @memcpy(atom_areas[0..n_atoms], result.atom_areas);
    total_area.* = result.total_area;

    return ZSASA_OK;
}

/// Calculate SASA using Shrake-Rupley bitmask algorithm with experimental correction.
export fn zsasa_calc_sr_bitmask_corrected(
    x: [*]const f64,
    y: [*]const f64,
    z: [*]const f64,
    radii: [*]const f64,
    n_atoms: usize,
    n_points: u32,
    probe_radius: f64,
    n_threads: usize,
    correction_coeff: f64,
    atom_areas: [*]f64,
    total_area: *f64,
) callconv(.c) c_int {
    if (n_atoms == 0 or n_points == 0 or !isPositiveFinite(f64, probe_radius) or
        correction_coeff < 0.0 or !std.math.isFinite(correction_coeff) or
        !validateAtomArraysF64(x, y, z, radii, n_atoms))
    {
        return ZSASA_ERROR_INVALID_INPUT;
    }

    if (!bitmask_lut.isSupportedNPoints(n_points)) {
        return ZSASA_ERROR_UNSUPPORTED_N_POINTS;
    }

    const r_copy = c_allocator.dupe(f64, radii[0..n_atoms]) catch {
        return ZSASA_ERROR_OUT_OF_MEMORY;
    };
    defer c_allocator.free(r_copy);

    const input = AtomInput{
        .x = x[0..n_atoms],
        .y = y[0..n_atoms],
        .z = z[0..n_atoms],
        .r = r_copy,
        .allocator = c_allocator,
    };

    const config = Config{
        .n_points = n_points,
        .probe_radius = probe_radius,
    };
    const correction = shrake_rupley_bitmask.BitmaskCorrectionGen(f64){
        .enabled = true,
        .coeff = correction_coeff,
    };

    const SRBitmask = shrake_rupley_bitmask.ShrakeRupleyBitmaskGen(f64);
    const result = if (n_threads == 1)
        SRBitmask.calculateSasaWithCorrection(c_allocator, input, config, correction) catch |err| {
            return calcErrorCode(err);
        }
    else
        SRBitmask.calculateSasaParallelWithCorrection(c_allocator, input, config, n_threads, correction) catch |err| {
            return calcErrorCode(err);
        };
    defer c_allocator.free(result.atom_areas);

    @memcpy(atom_areas[0..n_atoms], result.atom_areas);
    total_area.* = result.total_area;

    return ZSASA_OK;
}

// =============================================================================
// Batch Processing Infrastructure
// =============================================================================

/// Algorithm type for batch processing
const BatchAlgorithm = enum {
    shrake_rupley,
    lee_richards,
};

/// Common batch processing arguments
const BatchWorkerArgs = struct {
    coordinates: [*]const f32,
    n_atoms: usize,
    n_frames: usize,
    radii_f64: []f64,
    param: u32, // n_points for SR, n_slices for LR
    probe_radius: f64,
    atom_areas: [*]f32,
    error_code: *std.atomic.Value(c_int),
    thread_id: usize,
    n_threads: usize,
    algorithm: BatchAlgorithm,
};

/// Worker function for batch processing (shared by SR and LR)
fn batchWorkerFn(args: BatchWorkerArgs) void {
    // Pre-allocate coordinate buffers once per thread (reused across frames)
    const x = c_allocator.alloc(f64, args.n_atoms) catch {
        args.error_code.store(ZSASA_ERROR_OUT_OF_MEMORY, .release);
        return;
    };
    defer c_allocator.free(x);

    const y = c_allocator.alloc(f64, args.n_atoms) catch {
        args.error_code.store(ZSASA_ERROR_OUT_OF_MEMORY, .release);
        return;
    };
    defer c_allocator.free(y);

    const z = c_allocator.alloc(f64, args.n_atoms) catch {
        args.error_code.store(ZSASA_ERROR_OUT_OF_MEMORY, .release);
        return;
    };
    defer c_allocator.free(z);

    // Each thread processes frames: thread_id, thread_id + n_threads, ...
    var frame_idx = args.thread_id;
    while (frame_idx < args.n_frames) : (frame_idx += args.n_threads) {
        // Skip if error already occurred
        if (args.error_code.load(.acquire) != ZSASA_OK) return;

        const frame_offset = frame_idx * args.n_atoms * 3;
        const output_offset = frame_idx * args.n_atoms;

        // Convert f32 coordinates to f64 and split into x, y, z (reuse buffers)
        for (0..args.n_atoms) |i| {
            x[i] = @floatCast(args.coordinates[frame_offset + i * 3]);
            y[i] = @floatCast(args.coordinates[frame_offset + i * 3 + 1]);
            z[i] = @floatCast(args.coordinates[frame_offset + i * 3 + 2]);
        }

        // Create AtomInput for this frame
        const input = AtomInput{
            .x = x,
            .y = y,
            .z = z,
            .r = args.radii_f64,
            .allocator = c_allocator,
        };

        // Calculate SASA using the specified algorithm
        const result = switch (args.algorithm) {
            .shrake_rupley => blk: {
                const config = Config{
                    .n_points = args.param,
                    .probe_radius = args.probe_radius,
                };
                break :blk shrake_rupley.calculateSasa(c_allocator, input, config) catch |err| {
                    args.error_code.store(calcErrorCode(err), .release);
                    return;
                };
            },
            .lee_richards => blk: {
                const config = lee_richards.LeeRichardsConfig{
                    .n_slices = args.param,
                    .probe_radius = args.probe_radius,
                };
                break :blk lee_richards.calculateSasa(c_allocator, input, config) catch |err| {
                    args.error_code.store(calcErrorCode(err), .release);
                    return;
                };
            },
        };
        defer c_allocator.free(result.atom_areas);

        // Copy results to output buffer (convert f64 to f32)
        for (0..args.n_atoms) |i| {
            args.atom_areas[output_offset + i] = @floatCast(result.atom_areas[i]);
        }
    }
}

/// Common batch calculation logic
fn calculateBatch(
    coordinates: [*]const f32,
    n_frames: usize,
    n_atoms: usize,
    radii: [*]const f32,
    param: u32, // n_points for SR, n_slices for LR
    probe_radius: f32,
    n_threads: usize,
    atom_areas: [*]f32,
    algorithm: BatchAlgorithm,
) c_int {
    // Validate input
    if (n_frames == 0 or n_atoms == 0 or param == 0 or !isPositiveFinite(f32, probe_radius) or
        !validateBatchArraysF32(coordinates, n_frames, n_atoms, radii))
    {
        return ZSASA_ERROR_INVALID_INPUT;
    }

    // Determine actual thread count
    const actual_threads = if (n_threads == 0)
        @as(usize, @intCast(std.Thread.getCpuCount() catch 1))
    else
        n_threads;

    // Convert radii from f32 to f64 (done once, reused for all frames)
    const radii_f64 = c_allocator.alloc(f64, n_atoms) catch {
        return ZSASA_ERROR_OUT_OF_MEMORY;
    };
    defer c_allocator.free(radii_f64);

    for (0..n_atoms) |i| {
        radii_f64[i] = @floatCast(radii[i]);
    }

    // First error code reported by a worker (ZSASA_OK while none failed)
    var error_code = std.atomic.Value(c_int).init(ZSASA_OK);

    // Spawn worker threads
    const thread_count = @min(actual_threads, n_frames);
    const threads = c_allocator.alloc(std.Thread, thread_count) catch {
        return ZSASA_ERROR_OUT_OF_MEMORY;
    };
    defer c_allocator.free(threads);

    for (0..thread_count) |i| {
        threads[i] = std.Thread.spawn(.{}, batchWorkerFn, .{BatchWorkerArgs{
            .coordinates = coordinates,
            .n_atoms = n_atoms,
            .n_frames = n_frames,
            .radii_f64 = radii_f64,
            .param = param,
            .probe_radius = @floatCast(probe_radius),
            .atom_areas = atom_areas,
            .error_code = &error_code,
            .thread_id = i,
            .n_threads = thread_count,
            .algorithm = algorithm,
        }}) catch {
            // If thread spawn fails, set the error code and wait for already-spawned threads
            error_code.store(ZSASA_ERROR_CALCULATION, .release);
            for (threads[0..i]) |thread| {
                thread.join();
            }
            return ZSASA_ERROR_CALCULATION;
        };
    }

    // Wait for all threads to complete
    for (threads) |thread| {
        thread.join();
    }

    const worker_error = error_code.load(.acquire);
    if (worker_error != ZSASA_OK) {
        return worker_error;
    }

    return ZSASA_OK;
}

/// Calculate SASA for multiple frames using Shrake-Rupley algorithm (batch processing).
///
/// This function is optimized for MD trajectory analysis where the same atoms
/// are processed across multiple frames. It reuses buffers and parallelizes
/// across frames for maximum performance.
///
/// Parameters:
///   coordinates: Atom coordinates as contiguous array (n_frames * n_atoms * 3).
///                Layout: [frame0_atom0_xyz, frame0_atom1_xyz, ..., frame1_atom0_xyz, ...]
///                Units: Angstroms (NOT nm)
///   n_frames: Number of frames
///   n_atoms: Number of atoms per frame
///   radii: Atom radii (array of n_atoms elements, reused for all frames)
///   n_points: Number of test points per atom (e.g., 100)
///   probe_radius: Water probe radius in Angstroms (e.g., 1.4)
///   n_threads: Number of threads (0 = auto-detect)
///   atom_areas: Output buffer for per-atom SASA (n_frames * n_atoms elements)
///               Layout: [frame0_atom0, frame0_atom1, ..., frame1_atom0, ...]
///
/// Returns:
///   ZSASA_OK (0) on success, negative error code on failure.
export fn zsasa_calc_sr_batch(
    coordinates: [*]const f32,
    n_frames: usize,
    n_atoms: usize,
    radii: [*]const f32,
    n_points: u32,
    probe_radius: f32,
    n_threads: usize,
    atom_areas: [*]f32,
) callconv(.c) c_int {
    return calculateBatch(
        coordinates,
        n_frames,
        n_atoms,
        radii,
        n_points,
        probe_radius,
        n_threads,
        atom_areas,
        .shrake_rupley,
    );
}

/// Calculate SASA for multiple frames using Lee-Richards algorithm (batch processing).
///
/// Similar to zsasa_calc_sr_batch but uses Lee-Richards algorithm.
/// Arc angles are computed with exact trigonometry (`lee_richards.TrigMode.exact`).
///
/// Parameters:
///   coordinates: Atom coordinates as contiguous array (n_frames * n_atoms * 3)
///   n_frames: Number of frames
///   n_atoms: Number of atoms per frame
///   radii: Atom radii (array of n_atoms elements)
///   n_slices: Number of slices per atom (e.g., 20)
///   probe_radius: Water probe radius in Angstroms (e.g., 1.4)
///   n_threads: Number of threads (0 = auto-detect)
///   atom_areas: Output buffer for per-atom SASA (n_frames * n_atoms elements)
///
/// Returns:
///   ZSASA_OK (0) on success, negative error code on failure.
export fn zsasa_calc_lr_batch(
    coordinates: [*]const f32,
    n_frames: usize,
    n_atoms: usize,
    radii: [*]const f32,
    n_slices: u32,
    probe_radius: f32,
    n_threads: usize,
    atom_areas: [*]f32,
) callconv(.c) c_int {
    return calculateBatch(
        coordinates,
        n_frames,
        n_atoms,
        radii,
        n_slices,
        probe_radius,
        n_threads,
        atom_areas,
        .lee_richards,
    );
}

// =============================================================================
// Batch Processing Infrastructure (Pure f32 precision)
// =============================================================================

/// Worker arguments for pure f32 batch processing
const BatchWorkerArgsF32 = struct {
    coordinates: [*]const f32,
    n_atoms: usize,
    n_frames: usize,
    radii_f32: []f32,
    param: u32, // n_points for SR, n_slices for LR
    probe_radius: f32,
    atom_areas: [*]f32,
    error_code: *std.atomic.Value(c_int),
    thread_id: usize,
    n_threads: usize,
    algorithm: BatchAlgorithm,
};

/// Worker function for pure f32 batch processing
fn batchWorkerFnF32(args: BatchWorkerArgsF32) void {
    // Pre-allocate coordinate buffers (f64 for AtomInput compatibility)
    const x = c_allocator.alloc(f64, args.n_atoms) catch {
        args.error_code.store(ZSASA_ERROR_OUT_OF_MEMORY, .release);
        return;
    };
    defer c_allocator.free(x);

    const y = c_allocator.alloc(f64, args.n_atoms) catch {
        args.error_code.store(ZSASA_ERROR_OUT_OF_MEMORY, .release);
        return;
    };
    defer c_allocator.free(y);

    const z = c_allocator.alloc(f64, args.n_atoms) catch {
        args.error_code.store(ZSASA_ERROR_OUT_OF_MEMORY, .release);
        return;
    };
    defer c_allocator.free(z);

    // f64 radii for AtomInput (will be cast to f32 internally)
    const radii_f64 = c_allocator.alloc(f64, args.n_atoms) catch {
        args.error_code.store(ZSASA_ERROR_OUT_OF_MEMORY, .release);
        return;
    };
    defer c_allocator.free(radii_f64);

    for (0..args.n_atoms) |i| {
        radii_f64[i] = @floatCast(args.radii_f32[i]);
    }

    // Each thread processes frames: thread_id, thread_id + n_threads, ...
    var frame_idx = args.thread_id;
    while (frame_idx < args.n_frames) : (frame_idx += args.n_threads) {
        if (args.error_code.load(.acquire) != ZSASA_OK) return;

        const frame_offset = frame_idx * args.n_atoms * 3;
        const output_offset = frame_idx * args.n_atoms;

        // Convert f32 coordinates to f64 for AtomInput
        for (0..args.n_atoms) |i| {
            x[i] = @floatCast(args.coordinates[frame_offset + i * 3]);
            y[i] = @floatCast(args.coordinates[frame_offset + i * 3 + 1]);
            z[i] = @floatCast(args.coordinates[frame_offset + i * 3 + 2]);
        }

        const input = AtomInput{
            .x = x,
            .y = y,
            .z = z,
            .r = radii_f64,
            .allocator = c_allocator,
        };

        // Calculate SASA using f32 precision algorithm
        const result = switch (args.algorithm) {
            .shrake_rupley => blk: {
                const config = types.ConfigGen(f32){
                    .n_points = args.param,
                    .probe_radius = args.probe_radius,
                };
                break :blk shrake_rupley.calculateSasaf32(c_allocator, input, config) catch |err| {
                    args.error_code.store(calcErrorCode(err), .release);
                    return;
                };
            },
            .lee_richards => blk: {
                const config = lee_richards.LeeRichardsConfigGen(f32){
                    .n_slices = args.param,
                    .probe_radius = args.probe_radius,
                };
                break :blk lee_richards.calculateSasaf32(c_allocator, input, config) catch |err| {
                    args.error_code.store(calcErrorCode(err), .release);
                    return;
                };
            },
        };
        defer c_allocator.free(result.atom_areas);

        // Copy f32 results directly (no conversion needed)
        for (0..args.n_atoms) |i| {
            args.atom_areas[output_offset + i] = result.atom_areas[i];
        }
    }
}

/// Common batch calculation logic for pure f32 precision
fn calculateBatchF32(
    coordinates: [*]const f32,
    n_frames: usize,
    n_atoms: usize,
    radii: [*]const f32,
    param: u32,
    probe_radius: f32,
    n_threads: usize,
    atom_areas: [*]f32,
    algorithm: BatchAlgorithm,
) c_int {
    if (n_frames == 0 or n_atoms == 0 or param == 0 or !isPositiveFinite(f32, probe_radius) or
        !validateBatchArraysF32(coordinates, n_frames, n_atoms, radii))
    {
        return ZSASA_ERROR_INVALID_INPUT;
    }

    const actual_threads = if (n_threads == 0)
        @as(usize, @intCast(std.Thread.getCpuCount() catch 1))
    else
        n_threads;

    // Copy radii (kept as f32)
    const radii_f32 = c_allocator.alloc(f32, n_atoms) catch {
        return ZSASA_ERROR_OUT_OF_MEMORY;
    };
    defer c_allocator.free(radii_f32);

    for (0..n_atoms) |i| {
        radii_f32[i] = radii[i];
    }

    var error_code = std.atomic.Value(c_int).init(ZSASA_OK);

    const thread_count = @min(actual_threads, n_frames);
    const threads = c_allocator.alloc(std.Thread, thread_count) catch {
        return ZSASA_ERROR_OUT_OF_MEMORY;
    };
    defer c_allocator.free(threads);

    for (0..thread_count) |i| {
        threads[i] = std.Thread.spawn(.{}, batchWorkerFnF32, .{BatchWorkerArgsF32{
            .coordinates = coordinates,
            .n_atoms = n_atoms,
            .n_frames = n_frames,
            .radii_f32 = radii_f32,
            .param = param,
            .probe_radius = probe_radius,
            .atom_areas = atom_areas,
            .error_code = &error_code,
            .thread_id = i,
            .n_threads = thread_count,
            .algorithm = algorithm,
        }}) catch {
            error_code.store(ZSASA_ERROR_CALCULATION, .release);
            for (threads[0..i]) |thread| {
                thread.join();
            }
            return ZSASA_ERROR_CALCULATION;
        };
    }

    for (threads) |thread| {
        thread.join();
    }

    const worker_error = error_code.load(.acquire);
    if (worker_error != ZSASA_OK) {
        return worker_error;
    }

    return ZSASA_OK;
}

/// Calculate SASA for multiple frames using Shrake-Rupley algorithm (pure f32 precision).
///
/// Same as zsasa_calc_sr_batch but uses f32 precision throughout the calculation
/// for consistency with other f32-based tools (e.g., RustSASA).
export fn zsasa_calc_sr_batch_f32(
    coordinates: [*]const f32,
    n_frames: usize,
    n_atoms: usize,
    radii: [*]const f32,
    n_points: u32,
    probe_radius: f32,
    n_threads: usize,
    atom_areas: [*]f32,
) callconv(.c) c_int {
    return calculateBatchF32(
        coordinates,
        n_frames,
        n_atoms,
        radii,
        n_points,
        probe_radius,
        n_threads,
        atom_areas,
        .shrake_rupley,
    );
}

/// Calculate SASA for multiple frames using Lee-Richards algorithm (pure f32 precision).
///
/// Same as zsasa_calc_lr_batch but uses f32 precision throughout the calculation.
/// Arc angles are computed with exact trigonometry (`lee_richards.TrigMode.exact`).
export fn zsasa_calc_lr_batch_f32(
    coordinates: [*]const f32,
    n_frames: usize,
    n_atoms: usize,
    radii: [*]const f32,
    n_slices: u32,
    probe_radius: f32,
    n_threads: usize,
    atom_areas: [*]f32,
) callconv(.c) c_int {
    return calculateBatchF32(
        coordinates,
        n_frames,
        n_atoms,
        radii,
        n_slices,
        probe_radius,
        n_threads,
        atom_areas,
        .lee_richards,
    );
}

// =============================================================================
// Batch Processing Infrastructure (Bitmask LUT)
// =============================================================================

/// Worker arguments for bitmask batch processing (f64 internal precision)
const BatchWorkerArgsBitmask = struct {
    coordinates: [*]const f32,
    n_atoms: usize,
    n_frames: usize,
    radii_f64: []f64,
    n_points: u32,
    probe_radius: f64,
    atom_areas: [*]f32,
    error_code: *std.atomic.Value(c_int),
    thread_id: usize,
    n_threads: usize,
    lut: *const bitmask_lut.BitmaskLut,
    correction: shrake_rupley_bitmask.BitmaskCorrectionGen(f64) = .{},
};

/// Worker function for bitmask batch processing (f64 internal precision)
fn batchWorkerFnBitmask(args: BatchWorkerArgsBitmask) void {
    const x = c_allocator.alloc(f64, args.n_atoms) catch {
        args.error_code.store(ZSASA_ERROR_OUT_OF_MEMORY, .release);
        return;
    };
    defer c_allocator.free(x);

    const y = c_allocator.alloc(f64, args.n_atoms) catch {
        args.error_code.store(ZSASA_ERROR_OUT_OF_MEMORY, .release);
        return;
    };
    defer c_allocator.free(y);

    const z = c_allocator.alloc(f64, args.n_atoms) catch {
        args.error_code.store(ZSASA_ERROR_OUT_OF_MEMORY, .release);
        return;
    };
    defer c_allocator.free(z);

    const SRBitmask = shrake_rupley_bitmask.ShrakeRupleyBitmaskGen(f64);

    var frame_idx = args.thread_id;
    while (frame_idx < args.n_frames) : (frame_idx += args.n_threads) {
        if (args.error_code.load(.acquire) != ZSASA_OK) return;

        const frame_offset = frame_idx * args.n_atoms * 3;
        const output_offset = frame_idx * args.n_atoms;

        for (0..args.n_atoms) |i| {
            x[i] = @floatCast(args.coordinates[frame_offset + i * 3]);
            y[i] = @floatCast(args.coordinates[frame_offset + i * 3 + 1]);
            z[i] = @floatCast(args.coordinates[frame_offset + i * 3 + 2]);
        }

        const input = AtomInput{
            .x = x,
            .y = y,
            .z = z,
            .r = args.radii_f64,
            .allocator = c_allocator,
        };

        const config = Config{
            .n_points = args.n_points,
            .probe_radius = args.probe_radius,
        };

        const result = SRBitmask.calculateSasaWithLutAndCorrection(c_allocator, input, config, args.lut, args.correction) catch |err| {
            args.error_code.store(calcErrorCode(err), .release);
            return;
        };
        defer c_allocator.free(result.atom_areas);

        for (0..args.n_atoms) |i| {
            args.atom_areas[output_offset + i] = @floatCast(result.atom_areas[i]);
        }
    }
}

/// Common bitmask batch calculation logic (f64 internal precision)
fn calculateBatchBitmask(
    coordinates: [*]const f32,
    n_frames: usize,
    n_atoms: usize,
    radii: [*]const f32,
    n_points: u32,
    probe_radius: f32,
    n_threads: usize,
    correction_enabled: bool,
    correction_coeff: f64,
    atom_areas: [*]f32,
) c_int {
    if (n_frames == 0 or n_atoms == 0 or n_points == 0 or !isPositiveFinite(f32, probe_radius) or
        correction_coeff < 0.0 or !std.math.isFinite(correction_coeff) or
        !validateBatchArraysF32(coordinates, n_frames, n_atoms, radii))
    {
        return ZSASA_ERROR_INVALID_INPUT;
    }

    if (!bitmask_lut.isSupportedNPoints(n_points)) {
        return ZSASA_ERROR_UNSUPPORTED_N_POINTS;
    }

    const actual_threads = if (n_threads == 0)
        @as(usize, @intCast(std.Thread.getCpuCount() catch 1))
    else
        n_threads;

    // Build LUT once (shared across all threads/frames)
    var lut = bitmask_lut.BitmaskLut.init(c_allocator, n_points) catch {
        return ZSASA_ERROR_OUT_OF_MEMORY;
    };
    defer lut.deinit();

    // Convert radii from f32 to f64
    const radii_f64 = c_allocator.alloc(f64, n_atoms) catch {
        return ZSASA_ERROR_OUT_OF_MEMORY;
    };
    defer c_allocator.free(radii_f64);

    for (0..n_atoms) |i| {
        radii_f64[i] = @floatCast(radii[i]);
    }

    var error_code = std.atomic.Value(c_int).init(ZSASA_OK);
    const correction = shrake_rupley_bitmask.BitmaskCorrectionGen(f64){
        .enabled = correction_enabled,
        .coeff = correction_coeff,
    };

    const thread_count = @min(actual_threads, n_frames);
    const threads = c_allocator.alloc(std.Thread, thread_count) catch {
        return ZSASA_ERROR_OUT_OF_MEMORY;
    };
    defer c_allocator.free(threads);

    for (0..thread_count) |i| {
        threads[i] = std.Thread.spawn(.{}, batchWorkerFnBitmask, .{BatchWorkerArgsBitmask{
            .coordinates = coordinates,
            .n_atoms = n_atoms,
            .n_frames = n_frames,
            .radii_f64 = radii_f64,
            .n_points = n_points,
            .probe_radius = @floatCast(probe_radius),
            .atom_areas = atom_areas,
            .error_code = &error_code,
            .thread_id = i,
            .n_threads = thread_count,
            .lut = &lut,
            .correction = correction,
        }}) catch {
            error_code.store(ZSASA_ERROR_CALCULATION, .release);
            for (threads[0..i]) |thread| {
                thread.join();
            }
            return ZSASA_ERROR_CALCULATION;
        };
    }

    for (threads) |thread| {
        thread.join();
    }

    const worker_error = error_code.load(.acquire);
    if (worker_error != ZSASA_OK) {
        return worker_error;
    }

    return ZSASA_OK;
}

/// Worker arguments for bitmask batch processing (f32 internal precision)
const BatchWorkerArgsBitmaskF32 = struct {
    coordinates: [*]const f32,
    n_atoms: usize,
    n_frames: usize,
    radii_f32: []f32,
    n_points: u32,
    probe_radius: f32,
    atom_areas: [*]f32,
    error_code: *std.atomic.Value(c_int),
    thread_id: usize,
    n_threads: usize,
    lut: *const bitmask_lut.BitmaskLutf32,
    correction: shrake_rupley_bitmask.BitmaskCorrectionGen(f32) = .{},
};

/// Worker function for bitmask batch processing (f32 internal precision)
fn batchWorkerFnBitmaskF32(args: BatchWorkerArgsBitmaskF32) void {
    const x = c_allocator.alloc(f64, args.n_atoms) catch {
        args.error_code.store(ZSASA_ERROR_OUT_OF_MEMORY, .release);
        return;
    };
    defer c_allocator.free(x);

    const y = c_allocator.alloc(f64, args.n_atoms) catch {
        args.error_code.store(ZSASA_ERROR_OUT_OF_MEMORY, .release);
        return;
    };
    defer c_allocator.free(y);

    const z = c_allocator.alloc(f64, args.n_atoms) catch {
        args.error_code.store(ZSASA_ERROR_OUT_OF_MEMORY, .release);
        return;
    };
    defer c_allocator.free(z);

    const radii_f64 = c_allocator.alloc(f64, args.n_atoms) catch {
        args.error_code.store(ZSASA_ERROR_OUT_OF_MEMORY, .release);
        return;
    };
    defer c_allocator.free(radii_f64);

    for (0..args.n_atoms) |i| {
        radii_f64[i] = @floatCast(args.radii_f32[i]);
    }

    const SRBitmask = shrake_rupley_bitmask.ShrakeRupleyBitmaskGen(f32);

    var frame_idx = args.thread_id;
    while (frame_idx < args.n_frames) : (frame_idx += args.n_threads) {
        if (args.error_code.load(.acquire) != ZSASA_OK) return;

        const frame_offset = frame_idx * args.n_atoms * 3;
        const output_offset = frame_idx * args.n_atoms;

        for (0..args.n_atoms) |i| {
            x[i] = @floatCast(args.coordinates[frame_offset + i * 3]);
            y[i] = @floatCast(args.coordinates[frame_offset + i * 3 + 1]);
            z[i] = @floatCast(args.coordinates[frame_offset + i * 3 + 2]);
        }

        const input = AtomInput{
            .x = x,
            .y = y,
            .z = z,
            .r = radii_f64,
            .allocator = c_allocator,
        };

        const config = types.ConfigGen(f32){
            .n_points = args.n_points,
            .probe_radius = args.probe_radius,
        };

        const result = SRBitmask.calculateSasaWithLutAndCorrection(c_allocator, input, config, args.lut, args.correction) catch |err| {
            args.error_code.store(calcErrorCode(err), .release);
            return;
        };
        defer c_allocator.free(result.atom_areas);

        for (0..args.n_atoms) |i| {
            args.atom_areas[output_offset + i] = result.atom_areas[i];
        }
    }
}

/// Common bitmask batch calculation logic (f32 internal precision)
fn calculateBatchBitmaskF32(
    coordinates: [*]const f32,
    n_frames: usize,
    n_atoms: usize,
    radii: [*]const f32,
    n_points: u32,
    probe_radius: f32,
    n_threads: usize,
    correction_enabled: bool,
    correction_coeff: f64,
    atom_areas: [*]f32,
) c_int {
    if (n_frames == 0 or n_atoms == 0 or n_points == 0 or !isPositiveFinite(f32, probe_radius) or
        correction_coeff < 0.0 or !std.math.isFinite(correction_coeff) or
        !validateBatchArraysF32(coordinates, n_frames, n_atoms, radii))
    {
        return ZSASA_ERROR_INVALID_INPUT;
    }

    if (!bitmask_lut.isSupportedNPoints(n_points)) {
        return ZSASA_ERROR_UNSUPPORTED_N_POINTS;
    }

    const actual_threads = if (n_threads == 0)
        @as(usize, @intCast(std.Thread.getCpuCount() catch 1))
    else
        n_threads;

    // Build LUT once (shared across all threads/frames)
    var lut = bitmask_lut.BitmaskLutf32.init(c_allocator, n_points) catch {
        return ZSASA_ERROR_OUT_OF_MEMORY;
    };
    defer lut.deinit();

    // Copy radii (kept as f32)
    const radii_f32 = c_allocator.alloc(f32, n_atoms) catch {
        return ZSASA_ERROR_OUT_OF_MEMORY;
    };
    defer c_allocator.free(radii_f32);

    for (0..n_atoms) |i| {
        radii_f32[i] = radii[i];
    }

    var error_code = std.atomic.Value(c_int).init(ZSASA_OK);
    const correction = shrake_rupley_bitmask.BitmaskCorrectionGen(f32){
        .enabled = correction_enabled,
        .coeff = @floatCast(correction_coeff),
    };

    const thread_count = @min(actual_threads, n_frames);
    const threads = c_allocator.alloc(std.Thread, thread_count) catch {
        return ZSASA_ERROR_OUT_OF_MEMORY;
    };
    defer c_allocator.free(threads);

    for (0..thread_count) |i| {
        threads[i] = std.Thread.spawn(.{}, batchWorkerFnBitmaskF32, .{BatchWorkerArgsBitmaskF32{
            .coordinates = coordinates,
            .n_atoms = n_atoms,
            .n_frames = n_frames,
            .radii_f32 = radii_f32,
            .n_points = n_points,
            .probe_radius = probe_radius,
            .atom_areas = atom_areas,
            .error_code = &error_code,
            .thread_id = i,
            .n_threads = thread_count,
            .lut = &lut,
            .correction = correction,
        }}) catch {
            error_code.store(ZSASA_ERROR_CALCULATION, .release);
            for (threads[0..i]) |thread| {
                thread.join();
            }
            return ZSASA_ERROR_CALCULATION;
        };
    }

    for (threads) |thread| {
        thread.join();
    }

    const worker_error = error_code.load(.acquire);
    if (worker_error != ZSASA_OK) {
        return worker_error;
    }

    return ZSASA_OK;
}

/// Calculate SASA for multiple frames using bitmask Shrake-Rupley algorithm.
///
/// Uses a bitmask lookup-table optimization. The LUT is built once and reused
/// across all frames for better performance.
///
/// Parameters:
///   coordinates: Atom coordinates as contiguous array (n_frames * n_atoms * 3)
///   n_frames: Number of frames
///   n_atoms: Number of atoms per frame
///   radii: Atom radii (array of n_atoms elements)
///   n_points: Number of test points per atom (must be 1..1024)
///   probe_radius: Water probe radius in Angstroms (e.g., 1.4)
///   n_threads: Number of threads (0 = auto-detect)
///   atom_areas: Output buffer for per-atom SASA (n_frames * n_atoms elements)
///
/// Returns:
///   ZSASA_OK (0) on success, negative error code on failure.
export fn zsasa_calc_sr_batch_bitmask(
    coordinates: [*]const f32,
    n_frames: usize,
    n_atoms: usize,
    radii: [*]const f32,
    n_points: u32,
    probe_radius: f32,
    n_threads: usize,
    atom_areas: [*]f32,
) callconv(.c) c_int {
    return calculateBatchBitmask(
        coordinates,
        n_frames,
        n_atoms,
        radii,
        n_points,
        probe_radius,
        n_threads,
        false,
        shrake_rupley_bitmask.default_bitmask_correction_coeff,
        atom_areas,
    );
}

/// Calculate SASA for multiple frames using corrected bitmask Shrake-Rupley algorithm.
export fn zsasa_calc_sr_batch_bitmask_corrected(
    coordinates: [*]const f32,
    n_frames: usize,
    n_atoms: usize,
    radii: [*]const f32,
    n_points: u32,
    probe_radius: f32,
    n_threads: usize,
    correction_coeff: f64,
    atom_areas: [*]f32,
) callconv(.c) c_int {
    return calculateBatchBitmask(
        coordinates,
        n_frames,
        n_atoms,
        radii,
        n_points,
        probe_radius,
        n_threads,
        true,
        correction_coeff,
        atom_areas,
    );
}

/// Calculate SASA for multiple frames using bitmask Shrake-Rupley algorithm (f32 precision).
///
/// Same as zsasa_calc_sr_batch_bitmask but uses f32 precision throughout.
export fn zsasa_calc_sr_batch_bitmask_f32(
    coordinates: [*]const f32,
    n_frames: usize,
    n_atoms: usize,
    radii: [*]const f32,
    n_points: u32,
    probe_radius: f32,
    n_threads: usize,
    atom_areas: [*]f32,
) callconv(.c) c_int {
    return calculateBatchBitmaskF32(
        coordinates,
        n_frames,
        n_atoms,
        radii,
        n_points,
        probe_radius,
        n_threads,
        false,
        shrake_rupley_bitmask.default_bitmask_correction_coeff,
        atom_areas,
    );
}

/// Calculate SASA for multiple frames using corrected bitmask Shrake-Rupley algorithm (f32 precision).
export fn zsasa_calc_sr_batch_bitmask_f32_corrected(
    coordinates: [*]const f32,
    n_frames: usize,
    n_atoms: usize,
    radii: [*]const f32,
    n_points: u32,
    probe_radius: f32,
    n_threads: usize,
    correction_coeff: f64,
    atom_areas: [*]f32,
) callconv(.c) c_int {
    return calculateBatchBitmaskF32(
        coordinates,
        n_frames,
        n_atoms,
        radii,
        n_points,
        probe_radius,
        n_threads,
        true,
        correction_coeff,
        atom_areas,
    );
}

/// Calculate SASA using Lee-Richards algorithm.
///
/// Arc angles are computed with exact trigonometry (`lee_richards.TrigMode.exact`).
/// The approximation that was used up to zsasa 0.9.1, `TrigMode.fast`, is a
/// command-line option only and is not available through the C API.
///
/// Parameters:
///   x, y, z: Atom coordinates (arrays of n_atoms elements)
///   radii: Atom radii (array of n_atoms elements)
///   n_atoms: Number of atoms
///   n_slices: Number of slices per atom (e.g., 20)
///   probe_radius: Water probe radius in Angstroms (e.g., 1.4)
///   n_threads: Number of threads (0 = auto-detect)
///   atom_areas: Output buffer for per-atom SASA (must be pre-allocated, n_atoms elements)
///   total_area: Output pointer for total SASA
///
/// Returns:
///   ZSASA_OK (0) on success, negative error code on failure.
export fn zsasa_calc_lr(
    x: [*]const f64,
    y: [*]const f64,
    z: [*]const f64,
    radii: [*]const f64,
    n_atoms: usize,
    n_slices: u32,
    probe_radius: f64,
    n_threads: usize,
    atom_areas: [*]f64,
    total_area: *f64,
) callconv(.c) c_int {
    // Validate input
    if (n_atoms == 0 or n_slices == 0 or !isPositiveFinite(f64, probe_radius) or
        !validateAtomArraysF64(x, y, z, radii, n_atoms))
    {
        return ZSASA_ERROR_INVALID_INPUT;
    }

    // Duplicate radii (AtomInput.r requires []f64 for classifier support)
    const r_copy = c_allocator.dupe(f64, radii[0..n_atoms]) catch {
        return ZSASA_ERROR_OUT_OF_MEMORY;
    };
    defer c_allocator.free(r_copy);

    // Create AtomInput from raw arrays
    const input = AtomInput{
        .x = x[0..n_atoms],
        .y = y[0..n_atoms],
        .z = z[0..n_atoms],
        .r = r_copy,
        .allocator = c_allocator,
    };

    const config = lee_richards.LeeRichardsConfig{
        .n_slices = n_slices,
        .probe_radius = probe_radius,
    };

    // Calculate SASA
    const result = if (n_threads == 1)
        lee_richards.calculateSasa(c_allocator, input, config) catch |err| {
            return calcErrorCode(err);
        }
    else
        lee_richards.calculateSasaParallel(c_allocator, input, config, n_threads) catch |err| {
            return calcErrorCode(err);
        };
    defer {
        // Free the result's internal allocation
        c_allocator.free(result.atom_areas);
    }

    // Copy results to output buffers
    @memcpy(atom_areas[0..n_atoms], result.atom_areas);
    total_area.* = result.total_area;

    return ZSASA_OK;
}

// =============================================================================
// Classifier Functions
// =============================================================================

// Internal helper: get radius by classifier type
//
// This API has no element argument. For NACCESS and OONS an atom outside the
// tables (hydrogens, ligands) gets a radius guessed from its residue and atom
// name; CCD and ProtOr return null.
fn getRadiusByClassifier(classifier_type: c_int, residue: []const u8, atom: []const u8) ?f64 {
    return switch (classifier_type) {
        ZSASA_CLASSIFIER_NACCESS => classifier_naccess.getRadius(residue, atom) orelse
            classifier.guessRadiusFromResidueAtom(residue, atom),
        ZSASA_CLASSIFIER_PROTOR, ZSASA_CLASSIFIER_CCD => blk: {
            var ccd = classifier_ccd.CcdClassifier.init(std.heap.page_allocator);
            defer ccd.deinit();
            break :blk ccd.getRadius(residue, atom);
        },
        ZSASA_CLASSIFIER_OONS => classifier_oons.getRadius(residue, atom) orelse
            classifier.guessRadiusFromResidueAtom(residue, atom),
        else => null,
    };
}

// Internal helper: get class by classifier type
fn getClassByClassifier(classifier_type: c_int, residue: []const u8, atom: []const u8) classifier.AtomClass {
    return switch (classifier_type) {
        ZSASA_CLASSIFIER_NACCESS => classifier_naccess.getClass(residue, atom),
        ZSASA_CLASSIFIER_PROTOR, ZSASA_CLASSIFIER_CCD => blk: {
            var ccd = classifier_ccd.CcdClassifier.init(std.heap.page_allocator);
            defer ccd.deinit();
            break :blk ccd.getClass(residue, atom);
        },
        ZSASA_CLASSIFIER_OONS => classifier_oons.getClass(residue, atom),
        else => .unknown,
    };
}

// Internal helper: convert AtomClass to C int
fn atomClassToInt(atom_class: classifier.AtomClass) c_int {
    return switch (atom_class) {
        .polar => ZSASA_ATOM_CLASS_POLAR,
        .apolar => ZSASA_ATOM_CLASS_APOLAR,
        .unknown => ZSASA_ATOM_CLASS_UNKNOWN,
    };
}

/// Get van der Waals radius for an atom using the specified classifier.
///
/// Parameters:
///   classifier_type: Classifier to use (ZSASA_CLASSIFIER_NACCESS, etc.)
///   residue: Residue name (e.g., "ALA", "GLY") - null-terminated
///   atom: Atom name (e.g., "CA", "CB") - null-terminated
///
/// Returns:
///   Radius in Angstroms, or NaN if atom is not found in classifier.
///
/// NACCESS and OONS have no entries for hydrogens or ligands. For those atoms
/// the element is guessed from the names, because no element can be passed
/// here: a name starting with H, C, N, O, P or S is that element ("HG" is
/// hydrogen, "NA" in "HEM" is nitrogen), and an ion is recognized by a
/// residue name equal to its atom name ("ZN" in "ZN"). Callers that know the
/// element should use zsasa_guess_radius for atoms whose class is
/// ZSASA_ATOM_CLASS_UNKNOWN.
export fn zsasa_classifier_get_radius(
    classifier_type: c_int,
    residue: [*:0]const u8,
    atom: [*:0]const u8,
) callconv(.c) f64 {
    const radius = getRadiusByClassifier(
        classifier_type,
        std.mem.span(residue),
        std.mem.span(atom),
    );
    return radius orelse std.math.nan(f64);
}

/// Get atom polarity class using the specified classifier.
///
/// Parameters:
///   classifier_type: Classifier to use (ZSASA_CLASSIFIER_NACCESS, etc.)
///   residue: Residue name (e.g., "ALA", "GLY") - null-terminated
///   atom: Atom name (e.g., "CA", "CB") - null-terminated
///
/// Returns:
///   Atom class (ZSASA_ATOM_CLASS_POLAR, ZSASA_ATOM_CLASS_APOLAR, or ZSASA_ATOM_CLASS_UNKNOWN)
export fn zsasa_classifier_get_class(
    classifier_type: c_int,
    residue: [*:0]const u8,
    atom: [*:0]const u8,
) callconv(.c) c_int {
    const atom_class = getClassByClassifier(
        classifier_type,
        std.mem.span(residue),
        std.mem.span(atom),
    );
    return atomClassToInt(atom_class);
}

/// Guess van der Waals radius from element symbol.
///
/// Parameters:
///   element: Element symbol (e.g., "C", "N", "FE") - null-terminated
///            Case-insensitive, whitespace is trimmed.
///
/// Returns:
///   Radius in Angstroms, or NaN if element is not recognized.
export fn zsasa_guess_radius(
    element: [*:0]const u8,
) callconv(.c) f64 {
    const element_slice = std.mem.span(element);
    return classifier.guessRadius(element_slice) orelse std.math.nan(f64);
}

/// Guess van der Waals radius from PDB-style atom name.
/// Extracts element symbol from atom name and returns corresponding radius.
///
/// Parameters:
///   atom_name: PDB-style atom name (e.g., " CA ", "FE  ") - null-terminated
///              Following PDB conventions:
///              - Leading space indicates single-char element (e.g., " CA " = Carbon alpha)
///              - No leading space may indicate 2-char element (e.g., "FE  " = Iron)
///              Pass the name with its column padding: a trimmed "CA" or "HG"
///              is read as calcium or mercury.
///
/// Returns:
///   Radius in Angstroms, or NaN if element cannot be determined.
export fn zsasa_guess_radius_from_atom_name(
    atom_name: [*:0]const u8,
) callconv(.c) f64 {
    const atom_slice = std.mem.span(atom_name);
    return classifier.guessRadiusFromAtomName(atom_slice) orelse std.math.nan(f64);
}

/// Classify multiple atoms at once (batch operation).
///
/// This is more efficient than calling zsasa_classifier_get_radius
/// for each atom individually.
///
/// Parameters:
///   classifier_type: Classifier to use (ZSASA_CLASSIFIER_NACCESS, etc.)
///   residues: Array of residue name pointers (null-terminated strings)
///   atoms: Array of atom name pointers (null-terminated strings)
///   n_atoms: Number of atoms
///   radii_out: Output buffer for radii (pre-allocated, n_atoms elements)
///              NaN is written for atoms not found in classifier
///              (see zsasa_classifier_get_radius for NACCESS and OONS)
///   classes_out: Output buffer for classes (pre-allocated, n_atoms elements)
///                Can be NULL if classes are not needed
///
/// Returns:
///   ZSASA_OK on success, negative error code on failure.
export fn zsasa_classify_atoms(
    classifier_type: c_int,
    residues: [*]const [*:0]const u8,
    atoms: [*]const [*:0]const u8,
    n_atoms: usize,
    radii_out: [*]f64,
    classes_out: ?[*]c_int,
) callconv(.c) c_int {
    if (n_atoms == 0) {
        return ZSASA_OK; // Nothing to do
    }

    // Validate classifier type
    if (classifier_type < ZSASA_CLASSIFIER_NACCESS or classifier_type > ZSASA_CLASSIFIER_CCD) {
        return ZSASA_ERROR_INVALID_INPUT;
    }

    for (0..n_atoms) |i| {
        const residue_slice = std.mem.span(residues[i]);
        const atom_slice = std.mem.span(atoms[i]);

        // Get radius
        const radius = getRadiusByClassifier(classifier_type, residue_slice, atom_slice);
        radii_out[i] = radius orelse std.math.nan(f64);

        // Get class if requested
        if (classes_out) |classes| {
            const atom_class = getClassByClassifier(classifier_type, residue_slice, atom_slice);
            classes[i] = atomClassToInt(atom_class);
        }
    }

    return ZSASA_OK;
}

// =============================================================================
// RSA (Relative Solvent Accessibility) Functions
// =============================================================================

/// Get maximum SASA value for a standard amino acid.
/// Values from Tien et al. (2013) "Maximum allowed solvent accessibilities
/// of residues in proteins".
///
/// Parameters:
///   residue_name: 3-letter residue code (e.g., "ALA", "GLY") - null-terminated
///
/// Returns:
///   Maximum SASA in Angstroms², or NaN if residue is not a standard amino acid.
export fn zsasa_get_max_sasa(
    residue_name: [*:0]const u8,
) callconv(.c) f64 {
    const residue_slice = std.mem.span(residue_name);
    return analysis.MaxSASA.get(residue_slice) orelse std.math.nan(f64);
}

/// Calculate RSA (Relative Solvent Accessibility) for a single residue.
/// RSA = SASA / MaxSASA
///
/// Parameters:
///   sasa: Observed SASA value in Angstroms²
///   residue_name: 3-letter residue code (e.g., "ALA", "GLY") - null-terminated
///
/// Returns:
///   RSA value (0.0-1.0+), or NaN if residue is not a standard amino acid.
///   Note: RSA > 1.0 is possible for exposed terminal residues.
export fn zsasa_calculate_rsa(
    sasa: f64,
    residue_name: [*:0]const u8,
) callconv(.c) f64 {
    const residue_slice = std.mem.span(residue_name);
    const max_sasa = analysis.MaxSASA.get(residue_slice) orelse return std.math.nan(f64);
    if (max_sasa <= 0.0) {
        return std.math.nan(f64);
    }
    return sasa / max_sasa;
}

/// Calculate RSA for multiple residues at once (batch operation).
///
/// Parameters:
///   sasas: Array of SASA values in Angstroms²
///   residue_names: Array of residue name pointers (null-terminated strings)
///   n_residues: Number of residues
///   rsa_out: Output buffer for RSA values (pre-allocated, n_residues elements)
///            NaN is written for non-standard amino acids
///
/// Returns:
///   ZSASA_OK on success.
export fn zsasa_calculate_rsa_batch(
    sasas: [*]const f64,
    residue_names: [*]const [*:0]const u8,
    n_residues: usize,
    rsa_out: [*]f64,
) callconv(.c) c_int {
    for (0..n_residues) |i| {
        const residue_slice = std.mem.span(residue_names[i]);
        const max_sasa = analysis.MaxSASA.get(residue_slice);
        if (max_sasa) |max| {
            if (max > 0.0) {
                rsa_out[i] = sasas[i] / max;
            } else {
                rsa_out[i] = std.math.nan(f64);
            }
        } else {
            rsa_out[i] = std.math.nan(f64);
        }
    }
    return ZSASA_OK;
}

// =============================================================================
// XTC Trajectory Reader Functions
// =============================================================================

/// Error code: End of file reached
pub const ZSASA_XTC_END_OF_FILE: c_int = 1;

/// XTC reader handle (opaque pointer to internal struct)
const XtcHandle = struct {
    reader: xtc.XtcReader,
};

/// Open an XTC trajectory file.
///
/// Parameters:
///   path: Path to XTC file (null-terminated string)
///   natoms_out: Output pointer for number of atoms (set on success)
///   error_code: Output pointer for error code (set on failure):
///     ZSASA_ERROR_INVALID_INPUT: the file does not exist
///     ZSASA_ERROR_INVALID_FORMAT: the file is not an XTC file (another format,
///       empty or truncated)
///     ZSASA_ERROR_OUT_OF_MEMORY: allocation failed
///
/// Returns:
///   Opaque handle on success, null on failure.
///   Caller must call zsasa_xtc_close() to free resources.
export fn zsasa_xtc_open(
    path: [*:0]const u8,
    natoms_out: *i32,
    error_code: *c_int,
) callconv(.c) ?*anyopaque {
    const path_slice = std.mem.span(path);

    const handle = c_allocator.create(XtcHandle) catch {
        error_code.* = ZSASA_ERROR_OUT_OF_MEMORY;
        return null;
    };

    handle.reader = xtc.XtcReader.open(cIo(), c_allocator, path_slice) catch |err| {
        c_allocator.destroy(handle);
        error_code.* = trajectoryErrorCode(err);
        return null;
    };

    natoms_out.* = @intCast(handle.reader.nAtoms());
    error_code.* = ZSASA_OK;
    return handle;
}

/// Close an XTC trajectory file and free resources.
///
/// Parameters:
///   handle: Handle returned by zsasa_xtc_open()
export fn zsasa_xtc_close(
    handle: ?*anyopaque,
) callconv(.c) void {
    if (handle) |h| {
        const xtc_handle: *XtcHandle = @ptrCast(@alignCast(h));
        xtc_handle.reader.deinit();
        c_allocator.destroy(xtc_handle);
    }
}

/// Read the next frame from an XTC trajectory.
///
/// Parameters:
///   handle: Handle returned by zsasa_xtc_open()
///   coords_out: Output buffer for coordinates (natoms * 3 floats, in nm)
///   step_out: Output pointer for step number
///   time_out: Output pointer for frame time
///   box_out: Output buffer for box matrix (9 floats, row-major 3x3)
///   precision_out: Output pointer for precision value
///
/// Returns:
///   ZSASA_OK (0) on success
///   ZSASA_XTC_END_OF_FILE (1) when no more frames
///   ZSASA_ERROR_INVALID_FORMAT when the frame is corrupt or truncated
///   Other negative error code on failure
export fn zsasa_xtc_read_frame(
    handle: ?*anyopaque,
    coords_out: [*]f32,
    step_out: *i32,
    time_out: *f32,
    box_out: [*]f32,
    precision_out: *f32,
) callconv(.c) c_int {
    if (handle == null) {
        return ZSASA_ERROR_INVALID_INPUT;
    }

    const xtc_handle: *XtcHandle = @ptrCast(@alignCast(handle.?));

    const frame = (xtc_handle.reader.next() catch |err| {
        return trajectoryErrorCode(err);
    }) orelse return ZSASA_XTC_END_OF_FILE;

    // Copy step, time, precision. ztraj's high-level XTC reader normalizes to Å
    // and does not expose per-frame precision, so preserve the historical C ABI
    // contract with the default XTC precision.
    step_out.* = frame.step;
    time_out.* = frame.time;
    precision_out.* = 1000.0;

    // Copy box (3x3 matrix, row-major) in nm for the public XTC C/Python API.
    if (frame.box_vectors) |box| {
        for (0..3) |i| {
            for (0..3) |j| {
                box_out[i * 3 + j] = box[i][j] / 10.0;
            }
        }
    } else {
        for (0..9) |i| {
            box_out[i] = 0.0;
        }
    }

    // Copy coordinates in nm for the public XTC C/Python API.
    const natoms: usize = @intCast(xtc_handle.reader.nAtoms());
    for (0..natoms) |i| {
        coords_out[i * 3 + 0] = frame.x[i] / 10.0;
        coords_out[i * 3 + 1] = frame.y[i] / 10.0;
        coords_out[i * 3 + 2] = frame.z[i] / 10.0;
    }

    return ZSASA_OK;
}

/// Get the number of atoms in an opened XTC file.
///
/// Parameters:
///   handle: Handle returned by zsasa_xtc_open()
///
/// Returns:
///   Number of atoms, or -1 if handle is invalid
export fn zsasa_xtc_get_natoms(
    handle: ?*anyopaque,
) callconv(.c) i32 {
    if (handle == null) {
        return -1;
    }
    const xtc_handle: *XtcHandle = @ptrCast(@alignCast(handle.?));
    return @intCast(xtc_handle.reader.nAtoms());
}

// Tests
test "zsasa_version returns valid string" {
    const version = zsasa_version();
    try std.testing.expect(version[0] != 0);
}

test "zsasa_calc_sr with empty input returns error" {
    var total_area: f64 = 0.0;
    const result = zsasa_calc_sr(
        undefined,
        undefined,
        undefined,
        undefined,
        0, // n_atoms = 0
        100,
        1.4,
        1,
        undefined,
        &total_area,
    );
    try std.testing.expectEqual(ZSASA_ERROR_INVALID_INPUT, result);
}

test "zsasa_calc_sr rejects non-finite inputs" {
    var x = [_]f64{0.0};
    var y = [_]f64{0.0};
    var z = [_]f64{0.0};
    var radii = [_]f64{1.5};
    var atom_areas = [_]f64{0.0};
    var total_area: f64 = 0.0;

    x[0] = std.math.nan(f64);
    try std.testing.expectEqual(
        ZSASA_ERROR_INVALID_INPUT,
        zsasa_calc_sr(&x, &y, &z, &radii, 1, 100, 1.4, 1, &atom_areas, &total_area),
    );

    x[0] = 0.0;
    y[0] = std.math.inf(f64);
    try std.testing.expectEqual(
        ZSASA_ERROR_INVALID_INPUT,
        zsasa_calc_sr(&x, &y, &z, &radii, 1, 100, 1.4, 1, &atom_areas, &total_area),
    );

    y[0] = 0.0;
    radii[0] = std.math.nan(f64);
    try std.testing.expectEqual(
        ZSASA_ERROR_INVALID_INPUT,
        zsasa_calc_sr(&x, &y, &z, &radii, 1, 100, 1.4, 1, &atom_areas, &total_area),
    );

    radii[0] = 1.5;
    try std.testing.expectEqual(
        ZSASA_ERROR_INVALID_INPUT,
        zsasa_calc_sr(&x, &y, &z, &radii, 1, 100, std.math.inf(f64), 1, &atom_areas, &total_area),
    );
}

test "zsasa_calc_sr single atom" {
    const x = [_]f64{0.0};
    const y = [_]f64{0.0};
    const z = [_]f64{0.0};
    const radii = [_]f64{1.5};
    var atom_areas = [_]f64{0.0};
    var total_area: f64 = 0.0;

    const result = zsasa_calc_sr(
        &x,
        &y,
        &z,
        &radii,
        1,
        100,
        1.4,
        1,
        &atom_areas,
        &total_area,
    );

    try std.testing.expectEqual(ZSASA_OK, result);
    // Expected: 4π * (1.5 + 1.4)² ≈ 105.68 Ų
    try std.testing.expect(total_area > 100.0 and total_area < 110.0);
    try std.testing.expectEqual(total_area, atom_areas[0]);
}

test "zsasa_calc_lr single atom" {
    const x = [_]f64{0.0};
    const y = [_]f64{0.0};
    const z = [_]f64{0.0};
    const radii = [_]f64{1.5};
    var atom_areas = [_]f64{0.0};
    var total_area: f64 = 0.0;

    const result = zsasa_calc_lr(
        &x,
        &y,
        &z,
        &radii,
        1,
        20,
        1.4,
        1,
        &atom_areas,
        &total_area,
    );

    try std.testing.expectEqual(ZSASA_OK, result);
    // Expected: 4π * (1.5 + 1.4)² ≈ 105.68 Ų
    try std.testing.expect(total_area > 100.0 and total_area < 110.0);
    try std.testing.expectEqual(total_area, atom_areas[0]);
}

test "zsasa_calc_lr rejects non-finite inputs" {
    var x = [_]f64{0.0};
    var y = [_]f64{0.0};
    var z = [_]f64{0.0};
    var radii = [_]f64{1.5};
    var atom_areas = [_]f64{0.0};
    var total_area: f64 = 0.0;

    z[0] = -std.math.inf(f64);
    try std.testing.expectEqual(
        ZSASA_ERROR_INVALID_INPUT,
        zsasa_calc_lr(&x, &y, &z, &radii, 1, 20, 1.4, 1, &atom_areas, &total_area),
    );

    z[0] = 0.0;
    radii[0] = std.math.inf(f64);
    try std.testing.expectEqual(
        ZSASA_ERROR_INVALID_INPUT,
        zsasa_calc_lr(&x, &y, &z, &radii, 1, 20, 1.4, 1, &atom_areas, &total_area),
    );

    radii[0] = 1.5;
    try std.testing.expectEqual(
        ZSASA_ERROR_INVALID_INPUT,
        zsasa_calc_lr(&x, &y, &z, &radii, 1, 20, std.math.nan(f64), 1, &atom_areas, &total_area),
    );
}

test "zsasa_calc_sr_bitmask single atom" {
    const x = [_]f64{0.0};
    const y = [_]f64{0.0};
    const z = [_]f64{0.0};
    const radii = [_]f64{1.5};
    var atom_areas = [_]f64{0.0};
    var total_area: f64 = 0.0;

    const result = zsasa_calc_sr_bitmask(
        &x,
        &y,
        &z,
        &radii,
        1,
        128,
        1.4,
        1,
        &atom_areas,
        &total_area,
    );

    try std.testing.expectEqual(ZSASA_OK, result);
    // Expected: 4π * (1.5 + 1.4)² ≈ 105.68 Å²
    try std.testing.expect(total_area > 100.0 and total_area < 110.0);
    try std.testing.expectEqual(total_area, atom_areas[0]);
}

test "zsasa_calc_sr_bitmask with empty input returns error" {
    var total_area: f64 = 0.0;
    const result = zsasa_calc_sr_bitmask(
        undefined,
        undefined,
        undefined,
        undefined,
        0, // n_atoms = 0
        128,
        1.4,
        1,
        undefined,
        &total_area,
    );
    try std.testing.expectEqual(ZSASA_ERROR_INVALID_INPUT, result);
}

test "zsasa_calc_sr_bitmask with unsupported n_points returns error" {
    const x = [_]f64{0.0};
    const y = [_]f64{0.0};
    const z = [_]f64{0.0};
    const radii = [_]f64{1.5};
    var atom_areas = [_]f64{0.0};
    var total_area: f64 = 0.0;

    const result = zsasa_calc_sr_bitmask(
        &x,
        &y,
        &z,
        &radii,
        1,
        2000, // Not supported (must be 1..1024)
        1.4,
        1,
        &atom_areas,
        &total_area,
    );
    try std.testing.expectEqual(ZSASA_ERROR_UNSUPPORTED_N_POINTS, result);
}

test "zsasa_calc_sr_batch rejects non-finite inputs" {
    var coordinates = [_]f32{ 0.0, 0.0, 0.0 };
    var radii = [_]f32{1.5};
    var atom_areas = [_]f32{0.0};

    coordinates[0] = std.math.nan(f32);
    try std.testing.expectEqual(
        ZSASA_ERROR_INVALID_INPUT,
        zsasa_calc_sr_batch(&coordinates, 1, 1, &radii, 100, 1.4, 1, &atom_areas),
    );

    coordinates[0] = 0.0;
    radii[0] = std.math.inf(f32);
    try std.testing.expectEqual(
        ZSASA_ERROR_INVALID_INPUT,
        zsasa_calc_sr_batch(&coordinates, 1, 1, &radii, 100, 1.4, 1, &atom_areas),
    );

    radii[0] = 1.5;
    try std.testing.expectEqual(
        ZSASA_ERROR_INVALID_INPUT,
        zsasa_calc_sr_batch(&coordinates, 1, 1, &radii, 100, std.math.inf(f32), 1, &atom_areas),
    );
}

test "zsasa_calc_sr_batch_bitmask basic" {
    // 2 frames, 1 atom each
    const coordinates = [_]f32{
        // Frame 0: atom at origin
        0.0, 0.0, 0.0,
        // Frame 1: atom at origin
        0.0, 0.0, 0.0,
    };
    const radii = [_]f32{1.5};
    var atom_areas: [2]f32 = .{ 0.0, 0.0 };

    const result = zsasa_calc_sr_batch_bitmask(
        &coordinates,
        2, // n_frames
        1, // n_atoms
        &radii,
        128,
        1.4,
        1,
        &atom_areas,
    );

    try std.testing.expectEqual(ZSASA_OK, result);
    // Both frames should have the same area (single isolated atom)
    try std.testing.expect(atom_areas[0] > 100.0 and atom_areas[0] < 110.0);
    try std.testing.expectApproxEqAbs(atom_areas[0], atom_areas[1], 0.01);
}

test "zsasa_calc_sr_batch_bitmask with unsupported n_points returns error" {
    const coordinates = [_]f32{ 0.0, 0.0, 0.0 };
    const radii = [_]f32{1.5};
    var atom_areas = [_]f32{0.0};

    const result = zsasa_calc_sr_batch_bitmask(
        &coordinates,
        1,
        1,
        &radii,
        2000, // Not supported
        1.4,
        1,
        &atom_areas,
    );
    try std.testing.expectEqual(ZSASA_ERROR_UNSUPPORTED_N_POINTS, result);
}

test "zsasa_calc_sr_batch_bitmask_f32 basic" {
    const coordinates = [_]f32{ 0.0, 0.0, 0.0 };
    const radii = [_]f32{1.5};
    var atom_areas = [_]f32{0.0};

    const result = zsasa_calc_sr_batch_bitmask_f32(
        &coordinates,
        1,
        1,
        &radii,
        128,
        1.4,
        1,
        &atom_areas,
    );

    try std.testing.expectEqual(ZSASA_OK, result);
    try std.testing.expect(atom_areas[0] > 100.0 and atom_areas[0] < 110.0);
}

test "calcErrorCode separates out-of-memory from calculation errors" {
    try std.testing.expectEqual(ZSASA_ERROR_OUT_OF_MEMORY, calcErrorCode(error.OutOfMemory));
    try std.testing.expectEqual(ZSASA_ERROR_INVALID_INPUT, calcErrorCode(error.CoordinateRangeTooLarge));
    try std.testing.expectEqual(ZSASA_ERROR_CALCULATION, calcErrorCode(error.InvalidRadius));
}

test "zsasa_calc far-apart atoms are isolated spheres" {
    // Issue #428: the neighbor grid of the first input had 2^22 x 2^21 x 2^21 cells, a
    // product that wrapped around to zero; the second needed hundreds of MB.
    const cases = [_]struct { position: [3]f64, radius: f64, n_threads: usize }{
        .{ .position = .{ (4194304.0 - 2.0) * 4.0, (2097152.0 - 2.0) * 4.0, (2097152.0 - 2.0) * 4.0 }, .radius = 0.6, .n_threads = 1 },
        .{ .position = .{ 2000.0, 2000.0, 2000.0 }, .radius = 1.7, .n_threads = 2 },
    };

    for (cases) |case| {
        const x = [_]f64{ 0.0, case.position[0] };
        const y = [_]f64{ 0.0, case.position[1] };
        const z = [_]f64{ 0.0, case.position[2] };
        const radii = [_]f64{ case.radius, case.radius };
        const expected = 4.0 * std.math.pi * (case.radius + 1.4) * (case.radius + 1.4);

        var atom_areas = [_]f64{ 0.0, 0.0 };
        var total_area: f64 = 0.0;

        try std.testing.expectEqual(ZSASA_OK, zsasa_calc_sr(&x, &y, &z, &radii, 2, 100, 1.4, case.n_threads, &atom_areas, &total_area));
        try std.testing.expectApproxEqRel(expected, atom_areas[0], 1e-12);
        try std.testing.expectApproxEqRel(expected, atom_areas[1], 1e-12);

        atom_areas = .{ 0.0, 0.0 };
        try std.testing.expectEqual(ZSASA_OK, zsasa_calc_sr_bitmask(&x, &y, &z, &radii, 2, 64, 1.4, case.n_threads, &atom_areas, &total_area));
        try std.testing.expectApproxEqRel(expected, atom_areas[0], 1e-12);
        try std.testing.expectApproxEqRel(expected, atom_areas[1], 1e-12);

        atom_areas = .{ 0.0, 0.0 };
        try std.testing.expectEqual(ZSASA_OK, zsasa_calc_lr(&x, &y, &z, &radii, 2, 20, 1.4, case.n_threads, &atom_areas, &total_area));
        try std.testing.expectApproxEqRel(expected, atom_areas[0], 1e-12);
        try std.testing.expectApproxEqRel(expected, atom_areas[1], 1e-12);
        try std.testing.expectApproxEqRel(2.0 * expected, total_area, 1e-12);
    }
}

test "zsasa_calc batch far-apart atoms are isolated spheres" {
    // One frame of two atoms about 3.4e10 Å apart: a dense grid of 6.2 Å cells would have
    // more cells than a usize can count
    const coordinates = [_]f32{ 0.0, 0.0, 0.0, 34359738352.0, 34359738352.0, 0.0 };
    const radii = [_]f32{ 1.7, 1.7 };
    const expected: f32 = 4.0 * std.math.pi * 3.1 * 3.1;

    const BatchFn = *const fn ([*]const f32, usize, usize, [*]const f32, u32, f32, usize, [*]f32) callconv(.c) c_int;
    const batch_fns = [_]struct { func: BatchFn, param: u32 }{
        .{ .func = &zsasa_calc_sr_batch, .param = 100 },
        .{ .func = &zsasa_calc_lr_batch, .param = 20 },
        .{ .func = &zsasa_calc_sr_batch_f32, .param = 100 },
        .{ .func = &zsasa_calc_lr_batch_f32, .param = 20 },
        .{ .func = &zsasa_calc_sr_batch_bitmask, .param = 128 },
        .{ .func = &zsasa_calc_sr_batch_bitmask_f32, .param = 128 },
    };

    for (batch_fns) |batch_fn| {
        var atom_areas: [2]f32 = .{ 0.0, 0.0 };
        try std.testing.expectEqual(ZSASA_OK, batch_fn.func(&coordinates, 1, 2, &radii, batch_fn.param, 1.4, 1, &atom_areas));
        for (atom_areas) |area| {
            try std.testing.expectApproxEqRel(expected, area, 1e-5);
        }
    }
}

test "zsasa_calc batch reports a coordinate range too wide for f32" {
    // Finite f32 coordinates whose difference overflows f32
    const max = std.math.floatMax(f32);
    const coordinates = [_]f32{ -max, 0.0, 0.0, max, 0.0, 0.0 };
    const radii = [_]f32{ 1.7, 1.7 };
    var atom_areas: [2]f32 = .{ 0.0, 0.0 };

    try std.testing.expectEqual(
        ZSASA_ERROR_INVALID_INPUT,
        zsasa_calc_sr_batch_f32(&coordinates, 1, 2, &radii, 100, 1.4, 1, &atom_areas),
    );
    try std.testing.expectEqual(
        ZSASA_ERROR_INVALID_INPUT,
        zsasa_calc_lr_batch_f32(&coordinates, 1, 2, &radii, 20, 1.4, 1, &atom_areas),
    );

    // The same frame fits when the calculation runs in f64
    const expected: f32 = 4.0 * std.math.pi * 3.1 * 3.1;
    try std.testing.expectEqual(ZSASA_OK, zsasa_calc_sr_batch(&coordinates, 1, 2, &radii, 100, 1.4, 1, &atom_areas));
    try std.testing.expectApproxEqRel(expected, atom_areas[0], 1e-5);
    try std.testing.expectApproxEqRel(expected, atom_areas[1], 1e-5);
}

// =============================================================================
// Classifier Tests
// =============================================================================

test "zsasa_classifier_get_radius NACCESS" {
    // Standard backbone atoms
    const ca_radius = zsasa_classifier_get_radius(ZSASA_CLASSIFIER_NACCESS, "ALA", "CA");
    try std.testing.expectApproxEqAbs(1.87, ca_radius, 0.01);

    const n_radius = zsasa_classifier_get_radius(ZSASA_CLASSIFIER_NACCESS, "ALA", "N");
    try std.testing.expectApproxEqAbs(1.65, n_radius, 0.01);

    const o_radius = zsasa_classifier_get_radius(ZSASA_CLASSIFIER_NACCESS, "ALA", "O");
    try std.testing.expectApproxEqAbs(1.40, o_radius, 0.01);

    // Unknown atom should return NaN
    const unknown = zsasa_classifier_get_radius(ZSASA_CLASSIFIER_NACCESS, "ALA", "XX");
    try std.testing.expect(std.math.isNan(unknown));
}

test "zsasa_classifier_get_radius ProtoR" {
    const ca_radius = zsasa_classifier_get_radius(ZSASA_CLASSIFIER_PROTOR, "ALA", "CA");
    try std.testing.expect(ca_radius > 1.0 and ca_radius < 3.0);
}

test "zsasa_classifier_get_radius OONS" {
    const ca_radius = zsasa_classifier_get_radius(ZSASA_CLASSIFIER_OONS, "ALA", "CA");
    try std.testing.expect(ca_radius > 1.0 and ca_radius < 3.0);
}

test "zsasa_classifier_get_radius NACCESS/OONS atoms outside the tables" {
    for ([_]c_int{ ZSASA_CLASSIFIER_NACCESS, ZSASA_CLASSIFIER_OONS }) |ct| {
        // Hydrogen, not mercury
        try std.testing.expectApproxEqAbs(1.10, zsasa_classifier_get_radius(ct, "ALA", "H"), 1e-9);
        try std.testing.expectApproxEqAbs(1.10, zsasa_classifier_get_radius(ct, "SER", "HG"), 1e-9);
        try std.testing.expectApproxEqAbs(1.10, zsasa_classifier_get_radius(ct, "PRO", "HG2"), 1e-9);
        try std.testing.expectApproxEqAbs(1.10, zsasa_classifier_get_radius(ct, "ILE", "HG12"), 1e-9);
        try std.testing.expectApproxEqAbs(1.10, zsasa_classifier_get_radius(ct, "VAL", "HG21"), 1e-9);
        // Nitrogen, not sodium
        try std.testing.expectApproxEqAbs(1.55, zsasa_classifier_get_radius(ct, "HEM", "NA"), 1e-9);
        // Phosphorus, not lead
        try std.testing.expectApproxEqAbs(1.80, zsasa_classifier_get_radius(ct, "ATP", "PB"), 1e-9);
        // Carbon, not cadmium
        try std.testing.expectApproxEqAbs(1.70, zsasa_classifier_get_radius(ct, "PCA", "CD"), 1e-9);
        try std.testing.expectApproxEqAbs(1.70, zsasa_classifier_get_radius(ct, "LIG", "CD1"), 1e-9);
        try std.testing.expectApproxEqAbs(1.70, zsasa_classifier_get_radius(ct, "LIG", "CD2"), 1e-9);

        // An ion is a residue named after its atom
        try std.testing.expectApproxEqAbs(1.39, zsasa_classifier_get_radius(ct, "ZN", "ZN"), 1e-9);
        try std.testing.expectApproxEqAbs(2.27, zsasa_classifier_get_radius(ct, "NA", "NA"), 1e-9);
        try std.testing.expectApproxEqAbs(1.55, zsasa_classifier_get_radius(ct, "HG", "HG"), 1e-9);
        try std.testing.expectApproxEqAbs(1.58, zsasa_classifier_get_radius(ct, "CD", "CD"), 1e-9);
        try std.testing.expectApproxEqAbs(1.26, zsasa_classifier_get_radius(ct, "HEM", "FE"), 1e-9);

        // Guessed radii have no polarity class
        try std.testing.expectEqual(ZSASA_ATOM_CLASS_UNKNOWN, zsasa_classifier_get_class(ct, "SER", "HG"));
    }

    // Table entries are not affected by the guess
    try std.testing.expectApproxEqAbs(1.87, zsasa_classifier_get_radius(ZSASA_CLASSIFIER_NACCESS, "ARG", "CD"), 1e-9);
    try std.testing.expectApproxEqAbs(1.76, zsasa_classifier_get_radius(ZSASA_CLASSIFIER_NACCESS, "PHE", "CD1"), 1e-9);
    try std.testing.expectApproxEqAbs(2.00, zsasa_classifier_get_radius(ZSASA_CLASSIFIER_OONS, "ARG", "CD"), 1e-9);
    try std.testing.expectApproxEqAbs(1.75, zsasa_classifier_get_radius(ZSASA_CLASSIFIER_OONS, "PHE", "CD1"), 1e-9);

    // CCD and ProtOr do not guess
    try std.testing.expect(std.math.isNan(zsasa_classifier_get_radius(ZSASA_CLASSIFIER_CCD, "LIG", "HG")));
    try std.testing.expect(std.math.isNan(zsasa_classifier_get_radius(ZSASA_CLASSIFIER_PROTOR, "LIG", "HG")));
}

test "zsasa_classifier_get_radius invalid classifier" {
    const radius = zsasa_classifier_get_radius(99, "ALA", "CA");
    try std.testing.expect(std.math.isNan(radius));
}

test "zsasa_classifier_get_class" {
    // Carbon atoms are apolar
    const ca_class = zsasa_classifier_get_class(ZSASA_CLASSIFIER_NACCESS, "ALA", "CA");
    try std.testing.expectEqual(ZSASA_ATOM_CLASS_APOLAR, ca_class);

    // Oxygen atoms are polar
    const o_class = zsasa_classifier_get_class(ZSASA_CLASSIFIER_NACCESS, "ALA", "O");
    try std.testing.expectEqual(ZSASA_ATOM_CLASS_POLAR, o_class);

    // Nitrogen atoms are polar
    const n_class = zsasa_classifier_get_class(ZSASA_CLASSIFIER_NACCESS, "ALA", "N");
    try std.testing.expectEqual(ZSASA_ATOM_CLASS_POLAR, n_class);

    // Unknown atoms
    const unknown_class = zsasa_classifier_get_class(ZSASA_CLASSIFIER_NACCESS, "ALA", "XX");
    try std.testing.expectEqual(ZSASA_ATOM_CLASS_UNKNOWN, unknown_class);
}

test "zsasa_guess_radius" {
    // Common elements
    try std.testing.expectApproxEqAbs(1.70, zsasa_guess_radius("C"), 0.01);
    try std.testing.expectApproxEqAbs(1.55, zsasa_guess_radius("N"), 0.01);
    try std.testing.expectApproxEqAbs(1.52, zsasa_guess_radius("O"), 0.01);
    try std.testing.expectApproxEqAbs(1.80, zsasa_guess_radius("S"), 0.01);

    // Case insensitive
    try std.testing.expectApproxEqAbs(1.70, zsasa_guess_radius("c"), 0.01);

    // Two-character elements
    try std.testing.expectApproxEqAbs(1.26, zsasa_guess_radius("FE"), 0.01);
    try std.testing.expectApproxEqAbs(1.39, zsasa_guess_radius("ZN"), 0.01);

    // Unknown element returns NaN
    const unknown = zsasa_guess_radius("XX");
    try std.testing.expect(std.math.isNan(unknown));
}

test "zsasa_guess_radius_from_atom_name" {
    // Standard PDB atom names (leading space = single-char element)
    try std.testing.expectApproxEqAbs(1.70, zsasa_guess_radius_from_atom_name(" CA "), 0.01);
    try std.testing.expectApproxEqAbs(1.55, zsasa_guess_radius_from_atom_name(" N  "), 0.01);
    try std.testing.expectApproxEqAbs(1.52, zsasa_guess_radius_from_atom_name(" O  "), 0.01);

    // Metal atoms (no leading space = 2-char element)
    try std.testing.expectApproxEqAbs(1.26, zsasa_guess_radius_from_atom_name("FE  "), 0.01);
    try std.testing.expectApproxEqAbs(1.39, zsasa_guess_radius_from_atom_name("ZN  "), 0.01);
}

test "zsasa_classify_atoms batch" {
    const residues = [_][*:0]const u8{ "ALA", "ALA", "GLY" };
    const atoms = [_][*:0]const u8{ "CA", "O", "N" };
    var radii: [3]f64 = undefined;
    var classes: [3]c_int = undefined;

    const result = zsasa_classify_atoms(
        ZSASA_CLASSIFIER_NACCESS,
        &residues,
        &atoms,
        3,
        &radii,
        &classes,
    );

    try std.testing.expectEqual(ZSASA_OK, result);

    // Check radii
    try std.testing.expectApproxEqAbs(1.87, radii[0], 0.01); // CA
    try std.testing.expectApproxEqAbs(1.40, radii[1], 0.01); // O
    try std.testing.expectApproxEqAbs(1.65, radii[2], 0.01); // N

    // Check classes
    try std.testing.expectEqual(ZSASA_ATOM_CLASS_APOLAR, classes[0]); // CA
    try std.testing.expectEqual(ZSASA_ATOM_CLASS_POLAR, classes[1]); // O
    try std.testing.expectEqual(ZSASA_ATOM_CLASS_POLAR, classes[2]); // N
}

test "zsasa_classify_atoms without classes" {
    const residues = [_][*:0]const u8{ "ALA", "ALA" };
    const atoms = [_][*:0]const u8{ "CA", "CB" };
    var radii: [2]f64 = undefined;

    const result = zsasa_classify_atoms(
        ZSASA_CLASSIFIER_NACCESS,
        &residues,
        &atoms,
        2,
        &radii,
        null, // classes_out is null
    );

    try std.testing.expectEqual(ZSASA_OK, result);
    try std.testing.expectApproxEqAbs(1.87, radii[0], 0.01);
    try std.testing.expectApproxEqAbs(1.87, radii[1], 0.01);
}

test "zsasa_classify_atoms empty" {
    var radii: [0]f64 = undefined;
    const residues: [*]const [*:0]const u8 = undefined;
    const atoms: [*]const [*:0]const u8 = undefined;

    const result = zsasa_classify_atoms(
        ZSASA_CLASSIFIER_NACCESS,
        residues,
        atoms,
        0,
        &radii,
        null,
    );

    try std.testing.expectEqual(ZSASA_OK, result);
}

test "zsasa_classify_atoms invalid classifier" {
    const residues = [_][*:0]const u8{"ALA"};
    const atoms = [_][*:0]const u8{"CA"};
    var radii: [1]f64 = undefined;

    const result = zsasa_classify_atoms(
        99, // Invalid classifier type
        &residues,
        &atoms,
        1,
        &radii,
        null,
    );

    try std.testing.expectEqual(ZSASA_ERROR_INVALID_INPUT, result);
}

// =============================================================================
// RSA Tests
// =============================================================================

test "zsasa_get_max_sasa standard amino acids" {
    // Test known amino acids (values from Tien et al. 2013)
    try std.testing.expectApproxEqAbs(129.0, zsasa_get_max_sasa("ALA"), 0.01);
    try std.testing.expectApproxEqAbs(104.0, zsasa_get_max_sasa("GLY"), 0.01);
    try std.testing.expectApproxEqAbs(285.0, zsasa_get_max_sasa("TRP"), 0.01);
    try std.testing.expectApproxEqAbs(274.0, zsasa_get_max_sasa("ARG"), 0.01);
    try std.testing.expectApproxEqAbs(193.0, zsasa_get_max_sasa("ASP"), 0.01);
}

test "zsasa_get_max_sasa unknown residue" {
    const unknown = zsasa_get_max_sasa("XXX");
    try std.testing.expect(std.math.isNan(unknown));

    const water = zsasa_get_max_sasa("HOH");
    try std.testing.expect(std.math.isNan(water));
}

test "zsasa_calculate_rsa" {
    // ALA: RSA = 64.5 / 129.0 = 0.5
    const rsa_ala = zsasa_calculate_rsa(64.5, "ALA");
    try std.testing.expectApproxEqAbs(0.5, rsa_ala, 0.001);

    // GLY: RSA = 52.0 / 104.0 = 0.5
    const rsa_gly = zsasa_calculate_rsa(52.0, "GLY");
    try std.testing.expectApproxEqAbs(0.5, rsa_gly, 0.001);

    // RSA > 1.0 is possible for exposed terminal residues
    const rsa_exposed = zsasa_calculate_rsa(150.0, "GLY");
    try std.testing.expect(rsa_exposed > 1.0);
    try std.testing.expectApproxEqAbs(150.0 / 104.0, rsa_exposed, 0.001);
}

test "zsasa_calculate_rsa unknown residue" {
    const rsa = zsasa_calculate_rsa(100.0, "XXX");
    try std.testing.expect(std.math.isNan(rsa));
}

test "zsasa_calculate_rsa_batch" {
    const sasas = [_]f64{ 64.5, 52.0, 142.5 };
    const residues = [_][*:0]const u8{ "ALA", "GLY", "TRP" };
    var rsa_out: [3]f64 = undefined;

    const result = zsasa_calculate_rsa_batch(&sasas, &residues, 3, &rsa_out);

    try std.testing.expectEqual(ZSASA_OK, result);
    try std.testing.expectApproxEqAbs(0.5, rsa_out[0], 0.001); // ALA: 64.5/129
    try std.testing.expectApproxEqAbs(0.5, rsa_out[1], 0.001); // GLY: 52/104
    try std.testing.expectApproxEqAbs(0.5, rsa_out[2], 0.001); // TRP: 142.5/285
}

test "zsasa_calculate_rsa_batch with unknown" {
    const sasas = [_]f64{ 64.5, 100.0 };
    const residues = [_][*:0]const u8{ "ALA", "HOH" };
    var rsa_out: [2]f64 = undefined;

    const result = zsasa_calculate_rsa_batch(&sasas, &residues, 2, &rsa_out);

    try std.testing.expectEqual(ZSASA_OK, result);
    try std.testing.expectApproxEqAbs(0.5, rsa_out[0], 0.001); // ALA: known
    try std.testing.expect(std.math.isNan(rsa_out[1])); // HOH: unknown
}

// =============================================================================
// XTC Reader Tests
// =============================================================================

test "zsasa_xtc_open and close" {
    var natoms: i32 = 0;
    var error_code: c_int = 0;

    const handle = zsasa_xtc_open("test_data/1l2y.xtc", &natoms, &error_code);
    try std.testing.expect(handle != null);
    try std.testing.expectEqual(ZSASA_OK, error_code);
    try std.testing.expectEqual(@as(i32, 304), natoms);

    // Get natoms via function
    try std.testing.expectEqual(@as(i32, 304), zsasa_xtc_get_natoms(handle));

    zsasa_xtc_close(handle);
}

test "zsasa_xtc_open invalid file" {
    var natoms: i32 = 0;
    var error_code: c_int = 0;

    const handle = zsasa_xtc_open("nonexistent.xtc", &natoms, &error_code);
    try std.testing.expect(handle == null);
    try std.testing.expectEqual(ZSASA_ERROR_INVALID_INPUT, error_code);
}

test "zsasa_xtc_read_frame" {
    var natoms: i32 = 0;
    var error_code: c_int = 0;

    const handle = zsasa_xtc_open("test_data/1l2y.xtc", &natoms, &error_code);
    try std.testing.expect(handle != null);
    defer zsasa_xtc_close(handle);

    // Allocate buffers
    const natoms_u: usize = @intCast(natoms);
    const coords = c_allocator.alloc(f32, natoms_u * 3) catch unreachable;
    defer c_allocator.free(coords);
    var box: [9]f32 = undefined;
    var step: i32 = 0;
    var time: f32 = 0;
    var precision: f32 = 0;

    // Read first frame
    const result = zsasa_xtc_read_frame(handle, coords.ptr, &step, &time, &box, &precision);
    try std.testing.expectEqual(ZSASA_OK, result);
    try std.testing.expectEqual(@as(i32, 1), step);

    // Check first atom coordinates (in nm)
    const tolerance: f32 = 0.0001;
    try std.testing.expectApproxEqAbs(@as(f32, -0.8901), coords[0], tolerance);
    try std.testing.expectApproxEqAbs(@as(f32, 0.4127), coords[1], tolerance);
    try std.testing.expectApproxEqAbs(@as(f32, -0.0555), coords[2], tolerance);
}

test "zsasa_xtc_read_all_frames" {
    var natoms: i32 = 0;
    var error_code: c_int = 0;

    const handle = zsasa_xtc_open("test_data/1l2y.xtc", &natoms, &error_code);
    try std.testing.expect(handle != null);
    defer zsasa_xtc_close(handle);

    const natoms_u: usize = @intCast(natoms);
    const coords = c_allocator.alloc(f32, natoms_u * 3) catch unreachable;
    defer c_allocator.free(coords);
    var box: [9]f32 = undefined;
    var step: i32 = 0;
    var time: f32 = 0;
    var precision: f32 = 0;

    var frame_count: usize = 0;
    while (true) {
        const result = zsasa_xtc_read_frame(handle, coords.ptr, &step, &time, &box, &precision);
        if (result == ZSASA_XTC_END_OF_FILE) break;
        try std.testing.expectEqual(ZSASA_OK, result);
        frame_count += 1;
    }

    // 1l2y.xtc has 38 frames
    try std.testing.expectEqual(@as(usize, 38), frame_count);
}

test "zsasa_xtc_get_natoms null handle" {
    try std.testing.expectEqual(@as(i32, -1), zsasa_xtc_get_natoms(null));
}

test "zsasa_xtc_read_frame null handle" {
    var coords: [3]f32 = undefined;
    var box: [9]f32 = undefined;
    var step: i32 = 0;
    var time: f32 = 0;
    var precision: f32 = 0;

    const result = zsasa_xtc_read_frame(null, &coords, &step, &time, &box, &precision);
    try std.testing.expectEqual(ZSASA_ERROR_INVALID_INPUT, result);
}

// =============================================================================
// DCD Trajectory Reader Functions
// =============================================================================

/// Error code: End of DCD file reached
pub const ZSASA_DCD_END_OF_FILE: c_int = 2;

/// DCD reader handle (opaque pointer to internal struct)
const DcdHandle = struct {
    reader: dcd.DcdReader,
};

/// Open a DCD trajectory file.
///
/// Parameters:
///   path: Path to DCD file (null-terminated string)
///   natoms_out: Output pointer for number of atoms (set on success)
///   error_code: Output pointer for error code (set on failure):
///     ZSASA_ERROR_INVALID_INPUT: the file does not exist
///     ZSASA_ERROR_INVALID_FORMAT: the file is not a (supported) DCD file
///       (another format, empty or truncated, or unreadable)
///     ZSASA_ERROR_OUT_OF_MEMORY: allocation failed
///
/// Returns:
///   Opaque handle on success, null on failure.
///   Caller must call zsasa_dcd_close() to free resources.
export fn zsasa_dcd_open(
    path: [*:0]const u8,
    natoms_out: *i32,
    error_code: *c_int,
) callconv(.c) ?*anyopaque {
    const path_slice = std.mem.span(path);

    const handle = c_allocator.create(DcdHandle) catch {
        error_code.* = ZSASA_ERROR_OUT_OF_MEMORY;
        return null;
    };

    handle.reader = dcd.DcdReader.open(cIo(), c_allocator, path_slice) catch |err| {
        c_allocator.destroy(handle);
        error_code.* = trajectoryErrorCode(err);
        return null;
    };

    natoms_out.* = @intCast(handle.reader.nAtoms());
    error_code.* = ZSASA_OK;
    return handle;
}

/// Close a DCD trajectory file and free resources.
///
/// Parameters:
///   handle: Handle returned by zsasa_dcd_open()
export fn zsasa_dcd_close(
    handle: ?*anyopaque,
) callconv(.c) void {
    if (handle) |h| {
        const dcd_handle: *DcdHandle = @ptrCast(@alignCast(h));
        dcd_handle.reader.deinit();
        c_allocator.destroy(dcd_handle);
    }
}

/// Read the next frame from a DCD trajectory.
///
/// Parameters:
///   handle: Handle returned by zsasa_dcd_open()
///   coords_out: Output buffer for coordinates (natoms * 3 floats, in Angstroms)
///   step_out: Output pointer for step number
///   time_out: Output pointer for frame time
///   unitcell_out: Output buffer for unitcell (6 doubles), or null if not needed
///
/// Returns:
///   ZSASA_OK (0) on success
///   ZSASA_DCD_END_OF_FILE (2) when no more frames
///   ZSASA_ERROR_INVALID_FORMAT when the frame is corrupt or truncated
///   Other negative error code on failure
export fn zsasa_dcd_read_frame(
    handle: ?*anyopaque,
    coords_out: [*]f32,
    step_out: *i32,
    time_out: *f32,
    unitcell_out: ?[*]f64,
) callconv(.c) c_int {
    if (handle == null) {
        return ZSASA_ERROR_INVALID_INPUT;
    }

    const dcd_handle: *DcdHandle = @ptrCast(@alignCast(handle.?));

    const frame = (dcd_handle.reader.next() catch |err| {
        return trajectoryErrorCode(err);
    }) orelse return ZSASA_DCD_END_OF_FILE;

    // Copy step, time
    step_out.* = frame.step;
    time_out.* = frame.time;

    // Copy unitcell if present and output buffer provided
    if (unitcell_out) |uc_out| {
        if (frame.box_vectors) |box| {
            uc_out[0] = box[0][0];
            uc_out[1] = 90.0;
            uc_out[2] = box[1][1];
            uc_out[3] = 90.0;
            uc_out[4] = 90.0;
            uc_out[5] = box[2][2];
        } else {
            for (0..6) |i| {
                uc_out[i] = 0.0;
            }
        }
    }

    // Copy coordinates
    const natoms: usize = @intCast(dcd_handle.reader.nAtoms());
    for (0..natoms) |i| {
        coords_out[i * 3 + 0] = frame.x[i];
        coords_out[i * 3 + 1] = frame.y[i];
        coords_out[i * 3 + 2] = frame.z[i];
    }

    return ZSASA_OK;
}

/// Get the number of atoms in an opened DCD file.
///
/// Parameters:
///   handle: Handle returned by zsasa_dcd_open()
///
/// Returns:
///   Number of atoms, or -1 if handle is invalid
export fn zsasa_dcd_get_natoms(
    handle: ?*anyopaque,
) callconv(.c) i32 {
    if (handle == null) {
        return -1;
    }
    const dcd_handle: *DcdHandle = @ptrCast(@alignCast(handle.?));
    return @intCast(dcd_handle.reader.nAtoms());
}

// =============================================================================
// DCD Reader Tests
// =============================================================================

test "zsasa_dcd_open and close" {
    var natoms: i32 = 0;
    var error_code: c_int = 0;

    const handle = zsasa_dcd_open("test_data/1l2y.dcd", &natoms, &error_code);
    if (handle == null and error_code == ZSASA_ERROR_INVALID_INPUT) return; // Skip if not available
    try std.testing.expect(handle != null);
    try std.testing.expectEqual(ZSASA_OK, error_code);
    try std.testing.expectEqual(@as(i32, 304), natoms);

    try std.testing.expectEqual(@as(i32, 304), zsasa_dcd_get_natoms(handle));

    zsasa_dcd_close(handle);
}

test "zsasa_dcd_open invalid file" {
    var natoms: i32 = 0;
    var error_code: c_int = 0;

    const handle = zsasa_dcd_open("nonexistent.dcd", &natoms, &error_code);
    try std.testing.expect(handle == null);
    try std.testing.expectEqual(ZSASA_ERROR_INVALID_INPUT, error_code);
}

test "zsasa_dcd_read_frame" {
    var natoms: i32 = 0;
    var error_code: c_int = 0;

    const handle = zsasa_dcd_open("test_data/1l2y.dcd", &natoms, &error_code);
    if (handle == null and error_code == ZSASA_ERROR_INVALID_INPUT) return;
    try std.testing.expect(handle != null);
    defer zsasa_dcd_close(handle);

    const natoms_u: usize = @intCast(natoms);
    const coords = c_allocator.alloc(f32, natoms_u * 3) catch unreachable;
    defer c_allocator.free(coords);
    var unitcell: [6]f64 = undefined;
    var step: i32 = 0;
    var time: f32 = 0;

    const result = zsasa_dcd_read_frame(handle, coords.ptr, &step, &time, &unitcell);
    try std.testing.expectEqual(ZSASA_OK, result);

    // DCD coordinates are in Angstroms
    const tolerance: f32 = 0.05;
    try std.testing.expectApproxEqAbs(@as(f32, -8.901), coords[0], tolerance);
    try std.testing.expectApproxEqAbs(@as(f32, 4.127), coords[1], tolerance);
    try std.testing.expectApproxEqAbs(@as(f32, -0.555), coords[2], tolerance);
}

test "zsasa_dcd_read_all_frames" {
    var natoms: i32 = 0;
    var error_code: c_int = 0;

    const handle = zsasa_dcd_open("test_data/1l2y.dcd", &natoms, &error_code);
    if (handle == null and error_code == ZSASA_ERROR_INVALID_INPUT) return;
    try std.testing.expect(handle != null);
    defer zsasa_dcd_close(handle);

    const natoms_u: usize = @intCast(natoms);
    const coords = c_allocator.alloc(f32, natoms_u * 3) catch unreachable;
    defer c_allocator.free(coords);
    var unitcell: [6]f64 = undefined;
    var step: i32 = 0;
    var time: f32 = 0;

    var frame_count: usize = 0;
    while (true) {
        const result = zsasa_dcd_read_frame(handle, coords.ptr, &step, &time, &unitcell);
        if (result == ZSASA_DCD_END_OF_FILE) break;
        try std.testing.expectEqual(ZSASA_OK, result);
        frame_count += 1;
    }

    try std.testing.expectEqual(@as(usize, 38), frame_count);
}

test "zsasa_dcd_get_natoms null handle" {
    try std.testing.expectEqual(@as(i32, -1), zsasa_dcd_get_natoms(null));
}

test "zsasa_dcd_read_frame null handle" {
    var coords: [3]f32 = undefined;
    var step: i32 = 0;
    var time: f32 = 0;

    const result = zsasa_dcd_read_frame(null, &coords, &step, &time, null);
    try std.testing.expectEqual(ZSASA_ERROR_INVALID_INPUT, result);
}

// =============================================================================
// Batch Directory Processing Functions
// =============================================================================

/// Internal handle for batch directory processing results.
/// Stores the BatchResult and C-compatible null-terminated filename copies.
const BatchDirHandle = struct {
    result: batch.BatchResult,
    c_filenames: [][*:0]u8,

    fn deinit(self: *BatchDirHandle) void {
        for (self.c_filenames) |c_str| {
            // Free the null-terminated copy (allocated as slice of len+1)
            const len = std.mem.len(c_str);
            c_allocator.free(c_str[0 .. len + 1]);
        }
        c_allocator.free(self.c_filenames);
        self.result.deinit();
        c_allocator.destroy(self);
    }
};

/// Map an error returned by a directory batch run to a C API error code.
/// `stage` is the step of the run that failed: the same filesystem error means
/// a bad input directory while scanning and a bad output directory while
/// creating it.
fn batchErrorCode(err: anyerror, stage: batch.BatchStage) c_int {
    return switch (err) {
        error.OutOfMemory => ZSASA_ERROR_OUT_OF_MEMORY,
        error.OutputNameCollision => ZSASA_ERROR_OUTPUT_NAME_COLLISION,
        else => switch (stage) {
            .create_output_dir => ZSASA_ERROR_OUTPUT_DIR,
            .scan_inputs, .process => ZSASA_ERROR_FILE_IO,
        },
    };
}

/// Process all structure files in a directory.
///
/// Scans the directory for supported files (.pdb, .cif, .mmcif, .bcif, .json,
/// .ent and compressed variants), calculates SASA for each, and returns results
/// via an opaque handle.
///
/// Parameters:
///   input_dir: Path to directory containing structure files (null-terminated)
///   output_dir: Path to output directory for results (null-terminated), or null for no file output
///   algorithm: ZSASA_ALGORITHM_SR (0) or ZSASA_ALGORITHM_LR (1). Lee-Richards
///     uses exact trigonometry for its arc angles.
///   n_points: Number of test points (SR) or slices (LR)
///   probe_radius: Water probe radius in Angstroms (e.g., 1.4)
///   n_threads: Number of threads (0 = auto-detect)
///   classifier_type: ZSASA_CLASSIFIER_* constant, or -1 to use radii from input files
///   include_hydrogens: 0 to exclude, 1 to include hydrogen atoms
///   include_hetatm: 0 to exclude, 1 to include HETATM records
///   error_code: Output pointer for error code
///
/// Returns:
///   Opaque handle on success, null on failure.
///   Partial success (some files fail) returns a valid handle -- check per-file status.
///   Empty directory returns a valid handle with total_files=0.
///   When output_dir is set and several inputs would be written to the same
///   output file (e.g. 1crn.pdb and 1crn.cif), nothing is processed and
///   error_code is set to ZSASA_ERROR_OUTPUT_NAME_COLLISION.
///   When output_dir cannot be created (for example a path below a regular
///   file), error_code is set to ZSASA_ERROR_OUTPUT_DIR; an input_dir that
///   cannot be read gives ZSASA_ERROR_FILE_IO.
///   Safe to call from several threads at once.
///   Caller must call zsasa_batch_dir_free() to release resources.
export fn zsasa_batch_dir_process(
    input_dir: ?[*:0]const u8,
    output_dir: ?[*:0]const u8,
    algorithm: c_int,
    n_points: u32,
    probe_radius: f64,
    n_threads: usize,
    classifier_type: c_int,
    include_hydrogens: c_int,
    include_hetatm: c_int,
    error_code: ?*c_int,
) callconv(.c) ?*anyopaque {
    // Helper to set error code safely
    const setError = struct {
        fn call(ec: ?*c_int, code: c_int) void {
            if (ec) |p| p.* = code;
        }
    }.call;

    // Validate required parameters
    if (input_dir == null) {
        setError(error_code, ZSASA_ERROR_INVALID_INPUT);
        return null;
    }
    if (n_points == 0 or !isPositiveFinite(f64, probe_radius)) {
        setError(error_code, ZSASA_ERROR_INVALID_INPUT);
        return null;
    }
    if (algorithm != ZSASA_ALGORITHM_SR and algorithm != ZSASA_ALGORITHM_LR) {
        setError(error_code, ZSASA_ERROR_INVALID_INPUT);
        return null;
    }

    // Build BatchConfig
    const algo: batch.Algorithm = if (algorithm == ZSASA_ALGORITHM_LR) .lr else .sr;

    const ct: ?classifier.ClassifierType = switch (classifier_type) {
        ZSASA_CLASSIFIER_NACCESS => .naccess,
        ZSASA_CLASSIFIER_PROTOR => .protor,
        ZSASA_CLASSIFIER_OONS => .oons,
        ZSASA_CLASSIFIER_CCD => .ccd,
        -1 => null,
        else => {
            setError(error_code, ZSASA_ERROR_INVALID_INPUT);
            return null;
        },
    };

    const config = batch.BatchConfig{
        .n_threads = n_threads,
        .algorithm = algo,
        .n_points = n_points,
        .n_slices = n_points, // same parameter used for LR slices
        .probe_radius = probe_radius,
        .quiet = true,
        .precision = .f64,
        .classifier_type = ct,
        .include_hydrogens = include_hydrogens != 0,
        .include_hetatm = include_hetatm != 0,
    };

    const input_dir_slice = std.mem.span(input_dir.?);
    const output_dir_slice: ?[]const u8 = if (output_dir) |od| std.mem.span(od) else null;

    // Batch processing spawns N worker threads, so we need a multi-threaded Io
    // rather than the global single-threaded Io returned by cIo().
    const batch_io = sharedThreadedIo();

    // Run batch processing
    var stage: batch.BatchStage = .scan_inputs;
    var batch_result = batch.runBatchReportingStage(c_allocator, batch_io, input_dir_slice, output_dir_slice, config, &stage) catch |err| {
        setError(error_code, batchErrorCode(err, stage));
        return null;
    };

    // Build C-compatible null-terminated filename copies
    const n_files = batch_result.file_results.len;
    const c_filenames = c_allocator.alloc([*:0]u8, n_files) catch {
        batch_result.deinit();
        setError(error_code, ZSASA_ERROR_OUT_OF_MEMORY);
        return null;
    };

    for (batch_result.file_results, 0..) |fr, i| {
        const name = fr.filename;
        const buf = c_allocator.alloc(u8, name.len + 1) catch {
            // Free already-allocated filenames
            for (c_filenames[0..i]) |prev| {
                const prev_len = std.mem.len(prev);
                c_allocator.free(prev[0 .. prev_len + 1]);
            }
            c_allocator.free(c_filenames);
            batch_result.deinit();
            setError(error_code, ZSASA_ERROR_OUT_OF_MEMORY);
            return null;
        };
        @memcpy(buf[0..name.len], name);
        buf[name.len] = 0;
        c_filenames[i] = buf[0..name.len :0];
    }

    // Allocate handle
    const handle = c_allocator.create(BatchDirHandle) catch {
        for (c_filenames) |c_str| {
            const len = std.mem.len(c_str);
            c_allocator.free(c_str[0 .. len + 1]);
        }
        c_allocator.free(c_filenames);
        batch_result.deinit();
        setError(error_code, ZSASA_ERROR_OUT_OF_MEMORY);
        return null;
    };

    handle.* = .{
        .result = batch_result,
        .c_filenames = c_filenames,
    };

    setError(error_code, ZSASA_OK);
    return handle;
}

/// Get the total number of files found in the directory.
///
/// Parameters:
///   handle: Handle returned by zsasa_batch_dir_process()
///
/// Returns:
///   Number of files, or 0 if handle is null.
export fn zsasa_batch_dir_get_total_files(
    handle: ?*anyopaque,
) callconv(.c) usize {
    if (handle == null) return 0;
    const h: *BatchDirHandle = @ptrCast(@alignCast(handle.?));
    return h.result.total_files;
}

/// Get the number of successfully processed files.
///
/// Parameters:
///   handle: Handle returned by zsasa_batch_dir_process()
///
/// Returns:
///   Number of successful files, or 0 if handle is null.
export fn zsasa_batch_dir_get_successful(
    handle: ?*anyopaque,
) callconv(.c) usize {
    if (handle == null) return 0;
    const h: *BatchDirHandle = @ptrCast(@alignCast(handle.?));
    return h.result.successful;
}

/// Get the number of failed files.
///
/// Parameters:
///   handle: Handle returned by zsasa_batch_dir_process()
///
/// Returns:
///   Number of failed files, or 0 if handle is null.
export fn zsasa_batch_dir_get_failed(
    handle: ?*anyopaque,
) callconv(.c) usize {
    if (handle == null) return 0;
    const h: *BatchDirHandle = @ptrCast(@alignCast(handle.?));
    return h.result.failed;
}

/// Get the filename of a processed file by index (0-based).
///
/// Parameters:
///   handle: Handle returned by zsasa_batch_dir_process()
///   index: 0-based file index
///
/// Returns:
///   Null-terminated filename string, or null if handle is null or index is out of bounds.
///   The returned string is owned by the handle and valid until zsasa_batch_dir_free().
export fn zsasa_batch_dir_get_filename(
    handle: ?*anyopaque,
    index: usize,
) callconv(.c) ?[*:0]const u8 {
    if (handle == null) return null;
    const h: *BatchDirHandle = @ptrCast(@alignCast(handle.?));
    if (index >= h.result.file_results.len) return null;
    return h.c_filenames[index];
}

/// Get the number of atoms in a processed file by index (0-based).
///
/// Parameters:
///   handle: Handle returned by zsasa_batch_dir_process()
///   index: 0-based file index
///
/// Returns:
///   Number of atoms, or 0 if handle is null or index is out of bounds.
export fn zsasa_batch_dir_get_n_atoms(
    handle: ?*anyopaque,
    index: usize,
) callconv(.c) usize {
    if (handle == null) return 0;
    const h: *BatchDirHandle = @ptrCast(@alignCast(handle.?));
    if (index >= h.result.file_results.len) return 0;
    return h.result.file_results[index].n_atoms;
}

/// Get the total SASA of a processed file by index (0-based).
///
/// Parameters:
///   handle: Handle returned by zsasa_batch_dir_process()
///   index: 0-based file index
///
/// Returns:
///   Total SASA in Angstroms², or NaN if handle is null or index is out of bounds.
export fn zsasa_batch_dir_get_total_sasa(
    handle: ?*anyopaque,
    index: usize,
) callconv(.c) f64 {
    if (handle == null) return std.math.nan(f64);
    const h: *BatchDirHandle = @ptrCast(@alignCast(handle.?));
    if (index >= h.result.file_results.len) return std.math.nan(f64);
    const fr = h.result.file_results[index];
    if (fr.status == .err) return std.math.nan(f64);
    return fr.total_sasa;
}

/// Get the processing status of a file by index (0-based).
///
/// Parameters:
///   handle: Handle returned by zsasa_batch_dir_process()
///   index: 0-based file index
///
/// Returns:
///   1 = ok, 0 = failed, -1 = out of bounds or null handle.
export fn zsasa_batch_dir_get_status(
    handle: ?*anyopaque,
    index: usize,
) callconv(.c) c_int {
    if (handle == null) return -1;
    const h: *BatchDirHandle = @ptrCast(@alignCast(handle.?));
    if (index >= h.result.file_results.len) return -1;
    return switch (h.result.file_results[index].status) {
        .ok => 1,
        .err => 0,
    };
}

/// Free resources associated with a batch directory handle.
///
/// Parameters:
///   handle: Handle returned by zsasa_batch_dir_process(), or null (no-op).
export fn zsasa_batch_dir_free(
    handle: ?*anyopaque,
) callconv(.c) void {
    if (handle) |h| {
        const dir_handle: *BatchDirHandle = @ptrCast(@alignCast(h));
        dir_handle.deinit();
    }
}

// =============================================================================
// Batch Directory Processing Tests
// =============================================================================

test "zsasa_batch_dir_process null handle safety" {
    // All accessors return safe sentinels for null handle
    try std.testing.expectEqual(@as(usize, 0), zsasa_batch_dir_get_total_files(null));
    try std.testing.expectEqual(@as(usize, 0), zsasa_batch_dir_get_successful(null));
    try std.testing.expectEqual(@as(usize, 0), zsasa_batch_dir_get_failed(null));
    try std.testing.expect(zsasa_batch_dir_get_filename(null, 0) == null);
    try std.testing.expectEqual(@as(usize, 0), zsasa_batch_dir_get_n_atoms(null, 0));
    try std.testing.expect(std.math.isNan(zsasa_batch_dir_get_total_sasa(null, 0)));
    try std.testing.expectEqual(@as(c_int, -1), zsasa_batch_dir_get_status(null, 0));

    // Free null handle is a no-op
    zsasa_batch_dir_free(null);
}

test "zsasa_batch_dir_process invalid parameters" {
    var error_code: c_int = ZSASA_OK;

    // Null input_dir
    const h1 = zsasa_batch_dir_process(null, null, ZSASA_ALGORITHM_SR, 100, 1.4, 0, -1, 0, 0, &error_code);
    try std.testing.expect(h1 == null);
    try std.testing.expectEqual(ZSASA_ERROR_INVALID_INPUT, error_code);

    // n_points = 0
    const h2 = zsasa_batch_dir_process(".", null, ZSASA_ALGORITHM_SR, 0, 1.4, 0, -1, 0, 0, &error_code);
    try std.testing.expect(h2 == null);
    try std.testing.expectEqual(ZSASA_ERROR_INVALID_INPUT, error_code);

    // probe_radius <= 0
    const h3 = zsasa_batch_dir_process(".", null, ZSASA_ALGORITHM_SR, 100, 0.0, 0, -1, 0, 0, &error_code);
    try std.testing.expect(h3 == null);
    try std.testing.expectEqual(ZSASA_ERROR_INVALID_INPUT, error_code);

    // Invalid algorithm
    const h4 = zsasa_batch_dir_process(".", null, 99, 100, 1.4, 0, -1, 0, 0, &error_code);
    try std.testing.expect(h4 == null);
    try std.testing.expectEqual(ZSASA_ERROR_INVALID_INPUT, error_code);

    // Invalid classifier type
    const h5 = zsasa_batch_dir_process(".", null, ZSASA_ALGORITHM_SR, 100, 1.4, 0, 99, 0, 0, &error_code);
    try std.testing.expect(h5 == null);
    try std.testing.expectEqual(ZSASA_ERROR_INVALID_INPUT, error_code);
}

test "zsasa_batch_dir_process nonexistent directory" {
    var error_code: c_int = ZSASA_OK;

    const handle = zsasa_batch_dir_process(
        "/nonexistent/path/that/does/not/exist",
        null,
        ZSASA_ALGORITHM_SR,
        100,
        1.4,
        1,
        -1,
        0,
        0,
        &error_code,
    );
    try std.testing.expect(handle == null);
    try std.testing.expectEqual(ZSASA_ERROR_FILE_IO, error_code);
}

test "zsasa_batch_dir_process null error_code" {
    // Should not crash when error_code is null
    const handle = zsasa_batch_dir_process(null, null, ZSASA_ALGORITHM_SR, 100, 1.4, 0, -1, 0, 0, null);
    try std.testing.expect(handle == null);
}

/// One file of a batch-directory test fixture.
const BatchFixtureFile = struct {
    name: []const u8,
    data: []const u8,
};

/// Write `files` into the temporary directory and return its absolute path as a
/// null-terminated string, ready to hand to `zsasa_batch_dir_process`.
/// The caller frees the result with `std.testing.allocator`.
fn writeBatchFixtureDir(tmp: *std.testing.TmpDir, files: []const BatchFixtureFile) ![:0]u8 {
    for (files) |file| {
        try tmp.dir.writeFile(std.testing.io, .{ .sub_path = file.name, .data = file.data });
    }
    var buf: [std.fs.max_path_bytes]u8 = undefined;
    const len = try tmp.dir.realPath(std.testing.io, &buf);
    return std.testing.allocator.dupeZ(u8, buf[0..len]);
}

/// Index of the result entry whose filename is `name`, or null.
fn findBatchDirFile(handle: ?*anyopaque, name: []const u8) ?usize {
    for (0..zsasa_batch_dir_get_total_files(handle)) |i| {
        const filename = zsasa_batch_dir_get_filename(handle, i) orelse continue;
        if (std.mem.eql(u8, std.mem.span(filename), name)) return i;
    }
    return null;
}

/// Assert that `name` is a successful result with `n_atoms` atoms and the given total SASA.
fn expectBatchDirFile(handle: ?*anyopaque, name: []const u8, n_atoms: usize, sasa: f64) !void {
    const i = findBatchDirFile(handle, name) orelse {
        std.debug.print("batch result has no entry named {s}\n", .{name});
        return error.TestUnexpectedResult;
    };
    try std.testing.expectEqual(@as(c_int, 1), zsasa_batch_dir_get_status(handle, i));
    try std.testing.expectEqual(n_atoms, zsasa_batch_dir_get_n_atoms(handle, i));
    try std.testing.expectApproxEqAbs(sasa, zsasa_batch_dir_get_total_sasa(handle, i), 1e-6);
}

// Small inputs for the batch directory tests. The reference SASA values below
// were computed once with `zsasa calc --classifier=naccess` (probe 1.4 A,
// 100 test points for SR, 20 slices for LR, hydrogens and HETATM excluded).
const batch_fixture_ala_pdb =
    "ATOM      1  N   ALA A   1      -0.966   0.493   1.500  1.00  0.00           N\n" ++
    "ATOM      2  CA  ALA A   1       0.257   0.418   0.692  1.00  0.00           C\n" ++
    "ATOM      3  C   ALA A   1      -0.094   0.017  -0.716  1.00  0.00           C\n" ++
    "ATOM      4  O   ALA A   1      -1.056  -0.682  -0.923  1.00  0.00           O\n" ++
    "ATOM      5  CB  ALA A   1       1.204  -0.620   1.296  1.00  0.00           C\n" ++
    "END\n";
const batch_fixture_gly_ent =
    "ATOM      1  N   GLY A   1      10.000  10.000  10.000  1.00  0.00           N\n" ++
    "ATOM      2  CA  GLY A   1      11.450  10.000  10.000  1.00  0.00           C\n" ++
    "ATOM      3  C   GLY A   1      11.980  11.420  10.000  1.00  0.00           C\n" ++
    "ATOM      4  O   GLY A   1      11.230  12.390  10.000  1.00  0.00           O\n" ++
    "END\n";
const batch_fixture_ethanol_sdf =
    "ethanol\n" ++
    "     zsasa   3D\n" ++
    "\n" ++
    "  9  8  0  0  0  0  0  0  0  0999 V2000\n" ++
    "    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" ++
    "    1.5200    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" ++
    "    2.0800    1.2124    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" ++
    "   -0.5200    0.9400    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n" ++
    "   -0.5200   -0.5100    0.8900 H   0  0  0  0  0  0  0  0  0  0  0  0\n" ++
    "   -0.5200   -0.5100   -0.8900 H   0  0  0  0  0  0  0  0  0  0  0  0\n" ++
    "    1.8800   -0.5100    0.8900 H   0  0  0  0  0  0  0  0  0  0  0  0\n" ++
    "    1.8800   -0.5100   -0.8900 H   0  0  0  0  0  0  0  0  0  0  0  0\n" ++
    "    2.9200    1.2124    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n" ++
    "  1  2  1  0  0  0  0\n" ++
    "  1  4  1  0  0  0  0\n" ++
    "  1  5  1  0  0  0  0\n" ++
    "  1  6  1  0  0  0  0\n" ++
    "  2  3  1  0  0  0  0\n" ++
    "  2  7  1  0  0  0  0\n" ++
    "  2  8  1  0  0  0  0\n" ++
    "  3  9  1  0  0  0  0\n" ++
    "M  END\n" ++
    "$$$$\n";

const batch_fixture_files = [_]BatchFixtureFile{
    .{ .name = "ala.pdb", .data = batch_fixture_ala_pdb },
    .{ .name = "gly.ent", .data = batch_fixture_gly_ent },
    .{ .name = "ethanol.sdf", .data = batch_fixture_ethanol_sdf },
    .{ .name = "notes.txt", .data = "not a structure file\n" },
};

test "zsasa_batch_dir_process with small directory" {
    var tmp_dir = std.testing.tmpDir(.{});
    defer tmp_dir.cleanup();
    const input_dir = try writeBatchFixtureDir(&tmp_dir, &batch_fixture_files);
    defer std.testing.allocator.free(input_dir);

    var error_code: c_int = -999;
    const handle = zsasa_batch_dir_process(
        input_dir,
        null, // no file output
        ZSASA_ALGORITHM_SR,
        100,
        1.4,
        1, // single thread for test determinism
        ZSASA_CLASSIFIER_PROTOR,
        0, // exclude hydrogens
        0, // exclude HETATM
        &error_code,
    );
    try std.testing.expect(handle != null);
    defer zsasa_batch_dir_free(handle);
    try std.testing.expectEqual(ZSASA_OK, error_code);

    // notes.txt is not a supported structure file and is skipped.
    const total = zsasa_batch_dir_get_total_files(handle);
    try std.testing.expectEqual(@as(usize, 3), total);
    try std.testing.expectEqual(@as(usize, 3), zsasa_batch_dir_get_successful(handle));
    try std.testing.expectEqual(@as(usize, 0), zsasa_batch_dir_get_failed(handle));
    try std.testing.expect(findBatchDirFile(handle, "notes.txt") == null);

    try expectBatchDirFile(handle, "ala.pdb", 5, 214.19937915667546);
    try expectBatchDirFile(handle, "gly.ent", 4, 185.4280955820519);
    try expectBatchDirFile(handle, "ethanol_ethanol", 3, 169.50525994296797);

    // Out-of-bounds access returns safe sentinels
    try std.testing.expect(zsasa_batch_dir_get_filename(handle, total) == null);
    try std.testing.expectEqual(@as(usize, 0), zsasa_batch_dir_get_n_atoms(handle, total));
    try std.testing.expect(std.math.isNan(zsasa_batch_dir_get_total_sasa(handle, total)));
    try std.testing.expectEqual(@as(c_int, -1), zsasa_batch_dir_get_status(handle, total));
}

test "zsasa_batch_dir_process with LR algorithm" {
    var tmp_dir = std.testing.tmpDir(.{});
    defer tmp_dir.cleanup();
    const input_dir = try writeBatchFixtureDir(&tmp_dir, batch_fixture_files[0..2]);
    defer std.testing.allocator.free(input_dir);

    var error_code: c_int = -999;
    const lr = zsasa_batch_dir_process(
        input_dir,
        null,
        ZSASA_ALGORITHM_LR,
        20, // n_slices for LR
        1.4,
        1,
        ZSASA_CLASSIFIER_NACCESS,
        0,
        0,
        &error_code,
    );
    try std.testing.expect(lr != null);
    defer zsasa_batch_dir_free(lr);
    try std.testing.expectEqual(ZSASA_OK, error_code);

    try std.testing.expectEqual(@as(usize, 2), zsasa_batch_dir_get_total_files(lr));
    try std.testing.expectEqual(@as(usize, 2), zsasa_batch_dir_get_successful(lr));
    try std.testing.expectEqual(@as(usize, 0), zsasa_batch_dir_get_failed(lr));

    // Reference values from `zsasa calc --algorithm=lr --classifier=naccess`
    // (exact arc angles, the default).
    try expectBatchDirFile(lr, "ala.pdb", 5, 215.07786201196504);
    try expectBatchDirFile(lr, "gly.ent", 4, 188.2086929717762);

    // Shrake-Rupley on the same directory gives a different (but close) area,
    // so the algorithm argument is really honored.
    const sr = zsasa_batch_dir_process(
        input_dir,
        null,
        ZSASA_ALGORITHM_SR,
        100,
        1.4,
        1,
        ZSASA_CLASSIFIER_NACCESS,
        0,
        0,
        &error_code,
    );
    try std.testing.expect(sr != null);
    defer zsasa_batch_dir_free(sr);
    try std.testing.expectEqual(ZSASA_OK, error_code);

    // Reference values from `zsasa calc --algorithm=sr --classifier=naccess`.
    try expectBatchDirFile(sr, "ala.pdb", 5, 215.8249020274959);
    try expectBatchDirFile(sr, "gly.ent", 4, 187.0015182803052);
    for ([_][]const u8{ "ala.pdb", "gly.ent" }) |name| {
        const lr_area = zsasa_batch_dir_get_total_sasa(lr, findBatchDirFile(lr, name).?);
        const sr_area = zsasa_batch_dir_get_total_sasa(sr, findBatchDirFile(sr, name).?);
        try std.testing.expect(@abs(lr_area - sr_area) > 0.1);
        try std.testing.expect(@abs(lr_area - sr_area) < 0.02 * sr_area);
    }
}

/// Total area of atoms that are far enough apart to be isolated spheres:
/// the sum of 4 pi (r + probe)^2.
fn isolatedSpheresArea(radii: []const f64, probe: f64) f64 {
    var total: f64 = 0;
    for (radii) |r| total += 4.0 * std.math.pi * (r + probe) * (r + probe);
    return total;
}

test "zsasa_batch_dir_process classifier_type -1 uses input radii" {
    // Atoms 100 A apart never touch, so each one exposes its whole sphere.
    // A JSON file carries explicit radii. In a PDB file the radii come from the
    // parser (the element's van der Waals radius), which differs from the
    // NACCESS values a classifier would assign to the same atoms.
    var tmp_dir = std.testing.tmpDir(.{});
    defer tmp_dir.cleanup();
    const input_dir = try writeBatchFixtureDir(&tmp_dir, &.{
        .{
            .name = "radii.json",
            .data = "{\"x\":[0.0,100.0],\"y\":[0.0,0.0],\"z\":[0.0,0.0],\"r\":[1.5,2.0]}",
        },
        .{
            .name = "isolated.pdb",
            .data = "ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00  0.00           N\n" ++
                "ATOM      2  C   ALA A   1     100.000   0.000   0.000  1.00  0.00           C\n" ++
                "ATOM      3  O   ALA A   1     200.000   0.000   0.000  1.00  0.00           O\n" ++
                "END\n",
        },
    });
    defer std.testing.allocator.free(input_dir);

    const probe: f64 = 1.4;
    const json_area = isolatedSpheresArea(&.{ 1.5, 2.0 }, probe);
    const input_radii_area = isolatedSpheresArea(&.{ 1.55, 1.70, 1.52 }, probe); // N, C, O van der Waals
    const naccess_area = isolatedSpheresArea(&.{ 1.65, 1.76, 1.40 }, probe); // ALA N, C, O

    var error_code: c_int = -999;
    const handle = zsasa_batch_dir_process(
        input_dir,
        null,
        ZSASA_ALGORITHM_SR,
        100,
        probe,
        1,
        -1, // use radii from input files
        0,
        0,
        &error_code,
    );
    try std.testing.expect(handle != null);
    defer zsasa_batch_dir_free(handle);
    try std.testing.expectEqual(ZSASA_OK, error_code);

    try std.testing.expectEqual(@as(usize, 2), zsasa_batch_dir_get_total_files(handle));
    try std.testing.expectEqual(@as(usize, 2), zsasa_batch_dir_get_successful(handle));
    try std.testing.expectEqual(@as(usize, 0), zsasa_batch_dir_get_failed(handle));
    try expectBatchDirFile(handle, "radii.json", 2, json_area);
    try expectBatchDirFile(handle, "isolated.pdb", 3, input_radii_area);

    // Guard against a vacuous comparison: the NACCESS classifier would give the
    // PDB file another area, so matching the input radii proves no classifier ran.
    try std.testing.expect(@abs(naccess_area - input_radii_area) > 1.0);
}

// =============================================================================
// ABI Version, Error Mapping and Shared Io Tests
// =============================================================================

test "zsasa_abi_version is a positive constant" {
    try std.testing.expect(zsasa_abi_version() >= 1);
    try std.testing.expectEqual(zsasa_abi_version(), zsasa_abi_version());
}

test "trajectoryErrorCode separates missing, malformed and out-of-memory" {
    try std.testing.expectEqual(ZSASA_ERROR_INVALID_INPUT, trajectoryErrorCode(error.FileNotFound));
    for ([_]anyerror{
        error.InvalidMagic,
        error.BadFormat,
        error.EndOfFile,
        error.ReadError,
        error.DecompressionError,
        error.FixedAtomsNotSupported,
    }) |err| {
        try std.testing.expectEqual(ZSASA_ERROR_INVALID_FORMAT, trajectoryErrorCode(err));
    }
    try std.testing.expectEqual(ZSASA_ERROR_OUT_OF_MEMORY, trajectoryErrorCode(error.OutOfMemory));
    try std.testing.expectEqual(ZSASA_ERROR_CALCULATION, trajectoryErrorCode(error.Unexpected));
}

test "batchErrorCode reads the same filesystem error by stage" {
    try std.testing.expectEqual(ZSASA_ERROR_FILE_IO, batchErrorCode(error.NotDir, .scan_inputs));
    try std.testing.expectEqual(ZSASA_ERROR_OUTPUT_DIR, batchErrorCode(error.NotDir, .create_output_dir));
    try std.testing.expectEqual(ZSASA_ERROR_FILE_IO, batchErrorCode(error.AccessDenied, .process));
    try std.testing.expectEqual(ZSASA_ERROR_OUT_OF_MEMORY, batchErrorCode(error.OutOfMemory, .create_output_dir));
    try std.testing.expectEqual(ZSASA_ERROR_OUTPUT_NAME_COLLISION, batchErrorCode(error.OutputNameCollision, .scan_inputs));
}

/// Open `path` with `open` (zsasa_xtc_open or zsasa_dcd_open) and return the error code.
/// Fails the test when the file unexpectedly opens.
fn trajectoryOpenError(
    comptime open: fn ([*:0]const u8, *i32, *c_int) callconv(.c) ?*anyopaque,
    comptime close: fn (?*anyopaque) callconv(.c) void,
    path: [*:0]const u8,
) !c_int {
    var natoms: i32 = 0;
    var error_code: c_int = ZSASA_OK;
    const handle = open(path, &natoms, &error_code);
    if (handle != null) {
        close(handle);
        return error.TestUnexpectedResult;
    }
    return error_code;
}

test "trajectory open reports a file of another format as invalid format" {
    try std.testing.expectEqual(
        ZSASA_ERROR_INVALID_FORMAT,
        try trajectoryOpenError(zsasa_xtc_open, zsasa_xtc_close, "test_data/1l2y.dcd"),
    );
    try std.testing.expectEqual(
        ZSASA_ERROR_INVALID_FORMAT,
        try trajectoryOpenError(zsasa_dcd_open, zsasa_dcd_close, "test_data/1l2y.xtc"),
    );
    try std.testing.expectEqual(
        ZSASA_ERROR_INVALID_FORMAT,
        try trajectoryOpenError(zsasa_xtc_open, zsasa_xtc_close, "test_data/1l2y.pdb"),
    );
    try std.testing.expectEqual(
        ZSASA_ERROR_INVALID_FORMAT,
        try trajectoryOpenError(zsasa_dcd_open, zsasa_dcd_close, "test_data/1l2y.pdb"),
    );
}

test "trajectory open reports an empty file as invalid format and a missing one as invalid input" {
    var tmp_dir = std.testing.tmpDir(.{});
    defer tmp_dir.cleanup();
    const dir = try writeBatchFixtureDir(&tmp_dir, &.{.{ .name = "empty.bin", .data = "" }});
    defer std.testing.allocator.free(dir);
    const empty_path = try std.fmt.allocPrintSentinel(std.testing.allocator, "{s}/empty.bin", .{dir}, 0);
    defer std.testing.allocator.free(empty_path);
    const missing_path = try std.fmt.allocPrintSentinel(std.testing.allocator, "{s}/missing.bin", .{dir}, 0);
    defer std.testing.allocator.free(missing_path);

    try std.testing.expectEqual(
        ZSASA_ERROR_INVALID_FORMAT,
        try trajectoryOpenError(zsasa_xtc_open, zsasa_xtc_close, empty_path),
    );
    try std.testing.expectEqual(
        ZSASA_ERROR_INVALID_FORMAT,
        try trajectoryOpenError(zsasa_dcd_open, zsasa_dcd_close, empty_path),
    );
    try std.testing.expectEqual(
        ZSASA_ERROR_INVALID_INPUT,
        try trajectoryOpenError(zsasa_xtc_open, zsasa_xtc_close, missing_path),
    );
    try std.testing.expectEqual(
        ZSASA_ERROR_INVALID_INPUT,
        try trajectoryOpenError(zsasa_dcd_open, zsasa_dcd_close, missing_path),
    );
}

test "zsasa_batch_dir_process distinguishes a bad output directory from a bad input directory" {
    var tmp_dir = std.testing.tmpDir(.{});
    defer tmp_dir.cleanup();
    const input_dir = try writeBatchFixtureDir(&tmp_dir, batch_fixture_files[0..1]);
    defer std.testing.allocator.free(input_dir);

    // A path below a regular file can never be created as a directory.
    const blocked_output = try std.fmt.allocPrintSentinel(std.testing.allocator, "{s}/ala.pdb/sub", .{input_dir}, 0);
    defer std.testing.allocator.free(blocked_output);

    // One thread and several: both runners create the output directory.
    for ([_]usize{ 1, 4 }) |n_threads| {
        var error_code: c_int = ZSASA_OK;
        const handle = zsasa_batch_dir_process(
            input_dir,
            blocked_output,
            ZSASA_ALGORITHM_SR,
            100,
            1.4,
            n_threads,
            ZSASA_CLASSIFIER_PROTOR,
            0,
            0,
            &error_code,
        );
        try std.testing.expect(handle == null);
        try std.testing.expectEqual(ZSASA_ERROR_OUTPUT_DIR, error_code);
    }

    // The input directory is the one that is missing, with a usable output directory.
    const output_dir = try std.fmt.allocPrintSentinel(std.testing.allocator, "{s}/out", .{input_dir}, 0);
    defer std.testing.allocator.free(output_dir);
    var error_code: c_int = ZSASA_OK;
    const handle = zsasa_batch_dir_process(
        "/nonexistent/path/that/does/not/exist",
        output_dir,
        ZSASA_ALGORITHM_SR,
        100,
        1.4,
        1,
        ZSASA_CLASSIFIER_PROTOR,
        0,
        0,
        &error_code,
    );
    try std.testing.expect(handle == null);
    try std.testing.expectEqual(ZSASA_ERROR_FILE_IO, error_code);

    // An input directory that is a regular file is also an input problem.
    const file_as_input = try std.fmt.allocPrintSentinel(std.testing.allocator, "{s}/ala.pdb", .{input_dir}, 0);
    defer std.testing.allocator.free(file_as_input);
    const handle2 = zsasa_batch_dir_process(
        file_as_input,
        output_dir,
        ZSASA_ALGORITHM_SR,
        100,
        1.4,
        1,
        ZSASA_CLASSIFIER_PROTOR,
        0,
        0,
        &error_code,
    );
    try std.testing.expect(handle2 == null);
    try std.testing.expectEqual(ZSASA_ERROR_FILE_IO, error_code);
}

const SharedIoProbe = struct {
    userdata: ?*anyopaque = null,

    fn run(self: *SharedIoProbe) void {
        self.userdata = sharedThreadedIo().userdata;
    }
};

test "sharedThreadedIo creates one instance for concurrent first calls" {
    var probes: [8]SharedIoProbe = @splat(.{});
    var threads: [probes.len]std.Thread = undefined;
    for (&threads, &probes) |*thread, *probe| {
        thread.* = try std.Thread.spawn(.{}, SharedIoProbe.run, .{probe});
    }
    for (threads) |thread| thread.join();

    for (probes) |probe| {
        try std.testing.expect(probe.userdata != null);
        try std.testing.expectEqual(probes[0].userdata, probe.userdata);
    }
    // Later calls get the same instance too: nothing deinitializes it.
    try std.testing.expectEqual(probes[0].userdata, sharedThreadedIo().userdata);
}

const ConcurrentBatchRun = struct {
    input_dir: [*:0]const u8,
    error_code: c_int = -999,
    successful: usize = 0,

    fn run(self: *ConcurrentBatchRun) void {
        const handle = zsasa_batch_dir_process(
            self.input_dir,
            null,
            ZSASA_ALGORITHM_SR,
            100,
            1.4,
            2,
            ZSASA_CLASSIFIER_PROTOR,
            0,
            0,
            &self.error_code,
        );
        defer zsasa_batch_dir_free(handle);
        self.successful = zsasa_batch_dir_get_successful(handle);
    }
};

fn currentSigPipeHandler() ?*const anyopaque {
    var current: std.posix.Sigaction = undefined;
    std.posix.sigaction(.PIPE, null, &current);
    return @ptrCast(current.handler.handler);
}

test "zsasa_batch_dir_process is safe to call from several threads at once" {
    var tmp_dir = std.testing.tmpDir(.{});
    defer tmp_dir.cleanup();
    const input_dir = try writeBatchFixtureDir(&tmp_dir, batch_fixture_files[0..2]);
    defer std.testing.allocator.free(input_dir);

    // Once the shared Io exists, the SIGPIPE handler Threaded installed must
    // stay the same through any number of overlapping calls: no call may restore
    // a handler it found installed by another one.
    const have_handlers = std.posix.Sigaction != void;
    _ = sharedThreadedIo();
    const installed = if (have_handlers) currentSigPipeHandler() else null;

    var runs: [6]ConcurrentBatchRun = @splat(.{ .input_dir = input_dir.ptr });
    var threads: [runs.len]std.Thread = undefined;
    for (&threads, &runs) |*thread, *run| {
        thread.* = try std.Thread.spawn(.{}, ConcurrentBatchRun.run, .{run});
    }
    for (threads) |thread| thread.join();

    for (runs) |run| {
        try std.testing.expectEqual(ZSASA_OK, run.error_code);
        try std.testing.expectEqual(@as(usize, 2), run.successful);
    }
    if (have_handlers) try std.testing.expectEqual(installed, currentSigPipeHandler());
}
