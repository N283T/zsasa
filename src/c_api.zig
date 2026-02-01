//! C ABI interface for zsasa library.
//!
//! This module provides C-compatible functions that can be called from
//! other languages (Python, C, etc.) via FFI/ctypes.

const std = @import("std");
const shrake_rupley = @import("shrake_rupley.zig");
const lee_richards = @import("lee_richards.zig");
const types = @import("types.zig");
const classifier = @import("classifier.zig");
const classifier_naccess = @import("classifier_naccess.zig");
const classifier_protor = @import("classifier_protor.zig");
const classifier_oons = @import("classifier_oons.zig");
const analysis = @import("analysis.zig");

const AtomInput = types.AtomInput;
const Config = types.Config;

/// No error - calculation completed successfully
pub const ZSASA_OK: c_int = 0;
/// Invalid input parameters (n_atoms=0, n_points=0, n_slices=0, or invalid probe_radius)
/// Note: Passing NULL pointers results in undefined behavior
pub const ZSASA_ERROR_INVALID_INPUT: c_int = -1;
/// Memory allocation failed during calculation
pub const ZSASA_ERROR_OUT_OF_MEMORY: c_int = -2;
/// Internal calculation error
pub const ZSASA_ERROR_CALCULATION: c_int = -3;

// =============================================================================
// Classifier Types
// =============================================================================

/// NACCESS-compatible radii (default)
pub const ZSASA_CLASSIFIER_NACCESS: c_int = 0;
/// ProtOr radii based on hybridization state
pub const ZSASA_CLASSIFIER_PROTOR: c_int = 1;
/// OONS radii (older FreeSASA default)
pub const ZSASA_CLASSIFIER_OONS: c_int = 2;

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
const VERSION = "0.1.0";

/// Thread-safe allocator for C API (uses C allocator for simplicity)
const c_allocator = std.heap.c_allocator;

/// Get library version string.
export fn zsasa_version() callconv(.c) [*:0]const u8 {
    return VERSION;
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
    if (n_atoms == 0 or n_points == 0 or probe_radius <= 0.0) {
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
        shrake_rupley.calculateSasa(c_allocator, input, config) catch {
            return ZSASA_ERROR_CALCULATION;
        }
    else
        shrake_rupley.calculateSasaParallel(c_allocator, input, config, n_threads) catch {
            return ZSASA_ERROR_CALCULATION;
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
    error_flag: *std.atomic.Value(bool),
    thread_id: usize,
    n_threads: usize,
    algorithm: BatchAlgorithm,
};

/// Worker function for batch processing (shared by SR and LR)
fn batchWorkerFn(args: BatchWorkerArgs) void {
    // Pre-allocate coordinate buffers once per thread (reused across frames)
    const x = c_allocator.alloc(f64, args.n_atoms) catch {
        args.error_flag.store(true, .release);
        return;
    };
    defer c_allocator.free(x);

    const y = c_allocator.alloc(f64, args.n_atoms) catch {
        args.error_flag.store(true, .release);
        return;
    };
    defer c_allocator.free(y);

    const z = c_allocator.alloc(f64, args.n_atoms) catch {
        args.error_flag.store(true, .release);
        return;
    };
    defer c_allocator.free(z);

    // Each thread processes frames: thread_id, thread_id + n_threads, ...
    var frame_idx = args.thread_id;
    while (frame_idx < args.n_frames) : (frame_idx += args.n_threads) {
        // Skip if error already occurred
        if (args.error_flag.load(.acquire)) return;

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
                break :blk shrake_rupley.calculateSasa(c_allocator, input, config) catch {
                    args.error_flag.store(true, .release);
                    return;
                };
            },
            .lee_richards => blk: {
                const config = lee_richards.LeeRichardsConfig{
                    .n_slices = args.param,
                    .probe_radius = args.probe_radius,
                };
                break :blk lee_richards.calculateSasa(c_allocator, input, config) catch {
                    args.error_flag.store(true, .release);
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
    if (n_frames == 0 or n_atoms == 0 or param == 0 or probe_radius <= 0.0) {
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

    // Error flag shared across threads
    var error_flag = std.atomic.Value(bool).init(false);

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
            .error_flag = &error_flag,
            .thread_id = i,
            .n_threads = thread_count,
            .algorithm = algorithm,
        }}) catch {
            // If thread spawn fails, set error flag and wait for already-spawned threads
            error_flag.store(true, .release);
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

    if (error_flag.load(.acquire)) {
        return ZSASA_ERROR_CALCULATION;
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
    error_flag: *std.atomic.Value(bool),
    thread_id: usize,
    n_threads: usize,
    algorithm: BatchAlgorithm,
};

/// Worker function for pure f32 batch processing
fn batchWorkerFnF32(args: BatchWorkerArgsF32) void {
    // Pre-allocate coordinate buffers (f64 for AtomInput compatibility)
    const x = c_allocator.alloc(f64, args.n_atoms) catch {
        args.error_flag.store(true, .release);
        return;
    };
    defer c_allocator.free(x);

    const y = c_allocator.alloc(f64, args.n_atoms) catch {
        args.error_flag.store(true, .release);
        return;
    };
    defer c_allocator.free(y);

    const z = c_allocator.alloc(f64, args.n_atoms) catch {
        args.error_flag.store(true, .release);
        return;
    };
    defer c_allocator.free(z);

    // f64 radii for AtomInput (will be cast to f32 internally)
    const radii_f64 = c_allocator.alloc(f64, args.n_atoms) catch {
        args.error_flag.store(true, .release);
        return;
    };
    defer c_allocator.free(radii_f64);

    for (0..args.n_atoms) |i| {
        radii_f64[i] = @floatCast(args.radii_f32[i]);
    }

    // Each thread processes frames: thread_id, thread_id + n_threads, ...
    var frame_idx = args.thread_id;
    while (frame_idx < args.n_frames) : (frame_idx += args.n_threads) {
        if (args.error_flag.load(.acquire)) return;

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
                break :blk shrake_rupley.calculateSasaf32(c_allocator, input, config) catch {
                    args.error_flag.store(true, .release);
                    return;
                };
            },
            .lee_richards => blk: {
                const config = lee_richards.LeeRichardsConfigGen(f32){
                    .n_slices = args.param,
                    .probe_radius = args.probe_radius,
                };
                break :blk lee_richards.calculateSasaf32(c_allocator, input, config) catch {
                    args.error_flag.store(true, .release);
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
    if (n_frames == 0 or n_atoms == 0 or param == 0 or probe_radius <= 0.0) {
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

    var error_flag = std.atomic.Value(bool).init(false);

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
            .error_flag = &error_flag,
            .thread_id = i,
            .n_threads = thread_count,
            .algorithm = algorithm,
        }}) catch {
            error_flag.store(true, .release);
            for (threads[0..i]) |thread| {
                thread.join();
            }
            return ZSASA_ERROR_CALCULATION;
        };
    }

    for (threads) |thread| {
        thread.join();
    }

    if (error_flag.load(.acquire)) {
        return ZSASA_ERROR_CALCULATION;
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

/// Calculate SASA using Lee-Richards algorithm.
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
    if (n_atoms == 0 or n_slices == 0 or probe_radius <= 0.0) {
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
        lee_richards.calculateSasa(c_allocator, input, config) catch {
            return ZSASA_ERROR_CALCULATION;
        }
    else
        lee_richards.calculateSasaParallel(c_allocator, input, config, n_threads) catch {
            return ZSASA_ERROR_CALCULATION;
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
fn getRadiusByClassifier(classifier_type: c_int, residue: []const u8, atom: []const u8) ?f64 {
    return switch (classifier_type) {
        ZSASA_CLASSIFIER_NACCESS => classifier_naccess.getRadius(residue, atom),
        ZSASA_CLASSIFIER_PROTOR => classifier_protor.getRadius(residue, atom),
        ZSASA_CLASSIFIER_OONS => classifier_oons.getRadius(residue, atom),
        else => null,
    };
}

// Internal helper: get class by classifier type
fn getClassByClassifier(classifier_type: c_int, residue: []const u8, atom: []const u8) classifier.AtomClass {
    return switch (classifier_type) {
        ZSASA_CLASSIFIER_NACCESS => classifier_naccess.getClass(residue, atom),
        ZSASA_CLASSIFIER_PROTOR => classifier_protor.getClass(residue, atom),
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
    if (classifier_type < ZSASA_CLASSIFIER_NACCESS or classifier_type > ZSASA_CLASSIFIER_OONS) {
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
