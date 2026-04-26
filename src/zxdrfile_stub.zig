//! Stub for zxdrfile when -Dxtc=false. All XTC code paths must be gated
//! with `if (build_options.xtc)` so this stub is never executed.
//! Re-enable XTC support by passing -Dxtc=true (after zxdrfile 0.16 lands).

const std = @import("std");

pub const xtc = struct {
    pub const XtcError = error{
        EndOfFile,
        FileNotFound,
        InvalidMagic,
        OutOfMemory,
        DecompressionError,
        Unsupported,
    };

    pub const XtcFrame = struct {
        step: i32 = 0,
        time: f32 = 0,
        precision: f32 = 0,
        coords: []f32 = &.{},
        box: [3][3]f32 = .{.{0} ** 3} ** 3,

        pub fn deinit(self: *XtcFrame, allocator: std.mem.Allocator) void {
            _ = self;
            _ = allocator;
        }
    };

    pub const XtcReader = struct {
        natoms: i32 = 0,

        pub fn open(allocator: std.mem.Allocator, path: []const u8) XtcError!XtcReader {
            _ = allocator;
            _ = path;
            return error.Unsupported;
        }

        pub fn close(self: *XtcReader) void {
            _ = self;
        }

        pub fn readFrame(self: *XtcReader) XtcError!XtcFrame {
            _ = self;
            return error.Unsupported;
        }

        pub fn getNumAtoms(self: *const XtcReader) i32 {
            _ = self;
            return 0;
        }
    };
};
