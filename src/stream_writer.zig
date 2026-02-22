const std = @import("std");
const batch = @import("batch.zig");

/// Output format for streaming batch results
pub const StreamFormat = enum {
    ndjson, // Newline-delimited JSON: one JSON object per line
    json_array, // JSON array: all results wrapped in [ ... ]
};

/// Result for a single file, used for streaming output
pub const StreamResult = struct {
    filename: []const u8,
    n_atoms: usize,
    total_sasa: f64,
    status: batch.FileResult.Status,
    time_ns: u64,
    atom_areas: ?[]const f64 = null,
    error_msg: ?[]const u8 = null,

    /// Create a StreamResult from a FileResult.
    pub fn fromFileResult(fr: batch.FileResult) StreamResult {
        return .{
            .filename = fr.filename,
            .n_atoms = fr.n_atoms,
            .total_sasa = fr.total_sasa,
            .status = fr.status,
            .time_ns = fr.sasa_time_ns,
            .error_msg = fr.error_msg,
        };
    }
};

/// Lifecycle phase for StreamWriter protocol validation.
const Phase = enum {
    created,
    begun,
    ended,
};

/// Thread-safe streaming writer for batch SASA results.
///
/// Writes results incrementally as each file is processed,
/// avoiding the need to accumulate all results in memory.
///
/// Lifecycle: init -> begin -> writeResult* -> end
pub const StreamWriter = struct {
    writer: std.io.AnyWriter,
    format: StreamFormat,
    mutex: std.Thread.Mutex,
    first_written: bool,
    include_atom_areas: bool,
    phase: Phase,
    write_failures: usize,

    /// Initialize a StreamWriter.
    ///
    /// Parameters:
    ///   writer            - destination writer (e.g. stdout or file writer)
    ///   format            - output format (.ndjson or .json_array)
    ///   include_atom_areas - when true, per-atom SASA values are emitted if
    ///                        StreamResult.atom_areas is non-null.
    ///                        Currently batch streaming always passes false.
    pub fn init(
        writer: std.io.AnyWriter,
        format: StreamFormat,
        include_atom_areas: bool,
    ) StreamWriter {
        return StreamWriter{
            .writer = writer,
            .format = format,
            .mutex = std.Thread.Mutex{},
            .first_written = false,
            .include_atom_areas = include_atom_areas,
            .phase = .created,
            .write_failures = 0,
        };
    }

    /// Write the opening token(s) for the chosen format.
    ///
    /// For .json_array this emits "[\n"; for .ndjson it is a no-op.
    /// Must be called once before any writeResult calls.
    pub fn begin(self: *StreamWriter) !void {
        std.debug.assert(self.phase == .created);
        switch (self.format) {
            .json_array => try self.writer.writeAll("[\n"),
            .ndjson => {},
        }
        self.phase = .begun;
    }

    /// Serialize one result to the stream in a thread-safe manner.
    ///
    /// The entire JSON object is built in a local buffer before writing
    /// to the underlying writer, ensuring atomicity: either the full
    /// result is written or nothing is (no partial/corrupt objects).
    ///
    /// JSON object fields emitted:
    ///   "file"       - filename (JSON-escaped)
    ///   "n_atoms"    - integer atom count
    ///   "total_sasa" - float, 6 decimal places
    ///   "status"     - "ok" or "error"
    ///   "time_ms"    - float, 3 decimal places (converted from time_ns)
    ///   "error"      - string, only when error_msg is non-null
    ///   "atom_areas" - float array, only when include_atom_areas is true
    ///                  and atom_areas is non-null
    pub fn writeResult(self: *StreamWriter, result: StreamResult) !void {
        self.mutex.lock();
        defer self.mutex.unlock();

        std.debug.assert(self.phase == .begun);

        // Build the complete JSON object in a buffer first to ensure atomicity.
        var buf = std.ArrayList(u8){};
        defer buf.deinit(std.heap.page_allocator);
        const bw = buf.writer(std.heap.page_allocator).any();

        switch (self.format) {
            .json_array => try bw.writeAll("  {"),
            .ndjson => try bw.writeAll("{"),
        }

        try bw.writeAll("\"file\":");
        try writeJsonString(bw, result.filename);

        try bw.print(",\"n_atoms\":{d}", .{result.n_atoms});

        try bw.print(",\"total_sasa\":{d:.6}", .{result.total_sasa});

        const status_str: []const u8 = switch (result.status) {
            .ok => "ok",
            .err => "error",
        };
        try bw.print(",\"status\":\"{s}\"", .{status_str});

        const time_ms = @as(f64, @floatFromInt(result.time_ns)) / 1_000_000.0;
        try bw.print(",\"time_ms\":{d:.3}", .{time_ms});

        if (result.error_msg) |msg| {
            try bw.writeAll(",\"error\":");
            try writeJsonString(bw, msg);
        }

        if (self.include_atom_areas) {
            if (result.atom_areas) |areas| {
                try bw.writeAll(",\"atom_areas\":[");
                for (areas, 0..) |area, i| {
                    if (i > 0) try bw.writeAll(",");
                    try bw.print("{d:.6}", .{area});
                }
                try bw.writeAll("]");
            }
        }

        try bw.writeAll("}");
        if (self.format == .ndjson) {
            try bw.writeAll("\n");
        }

        // Write the buffered object atomically to the real writer.
        // first_written is only updated after a successful write.
        const w = self.writer;
        if (self.format == .json_array and self.first_written) {
            try w.writeAll(",\n");
        }
        try w.writeAll(buf.items);
        self.first_written = true;
    }

    /// Write the closing token(s) for the chosen format.
    ///
    /// For .json_array: emits "\n]\n" when entries were written, or "]\n"
    /// for an empty array. For .ndjson it is a no-op.
    /// Must be called once after all writeResult calls.
    pub fn end(self: *StreamWriter) !void {
        std.debug.assert(self.phase == .begun);
        switch (self.format) {
            .json_array => {
                if (self.first_written) {
                    try self.writer.writeAll("\n]\n");
                } else {
                    try self.writer.writeAll("]\n");
                }
            },
            .ndjson => {},
        }
        self.phase = .ended;
    }
};

/// Write a JSON-escaped string surrounded by double quotes.
///
/// Handles escape sequences required by RFC 8259:
///   \" \\ \n \r \t           (named escapes)
///   \u00XX                    (all other control characters U+0000..U+001F)
/// All other bytes are written as-is (assumes valid UTF-8 input).
pub fn writeJsonString(writer: std.io.AnyWriter, s: []const u8) !void {
    try writer.writeAll("\"");
    for (s) |c| {
        switch (c) {
            '"' => try writer.writeAll("\\\""),
            '\\' => try writer.writeAll("\\\\"),
            '\n' => try writer.writeAll("\\n"),
            '\r' => try writer.writeAll("\\r"),
            '\t' => try writer.writeAll("\\t"),
            0x00...0x08, 0x0B, 0x0C, 0x0E...0x1F => {
                try writer.print("\\u{x:0>4}", .{@as(u16, c)});
            },
            else => try writer.writeByte(c),
        }
    }
    try writer.writeAll("\"");
}

// ---- Tests ---------------------------------------------------------------

test "NDJSON writes one result per line" {
    var buf = std.ArrayList(u8){};
    defer buf.deinit(std.testing.allocator);

    var sw = StreamWriter.init(buf.writer(std.testing.allocator).any(), .ndjson, false);
    try sw.begin();
    try sw.writeResult(.{
        .filename = "test.pdb",
        .n_atoms = 42,
        .total_sasa = 1234.567890,
        .status = .ok,
        .time_ns = 500_000, // 0.5 ms
    });
    try sw.end();

    const output = buf.items;

    // Must end with newline
    try std.testing.expect(output.len > 0);
    try std.testing.expectEqual('\n', output[output.len - 1]);

    // Single line: only one newline at the very end
    const newline_count = std.mem.count(u8, output, "\n");
    try std.testing.expectEqual(@as(usize, 1), newline_count);

    // NDJSON: no leading whitespace
    try std.testing.expectEqual('{', output[0]);

    // Must contain expected field values
    try std.testing.expect(std.mem.indexOf(u8, output, "\"file\":\"test.pdb\"") != null);
    try std.testing.expect(std.mem.indexOf(u8, output, "\"n_atoms\":42") != null);
    try std.testing.expect(std.mem.indexOf(u8, output, "\"total_sasa\":1234.567890") != null);
    try std.testing.expect(std.mem.indexOf(u8, output, "\"status\":\"ok\"") != null);
    try std.testing.expect(std.mem.indexOf(u8, output, "\"time_ms\":0.500") != null);
}

test "JSON array wraps results in brackets" {
    var buf = std.ArrayList(u8){};
    defer buf.deinit(std.testing.allocator);

    var sw = StreamWriter.init(buf.writer(std.testing.allocator).any(), .json_array, false);
    try sw.begin();
    try sw.writeResult(.{
        .filename = "a.pdb",
        .n_atoms = 10,
        .total_sasa = 100.0,
        .status = .ok,
        .time_ns = 1_000_000,
    });
    try sw.writeResult(.{
        .filename = "b.pdb",
        .n_atoms = 20,
        .total_sasa = 200.0,
        .status = .ok,
        .time_ns = 2_000_000,
    });
    try sw.end();

    const output = buf.items;

    // Must start with '[' and end with ']'
    try std.testing.expect(std.mem.startsWith(u8, output, "["));
    try std.testing.expect(std.mem.endsWith(u8, std.mem.trimRight(u8, output, "\n"), "]"));

    // Both files must be present
    try std.testing.expect(std.mem.indexOf(u8, output, "\"file\":\"a.pdb\"") != null);
    try std.testing.expect(std.mem.indexOf(u8, output, "\"file\":\"b.pdb\"") != null);

    // The two entries must be separated by a comma
    try std.testing.expect(std.mem.indexOf(u8, output, "},\n") != null);
}

test "JSON array with single result has no trailing comma" {
    var buf = std.ArrayList(u8){};
    defer buf.deinit(std.testing.allocator);

    var sw = StreamWriter.init(buf.writer(std.testing.allocator).any(), .json_array, false);
    try sw.begin();
    try sw.writeResult(.{
        .filename = "only.pdb",
        .n_atoms = 5,
        .total_sasa = 50.0,
        .status = .ok,
        .time_ns = 100_000,
    });
    try sw.end();

    const output = buf.items;

    // Must be valid JSON array
    try std.testing.expect(std.mem.startsWith(u8, output, "["));
    try std.testing.expect(std.mem.endsWith(u8, std.mem.trimRight(u8, output, "\n"), "]"));

    // No comma in output (single result)
    try std.testing.expect(std.mem.indexOf(u8, output, "},\n") == null);
    try std.testing.expect(std.mem.indexOf(u8, output, "\"file\":\"only.pdb\"") != null);
}

test "NDJSON with atom_areas when full mode enabled" {
    var buf = std.ArrayList(u8){};
    defer buf.deinit(std.testing.allocator);

    const areas = [_]f64{ 10.5, 20.3, 5.0 };

    var sw = StreamWriter.init(buf.writer(std.testing.allocator).any(), .ndjson, true);
    try sw.begin();
    try sw.writeResult(.{
        .filename = "full.pdb",
        .n_atoms = 3,
        .total_sasa = 35.8,
        .status = .ok,
        .time_ns = 1_000_000,
        .atom_areas = &areas,
    });
    try sw.end();

    const output = buf.items;

    // atom_areas field must be present
    try std.testing.expect(std.mem.indexOf(u8, output, "\"atom_areas\":[10.500000,20.300000,5.000000]") != null);
}

test "atom_areas omitted when include_atom_areas is false" {
    var buf = std.ArrayList(u8){};
    defer buf.deinit(std.testing.allocator);

    const areas = [_]f64{ 1.0, 2.0 };

    var sw = StreamWriter.init(buf.writer(std.testing.allocator).any(), .ndjson, false);
    try sw.begin();
    try sw.writeResult(.{
        .filename = "no_areas.pdb",
        .n_atoms = 2,
        .total_sasa = 3.0,
        .status = .ok,
        .time_ns = 100_000,
        .atom_areas = &areas,
    });
    try sw.end();

    const output = buf.items;

    // atom_areas must NOT be present
    try std.testing.expect(std.mem.indexOf(u8, output, "atom_areas") == null);
}

test "error result includes error message" {
    var buf = std.ArrayList(u8){};
    defer buf.deinit(std.testing.allocator);

    var sw = StreamWriter.init(buf.writer(std.testing.allocator).any(), .ndjson, false);
    try sw.begin();
    try sw.writeResult(.{
        .filename = "bad.pdb",
        .n_atoms = 0,
        .total_sasa = 0.0,
        .status = .err,
        .time_ns = 0,
        .error_msg = "file not found",
    });
    try sw.end();

    const output = buf.items;

    try std.testing.expect(std.mem.indexOf(u8, output, "\"status\":\"error\"") != null);
    try std.testing.expect(std.mem.indexOf(u8, output, "\"error\":\"file not found\"") != null);
}

test "empty batch produces valid output" {
    // NDJSON: empty output
    {
        var buf = std.ArrayList(u8){};
        defer buf.deinit(std.testing.allocator);

        var sw = StreamWriter.init(buf.writer(std.testing.allocator).any(), .ndjson, false);
        try sw.begin();
        try sw.end();

        try std.testing.expectEqualStrings("", buf.items);
    }

    // JSON array: minimal valid array
    {
        var buf = std.ArrayList(u8){};
        defer buf.deinit(std.testing.allocator);

        var sw = StreamWriter.init(buf.writer(std.testing.allocator).any(), .json_array, false);
        try sw.begin();
        try sw.end();

        try std.testing.expectEqualStrings("[\n]\n", buf.items);
    }
}

test "writeJsonString escapes special characters" {
    var buf = std.ArrayList(u8){};
    defer buf.deinit(std.testing.allocator);

    const w = buf.writer(std.testing.allocator).any();

    // Test backslash and double quote
    try writeJsonString(w, "path\\to\"file");
    try std.testing.expectEqualStrings("\"path\\\\to\\\"file\"", buf.items);

    // Reset and test control characters
    buf.clearRetainingCapacity();
    try writeJsonString(w, "line1\nline2\r\ttab");
    try std.testing.expectEqualStrings("\"line1\\nline2\\r\\ttab\"", buf.items);

    // Reset and test low control characters (0x01, 0x1F)
    buf.clearRetainingCapacity();
    try writeJsonString(w, &[_]u8{ 'a', 0x01, 'b', 0x1F, 'c' });
    try std.testing.expectEqualStrings("\"a\\u0001b\\u001fc\"", buf.items);
}

test "fromFileResult converts correctly" {
    const fr = batch.FileResult{
        .filename = "test.pdb",
        .n_atoms = 100,
        .sasa_time_ns = 5_000_000,
        .total_sasa = 1234.5,
        .status = .ok,
        .error_msg = null,
    };

    const sr = StreamResult.fromFileResult(fr);

    try std.testing.expectEqualStrings("test.pdb", sr.filename);
    try std.testing.expectEqual(@as(usize, 100), sr.n_atoms);
    try std.testing.expectEqual(@as(u64, 5_000_000), sr.time_ns);
    try std.testing.expectEqual(@as(f64, 1234.5), sr.total_sasa);
    try std.testing.expectEqual(batch.FileResult.Status.ok, sr.status);
    try std.testing.expectEqual(@as(?[]const u8, null), sr.error_msg);
    try std.testing.expectEqual(@as(?[]const f64, null), sr.atom_areas);
}
