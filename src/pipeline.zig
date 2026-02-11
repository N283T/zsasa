const std = @import("std");
const format_detect = @import("format_detect.zig");
const json_parser = @import("json_parser.zig");
const mmcif_parser = @import("mmcif_parser.zig");
const pdb_parser = @import("pdb_parser.zig");
const types = @import("types.zig");

const Allocator = std.mem.Allocator;
const AtomInput = types.AtomInput;

/// A prefetched file ready for processing
pub const PrefetchedFile = struct {
    filename: []const u8,
    input: AtomInput,
    arena: *std.heap.ArenaAllocator,

    pub fn deinit(self: *PrefetchedFile, allocator: Allocator) void {
        self.input.deinit();
        self.arena.deinit();
        allocator.destroy(self.arena);
    }
};

/// Bounded queue for prefetched files
pub const PrefetchQueue = struct {
    /// Number of prefetch slots. Sized to keep I/O thread busy while processing,
    /// but small enough to limit memory usage (each slot holds a fully parsed file).
    /// 4 slots provides good overlap between I/O and computation without excessive memory.
    const QUEUE_SIZE = 4;

    slots: [QUEUE_SIZE]?PrefetchedFile,
    head: usize, // Consumer reads from here
    tail: usize, // Producer writes here
    count: usize,
    done: bool, // I/O thread finished
    mutex: std.Thread.Mutex,
    not_empty: std.Thread.Condition,
    not_full: std.Thread.Condition,
    allocator: Allocator,

    pub fn init(allocator: Allocator) PrefetchQueue {
        return PrefetchQueue{
            .slots = [_]?PrefetchedFile{null} ** QUEUE_SIZE,
            .head = 0,
            .tail = 0,
            .count = 0,
            .done = false,
            .mutex = .{},
            .not_empty = .{},
            .not_full = .{},
            .allocator = allocator,
        };
    }

    pub fn deinit(self: *PrefetchQueue) void {
        // Clean up any remaining items
        self.mutex.lock();
        defer self.mutex.unlock();

        while (self.count > 0) {
            if (self.slots[self.head]) |*item| {
                item.deinit(self.allocator);
            }
            self.slots[self.head] = null;
            self.head = (self.head + 1) % QUEUE_SIZE;
            self.count -= 1;
        }
    }

    /// Producer: Add a prefetched file to the queue (blocks if full)
    pub fn push(self: *PrefetchQueue, file: PrefetchedFile) void {
        self.mutex.lock();
        defer self.mutex.unlock();

        // Wait until there's space
        while (self.count >= QUEUE_SIZE) {
            self.not_full.wait(&self.mutex);
        }

        self.slots[self.tail] = file;
        self.tail = (self.tail + 1) % QUEUE_SIZE;
        self.count += 1;

        // Signal consumer
        self.not_empty.signal();
    }

    /// Consumer: Get next file (blocks if empty, returns null if done)
    pub fn pop(self: *PrefetchQueue) ?PrefetchedFile {
        self.mutex.lock();
        defer self.mutex.unlock();

        // Wait until there's an item or we're done
        while (self.count == 0 and !self.done) {
            self.not_empty.wait(&self.mutex);
        }

        if (self.count == 0) {
            return null; // Queue empty and done
        }

        const item = self.slots[self.head];
        self.slots[self.head] = null;
        self.head = (self.head + 1) % QUEUE_SIZE;
        self.count -= 1;

        // Signal producer
        self.not_full.signal();

        return item;
    }

    /// Mark I/O as complete
    pub fn markDone(self: *PrefetchQueue) void {
        self.mutex.lock();
        defer self.mutex.unlock();

        self.done = true;
        self.not_empty.broadcast(); // Wake up any waiting consumers
    }
};

/// Record of a failed file during I/O
pub const FailedFile = struct {
    filename: []const u8,
    reason: FailureReason,

    pub const FailureReason = enum {
        allocation_failed,
        path_join_failed,
        read_parse_failed,
    };
};

/// I/O thread context
pub const IoContext = struct {
    queue: *PrefetchQueue,
    files: []const []const u8,
    input_dir: []const u8,
    allocator: Allocator,
    // Thread-safe tracking of failed files
    failed_files: std.ArrayListUnmanaged(FailedFile),
    failed_mutex: std.Thread.Mutex,

    /// Filter options for parsers
    skip_hydrogens: bool,
    atom_only: bool,

    pub fn init(
        allocator: Allocator,
        queue: *PrefetchQueue,
        files: []const []const u8,
        input_dir: []const u8,
    ) IoContext {
        return IoContext{
            .queue = queue,
            .files = files,
            .input_dir = input_dir,
            .allocator = allocator,
            .failed_files = .{},
            .failed_mutex = .{},
            .skip_hydrogens = true,
            .atom_only = true,
        };
    }

    pub fn deinit(self: *IoContext) void {
        self.failed_files.deinit(self.allocator);
    }

    /// Record a failed file (thread-safe)
    fn recordFailure(self: *IoContext, filename: []const u8, reason: FailedFile.FailureReason) void {
        self.failed_mutex.lock();
        defer self.failed_mutex.unlock();
        self.failed_files.append(self.allocator, .{ .filename = filename, .reason = reason }) catch {};
    }

    /// Get the number of failed files
    pub fn failedCount(self: *IoContext) usize {
        self.failed_mutex.lock();
        defer self.failed_mutex.unlock();
        return self.failed_files.items.len;
    }

    /// I/O thread function: reads files and pushes to queue
    pub fn run(self: *IoContext) void {
        for (self.files) |filename| {
            // Allocate arena for this file
            const arena = self.allocator.create(std.heap.ArenaAllocator) catch {
                self.recordFailure(filename, .allocation_failed);
                continue;
            };
            arena.* = std.heap.ArenaAllocator.init(std.heap.page_allocator);

            // Build input path
            const input_path = std.fs.path.join(arena.allocator(), &.{ self.input_dir, filename }) catch {
                arena.deinit();
                self.allocator.destroy(arena);
                self.recordFailure(filename, .path_join_failed);
                continue;
            };

            // Read and parse input (auto-detect format)
            const format = format_detect.detectInputFormat(input_path);
            const input = switch (format) {
                .json => json_parser.readAtomInputFromFile(arena.allocator(), input_path),
                .mmcif => blk: {
                    var parser = mmcif_parser.MmcifParser.init(arena.allocator());
                    parser.skip_hydrogens = self.skip_hydrogens;
                    parser.atom_only = self.atom_only;
                    break :blk parser.parseFile(input_path);
                },
                .pdb => blk: {
                    var parser = pdb_parser.PdbParser.init(arena.allocator());
                    parser.skip_hydrogens = self.skip_hydrogens;
                    parser.atom_only = self.atom_only;
                    break :blk parser.parseFile(input_path);
                },
            } catch {
                arena.deinit();
                self.allocator.destroy(arena);
                self.recordFailure(filename, .read_parse_failed);
                continue;
            };

            // Push to queue
            self.queue.push(.{
                .filename = filename,
                .input = input,
                .arena = arena,
            });
        }

        // Signal completion
        self.queue.markDone();
    }
};

test "PrefetchQueue basic operations" {
    const allocator = std.testing.allocator;
    var queue = PrefetchQueue.init(allocator);
    defer queue.deinit();

    // Test empty queue after marking done
    queue.markDone();
    const result = queue.pop();
    try std.testing.expectEqual(@as(?PrefetchedFile, null), result);
}

test "PrefetchQueue count tracking" {
    const allocator = std.testing.allocator;
    var queue = PrefetchQueue.init(allocator);
    defer queue.deinit();

    // Initial state
    try std.testing.expectEqual(@as(usize, 0), queue.count);
    try std.testing.expectEqual(false, queue.done);

    // Mark done
    queue.markDone();
    try std.testing.expectEqual(true, queue.done);
}

test "IoContext failure tracking" {
    const allocator = std.testing.allocator;
    var queue = PrefetchQueue.init(allocator);
    defer queue.deinit();

    const files = [_][]const u8{};
    var ctx = IoContext.init(allocator, &queue, &files, "/nonexistent");
    defer ctx.deinit();

    // Initially no failures
    try std.testing.expectEqual(@as(usize, 0), ctx.failedCount());

    // Record some failures
    ctx.recordFailure("file1.json", .read_parse_failed);
    try std.testing.expectEqual(@as(usize, 1), ctx.failedCount());

    ctx.recordFailure("file2.json", .allocation_failed);
    try std.testing.expectEqual(@as(usize, 2), ctx.failedCount());

    // Verify failure details
    try std.testing.expectEqual(@as(usize, 2), ctx.failed_files.items.len);
    try std.testing.expectEqualStrings("file1.json", ctx.failed_files.items[0].filename);
    try std.testing.expectEqual(FailedFile.FailureReason.read_parse_failed, ctx.failed_files.items[0].reason);
}
