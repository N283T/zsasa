const std = @import("std");
const json_parser = @import("json_parser.zig");
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
    const QUEUE_SIZE = 4; // Number of prefetch slots

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

/// I/O thread context
pub const IoContext = struct {
    queue: *PrefetchQueue,
    files: []const []const u8,
    input_dir: []const u8,
    allocator: Allocator,

    /// I/O thread function: reads files and pushes to queue
    pub fn run(self: *IoContext) void {
        for (self.files) |filename| {
            // Allocate arena for this file
            const arena = self.allocator.create(std.heap.ArenaAllocator) catch {
                continue; // Skip on allocation failure
            };
            arena.* = std.heap.ArenaAllocator.init(std.heap.page_allocator);

            // Build input path
            const input_path = std.fs.path.join(arena.allocator(), &.{ self.input_dir, filename }) catch {
                arena.deinit();
                self.allocator.destroy(arena);
                continue;
            };

            // Read and parse input
            const input = json_parser.readAtomInputFromFile(arena.allocator(), input_path) catch {
                arena.deinit();
                self.allocator.destroy(arena);
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
