const std = @import("std");
const Allocator = std.mem.Allocator;

/// A simple thread pool for parallel work distribution.
/// Designed for batch processing where all tasks are known upfront.
pub fn ThreadPool(comptime Context: type, comptime Result: type) type {
    return struct {
        const Self = @This();

        /// Work function type: fn(context: Context, chunk_start: usize, chunk_end: usize) Result
        pub const WorkFn = *const fn (Context, usize, usize) Result;

        allocator: Allocator,
        threads: []std.Thread,
        results: []Result,
        work_fn: WorkFn,
        context: Context,
        total_items: usize,
        chunk_size: usize,
        next_chunk: std.atomic.Value(usize),
        total_chunks: usize,

        /// Initialize thread pool with specified number of worker threads.
        /// Returns error.InvalidChunkSize if chunk_size is 0.
        /// Returns error.InvalidThreadCount if n_threads is 0.
        pub fn init(
            allocator: Allocator,
            n_threads: usize,
            work_fn: WorkFn,
            context: Context,
            total_items: usize,
            chunk_size: usize,
        ) !Self {
            if (chunk_size == 0) return error.InvalidChunkSize;
            if (n_threads == 0) return error.InvalidThreadCount;

            const cpu_count = try std.Thread.getCpuCount();
            const actual_threads = @min(n_threads, cpu_count);
            const total_chunks = (total_items + chunk_size - 1) / chunk_size;

            const threads = try allocator.alloc(std.Thread, actual_threads);
            errdefer allocator.free(threads);

            const results = try allocator.alloc(Result, total_chunks);
            errdefer allocator.free(results);

            return Self{
                .allocator = allocator,
                .threads = threads,
                .results = results,
                .work_fn = work_fn,
                .context = context,
                .total_items = total_items,
                .chunk_size = chunk_size,
                .next_chunk = std.atomic.Value(usize).init(0),
                .total_chunks = total_chunks,
            };
        }

        /// Start all worker threads and wait for completion.
        pub fn run(self: *Self) !void {
            // Spawn worker threads
            for (self.threads, 0..) |*thread, i| {
                thread.* = try std.Thread.spawn(.{}, workerLoop, .{ self, i });
            }

            // Wait for all threads to complete
            for (self.threads) |thread| {
                thread.join();
            }
        }

        /// Worker thread main loop - grab chunks and process them.
        fn workerLoop(self: *Self, thread_id: usize) void {
            _ = thread_id;

            while (true) {
                // .monotonic: each chunk_idx is unique; results are read after join()
                const chunk_idx = self.next_chunk.fetchAdd(1, .monotonic);

                if (chunk_idx >= self.total_chunks) {
                    break; // No more work
                }

                // Calculate chunk bounds
                const chunk_start = chunk_idx * self.chunk_size;
                const chunk_end = @min(chunk_start + self.chunk_size, self.total_items);

                // Execute work function
                const result = self.work_fn(self.context, chunk_start, chunk_end);
                self.results[chunk_idx] = result;
            }
        }

        /// Get all results after run() completes.
        pub fn getResults(self: *const Self) []const Result {
            return self.results[0..self.total_chunks];
        }

        /// Deinitialize and free resources.
        pub fn deinit(self: *Self) void {
            self.allocator.free(self.threads);
            self.allocator.free(self.results);
        }
    };
}

/// Simplified parallel for loop - executes work_fn for each chunk in parallel.
/// Returns the aggregated result using the provided reduce function.
pub fn parallelFor(
    comptime Context: type,
    comptime Result: type,
    allocator: Allocator,
    n_threads: usize,
    work_fn: *const fn (Context, usize, usize) Result,
    context: Context,
    total_items: usize,
    chunk_size: usize,
    reduce_fn: *const fn ([]const Result) Result,
) !Result {
    if (total_items == 0) {
        const empty: []const Result = &.{};
        return reduce_fn(empty);
    }

    // For single-threaded or small workloads, run directly
    if (n_threads <= 1 or total_items <= chunk_size) {
        const result = work_fn(context, 0, total_items);
        const single: []const Result = &.{result};
        return reduce_fn(single);
    }

    var pool = try ThreadPool(Context, Result).init(
        allocator,
        n_threads,
        work_fn,
        context,
        total_items,
        chunk_size,
    );
    defer pool.deinit();

    try pool.run();

    return reduce_fn(pool.getResults());
}

// Tests

test "ThreadPool - basic functionality" {
    const allocator = std.testing.allocator;

    const Context = struct {
        data: []const i32,
    };

    const work_fn = struct {
        fn call(ctx: Context, start: usize, end: usize) i64 {
            var sum: i64 = 0;
            for (ctx.data[start..end]) |val| {
                sum += val;
            }
            return sum;
        }
    }.call;

    const data = [_]i32{ 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 };
    const context = Context{ .data = &data };

    var pool = try ThreadPool(Context, i64).init(
        allocator,
        4,
        work_fn,
        context,
        data.len,
        3, // chunk size
    );
    defer pool.deinit();

    try pool.run();

    const results = pool.getResults();
    var total: i64 = 0;
    for (results) |r| {
        total += r;
    }

    // 1+2+3+4+5+6+7+8+9+10 = 55
    try std.testing.expectEqual(@as(i64, 55), total);
}

test "ThreadPool - single item" {
    const allocator = std.testing.allocator;

    const Context = struct {
        value: i32,
    };

    const work_fn = struct {
        fn call(ctx: Context, start: usize, end: usize) i64 {
            _ = start;
            _ = end;
            return ctx.value;
        }
    }.call;

    const context = Context{ .value = 42 };

    var pool = try ThreadPool(Context, i64).init(
        allocator,
        4,
        work_fn,
        context,
        1,
        1,
    );
    defer pool.deinit();

    try pool.run();

    const results = pool.getResults();
    try std.testing.expectEqual(@as(usize, 1), results.len);
    try std.testing.expectEqual(@as(i64, 42), results[0]);
}

test "parallelFor - sum reduction" {
    const allocator = std.testing.allocator;

    const Context = struct {
        data: []const i32,
    };

    const work_fn = struct {
        fn call(ctx: Context, start: usize, end: usize) i64 {
            var sum: i64 = 0;
            for (ctx.data[start..end]) |val| {
                sum += val;
            }
            return sum;
        }
    }.call;

    const reduce_fn = struct {
        fn call(results: []const i64) i64 {
            var total: i64 = 0;
            for (results) |r| {
                total += r;
            }
            return total;
        }
    }.call;

    const data = [_]i32{ 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 };
    const context = Context{ .data = &data };

    const result = try parallelFor(
        Context,
        i64,
        allocator,
        4,
        work_fn,
        context,
        data.len,
        3,
        reduce_fn,
    );

    try std.testing.expectEqual(@as(i64, 55), result);
}

test "parallelFor - empty input" {
    const allocator = std.testing.allocator;

    const Context = struct {};

    const work_fn = struct {
        fn call(_: Context, _: usize, _: usize) i64 {
            return 0;
        }
    }.call;

    const reduce_fn = struct {
        fn call(results: []const i64) i64 {
            var total: i64 = 0;
            for (results) |r| {
                total += r;
            }
            return total;
        }
    }.call;

    const result = try parallelFor(
        Context,
        i64,
        allocator,
        4,
        work_fn,
        Context{},
        0,
        10,
        reduce_fn,
    );

    try std.testing.expectEqual(@as(i64, 0), result);
}

test "parallelFor - single thread fallback" {
    const allocator = std.testing.allocator;

    const Context = struct {
        data: []const i32,
    };

    const work_fn = struct {
        fn call(ctx: Context, start: usize, end: usize) i64 {
            var sum: i64 = 0;
            for (ctx.data[start..end]) |val| {
                sum += val;
            }
            return sum;
        }
    }.call;

    const reduce_fn = struct {
        fn call(results: []const i64) i64 {
            var total: i64 = 0;
            for (results) |r| {
                total += r;
            }
            return total;
        }
    }.call;

    const data = [_]i32{ 1, 2, 3, 4, 5 };
    const context = Context{ .data = &data };

    // Single thread
    const result = try parallelFor(
        Context,
        i64,
        allocator,
        1,
        work_fn,
        context,
        data.len,
        10,
        reduce_fn,
    );

    try std.testing.expectEqual(@as(i64, 15), result);
}

test "ThreadPool - zero chunk size returns error" {
    const allocator = std.testing.allocator;

    const Context = struct {};
    const work_fn = struct {
        fn call(_: Context, _: usize, _: usize) i64 {
            return 0;
        }
    }.call;

    const result = ThreadPool(Context, i64).init(
        allocator,
        4,
        work_fn,
        Context{},
        10,
        0, // Invalid: zero chunk size
    );

    try std.testing.expectError(error.InvalidChunkSize, result);
}

test "ThreadPool - zero threads returns error" {
    const allocator = std.testing.allocator;

    const Context = struct {};
    const work_fn = struct {
        fn call(_: Context, _: usize, _: usize) i64 {
            return 0;
        }
    }.call;

    const result = ThreadPool(Context, i64).init(
        allocator,
        0, // Invalid: zero threads
        work_fn,
        Context{},
        10,
        5,
    );

    try std.testing.expectError(error.InvalidThreadCount, result);
}
