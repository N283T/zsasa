const std = @import("std");
const builtin = @import("builtin");
const Allocator = std.mem.Allocator;

/// Fault injection and bookkeeping for tests. Outside test builds this is an
/// empty namespace and every use of it below is compiled out.
pub const testing = if (builtin.is_test) struct {
    /// Number of worker spawns that still succeed before the next one reports
    /// `error.ThreadQuotaExceeded`. `null` disables the injection.
    pub var spawns_until_failure: ?usize = null;
    /// Worker threads that were started and have not returned yet.
    pub var live_workers: std.atomic.Value(usize) = .init(0);
} else struct {};

/// Number of chunks `total_items` is split into: the number of `work_fn` calls
/// of a pool run, each starting at a multiple of `chunk_size`.
pub fn chunkCount(total_items: usize, chunk_size: usize) usize {
    return (total_items + chunk_size - 1) / chunk_size;
}

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
            const total_chunks = chunkCount(total_items, chunk_size);

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

        /// Start the worker threads and wait until every chunk is processed.
        ///
        /// A worker that cannot be spawned (thread quota, `pids.max`, memory) is
        /// not an error: the chunk queue is shared, so the workers that did start
        /// drain it, and the calling thread takes the place of the missing ones.
        /// Every started worker is joined before this returns, so nothing still
        /// references the pool, the context or the caller's buffers afterwards.
        pub fn run(self: *Self) void {
            var spawned: usize = 0;
            for (self.threads) |*thread| {
                thread.* = self.spawnWorker() catch break;
                spawned += 1;
            }

            if (spawned < self.threads.len) self.workerLoop();

            for (self.threads[0..spawned]) |thread| {
                thread.join();
            }
        }

        fn spawnWorker(self: *Self) std.Thread.SpawnError!std.Thread {
            if (builtin.is_test) {
                if (testing.spawns_until_failure) |*remaining| {
                    if (remaining.* == 0) return error.ThreadQuotaExceeded;
                    remaining.* -= 1;
                }
                _ = testing.live_workers.fetchAdd(1, .monotonic);
            }
            errdefer if (builtin.is_test) {
                _ = testing.live_workers.fetchSub(1, .monotonic);
            };
            return std.Thread.spawn(.{}, workerMain, .{self});
        }

        /// Entry point of a spawned worker thread.
        fn workerMain(self: *Self) void {
            self.workerLoop();
            if (builtin.is_test) _ = testing.live_workers.fetchSub(1, .monotonic);
        }

        /// Grab chunks and process them until none are left. Runs on every
        /// worker thread, and on the calling thread when a spawn failed.
        fn workerLoop(self: *Self) void {
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
///
/// The only errors are those of `ThreadPool.init`, which is reached before any
/// thread starts. Once `work_fn` has been called, every chunk is processed and
/// no worker is left running, even if some threads could not be spawned, so the
/// caller may free what `context` points to as soon as this returns.
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

    pool.run();

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

    pool.run();

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

    pool.run();

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

// Spawn-failure tests. `testing.spawns_until_failure` makes the pool behave as
// if the process hit its thread limit after that many workers were started.

/// Keeps a chunk busy for a moment, so that the workers which did start are
/// still running when the failing spawn is reached.
fn testStall() void {
    for (0..50) |_| std.Thread.yield() catch {};
}

test "ThreadPool - spawn failure after K workers still processes every chunk" {
    const allocator = std.testing.allocator;

    const n_items = 403;
    const chunk_size = 10;
    const n_threads = 4;

    const Context = struct {
        visits: []std.atomic.Value(u32),
    };

    const work_fn = struct {
        fn call(ctx: Context, start: usize, end: usize) usize {
            testStall();
            var sum: usize = 0;
            for (start..end) |i| {
                _ = ctx.visits[i].fetchAdd(1, .monotonic);
                sum += i;
            }
            return sum;
        }
    }.call;

    // K == n_threads is the run without a failure.
    for (0..n_threads + 1) |k| {
        var visits: [n_items]std.atomic.Value(u32) = @splat(.init(0));

        var pool = try ThreadPool(Context, usize).init(
            allocator,
            n_threads,
            work_fn,
            Context{ .visits = &visits },
            n_items,
            chunk_size,
        );
        defer pool.deinit();

        testing.spawns_until_failure = k;
        defer testing.spawns_until_failure = null;

        pool.run();

        // Nothing may still be running once run() has returned: the pool and
        // everything the context points to are about to be released.
        try std.testing.expectEqual(@as(usize, 0), testing.live_workers.load(.monotonic));

        // Every item was processed exactly once, whoever processed it.
        for (&visits) |*v| {
            try std.testing.expectEqual(@as(u32, 1), v.load(.monotonic));
        }

        // Every chunk has its own result.
        const results = pool.getResults();
        try std.testing.expectEqual(@as(usize, 41), results.len);
        for (results, 0..) |result, chunk_idx| {
            const start = chunk_idx * chunk_size;
            const end = @min(start + chunk_size, n_items);
            var expected: usize = 0;
            for (start..end) |i| expected += i;
            try std.testing.expectEqual(expected, result);
        }
    }
}

test "ThreadPool - runs on the calling thread when no worker can be spawned" {
    const allocator = std.testing.allocator;

    const Context = struct {
        caller: std.Thread.Id,
        chunks_on_caller: *std.atomic.Value(usize),
    };

    const work_fn = struct {
        fn call(ctx: Context, _: usize, _: usize) void {
            if (std.Thread.getCurrentId() == ctx.caller) {
                _ = ctx.chunks_on_caller.fetchAdd(1, .monotonic);
            }
        }
    }.call;

    var chunks_on_caller = std.atomic.Value(usize).init(0);

    var pool = try ThreadPool(Context, void).init(
        allocator,
        4,
        work_fn,
        Context{ .caller = std.Thread.getCurrentId(), .chunks_on_caller = &chunks_on_caller },
        100,
        10,
    );
    defer pool.deinit();

    testing.spawns_until_failure = 0;
    defer testing.spawns_until_failure = null;

    pool.run();

    try std.testing.expectEqual(@as(usize, 0), testing.live_workers.load(.monotonic));
    try std.testing.expectEqual(@as(usize, 10), chunks_on_caller.load(.monotonic));
}

test "parallelFor - spawn failure after K workers returns the complete result" {
    const allocator = std.testing.allocator;

    const Context = struct {
        data: []const i32,
    };

    const work_fn = struct {
        fn call(ctx: Context, start: usize, end: usize) i64 {
            testStall();
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

    const n_threads = 4;

    for (0..n_threads + 1) |k| {
        // Heap-allocated so that a worker outliving the call would touch freed
        // memory, as it did with the buffers of the SASA callers.
        const data = try allocator.alloc(i32, 500);
        defer allocator.free(data);
        for (data, 0..) |*d, i| d.* = @intCast(i);

        testing.spawns_until_failure = k;
        defer testing.spawns_until_failure = null;

        const result = try parallelFor(
            Context,
            i64,
            allocator,
            n_threads,
            work_fn,
            Context{ .data = data },
            data.len,
            7,
            reduce_fn,
        );

        try std.testing.expectEqual(@as(usize, 0), testing.live_workers.load(.monotonic));
        // 0 + 1 + ... + 499
        try std.testing.expectEqual(@as(i64, 124750), result);
    }
}
