# Phase 4: Multi-threading Optimization

## Goal

Parallelize SASA calculations across multiple CPU cores using Zig's `std.Thread`.

**Current**: ~11ms (SIMD optimized, single-threaded)
**Target**: ~3-5ms (2-4x improvement with multi-threading)

## Sub-Phases

### Phase 4.1: Thread Pool Infrastructure

**Goal**: Create reusable thread pool for parallel work distribution.

**Tasks**:
- [x] Create `src/thread_pool.zig` with basic thread pool implementation
- [x] Implement work queue with thread-safe task distribution
- [x] Add graceful shutdown and error handling
- [x] Unit tests for thread pool

**Files**:
| File | Action |
|------|--------|
| `src/thread_pool.zig` | CREATE |

---
- [x] **DONE** - Sub-phase 4.1 complete

---

### Phase 4.2: Atom-level Parallelization

**Goal**: Parallelize per-atom SASA calculations.

**Tasks**:
- [ ] Add `calculateSasaParallel()` function in `shrake_rupley.zig`
- [ ] Distribute atoms across worker threads
- [ ] Aggregate results from all threads
- [ ] Ensure thread-safe accumulation of total_area
- [ ] Integration tests comparing parallel vs sequential results

**Key Insight**: Each atom's SASA is independent - embarrassingly parallel.

```zig
// Pseudo-code
pub fn calculateSasaParallel(allocator, input, config, n_threads) !SasaResult {
    // 1. Pre-compute shared data (positions, radii_sq, neighbor_list)
    // 2. Divide atoms into chunks
    // 3. Spawn worker threads, each processes a chunk
    // 4. Collect results and sum total_area
}
```

**Files**:
| File | Action |
|------|--------|
| `src/shrake_rupley.zig` | MODIFY |

---
- [ ] **DONE** - Sub-phase 4.2 complete

---

### Phase 4.3: Chunk-based Processing

**Goal**: Optimize cache efficiency with proper chunk sizing.

**Tasks**:
- [ ] Implement chunk-based work distribution
- [ ] Experiment with chunk sizes (64, 128, 256, 512 atoms)
- [ ] Consider cache line alignment for shared data
- [ ] Minimize false sharing between threads

**Considerations**:
- Too small chunks → thread overhead dominates
- Too large chunks → load imbalance
- Sweet spot: ~100-500 atoms per chunk (depends on neighbor count)

**Files**:
| File | Action |
|------|--------|
| `src/shrake_rupley.zig` | MODIFY |

---
- [ ] **DONE** - Sub-phase 4.3 complete

---

### Phase 4.4: Benchmark & Tuning

**Goal**: Find optimal configuration and validate performance.

**Tasks**:
- [ ] Benchmark with varying thread counts (1, 2, 4, 8, N_CPU)
- [ ] Benchmark with varying chunk sizes
- [ ] Compare against single-threaded baseline
- [ ] Document optimal settings for different protein sizes
- [ ] Update README with multi-threading usage

**Expected Results Table**:
| Threads | Time (ms) | Speedup |
|---------|-----------|---------|
| 1 | ~11 | 1.0x |
| 2 | ~6 | ~1.8x |
| 4 | ~3-4 | ~3x |
| 8 | ~2-3 | ~4x |

**Files**:
| File | Action |
|------|--------|
| `scripts/benchmark.py` | MODIFY (add thread count option) |
| `README.md` | MODIFY |

---
- [ ] **DONE** - Sub-phase 4.4 complete

---

## Architecture

```
┌─────────────────────────────────────────────────────────┐
│                    Main Thread                          │
│  - Build neighbor list (once)                           │
│  - Pre-compute radii_with_probe_sq                      │
│  - Allocate result arrays                               │
└─────────────────────┬───────────────────────────────────┘
                      │
         ┌────────────┼────────────┐
         ▼            ▼            ▼
    ┌─────────┐  ┌─────────┐  ┌─────────┐
    │ Worker 1│  │ Worker 2│  │ Worker 3│  ...
    │ Atoms   │  │ Atoms   │  │ Atoms   │
    │ 0-999   │  │ 1000-   │  │ 2000-   │
    │         │  │ 1999    │  │ 2999    │
    └────┬────┘  └────┬────┘  └────┬────┘
         │            │            │
         └────────────┼────────────┘
                      ▼
              ┌───────────────┐
              │ Aggregate     │
              │ Results       │
              └───────────────┘
```

## Thread Safety Analysis

| Data | Access Pattern | Thread Safety |
|------|----------------|---------------|
| `positions` | Read-only | Safe (shared) |
| `radii_with_probe_sq` | Read-only | Safe (shared) |
| `neighbor_list` | Read-only | Safe (shared) |
| `test_points` | Read-only | Safe (shared) |
| `atom_areas[i]` | Write (unique i) | Safe (no overlap) |
| `total_area` | Write (accumulate) | Need atomic or per-thread sum |

## Success Criteria

- [ ] All existing tests pass
- [ ] SASA values identical to single-threaded (< 0.0001% difference)
- [ ] Performance improvement with multi-threading
- [ ] No data races (verify with ThreadSanitizer if available)
- [ ] Graceful fallback to single-threaded if n_threads=1

## Verification Commands

```bash
# Build and test
zig build test

# Benchmark single-threaded
time ./zig-out/bin/freesasa_zig examples/input_1a0q.json /tmp/out.json

# Benchmark multi-threaded (after implementation)
time ./zig-out/bin/freesasa_zig examples/input_1a0q.json /tmp/out.json --threads 4

# Full benchmark
uv run scripts/benchmark.py examples/1A0Q.cif.gz --runs 5
```

## Notes

- Zig's `std.Thread.Pool` might be usable, or we can implement a simple one
- Consider using `@atomicRmw` for thread-safe accumulation
- Start with simple implementation, optimize later
- Amdahl's Law: speedup limited by sequential portions (neighbor list build)

---
- [ ] **DONE** - Phase 4 complete
