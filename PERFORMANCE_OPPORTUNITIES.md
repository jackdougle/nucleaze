# Nucleaze Performance Improvement Opportunities

Analysis of the current architecture with proposed optimizations that maintain full correctness and existing functionality.

---

## Current Architecture Summary

- **Bloom mode**: Cache-line blocked bloom filter (512-bit blocks), depth-16 prefetch ring buffer, double-hashing with 9-bit probe extraction
- **Threading**: Rayon work-stealing pool + dedicated reader thread + crossbeam backpressure channels
- **K-mer processing**: ~17 cycles/iteration in the hot loop, every k-mer in every read is queried against the bloom filter
- **I/O**: needletail parser, single-threaded decompression, 10k-read arena chunks

---

## 1. Compiler Profile Tuning (Estimated: 5-15% overall speedup)

**What**: Add release profile optimizations to `Cargo.toml`.

```toml
[profile.release]
lto = true
codegen-units = 1
target-cpu = "native"    # or set via RUSTFLAGS
```

**Why it helps**: Currently no LTO or codegen-units settings exist. With `codegen-units = 1`, LLVM sees the entire crate as one compilation unit, enabling cross-function inlining and better register allocation across the hot loop in `process_read_bloom`. LTO extends this across crate boundaries (e.g., inlining needletail's record parsing). `target-cpu=native` unlocks architecture-specific instruction scheduling and wider SIMD registers.

**Risk**: Zero. Purely a build-time change. Increases compile time but produces identical behavior. Can be gated behind a `--release` build flag.

---

## 2. Minimizers for Index Size Reduction (Estimated: 2-4x smaller index, lower FPR)

**What**: Instead of inserting *every* k-mer from the reference into the bloom filter, insert only minimizer k-mers (the lexicographically or hash-smallest k-mer in each window of w consecutive k-mers).

**How it works**:
- During indexing: slide a window of size w over the reference, select the minimizer (smallest hash) from each window, insert only those into the bloom filter
- During querying: same window + selection logic on the read, query only minimizer k-mers
- Both paths use the same minimizer scheme, so any read that truly comes from the reference will select the same minimizers as the reference did

**Benefits**:
- **Index size reduction**: ~w-fold fewer k-mers inserted. For w=10, ~10x fewer entries. This means either:
  - Same memory, much lower FPR (fewer bits set per block)
  - Same FPR, much less memory
- **Fewer bloom lookups per read**: ~(read_len - k - w + 2) / w queries instead of (read_len - k + 1). For 150bp reads with k=31 and w=10, that's ~13 queries instead of ~120. Each query is a cache-line fetch from DRAM, so this directly reduces memory bandwidth pressure.
- **No false negatives**: If a contiguous region of a read matches the reference, the minimizer in that region will be identical. The bloom filter still has no false negatives for inserted k-mers.

**Caveats**:
- Adds per-base overhead for minimizer window maintenance (monotone deque, ~3-5 extra ops per base). The current hot loop is ~17 cycles/k-mer. With minimizers, you do ~17 cycles of k-mer extraction + ~5 cycles of deque maintenance per base, but only ~1/w of the bloom lookups. Net throughput depends on whether the workload is compute-bound (k-mer extraction) or memory-bound (bloom lookups).
- **When it helps most**: Large references where the bloom filter is large and memory bandwidth is the bottleneck. For small references (the 30kb ribodepletion case), the entire bloom filter fits in L2 cache and lookups are ~4 cycles, so minimizers would add overhead with minimal benefit.
- **When it doesn't help**: If reads are already very short relative to k+w, fewer minimizers means fewer chances to hit the threshold.

**Recommendation**: Implement as an opt-in flag (`--minimizer-window` or `--window`). Default to off (current behavior) for backward compatibility. For large-reference use cases (decontamination against human genome, etc.), this could be a significant win.

---

## 3. Parallel Gzip Decompression (Estimated: 1.5-3x I/O throughput for .gz inputs)

**What**: Replace single-threaded flate2/miniz decompression with a parallel decompressor like `gzp` (Rust wrapper around libdeflate with thread pool) or pigz-style block-parallel decompression.

**Why it helps**: The current architecture has a single reader thread that parses and decompresses sequentially. For gzipped FASTQ (the common case in production), decompression can be the bottleneck -- a single core decompressing at ~500 MB/s while 8+ cores sit idle waiting for data on the crossbeam channel. The bounded channel capacity of 16 chunks (10k reads each) provides backpressure, but if the reader can't fill it fast enough, compute threads starve.

**Implementation**:
- Use `niffler` or `gzp` crate for transparent parallel decompression
- Or: pre-decompress into a pipe with pigz and read from stdin (no code change needed, just a usage pattern)
- Or: use `mgzip` block-parallel format if users can re-compress their inputs

**Risk**: Very low. The decompression output is byte-identical. Only the decompression *speed* changes.

---

## 4. SIMD-Accelerated Base Encoding (Estimated: 2-4x encoding throughput)

**What**: Replace the scalar table-lookup `encode()` function with SIMD vectorized encoding that processes 16 (SSE2/NEON) or 32 (AVX2) bases simultaneously.

**How**:
- Load 16 bytes of sequence into a SIMD register
- Use `_mm_shuffle_epi8` (SSSE3) or equivalent NEON intrinsic as a parallel lookup table
- Compare against 0xFF sentinel to create a valid-bases mask
- Extract valid 2-bit encodings and pack into k-mer integers

**Why it helps**: The current `encode()` is called once per base in the innermost loop. While the table lookup is fast (~4 cycles), it's scalar. SIMD encoding would process 16 bases per instruction, reducing encoding from ~4 cycles/base to ~0.25 cycles/base. This matters because encoding + k-mer shifting is the dominant cost when the bloom filter fits in cache (small references).

**Caveats**:
- Adds platform-specific code paths (SSE2, SSSE3, AVX2, NEON)
- The benefit is largest when the bloom filter is cache-resident (small refs) and the bottleneck is CPU, not memory
- Needs careful handling of ambiguous bases (N) within SIMD lanes

**Risk**: Low. Can be implemented behind `#[cfg(target_arch)]` gates with scalar fallback. Output is identical.

---

## 5. Adaptive Chunk Sizing (Estimated: 5-10% throughput improvement)

**What**: Instead of a fixed 10,000-read chunk size, dynamically tune chunk size based on read length and thread count.

**Why**:
- For short reads (50bp), 10,000 reads = 500KB of sequence data -- too small, scheduling overhead dominates
- For long reads (10,000bp), 10,000 reads = 100MB of sequence data -- too large, causes memory spikes and load imbalance at the end of the file
- Optimal chunk size balances: scheduling overhead, memory usage, and tail-effect load balancing

**Implementation**: Set chunk size to `max(1000, min(50000, target_bytes / avg_read_len))` where `target_bytes` is ~2-5MB (enough work to amortize Rayon scheduling, small enough for good load balancing).

**Risk**: Zero. Only affects work granularity, not results.

---

## 6. Early Exit Optimization for High-Hit Reads (Estimated: 0-20% depending on hit rate)

**What**: In `process_read_bloom`, once `hits >= threshold`, the function returns `true` immediately (line 767). This already exists. However, the prefetch ring buffer has up to 15 pending lookups that will never be used. The drain loops at lines 780-789 and 800-809 also have early exit. This is already well-optimized.

**What could be added**: For the *inverse* case (keeping non-hits, which is the ribodepletion scenario), consider a "likely miss" fast path: if the first N k-mers all miss, skip the rest of the read. This is speculative and could introduce false negatives, so it's **not recommended** as a default. Mentioned only for completeness.

---

## 7. Memory-Mapped Bloom Filter Loading (Estimated: faster startup for pre-built indices)

**What**: Instead of reading the `.nkb` file into a `Vec<AtomicU64>` via `read_to_end()` + parse, use `mmap` to map the bloom filter data directly from disk.

**Why it helps**: For large bloom filters (e.g., human genome at 1% FPR = ~470MB), the current approach allocates memory, reads the file, and copies data. `mmap` skips the copy -- pages are loaded on-demand as the bloom filter is accessed during read processing. First-touch overhead is amortized across the processing run.

**Caveats**:
- Requires the on-disk format to match the in-memory layout (currently it does -- little-endian u64 words)
- `AtomicU64` has the same layout as `u64` on all supported platforms, so transmuting mmap'd data works
- During indexing (concurrent writes), mmap with `MAP_SHARED` + `msync` could work but adds complexity. Recommend mmap only for the read-processing (query-only) path.

**Risk**: Low. The file format already stores raw u64 words in the correct order.

---

## 8. Bloom Filter: Wider Blocks (Experimental)

**What**: Use 1024-bit blocks (16 x u64, 2 cache lines) instead of 512-bit.

**Why it might help**: Wider blocks reduce the collision rate within a block, lowering FPR for the same total memory. The 2x correction factor in `bloom_params` compensates for the blocked design's higher FPR vs. standard bloom -- wider blocks would reduce that correction.

**Why it might not help**: Each query now touches 2 cache lines instead of 1. The prefetch ring buffer would need to prefetch 128 bytes instead of 64. On modern CPUs with 64-byte cache lines, this doubles the DRAM bandwidth per query. Only worthwhile if FPR reduction matters more than raw throughput.

**Recommendation**: Not recommended unless FPR precision is more important than throughput. The current 512-bit design is the sweet spot.

---

## Priority Ranking

| # | Optimization | Impact | Effort | Risk |
|---|-------------|--------|--------|------|
| 1 | Compiler profile tuning | 5-15% | Trivial | None |
| 2 | Minimizers (opt-in) | 2-4x index size, fewer lookups | Medium | None (opt-in) |
| 3 | Parallel gzip decompression | 1.5-3x I/O for .gz | Low | None |
| 4 | SIMD base encoding | 2-4x encoding speed | Medium | None (fallback) |
| 5 | Adaptive chunk sizing | 5-10% | Low | None |
| 6 | mmap bloom loading | Faster startup | Low | None |
| 7 | Wider bloom blocks | Lower FPR | Medium | Tradeoff |

**Recommended implementation order**: 1, 3, 5, 2, 6, 4 (prioritizing low-effort, high-impact, zero-risk items first).

---

## What NOT to Change

- **The prefetch ring buffer**: Already excellent. Depth-16 perfectly hides DRAM latency.
- **Cache-line blocking**: 512-bit blocks are the right granularity for modern CPUs.
- **The hash function**: Simple, fast, good entropy distribution. No need for more expensive hashes.
- **Atomic relaxed ordering**: Correct for the indexing use case (no ordering constraints needed between k-mer insertions).
- **Arena-based chunking**: The double-buffer swap pattern is efficient and avoids allocation pressure.
