# MiniBloom Changelog

## Multi-hash bloom filter with double hashing

**Problem**: Single hash function gave ~10% per-k-mer FPR regardless of `--fpr` setting, because the sizing formula assumes optimal k hash functions but only 1 was used. With ~130 k-mers per 150bp read and `min_hits=1`, per-read FPR reached ~100%.

**Fix**: Added double hashing (Kirsch-Mitzenmacher scheme) to derive k independent bit positions from two base hashes:

- `h1`: murmur3 64-bit finalizer
- `h2`: Stafford Mix variant (golden-ratio seeded)
- Position i: `h1 + i * h2`

Optimal k is computed automatically: `k = (m/n) * ln(2)`. For 100M k-mers at 1% FPR this gives 7 hash functions.

`contains()` short-circuits on the first unset bit, so non-matching k-mers (the common case) typically bail after 1-2 hash checks.

## Power-of-2 filter sizing

**Problem**: `num_bits` was not a power of 2, so `kmer & mask` did not distribute uniformly across all bit positions.

**Fix**: `num_bits` is now rounded up to the next power of 2 via `next_power_of_two()`. Wastes a small amount of memory but guarantees uniform distribution with a simple bitmask.

## Parallel `absorb_shards` and `from_serialized`

**Problem**: After the parallel reference indexing phase, `absorb_shards` inserted 100M k-mers sequentially. With 7 hash functions, this meant 700M bit-set operations on a single thread, taking ~3s.

**Fix**: Both `absorb_shards` and `from_serialized` now temporarily reinterpret `Vec<u64>` as `&[AtomicU64]` and use `fetch_or` with `Ordering::Relaxed` to insert across rayon threads in parallel. Atomics are only used during bulk-load; normal `insert`/`contains` remain non-atomic.

## `contains` correctness fix (prior session)

**Problem**: `contains()` used `fetch_and` (a mutating atomic op) instead of a read-only check. This destroyed bits in the filter during queries, causing subsequent lookups to return false negatives.

**Fix**: Replaced with `load()` + bitwise AND check (read-only).

## Buffered `stage` (prior session)

**Problem**: `MiniBloom::stage` inserted directly into the bloom filter via per-k-mer atomics during the parallel reference indexing phase. This caused heavy cross-thread cache-line contention (~200M atomic ops) and double-inserted every k-mer since `absorb_shards` ran afterward.

**Fix**: `stage` now buffers k-mers into the subindex (same as `HashShards`). Bloom insertion happens once during `absorb_shards`. `needs_absorb()` returns `true`.

---

## Blocked bloom filter FPR correction

**Problem**: `bloom_params` used the standard (non-blocked) bloom filter sizing
formula `m = -(n * ln(fpr)) / (ln(2))^2`. In a cache-line blocked filter where
all probes for one item target the same 512-bit block, items are Poisson-distributed
across blocks. Overloaded blocks dominate the aggregate FPR (Jensen's inequality),
causing 3-4x higher FPR than the formula predicts.

**Fix**: Applied a 2x correction factor to the bit count before computing block
count and hash count. This is the standard approach for cache-line blocked bloom
filters (Putze, Sanders, Singler 2007). Memory usage doubles but remains small
relative to the input data.

---

## Note on --fpr semantics

The `--fpr` flag controls the **per-read** false positive rate target. Internally,
this is converted to a per-k-mer FPR using:

```
per_kmer_fpr = 1 - (1 - per_read_fpr)^(1 / num_kmers_per_read)
```

where `num_kmers_per_read` is estimated as `150 - k + 1` for a typical 150bp read.

For low hit-rate scenarios with `--minhits 1`, even a small per-read FPR can cause
poor precision (many false positives among reported hits). Mitigations:
- Increase `--minhits` (e.g., 2-3) so a single FP k-mer doesn't trigger a match
- Lower `--fpr` (increases memory)
