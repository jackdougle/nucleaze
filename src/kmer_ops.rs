//! Operations for extracting and measuring k-mers from FASTX records
use needletail::bitkmer::canonical;
use rayon::prelude::*;
use rustc_hash::FxHashSet;
use std::fs::File;
use std::io::{BufReader, BufWriter, Read as IoRead, Write as IoWrite};
use std::sync::atomic::{AtomicU64, Ordering as AtomicOrdering};

pub trait KmerStore: Send + Sync {
    fn new(size: usize, fpr: f64) -> Self;
    fn map(&self, kmer: &u64) -> usize;
    fn stage(&self, kmer: u64, subindex: &mut Vec<Vec<u64>>);
    fn insert(&mut self, kmer: &u64);
    fn contains(&self, kmer: &u64) -> bool;
    fn count(&self) -> usize;
    fn clear(&mut self);
    fn from_serialized(&mut self, data: Vec<Vec<u64>>);
    /// Bulk-load pre-sharded k-mer sets from the reference indexing accumulator.
    fn absorb_shards(&mut self, shards: Vec<FxHashSet<u64>>);
    /// Whether `stage` buffers into the subindex (true) or inserts directly (false).
    /// When false, `absorb_shards` is skipped after reference indexing.
    fn needs_absorb(&self) -> bool;
}

pub struct HashShards {
    index: Vec<FxHashSet<u64>>,
    mask: usize,
}

impl KmerStore for HashShards {
    fn new(_: usize, _: f64) -> Self {
        let num_idx = 1024;
        HashShards {
            index: vec![FxHashSet::default(); num_idx],
            mask: num_idx - 1,
        }
    }

    /// Insert k-mer into respective shard.
    #[inline(always)]
    fn insert(&mut self, kmer: &u64) {
        let idx = self.map(kmer);
        self.index[idx].insert(*kmer);
    }

    /// Hash k-mer and return shard index.
    #[inline(always)]
    fn map(&self, kmer: &u64) -> usize {
        let hashed = kmer ^ (kmer >> 12);
        hashed as usize & self.mask
    }

    /// Add a k-mer to a temporary & external subindex.
    #[inline(always)]
    fn stage(&self, kmer: u64, subindex: &mut Vec<Vec<u64>>) {
        let idx = self.map(&kmer);
        subindex[idx].push(kmer);
    }

    /// Check shard for presence for k-mer.
    #[inline(always)]
    fn contains(&self, kmer: &u64) -> bool {
        let idx = self.map(kmer);
        unsafe { self.index.get_unchecked(idx).contains(kmer) }
    }

    /// Return number of k-mers.
    fn count(&self) -> usize {
        self.index.iter().map(|i| i.len()).sum()
    }

    /// Remove k-mers from index.
    fn clear(&mut self) {
        self.index = vec![FxHashSet::default(); self.mask + 1];
    }

    /// Absorb pre-sharded sets directly into the index.
    fn absorb_shards(&mut self, shards: Vec<FxHashSet<u64>>) {
        self.index = shards;
    }

    fn needs_absorb(&self) -> bool {
        true
    }

    /// Add serialized k-mers to reference index.
    fn from_serialized(&mut self, data: Vec<Vec<u64>>) {
        self.index
            .par_iter_mut()
            .zip(data.into_par_iter())
            .for_each(|(set, vec)| {
                // Reserve space to avoid reallocations
                set.reserve(vec.len());
                for kmer in vec {
                    set.insert(kmer);
                }
            });
    }
}

impl HashShards {
    /// Tranpose reference index to a serializable format.
    pub fn serialize(&self) -> Vec<Vec<u64>> {
        self.index
            .iter()
            .map(|set| set.iter().cloned().collect())
            .collect()
    }
}

/// Compute optimal blocked bloom filter parameters for `n` items at target `fpr`.
/// Returns (num_blocks, n_hashes). Each block is 512 bits (1 cache line).
pub fn bloom_params(n: usize, fpr: f64) -> (usize, u8) {
    let n = n.max(64) as f64;
    let m_standard = (-(n * fpr.ln()) / (2.0f64.ln().powi(2))).ceil() as u64;
    // 2x correction for cache-line blocked bloom filter: all probes for one item
    // target the same 512-bit block, so block occupancy follows a Poisson distribution.
    // By Jensen's inequality the average FPR exceeds the FPR at mean occupancy;
    // doubling the bit count compensates (Putze, Sanders, Singler 2007).
    let m_corrected = m_standard * 2;
    let k_h = ((m_corrected as f64 / n) * 2.0f64.ln()).round() as u8;
    let k_h = k_h.max(1).min(8);
    let num_blocks_raw = ((m_corrected + 511) / 512).max(1);
    let num_blocks = num_blocks_raw.next_power_of_two() as usize;
    (num_blocks, k_h)
}

/// Compute per-k-mer FPR from per-read FPR, assuming ~150bp reads.
pub fn per_kmer_fpr(read_fpr: f64, k: usize) -> f64 {
    let queries_per_read = 150usize.saturating_sub(k).max(1) + 1;
    1.0 - (1.0 - read_fpr).powf(1.0 / queries_per_read as f64)
}

/// Two independent hash values for bloom filter probing.
/// Uses multiply + add + xor-shift to ensure full-width entropy
/// (needed for bit-extraction probing where every 9-bit window matters).
#[inline(always)]
fn bloom_hash(item: u64) -> (u64, u64) {
    let mut h1 = item
        .wrapping_mul(0x9E3779B97F4A7C15)
        .wrapping_add(0x9E3779B97F4A7C15);
    h1 ^= h1 >> 17;
    let mut h2 = item
        .wrapping_mul(0xBF58476D1CE4E5B9)
        .wrapping_add(0xBF58476D1CE4E5B9);
    h2 ^= h2 >> 17;
    (h1, h2)
}

/// Issue a prefetch for the given address into L1 cache.
#[inline(always)]
fn prefetch_read(ptr: *const u8) {
    #[cfg(target_arch = "aarch64")]
    unsafe {
        std::arch::asm!("prfm pldl1keep, [{ptr}]", ptr = in(reg) ptr);
    }
    #[cfg(target_arch = "x86_64")]
    unsafe {
        std::arch::x86_64::_mm_prefetch::<{ std::arch::x86_64::_MM_HINT_T0 }>(ptr as *const i8);
    }
}

/// Cache-line blocked bloom filter.
///
/// Each block is 512 bits (8 × u64 = 1 cache line). Queries touch exactly
/// one block, so each lookup is a single cache-line fetch from DRAM.
pub struct MiniBloom {
    data: Vec<AtomicU64>,
    num_blocks_mask: u64,
    pub num_hashes: u32,
    pub num_bits: u64,
    count: usize,
}

unsafe impl Send for MiniBloom {}
unsafe impl Sync for MiniBloom {}

impl MiniBloom {
    /// Extract a 9-bit probe position from pre-mixed hash values.
    /// Probes 0–6 use non-overlapping 9-bit windows of h2 (63 bits).
    /// Probes 7+ use upper bits of h1 (above block-selection bits).
    #[inline(always)]
    fn extract_probe(h2: u64, h1_upper: u64, i: u32) -> u64 {
        if i < 7 {
            (h2 >> (i * 9)) & 511
        } else {
            (h1_upper >> ((i - 7) * 9)) & 511
        }
    }

    /// Number of 512-bit blocks in this filter.
    #[inline(always)]
    pub fn num_blocks(&self) -> usize {
        (self.num_blocks_mask + 1) as usize
    }

    /// Compute hash, issue prefetch for the target block, return prehashed state.
    #[inline(always)]
    pub fn hash_and_prefetch(&self, item: u64) -> (usize, u64, u64) {
        let (h1, h2) = bloom_hash(item);
        let block_start = ((h1 & self.num_blocks_mask) as usize) << 3;
        let h1_upper = h1 >> (self.num_blocks_mask + 1).trailing_zeros();
        let ptr = unsafe { self.data.as_ptr().add(block_start) };
        prefetch_read(ptr as *const u8);
        (block_start, h2, h1_upper)
    }

    /// Insert an item using atomic operations (safe to call from &self).
    #[inline(always)]
    pub fn insert_prehashed_item(&self, item: u64) {
        let (h1, h2) = bloom_hash(item);
        let block_start = ((h1 & self.num_blocks_mask) as usize) << 3;
        let h1_upper = h1 >> (self.num_blocks_mask + 1).trailing_zeros();
        self.insert_prehashed(block_start, h2, h1_upper);
    }

    /// Insert with compile-time-known hash count for full loop unrolling.
    #[inline(always)]
    fn insert_prehashed_unrolled<const HASHES: u32>(
        &self,
        block_start: usize,
        h2: u64,
        h1_upper: u64,
    ) {
        for i in 0..HASHES {
            let bit_pos = Self::extract_probe(h2, h1_upper, i);
            let word_idx = (bit_pos >> 6) as usize;
            let bit_mask = 1u64 << (bit_pos & 63);
            unsafe { self.data.get_unchecked(block_start + word_idx) }
                .fetch_or(bit_mask, AtomicOrdering::Relaxed);
        }
    }

    /// Query with compile-time-known hash count for full loop unrolling.
    #[inline(always)]
    fn contains_prehashed_unrolled<const HASHES: u32>(
        &self,
        block_start: usize,
        h2: u64,
        h1_upper: u64,
    ) -> bool {
        for i in 0..HASHES {
            let bit_pos = Self::extract_probe(h2, h1_upper, i);
            let word_idx = (bit_pos >> 6) as usize;
            let bit_mask = 1u64 << (bit_pos & 63);
            if unsafe {
                self.data
                    .get_unchecked(block_start + word_idx)
                    .load(AtomicOrdering::Relaxed)
            } & bit_mask
                == 0
            {
                return false;
            }
        }
        true
    }

    /// Insert using pre-computed hash values (block already prefetched).
    /// Dispatches to a const-generic monomorphized version so LLVM can fully
    /// unroll the probe loop and resolve all shift amounts at compile time.
    #[inline(always)]
    pub fn insert_prehashed(&self, block_start: usize, h2: u64, h1_upper: u64) {
        match self.num_hashes {
            1 => self.insert_prehashed_unrolled::<1>(block_start, h2, h1_upper),
            2 => self.insert_prehashed_unrolled::<2>(block_start, h2, h1_upper),
            3 => self.insert_prehashed_unrolled::<3>(block_start, h2, h1_upper),
            4 => self.insert_prehashed_unrolled::<4>(block_start, h2, h1_upper),
            5 => self.insert_prehashed_unrolled::<5>(block_start, h2, h1_upper),
            6 => self.insert_prehashed_unrolled::<6>(block_start, h2, h1_upper),
            7 => self.insert_prehashed_unrolled::<7>(block_start, h2, h1_upper),
            _ => self.insert_prehashed_unrolled::<8>(block_start, h2, h1_upper),
        }
    }

    /// Query using pre-computed hash values (block already prefetched).
    /// Dispatches to a const-generic monomorphized version so LLVM can fully
    /// unroll the probe loop and resolve all shift amounts at compile time.
    #[inline(always)]
    pub fn contains_prehashed(&self, block_start: usize, h2: u64, h1_upper: u64) -> bool {
        match self.num_hashes {
            1 => self.contains_prehashed_unrolled::<1>(block_start, h2, h1_upper),
            2 => self.contains_prehashed_unrolled::<2>(block_start, h2, h1_upper),
            3 => self.contains_prehashed_unrolled::<3>(block_start, h2, h1_upper),
            4 => self.contains_prehashed_unrolled::<4>(block_start, h2, h1_upper),
            5 => self.contains_prehashed_unrolled::<5>(block_start, h2, h1_upper),
            6 => self.contains_prehashed_unrolled::<6>(block_start, h2, h1_upper),
            7 => self.contains_prehashed_unrolled::<7>(block_start, h2, h1_upper),
            _ => self.contains_prehashed_unrolled::<8>(block_start, h2, h1_upper),
        }
    }

    /// Compute per-block fill statistics: (avg_fill, min_fill, max_fill).
    pub fn fill_stats(&self) -> (f64, f64, f64) {
        let num_blocks = self.num_blocks();
        let mut total_bits_set = 0u64;
        let mut min_block = 512u32;
        let mut max_block = 0u32;
        for b in 0..num_blocks {
            let base = b * 8;
            let mut block_bits = 0u32;
            for w in 0..8 {
                block_bits += self.data[base + w]
                    .load(AtomicOrdering::Relaxed)
                    .count_ones();
            }
            total_bits_set += block_bits as u64;
            min_block = min_block.min(block_bits);
            max_block = max_block.max(block_bits);
        }
        let avg = total_bits_set as f64 / (num_blocks as f64 * 512.0);
        let min_f = min_block as f64 / 512.0;
        let max_f = max_block as f64 / 512.0;
        (avg, min_f, max_f)
    }

    /// Persist bloom filter to disk in nucleaze binary format.
    pub fn save(&self, path: &str, k: u64, kmer_count: u64) -> std::io::Result<()> {
        let mut w = BufWriter::new(File::create(path)?);
        w.write_all(b"NBLM")?;
        w.write_all(&2u32.to_le_bytes())?; // version 2 = blocked
        w.write_all(&(self.num_blocks() as u64).to_le_bytes())?;
        w.write_all(&self.num_hashes.to_le_bytes())?;
        w.write_all(&k.to_le_bytes())?;
        w.write_all(&kmer_count.to_le_bytes())?;
        for word in &self.data {
            w.write_all(&word.load(AtomicOrdering::Relaxed).to_le_bytes())?;
        }
        w.flush()
    }

    /// Load bloom filter from nucleaze binary format.
    pub fn load(path: &str) -> std::io::Result<(Self, u64, u64)> {
        let mut r = BufReader::new(File::open(path)?);
        let mut magic = [0u8; 4];
        r.read_exact(&mut magic)?;
        if &magic != b"NBLM" {
            return Err(std::io::Error::new(
                std::io::ErrorKind::InvalidData,
                "not a nucleaze bloom file",
            ));
        }
        let mut buf4 = [0u8; 4];
        let mut buf8 = [0u8; 8];
        r.read_exact(&mut buf4)?; // version
        r.read_exact(&mut buf8)?;
        let num_blocks = u64::from_le_bytes(buf8) as usize;
        r.read_exact(&mut buf4)?;
        let num_hashes = u32::from_le_bytes(buf4);
        r.read_exact(&mut buf8)?;
        let k = u64::from_le_bytes(buf8);
        r.read_exact(&mut buf8)?;
        let kmer_count = u64::from_le_bytes(buf8);

        let total_words = num_blocks * 8;
        let mut data = Vec::with_capacity(total_words);
        for _ in 0..total_words {
            r.read_exact(&mut buf8)?;
            data.push(AtomicU64::new(u64::from_le_bytes(buf8)));
        }

        Ok((
            MiniBloom {
                data,
                num_blocks_mask: num_blocks as u64 - 1,
                num_hashes,
                num_bits: (num_blocks * 512) as u64,
                count: kmer_count as usize,
            },
            k,
            kmer_count,
        ))
    }
}

impl KmerStore for MiniBloom {
    fn new(size: usize, fpr: f64) -> Self {
        let (num_blocks, num_hashes) = bloom_params(size, fpr);
        let total_words = num_blocks * 8;
        let num_bits = (num_blocks * 512) as u64;

        println!(
            "Bloom filter: {} bits ({:.1} MB), {} hash functions, {} blocks",
            num_bits,
            total_words as f64 * 8.0 / 1_000_000.0,
            num_hashes,
            num_blocks,
        );

        MiniBloom {
            data: {
                let raw: Vec<u64> = vec![0u64; total_words];
                // SAFETY: AtomicU64 has the same size, alignment, and bit validity as u64.
                // Zero bits is a valid AtomicU64 representing the value 0.
                unsafe { std::mem::transmute::<Vec<u64>, Vec<AtomicU64>>(raw) }
            },
            num_blocks_mask: num_blocks as u64 - 1,
            num_hashes: num_hashes as u32,
            num_bits,
            count: 0,
        }
    }

    /// Shard index for subindex buffering (not used for bloom bit positions).
    fn map(&self, kmer: &u64) -> usize {
        let hashed = kmer ^ (kmer >> 12);
        hashed as usize & 1023
    }

    /// Insert k-mer by setting bits in a single cache-line block.
    fn insert(&mut self, kmer: &u64) {
        let (h1, h2) = bloom_hash(*kmer);
        let block_start = ((h1 & self.num_blocks_mask) as usize) << 3;
        let h1_upper = h1 >> (self.num_blocks_mask + 1).trailing_zeros();
        for i in 0..self.num_hashes {
            let bit_pos = Self::extract_probe(h2, h1_upper, i);
            let word_idx = (bit_pos >> 6) as usize;
            let bit_mask = 1u64 << (bit_pos & 63);
            self.data[block_start + word_idx].fetch_or(bit_mask, AtomicOrdering::Relaxed);
        }
        self.count += 1;
    }

    /// Buffer k-mer into subindex for later batch insertion.
    fn stage(&self, kmer: u64, subindex: &mut Vec<Vec<u64>>) {
        let idx = self.map(&kmer);
        subindex[idx].push(kmer);
    }

    /// Query bloom filter: all probes hit the same cache-line block.
    /// Early-exits on first missing bit.
    fn contains(&self, kmer: &u64) -> bool {
        let (h1, h2) = bloom_hash(*kmer);
        let block_start = ((h1 & self.num_blocks_mask) as usize) << 3;
        let h1_upper = h1 >> (self.num_blocks_mask + 1).trailing_zeros();
        for i in 0..self.num_hashes {
            let bit_pos = Self::extract_probe(h2, h1_upper, i);
            let word_idx = (bit_pos >> 6) as usize;
            let bit_mask = 1u64 << (bit_pos & 63);
            if unsafe {
                self.data
                    .get_unchecked(block_start + word_idx)
                    .load(AtomicOrdering::Relaxed)
            } & bit_mask
                == 0
            {
                return false;
            }
        }
        true
    }

    /// Clear Bloom filter.
    fn clear(&mut self) {
        for word in &self.data {
            word.store(0, AtomicOrdering::Relaxed);
        }
        self.count = 0;
    }

    /// Return k-mer counter.
    fn count(&self) -> usize {
        self.count
    }

    /// Add serialized k-mers to Bloom filter (parallel).
    fn from_serialized(&mut self, data: Vec<Vec<u64>>) {
        let num_blocks_mask = self.num_blocks_mask;
        let num_hashes = self.num_hashes;
        let bloom_data = &self.data;

        let total: usize = data
            .par_iter()
            .map(|shard| {
                for kmer in shard {
                    let (h1, h2) = bloom_hash(*kmer);
                    let block_start = ((h1 & num_blocks_mask) as usize) << 3;
                    let h1_upper = h1 >> (num_blocks_mask + 1).trailing_zeros();
                    for i in 0..num_hashes {
                        let bit_pos = MiniBloom::extract_probe(h2, h1_upper, i);
                        let word_idx = (bit_pos >> 6) as usize;
                        let bit_mask = 1u64 << (bit_pos & 63);
                        unsafe { bloom_data.get_unchecked(block_start + word_idx) }
                            .fetch_or(bit_mask, AtomicOrdering::Relaxed);
                    }
                }
                shard.len()
            })
            .sum();
        self.count += total;
    }

    /// Insert all buffered k-mers from the sharded accumulator into the Bloom filter (parallel).
    fn absorb_shards(&mut self, shards: Vec<FxHashSet<u64>>) {
        let num_blocks_mask = self.num_blocks_mask;
        let num_hashes = self.num_hashes;
        let bloom_data = &self.data;

        let total: usize = shards
            .par_iter()
            .map(|set| {
                for &kmer in set {
                    let (h1, h2) = bloom_hash(kmer);
                    let block_start = ((h1 & num_blocks_mask) as usize) << 3;
                    let h1_upper = h1 >> (num_blocks_mask + 1).trailing_zeros();
                    for i in 0..num_hashes {
                        let bit_pos = MiniBloom::extract_probe(h2, h1_upper, i);
                        let word_idx = (bit_pos >> 6) as usize;
                        let bit_mask = 1u64 << (bit_pos & 63);
                        unsafe { bloom_data.get_unchecked(block_start + word_idx) }
                            .fetch_or(bit_mask, AtomicOrdering::Relaxed);
                    }
                }
                set.len()
            })
            .sum();
        self.count += total;
    }

    fn needs_absorb(&self) -> bool {
        true
    }
}

/// K-mer encoding, storage, and sequence processing operations
pub struct KmerProcessor<S: KmerStore> {
    pub k: usize,
    pub k_cap: u64,
    pub threshold: u8,
    pub use_canonical: bool,
    pub ref_kmers: S,
}

impl<S: KmerStore> KmerProcessor<S> {
    pub fn new(k: usize, threshold: u8, use_canonical: bool, ref_kmers: S) -> Self {
        KmerProcessor {
            k,
            k_cap: if k >= 32 {
                u64::MAX
            } else {
                (1 << (k * 2)) - 1
            },
            threshold,
            use_canonical,
            ref_kmers,
        }
    }

    /// Add k-mers from sequence to reference index.
    pub fn process_ref(&self, seq: &[u8], subindex: &mut Vec<Vec<u64>>) {
        let mut kmer = 0u64;
        let mut valid = 0usize;

        for &base in seq {
            if let Some(bits) = encode(base) {
                kmer = ((kmer << 2) | bits) & self.k_cap;
                valid += 1;

                if valid >= self.k {
                    let kmer_to_store = if self.use_canonical {
                        canonical((kmer, self.k as u8)).0.0
                    } else {
                        kmer
                    };

                    self.stage_kmer(kmer_to_store, subindex);
                }
            } else {
                kmer = 0;
                valid = 0;
            }
        }
    }

    /// Check if a read has enough matching k-mers against the reference.
    pub fn process_read(&self, seq: &[u8]) -> bool {
        if seq.len() < self.k {
            return false;
        }

        let mut hits = 0;
        let mut kmer = 0;
        let mut rc_kmer = 0;
        let mut valid = 0;

        for &base in seq {
            match encode(base) {
                Some(base) => {
                    kmer = ((kmer << 2) | base) & self.k_cap;

                    // Update rc kmer: (rc >> 2) | (complement << shift)
                    let rc_base = base ^ 0b11; // Inverse of bits (A -> T)
                    rc_kmer = (rc_kmer >> 2) | (rc_base << (2 * (self.k - 1)));

                    valid += 1;

                    if valid >= self.k {
                        let is_hit = if self.use_canonical {
                            self.contains_kmer(&std::cmp::min(kmer, rc_kmer))
                        } else {
                            self.contains_kmer(&kmer)
                        };

                        if is_hit {
                            hits += 1;
                            if hits >= self.threshold {
                                return true;
                            }
                        }
                    }
                }
                None => {
                    valid = 0;
                    kmer = 0;
                    rc_kmer = 0;
                }
            }
        }

        false
    }

    // KmerProcessor wrapper for KmerStore::stage().
    #[inline(always)]
    pub fn stage_kmer(&self, kmer: u64, subindex: &mut Vec<Vec<u64>>) {
        self.ref_kmers.stage(kmer, subindex);
    }

    // KmerProcessor wrapper for KmerStore::insert().
    #[inline(always)]
    pub fn insert_kmer(&mut self, kmer: &u64) {
        self.ref_kmers.insert(kmer);
    }

    // KmerProcessor wrapper for KmerStore::contains().
    #[inline(always)]
    pub fn contains_kmer(&self, kmer: &u64) -> bool {
        self.ref_kmers.contains(kmer)
    }

    // KmerProcessor wrapper for KmerStore::count().
    pub fn num_kmers(&self) -> usize {
        self.ref_kmers.count()
    }

    pub fn needs_absorb(&self) -> bool {
        self.ref_kmers.needs_absorb()
    }
}

impl KmerProcessor<MiniBloom> {
    /// Extract k-mers from a sequence and insert them into the bloom filter.
    /// Uses a depth-16 ring buffer to prefetch cache lines far enough ahead
    /// to hide DRAM latency (~250 cycles) behind useful computation.
    /// Each loop iteration is ~17 cycles, so depth 16 gives ~272 cycles of
    /// lead time — enough to fully cover the round-trip to main memory.
    #[inline(always)]
    pub fn insert_kmers_bloom(&self, seq: &[u8]) -> u64 {
        const DEPTH: usize = 16;

        let bloom = &self.ref_kmers;
        let mut kmer = 0u64;
        let mut rc_kmer = 0u64;
        let mut valid_bases = 0usize;
        let mut count = 0u64;
        let rc_shift = 2 * (self.k - 1);

        let mut queue = [(0usize, 0u64, 0u64); DEPTH];
        let mut q_head: usize = 0;
        let mut q_len: usize = 0;

        for &base in seq {
            if let Some(bits) = encode(base) {
                kmer = ((kmer << 2) | bits) & self.k_cap;
                rc_kmer = (rc_kmer >> 2) | ((bits ^ 3) << rc_shift);
                valid_bases += 1;

                if valid_bases >= self.k {
                    let canonical = if self.use_canonical {
                        std::cmp::min(kmer, rc_kmer)
                    } else {
                        kmer
                    };
                    let pending = bloom.hash_and_prefetch(canonical);

                    // When full, the oldest prefetch has had DEPTH iterations to arrive in L1
                    if q_len == DEPTH {
                        let (bs, h2, h1u) = queue[q_head];
                        bloom.insert_prehashed(bs, h2, h1u);
                        count += 1;
                    } else {
                        q_len += 1;
                    }

                    queue[q_head] = pending;
                    q_head = (q_head + 1) & (DEPTH - 1);
                }
            } else {
                // Drain queue on ambiguous base
                while q_len > 0 {
                    let idx = (q_head + DEPTH - q_len) & (DEPTH - 1);
                    let (bs, h2, h1u) = queue[idx];
                    bloom.insert_prehashed(bs, h2, h1u);
                    count += 1;
                    q_len -= 1;
                }
                kmer = 0;
                rc_kmer = 0;
                valid_bases = 0;
                q_head = 0;
            }
        }

        // Drain remaining
        while q_len > 0 {
            let idx = (q_head + DEPTH - q_len) & (DEPTH - 1);
            let (bs, h2, h1u) = queue[idx];
            bloom.insert_prehashed(bs, h2, h1u);
            count += 1;
            q_len -= 1;
        }
        count
    }

    /// Check if a read matches the reference bloom filter.
    /// Uses a depth-16 prefetch ring buffer so the cache line for each k-mer's
    /// block is fetched from DRAM ~16 iterations (~272 cycles) before it is
    /// read, fully hiding main memory latency behind k-mer extraction work.
    #[inline(always)]
    pub fn process_read_bloom(&self, seq: &[u8]) -> bool {
        const DEPTH: usize = 16;

        if seq.len() < self.k {
            return false;
        }
        let bloom = &self.ref_kmers;
        let mut hits: u8 = 0;
        let mut kmer: u64 = 0;
        let mut rc_kmer: u64 = 0;
        let mut valid_bases: usize = 0;
        let rc_shift = 2 * (self.k - 1);

        let mut queue = [(0usize, 0u64, 0u64); DEPTH];
        let mut q_head: usize = 0;
        let mut q_len: usize = 0;

        for &base in seq {
            match encode(base) {
                Some(bits) => {
                    kmer = ((kmer << 2) | bits) & self.k_cap;
                    rc_kmer = (rc_kmer >> 2) | ((bits ^ 3) << rc_shift);
                    valid_bases += 1;

                    if valid_bases >= self.k {
                        let canonical = if self.use_canonical {
                            std::cmp::min(kmer, rc_kmer)
                        } else {
                            kmer
                        };
                        let pending = bloom.hash_and_prefetch(canonical);

                        if q_len == DEPTH {
                            let (bs, h2, h1u) = queue[q_head];
                            if bloom.contains_prehashed(bs, h2, h1u) {
                                hits += 1;
                                if hits >= self.threshold {
                                    return true;
                                }
                            }
                        } else {
                            q_len += 1;
                        }

                        queue[q_head] = pending;
                        q_head = (q_head + 1) & (DEPTH - 1);
                    }
                }
                None => {
                    // Drain queue on ambiguous base
                    while q_len > 0 {
                        let idx = (q_head + DEPTH - q_len) & (DEPTH - 1);
                        let (bs, h2, h1u) = queue[idx];
                        if bloom.contains_prehashed(bs, h2, h1u) {
                            hits += 1;
                            if hits >= self.threshold {
                                return true;
                            }
                        }
                        q_len -= 1;
                    }
                    valid_bases = 0;
                    kmer = 0;
                    rc_kmer = 0;
                    q_head = 0;
                }
            }
        }

        // Drain remaining
        while q_len > 0 {
            let idx = (q_head + DEPTH - q_len) & (DEPTH - 1);
            let (bs, h2, h1u) = queue[idx];
            if bloom.contains_prehashed(bs, h2, h1u) {
                hits += 1;
                if hits >= self.threshold {
                    return true;
                }
            }
            q_len -= 1;
        }
        false
    }
}

impl KmerProcessor<HashShards> {
    pub fn serialize_kmers(&self) -> Vec<Vec<u64>> {
        self.ref_kmers.serialize()
    }
}

#[inline(always)]
/// Encodes UTF8 bases to 2 bits.
pub fn encode(b: u8) -> Option<u64> {
    // A=00, C=01, G=10, T/U=11
    static BASE_TABLE: [u8; 256] = {
        let mut bases = [0xFF; 256];
        bases[b'A' as usize] = 0;
        bases[b'a' as usize] = 0;
        bases[b'C' as usize] = 1;
        bases[b'c' as usize] = 1;
        bases[b'G' as usize] = 2;
        bases[b'g' as usize] = 2;
        bases[b'T' as usize] = 3;
        bases[b't' as usize] = 3;
        bases[b'U' as usize] = 3;
        bases[b'u' as usize] = 3;
        bases
    };

    let v = unsafe { *BASE_TABLE.get_unchecked(b as usize) };
    if v == 0xFF { None } else { Some(v as u64) }
}
