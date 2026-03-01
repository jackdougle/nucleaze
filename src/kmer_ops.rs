//! Operations for extracting and measuring k-mers from FASTX records
use needletail::bitkmer::canonical;
use rayon::prelude::*;
use rustc_hash::FxHashSet;

pub trait KmerStore: Send + Sync {
    fn new(size: usize, fpr: f64) -> Self;
    fn map(&self, kmer: &u64) -> usize;
    fn insert(&mut self, kmer: &u64);
    fn contains(&self, kmer: &u64) -> bool;
    fn count(&self) -> usize;
    fn clear(&mut self);
    fn from_serialized(&mut self, data: Vec<Vec<u64>>);
    /// Bulk-load pre-sharded k-mer sets from the reference indexing accumulator.
    fn absorb_shards(&mut self, shards: Vec<FxHashSet<u64>>);
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

    /// Check which shard would contain the k-mer.
    #[inline(always)]
    fn map(&self, kmer: &u64) -> usize {
        let hashed = kmer ^ (kmer >> 12);
        hashed as usize & self.mask
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

pub struct MiniBloom {
    index: Vec<u64>,
    mask: u64,
    count: usize,
}

impl KmerStore for MiniBloom {
    // Create a new Bloom filter sized for the number of items and fpr.
    fn new(size: usize, fpr: f64) -> Self {
        // Compute number of u64 words needed for optimal Bloom filter size.
        // m_bits = -(n * ln(p)) / (ln2)^2, then divide by 64 for u64 count.
        fn compute_size(n: usize, fpp: f64) -> usize {
            use std::f64::consts::LN_2;
            let ln2_2 = (LN_2 as f64) * (LN_2 as f64);
            let m_bits = -((n as f64) * f64::ln(fpp)) / ln2_2;
            (m_bits / 64.0).ceil().max(1.0) as usize
        }
        let num_idx = compute_size(size, fpr);
        println!("Number of indices in Bloom filter: {}", num_idx);

        MiniBloom {
            index: vec![0u64; num_idx],
            mask: (num_idx as u64 * 64 - 1),
            count: 0,
        }
    }

    /// O(1) Bloom filter insert.
    fn insert(&mut self, kmer: &u64) {
        let pos = kmer & self.mask;
        let i = (pos / 64) as usize;
        self.index[i] |= 1 << (pos % 64);
        self.count += 1;
    }

    /// Find k-mer's Bloom index.
    fn map(&self, kmer: &u64) -> usize {
        (kmer & self.mask) as usize
    }

    /// O(1) Bloom filter query.
    fn contains(&self, kmer: &u64) -> bool {
        let pos = self.map(kmer);
        let idx = (pos / 64) as usize;
        if self.index[idx] & 1 << (pos % 64) == 0 {
            return false;
        }
        true
    }

    /// Clear Bloom filter.
    fn clear(&mut self) {
        self.index.fill(0);
        self.count = 0;
    }

    /// Return k-mer counter.
    fn count(&self) -> usize {
        self.count
    }

    /// Add serialized k-mers to Bloom filter.
    fn from_serialized(&mut self, data: Vec<Vec<u64>>) {
        use std::sync::atomic::{AtomicU64, Ordering};
        let idx: &[AtomicU64] = unsafe {
            std::slice::from_raw_parts(self.index.as_ptr() as *const AtomicU64, self.index.len())
        };
        let total: usize = data
            .par_iter()
            .map(|shard| {
                for kmer in shard {
                    let pos = self.map(kmer);
                    let i = (pos / 64) as usize;
                    // Atomic insertion process
                    idx[i].fetch_or(1 << (pos % 64), Ordering::Relaxed);
                }
                shard.len()
            })
            .sum();
        self.count += total;
    }

    /// Absorb pre-sharded sets by inserting each k-mer into the Bloom filter.
    fn absorb_shards(&mut self, shards: Vec<FxHashSet<u64>>) {
        for set in shards {
            for kmer in set {
                self.insert(&kmer);
            }
        }
    }
}

/// K-mer encoding, storage, and sequence processing operations
#[derive(Clone)]
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
    pub fn process_ref(&self, seq: &[u8], idx: &mut Vec<Vec<u64>>) {
        let mut kmer = 0u64;
        let mut valid = 0usize;

        for &base in seq {
            if let Some(bits) = encode(base) {
                kmer = ((kmer << 2) | bits) & self.k_cap;
                valid += 1;

                if valid >= self.k {
                    let final_kmer = if self.use_canonical {
                        canonical((kmer, self.k as u8)).0.0
                    } else {
                        kmer
                    };

                    let sid = self.map_kmer(&final_kmer);
                    idx[sid].push(final_kmer);
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

    fn map_kmer(&self, kmer: &u64) -> usize {
        self.ref_kmers.map(kmer)
    }

    #[inline(always)]
    pub fn insert_kmer(&mut self, kmer: &u64) {
        self.ref_kmers.insert(kmer);
    }

    #[inline(always)]
    pub fn contains_kmer(&self, kmer: &u64) -> bool {
        self.ref_kmers.contains(kmer)
    }

    pub fn num_kmers(&self) -> usize {
        self.ref_kmers.count()
    }

    pub fn num_shards(&self) -> usize {
        1024
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
