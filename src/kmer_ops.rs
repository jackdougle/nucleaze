//! Operations for extracting and measuring k-mers from FASTX records
use needletail::bitkmer::canonical;
use std::{cmp::min, u64};

#[derive(Clone)]
struct QuickBloom {
    bits: Vec<u64>,
    fpr: f64,
    mask: u64,
    count: usize,
}

impl QuickBloom {
    // Create a new Bloom filter sized for the number of items and fpr
    fn new(size: Option<usize>, fpr: f64) -> Self {
        match size {
            Some(s) => {
                return QuickBloom {
                    bits: vec![0u64; s / 8],
                    fpr,
                    mask: s as u64 * 8 - 1,
                    count: 0,
                };
            }
            None => {
                let n = 5_000_000;
                return QuickBloom {
                    bits: vec![0u64; n],
                    fpr,
                    mask: n as u64 * 8 - 1,
                    count: 0,
                };
            }
        }
    }

    /// O(1) Bloom filter insert
    fn insert(&mut self, kmer: &u64) {
        let pos = kmer & self.mask;
        let idx = (pos / 64) as usize;
        self.bits[idx] |= 1 << (pos % 64);
        self.count += 1;
    }

    /// O(1) Bloom filter query
    fn contains(&self, kmer: &u64) -> bool {
        let pos = kmer & self.mask;
        let idx = (pos / 64) as usize;
        if self.bits[idx] & 1 << (pos % 64) == 0 {
            return false;
        }
        true
    }

    fn clear(&mut self) {
        self.bits.fill(0);
        self.count = 0;
    }

    fn len(&self) -> usize {
        self.count
    }
}

/// K-mer encoding, storage, and sequence processing operations
#[derive(Clone)]
pub struct KmerProcessor {
    pub k: usize,
    pub k_cap: u64,
    pub threshold: u8,
    pub use_canonical: bool,
    ref_kmers: QuickBloom,
    pub idx_mask: usize,
}

impl KmerProcessor {
    pub fn new(
        k: usize,
        threshold: u8,
        use_canonical: bool,
        size: Option<usize>,
        fpr: f64,
    ) -> Self {
        let num_idx = 1024;
        let _ = if k > 12 { 10 } else { (k / 2).max(1) }; // minimer length
        KmerProcessor {
            k,
            k_cap: if k >= 32 {
                u64::MAX
            } else {
                (1 << (k * 2)) - 1
            },
            threshold,
            use_canonical,
            ref_kmers: QuickBloom::new(size, fpr),
            idx_mask: num_idx - 1,
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
        let mut valids = 0;

        for &base in seq {
            match encode(base) {
                Some(base) => {
                    kmer = ((kmer << 2) | base) & self.k_cap;

                    // Update rc kmer: (rc >> 2) | (complement << shift)
                    let rc_base = base ^ 0b11; // Inverse of bits (A -> T)
                    rc_kmer = (rc_kmer >> 2) | (rc_base << (2 * (self.k - 1)));

                    valids += 1;

                    if valids >= self.k {
                        let is_hit = if self.use_canonical {
                            self.contains_kmer(&min(kmer, rc_kmer))
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
                    valids = 0;
                    kmer = 0;
                    rc_kmer = 0;
                }
            }
        }

        false
    }

    #[inline(always)]
    /// Map k-mer to shard using middle bases.
    fn map_kmer(&self, kmer: &u64) -> usize {
        let kmer = kmer ^ (kmer >> 12); // Spread entropy
        kmer as usize & self.idx_mask
    }

    pub fn contains_kmer(&self, kmer: &u64) -> bool {
        self.ref_kmers.contains(kmer)
    }

    /// Insert a k-mer into the bloom filter.
    pub fn insert_kmer(&mut self, kmer: &u64) {
        self.ref_kmers.insert(kmer);
    }

    /// Clear all k-mers from the bloom filter.
    pub fn clear_kmers(&mut self) {
        self.ref_kmers.clear();
    }

    /// Return the number of k-mers inserted.
    pub fn num_kmers(&self) -> usize {
        self.ref_kmers.len()
    }

    /// Return the number of collection shards.
    pub fn num_shards(&self) -> usize {
        self.idx_mask + 1
    }

    /// Tranpose reference index to a serializable format.
    pub fn to_serializable(&self) -> Vec<Vec<u64>> {
        unimplemented!()
        // self.ref_kmers
        //     .iter()
        //     .map(|set| set.iter().cloned().collect())
        //     .collect()
    }

    /// Add serialized k-mers to reference index.
    pub fn add_serializable_kmers(&mut self, _data: Vec<Vec<u64>>) {
        unimplemented!()
        // self.ref_kmers
        //     .par_iter_mut()
        //     .zip(data.into_par_iter())
        //     .for_each(|(set, vec)| {
        //         // Reserve space to avoid reallocations
        //         set.reserve(vec.len());
        //         for kmer in vec {
        //             set.insert(kmer);
        //         }
        //     });
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
