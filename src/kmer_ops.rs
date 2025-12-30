use std::cmp::min;

/// K-mer processor for building reference indices and filtering reads
#[derive(Clone)]
pub struct KmerProcessor {
    pub k: usize,                   
    pub threshold: u8,              
    pub use_canonical: bool,        
    pub ref_kmers: Vec<Vec<u64>>,   
    pub bit_cap: u64,               
    pub kmer_id_length: usize,      
}

impl KmerProcessor {
    pub fn new(k: usize, threshold: u8, use_canonical: bool) -> Self {
        // Initialize with empty vectors to avoid index out of bounds during parallel build
        let kmer_id_length = 12;
        let num_partitions = 1 << kmer_id_length;
        let mut ref_kmers = vec![Vec::new(); num_partitions];
        
        // Store metadata in the first bucket (optional, kept for compatibility)
        let metadata = u64::MAX ^ k as u64;
        ref_kmers[0].push(metadata);

        KmerProcessor {
            k,
            threshold,
            use_canonical,
            ref_kmers,
            // Safe bit shifting for u64
            bit_cap: if k >= 32 { u64::MAX } else { (1u64 << (k * 2)) - 1 },
            kmer_id_length,
        }
    }

    /// Build reference k-mer index from a sequence
    pub fn process_ref(&self, seq: &[u8], kmer_index: &mut Vec<Vec<u64>>) {
        let mut kmer = 0u64;
        let mut valid_bases = 0usize;

        for &base in seq {
            if let Some(bits) = encode(base) {
                kmer = ((kmer << 2) | bits) & self.bit_cap;
                valid_bases += 1;

                if valid_bases >= self.k {
                    let kmer_to_store = if self.use_canonical {
                        needletail::bitkmer::canonical((kmer, 0)).0.0
                    } else {
                        kmer
                    };
                    
                    let sid = self.map_kmer(&kmer_to_store);
                    // Just push; sorting happens during merge phase
                    kmer_index[sid].push(kmer_to_store);
                }
            } else {
                kmer = 0;
                valid_bases = 0;
            }
        }
    }

    /// Check if a read has enough matching k-mers against the reference
    pub fn process_read(&self, seq: &[u8]) -> bool {
        if seq.len() < self.k {
            return false;
        }

        let mut hits: u8 = 0;
        let mut kmer = 0u64;
        let mut rc_kmer = 0u64;
        let mut valid_bases = 0;

        for &base in seq {
            match encode(base) {
                Some(encoded_base) => {
                    kmer = ((kmer << 2) | encoded_base) & self.bit_cap;
                    
                    // Update rc kmer: (rc >> 2) | (complement << shift)
                    let encoded_rc_base = encoded_base ^ 0b11; // Inverse of bits (A=00->T=11)
                    rc_kmer = (rc_kmer >> 2) | (encoded_rc_base << (2 * (self.k - 1)));
                    
                    valid_bases += 1;

                    if valid_bases >= self.k {
                        let is_hit = if self.use_canonical {
                            let canonical_kmer = min(kmer, rc_kmer);
                            self.contains_kmer(&canonical_kmer)
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
                    valid_bases = 0;
                    kmer = 0;
                    rc_kmer = 0;
                }
            }
        }

        false
    }

    #[inline(always)]
    fn contains_kmer(&self, kmer: &u64) -> bool {
        let kmer_id = self.map_kmer(kmer);
        self.ref_kmers[kmer_id].binary_search(kmer).is_ok()
    }

    #[inline(always)]
    fn map_kmer(&self, kmer: &u64) -> usize {
        (kmer >> (2 * self.k - self.kmer_id_length)) as usize
    }

}

// Optimized base lookup using const array (avoids stack init overhead)
#[inline(always)]
fn encode(b: u8) -> Option<u64> {
    const INVALID: u8 = 0xFF;
    // A=00, C=01, G=10, T/U=11
    static BASE_TABLE: [u8; 256] = {
        let mut bases = [INVALID; 256];
        bases[b'A' as usize] = 0; bases[b'a' as usize] = 0;
        bases[b'C' as usize] = 1; bases[b'c' as usize] = 1;
        bases[b'G' as usize] = 2; bases[b'g' as usize] = 2;
        bases[b'T' as usize] = 3; bases[b't' as usize] = 3;
        bases[b'U' as usize] = 3; bases[b'u' as usize] = 3;
        bases
    };
    
    let v = unsafe { *BASE_TABLE.get_unchecked(b as usize) };
    if v == INVALID { None } else { Some(v as u64) }
}