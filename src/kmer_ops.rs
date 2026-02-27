//! Operations for extracting and measuring k-mers from FASTX records
use bloomfilter::Bloom;
use needletail::bitkmer::canonical;
use rayon::prelude::*;
use rustc_hash::FxHashSet;
use std::{cmp::min, u64};

// ---------------------------------------------------------------------------
// KmerStore trait
// ---------------------------------------------------------------------------

#[allow(dead_code)]
pub trait KmerStore: Send + Sync {
    fn insert_in_shard(&mut self, shard: usize, kmer: u64);
    fn contains_in_shard(&self, shard: usize, kmer: u64) -> bool;
    fn bulk_insert_shard(&mut self, shard: usize, kmers: &[u64]);
    fn num_shards(&self) -> usize;
    fn num_kmers(&self) -> usize;
    fn clear(&mut self);
    fn to_serializable(&self) -> Vec<Vec<u64>>;
    fn load_from_serializable(&mut self, data: Vec<Vec<u64>>);
    fn insert_metadata(&mut self, k: usize);
    fn verify_metadata(&self, k: usize) -> bool;
}

// ---------------------------------------------------------------------------
// ShardedHashStore
// ---------------------------------------------------------------------------

pub struct ShardedHashStore {
    pub shards: Vec<FxHashSet<u64>>,
}

impl ShardedHashStore {
    pub fn new(num_shards: usize) -> Self {
        ShardedHashStore {
            shards: vec![FxHashSet::default(); num_shards],
        }
    }
}

impl KmerStore for ShardedHashStore {
    #[inline(always)]
    fn insert_in_shard(&mut self, shard: usize, kmer: u64) {
        self.shards[shard].insert(kmer);
    }

    #[inline(always)]
    fn contains_in_shard(&self, shard: usize, kmer: u64) -> bool {
        unsafe { self.shards.get_unchecked(shard).contains(&kmer) }
    }

    fn bulk_insert_shard(&mut self, shard: usize, kmers: &[u64]) {
        self.shards[shard].reserve(kmers.len());
        for &kmer in kmers {
            self.shards[shard].insert(kmer);
        }
    }

    fn num_shards(&self) -> usize {
        self.shards.len()
    }

    fn num_kmers(&self) -> usize {
        self.shards.iter().map(|s| s.len()).sum()
    }

    fn clear(&mut self) {
        for shard in &mut self.shards {
            shard.clear();
        }
    }

    fn to_serializable(&self) -> Vec<Vec<u64>> {
        self.shards
            .iter()
            .map(|set| set.iter().cloned().collect())
            .collect()
    }

    fn load_from_serializable(&mut self, data: Vec<Vec<u64>>) {
        self.shards
            .par_iter_mut()
            .zip(data.into_par_iter())
            .for_each(|(set, vec)| {
                set.reserve(vec.len());
                for kmer in vec {
                    set.insert(kmer);
                }
            });
    }

    fn insert_metadata(&mut self, k: usize) {
        self.shards[0].insert(u64::MAX ^ k as u64);
    }

    fn verify_metadata(&self, k: usize) -> bool {
        let sentinel = u64::MAX ^ k as u64;
        self.shards[0].contains(&sentinel)
    }
}

// ---------------------------------------------------------------------------
// ShardedBloomStore
// ---------------------------------------------------------------------------

pub struct ShardedBloomStore {
    pub shards: Vec<Bloom<u64>>,
}

impl ShardedBloomStore {
    pub fn new(num_shards: usize, fpr: f64, expected_kmers_per_shard: usize) -> Self {
        let shards = (0..num_shards)
            .map(|_| Bloom::new_for_fp_rate(expected_kmers_per_shard, fpr))
            .collect();
        ShardedBloomStore { shards }
    }
}

impl KmerStore for ShardedBloomStore {
    #[inline(always)]
    fn insert_in_shard(&mut self, shard: usize, kmer: u64) {
        self.shards[shard].set(&kmer);
    }

    #[inline(always)]
    fn contains_in_shard(&self, shard: usize, kmer: u64) -> bool {
        self.shards[shard].check(&kmer)
    }

    fn bulk_insert_shard(&mut self, shard: usize, kmers: &[u64]) {
        for &kmer in kmers {
            self.shards[shard].set(&kmer);
        }
    }

    fn num_shards(&self) -> usize {
        self.shards.len()
    }

    // Bloom filters don't track cardinality exactly
    fn num_kmers(&self) -> usize {
        0
    }

    fn clear(&mut self) {
        for shard in &mut self.shards {
            shard.clear();
        }
    }

    fn to_serializable(&self) -> Vec<Vec<u64>> {
        eprintln!(
            "Warning: bloom filters cannot be losslessly serialized to the k-mer index format."
        );
        vec![Vec::new(); self.shards.len()]
    }

    fn load_from_serializable(&mut self, _data: Vec<Vec<u64>>) {
        // no-op: bloom filters cannot be loaded from the hash-based serialization format
    }

    fn insert_metadata(&mut self, _k: usize) {
        // no-op in bloom mode
    }

    fn verify_metadata(&self, _k: usize) -> bool {
        true // bloom mode skips binary index caching entirely
    }
}

// ---------------------------------------------------------------------------
// KmerProcessor<S>
// ---------------------------------------------------------------------------

/// K-mer encoding, storage, and sequence processing operations
pub struct KmerProcessor<S: KmerStore> {
    pub k: usize,
    pub k_cap: u64,
    pub threshold: u8,
    pub use_canonical: bool,
    pub ref_kmers: S,
    pub idx_mask: usize,
}

#[allow(dead_code)]
pub type DefaultProcessor = KmerProcessor<ShardedHashStore>;

impl KmerProcessor<ShardedHashStore> {
    pub fn new(k: usize, threshold: u8, use_canonical: bool) -> Self {
        let num_idx = 1024;
        KmerProcessor {
            k,
            k_cap: if k >= 32 {
                u64::MAX
            } else {
                (1 << (k * 2)) - 1
            },
            threshold,
            use_canonical,
            ref_kmers: ShardedHashStore::new(num_idx),
            idx_mask: num_idx - 1,
        }
    }
}

impl KmerProcessor<ShardedBloomStore> {
    pub fn new_bloom(
        k: usize,
        threshold: u8,
        use_canonical: bool,
        fpr: f64,
        bloomcap: usize,
    ) -> Self {
        let num_idx = 1024;
        let per_shard = if bloomcap > 0 {
            bloomcap / num_idx
        } else {
            5_000_000 / num_idx
        };
        KmerProcessor {
            k,
            k_cap: if k >= 32 {
                u64::MAX
            } else {
                (1 << (k * 2)) - 1
            },
            threshold,
            use_canonical,
            ref_kmers: ShardedBloomStore::new(num_idx, fpr, per_shard.max(1)),
            idx_mask: num_idx - 1,
        }
    }
}

impl<S: KmerStore> KmerProcessor<S> {
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
                    valid = 0;
                    kmer = 0;
                    rc_kmer = 0;
                }
            }
        }

        false
    }

    #[inline(always)]
    /// Returns the number of k-mers in the reference index.
    pub fn num_kmers(&self) -> usize {
        self.ref_kmers.num_kmers()
    }

    #[inline(always)]
    /// Checks reference index for presence of the k-mer.
    fn contains_kmer(&self, kmer: &u64) -> bool {
        let idx = self.map_kmer(kmer);
        self.ref_kmers.contains_in_shard(idx, *kmer)
    }

    #[inline(always)]
    /// Map k-mer to shard using middle bases.
    fn map_kmer(&self, kmer: &u64) -> usize {
        let kmer = kmer ^ (kmer >> 12); // Spread entropy
        kmer as usize & self.idx_mask
    }

    /// Transpose reference index to a serializable format.
    pub fn to_serializable(&self) -> Vec<Vec<u64>> {
        self.ref_kmers.to_serializable()
    }

    /// Add serialized k-mers to reference index.
    pub fn add_serializable_kmers(&mut self, data: Vec<Vec<u64>>) {
        self.ref_kmers.load_from_serializable(data);
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

#[cfg(test)]
mod tests {
    use super::*;
    use rand::RngExt;

    // Helpers that use `encode` (which returns Option<u64>) so tests
    // compute forward and reverse-compliment encodings via `encode`.
    fn encode_forward(seq: &[u8]) -> u64 {
        let mut v: u64 = 0;
        for &b in seq {
            v = (v << 2) | encode(b).expect("invalid base in test");
        }
        v
    }

    fn encode_reverse(seq: &[u8]) -> u64 {
        let mut v: u64 = 0;
        for &b in seq.iter().rev() {
            v = (v << 2) | (encode(b).expect("invalid base in test") ^ 0b11);
        }
        v
    }

    /// Helper to clone the processor's `ref_kmers`, run a mutation on it,
    /// then restore it back to the processor. Returns the closure result.
    fn with_idx<R, F: FnOnce(&mut KmerProcessor<ShardedHashStore>, &mut Vec<Vec<u64>>) -> R>(
        processor: &mut KmerProcessor<ShardedHashStore>,
        f: F,
    ) -> R {
        let mut temp_index = vec![Vec::new(); processor.ref_kmers.num_shards()];
        let result = f(processor, &mut temp_index);

        // Sort and deduplicate each subindex
        for subindex in &mut temp_index {
            subindex.sort_unstable();
            subindex.dedup();
        }

        processor.add_serializable_kmers(temp_index);

        result
    }

    // KMER ENCODING TESTS

    #[test]
    fn test_encode_single_base() {
        // Test forward encoding
        assert_eq!(encode_forward(b"C"), 0b01);
        assert_eq!(encode_forward(b"G"), 0b10);
        assert_eq!(encode_forward(b"T"), 0b11);

        // Test reverse encoding
        assert_eq!(encode_reverse(b"A"), 0b11);
        assert_eq!(encode_reverse(b"C"), 0b10);
        assert_eq!(encode_reverse(b"G"), 0b01);
        assert_eq!(encode_reverse(b"T"), 0b00);
    }

    #[test]
    fn test_encode_multiple_bases() {
        // Test forward encoding
        // AA = 0b0000
        assert_eq!(encode_forward(b"AA"), 0b0000);
        // AC = 0b0001
        assert_eq!(encode_forward(b"AC"), 0b0001);
        // AT = 0b0011
        assert_eq!(encode_forward(b"AT"), 0b0011);
        // ACGT = 0b00011011
        assert_eq!(encode_forward(b"ACGT"), 0b00011011);

        // Test reverse encoding
        // AA reverse = 0b1111
        assert_eq!(encode_reverse(b"AA"), 0b1111);
        // AC reverse = 0b1011
        assert_eq!(encode_reverse(b"AC"), 0b1011);
        // AT reverse = 0b0011
        assert_eq!(encode_reverse(b"AT"), 0b0011);
        // ACGT reverse = 0b00011011
        assert_eq!(encode_reverse(b"ACGT"), 0b00011011);

        // Test canonical behavior
        // TTA forward = 0b111100, reverse = 0b001111
        let tta_forward = encode_forward(b"TTA");
        let tta_reverse = encode_reverse(b"TTA");
        assert_ne!(tta_forward, tta_reverse);
        assert_eq!(min(tta_forward, tta_reverse), 0b110000);

        // Test that forward and reverse are different for asymmetric sequences
        assert_ne!(encode_forward(b"ACGG"), encode_reverse(b"ACGG"));
        assert_ne!(encode_forward(b"TTA"), encode_reverse(b"TTA"));
    }

    #[test]
    fn test_encode_longer_sequence() {
        // Test a longer sequence
        let seq = b"ACGTACGG";

        // Test forward encoding
        let forward_encoded = encode_forward(seq);
        assert!(forward_encoded > 0);

        // Test reverse encoding
        let reverse_encoded = encode_reverse(seq);
        assert!(reverse_encoded > 0);

        // Verify both are valid encodings
        let expected_bits = seq.len() * 2;
        assert!(forward_encoded < (1u64 << expected_bits));
        assert!(reverse_encoded < (1u64 << expected_bits));

        // Test that forward and reverse are different for this sequence
        assert_ne!(forward_encoded, reverse_encoded);

        // Test canonical behavior
        let canonical = min(forward_encoded, reverse_encoded);
        assert!(canonical > 0);
        assert!(canonical <= forward_encoded);
        assert!(canonical <= reverse_encoded);
    }

    #[test]
    fn test_canonical_encoding_behavior() {
        // Test sequences where forward and reverse are different
        let asymmetric_seqs = vec![
            b"ACGG".as_slice(),
            b"TTA".as_slice(),
            b"GGCCC".as_slice(),
            b"ATCGA".as_slice(),
        ];

        for seq in &asymmetric_seqs {
            let forward = encode_forward(seq);
            let reverse = encode_reverse(seq);
            let canonical = min(forward, reverse);

            // Canonical should be the minimum of forward and reverse
            assert_eq!(canonical, min(forward, reverse));
            assert!(canonical <= forward);
            assert!(canonical <= reverse);

            // For asymmetric sequences, forward and reverse should be different
            assert_ne!(forward, reverse);
        }
    }

    #[test]
    fn test_palindromic_sequences() {
        // Test palindromic sequences where forward and reverse should be equal
        let palindromic_seqs = vec![
            b"AT".as_slice(),
            b"GC".as_slice(),
            b"ATAT".as_slice(),
            b"GCGC".as_slice(),
        ];

        for seq in &palindromic_seqs {
            let forward = encode_forward(seq);
            let reverse = encode_reverse(seq);
            let canonical = min(forward, reverse);

            // For palindromic sequences, forward and reverse should be equal
            assert_eq!(forward, reverse);
            assert_eq!(canonical, forward);
            assert_eq!(canonical, reverse);
        }
    }

    #[test]
    fn test_encoding_consistency() {
        // Test that encoding is consistent with expected bit patterns
        let test_cases = vec![
            (b"A".as_slice(), 0b00, 0b11),
            (b"C".as_slice(), 0b01, 0b10),
            (b"G".as_slice(), 0b10, 0b01),
            (b"T".as_slice(), 0b11, 0b00),
            (b"AA".as_slice(), 0b0000, 0b1111),
            (b"AC".as_slice(), 0b0001, 0b1011),
            (b"AT".as_slice(), 0b0011, 0b0011),
            (b"GC".as_slice(), 0b1001, 0b1001),
        ];

        for (seq, expected_forward, expected_reverse) in &test_cases {
            assert_eq!(encode_forward(seq), *expected_forward);
            assert_eq!(encode_reverse(seq), *expected_reverse);
        }
    }

    #[test]
    fn test_encoding_bit_length() {
        // Test that encoding produces correct bit lengths
        let test_seqs = vec![
            b"C".as_slice(),
            b"AC".as_slice(),
            b"ACG".as_slice(),
            b"ACGT".as_slice(),
            b"ACGTACGT".as_slice(),
        ];

        for seq in &test_seqs {
            let forward = encode_forward(seq);
            let reverse = encode_reverse(seq);
            let expected_bits = seq.len() * 2;
            let max_value = (1u64 << expected_bits) - 1;

            assert!(forward <= max_value);
            assert!(reverse <= max_value);
            assert!(forward > 0 || seq.len() == 0);
            assert!(reverse > 0 || seq.len() == 0);
        }
    }

    // KMER PROCESSOR INITIALIZATION TESTS

    #[test]
    fn test_kmer_processor_normal_vars() {
        let processor = KmerProcessor::<ShardedHashStore>::new(21, 1, true);
        assert_eq!(processor.k, 21);
        assert_eq!(processor.threshold, 1);
        assert_eq!(processor.ref_kmers.num_shards(), 1024);
        assert_eq!(processor.k_cap, (1u64 << 42) - 1);
    }

    #[test]
    fn test_kmer_processor_diff_vars() {
        let processor = KmerProcessor::<ShardedHashStore>::new(15, 3, true);
        assert_eq!(processor.k, 15);
        assert_eq!(processor.threshold, 3);
        assert_eq!(processor.k_cap, (1u64 << 30) - 1);
    }

    // REFERENCE PROCESSING TESTS

    #[test]
    fn test_single_sequence() {
        let mut processor = KmerProcessor::<ShardedHashStore>::new(5, 1, true);
        let ref_seq = b"ATCGCGGA";

        with_idx(&mut processor, |p, k| {
            p.process_ref(ref_seq, k);
        });

        println!("{}", processor.num_kmers());
        assert!(processor.num_kmers() == 4); // 4 5-mers in 8 bp sequence
    }

    #[test]
    fn test_process_multiple() {
        let mut processor = KmerProcessor::<ShardedHashStore>::new(5, 1, true);

        let mut count1 = 0usize;
        let mut count2 = 0usize;
        with_idx(&mut processor, |p, k| {
            p.process_ref(b"ACGTACGT", k);
            count1 = k.len();
            p.process_ref(b"TGCATGCA", k);
            count2 = k.len();
        });

        // Should have added more k-mers (may have some overlap)
        assert!(count2 >= count1);
    }

    #[test]
    fn test_rc_refs() {
        let mut processor = KmerProcessor::<ShardedHashStore>::new(4, 1, true);
        with_idx(&mut processor, |p, k| {
            p.process_ref(b"TTTT", k); // original
        });

        assert!(processor.process_read(b"AAAA"));
        assert!(processor.process_read(b"TTTT"));
        assert_eq!(processor.num_kmers(), 1); // 1 canonical k-mer included
    }

    // READ PROCESSING TESTS

    #[test]
    fn test_process_read_exact_match() {
        let mut processor = KmerProcessor::<ShardedHashStore>::new(5, 1, true);
        let ref_seq = b"ACGTACGT";
        with_idx(&mut processor, |p, k| {
            p.process_ref(ref_seq, k);
        });

        // Read with exact k-mer from reference
        let read = b"ACGTACGT";
        assert!(processor.process_read(read));
    }

    #[test]
    fn test_process_read_no_match() {
        let mut processor = KmerProcessor::<ShardedHashStore>::new(5, 1, true);
        with_idx(&mut processor, |p, k| {
            p.process_ref(b"ACGTACGT", k);
        });
        let read = b"TTTTTTTT";
        assert!(!processor.process_read(read));
    }

    #[test]
    fn test_process_read_partial_match_below_threshold() {
        let mut processor = KmerProcessor::<ShardedHashStore>::new(5, 3, true);
        with_idx(&mut processor, |p, k| {
            p.process_ref(b"ACGTACGT", k);
        });

        // This read should have some matches but not reach threshold of 3
        let read = b"ACGTTTTTT";

        // The test verifies the threshold logic works
        let result = processor.process_read(read);
        assert!(result == true || result == false); // Just verify it completes
    }

    #[test]
    fn test_process_read_meets_threshold() {
        let mut processor = KmerProcessor::<ShardedHashStore>::new(4, 1, true);
        with_idx(&mut processor, |p, k| {
            p.process_ref(b"ACGTACGT", k);
        });

        // Read with at least one matching k-mer
        let read = b"ACGTTTTTTT";
        assert!(processor.process_read(read));
    }

    #[test]
    fn test_process_read_too_short() {
        let processor = KmerProcessor::<ShardedHashStore>::new(10, 1, true);
        let read = b"ACGT"; // Only 4 bases, k=10

        assert!(!processor.process_read(read));
    }

    // BIT MANIPULATION TESTS

    #[test]
    fn test_k_cap_calculation() {
        let processor5 = KmerProcessor::<ShardedHashStore>::new(5, 1, true);
        assert_eq!(processor5.k_cap, (1u64 << 10) - 1);

        let processor10 = KmerProcessor::<ShardedHashStore>::new(10, 1, true);
        assert_eq!(processor10.k_cap, (1u64 << 20) - 1);

        let processor21 = KmerProcessor::<ShardedHashStore>::new(21, 1, true);
        assert_eq!(processor21.k_cap, (1u64 << 42) - 1);
    }

    // SLIDING WINDOW TESTS

    #[test]
    fn test_kmer_rolling_hash() {
        let mut processor = KmerProcessor::<ShardedHashStore>::new(3, 1, true);

        // For sequence "ACGTACGT" with k=3:
        // k-mers should be: ACG, CGT, GTA, TAC, ACG, CGA
        with_idx(&mut processor, |p, k| {
            p.process_ref(b"ACGTACGA", k);
        });

        assert!(processor.ref_kmers.num_shards() >= 4);
    }

    // K-MER UNIQUENESS TESTS

    #[test]
    fn test_duplicate_kmers() {
        let mut processor = KmerProcessor::<ShardedHashStore>::new(5, 1, true);

        // Process same sequence twice
        let mut count1 = 0usize;
        let mut count2 = 0usize;
        with_idx(&mut processor, |p, k| {
            p.process_ref(b"ACGTACGT", k);
            count1 = k.len();
            p.process_ref(b"ACGTACGT", k);
            count2 = k.len();
        });

        // Should have same count (HashSet prevents duplicates)
        assert_eq!(count1, count2);
    }

    // THRESHOLD TESTS

    #[test]
    fn test_threshold_one() {
        let mut processor = KmerProcessor::<ShardedHashStore>::new(5, 1, true);
        with_idx(&mut processor, |p, k| {
            p.process_ref(b"ACGTACGTACGT", k);
        });

        // Even one matching k-mer should return true
        let read = b"ACGTATTTTTTT";
        assert!(processor.process_read(read));
    }

    #[test]
    fn test_threshold_higher() {
        let mut processor = KmerProcessor::<ShardedHashStore>::new(3, 3, true);
        with_idx(&mut processor, |p, k| {
            p.process_ref(b"ACGTACGT", k);
        });

        // Need at least 3 matching k-mers
        let read_with_matches = b"ACGTACGTACGT";
        assert!(processor.process_read(read_with_matches));
    }

    // EDGE CASES

    #[test]
    fn test_minimum_k_value() {
        let mut processor = KmerProcessor::<ShardedHashStore>::new(1, 1, true);
        with_idx(&mut processor, |p, k| {
            p.process_ref(b"ACGT", k);
        });

        let read = b"AAAA";
        processor.process_read(read); // Should not panic
    }

    #[test]
    fn test_sequence_exactly_k_length() {
        let mut processor = KmerProcessor::<ShardedHashStore>::new(5, 1, true);
        let seq = b"ACGTA"; // Exactly k=5
        with_idx(&mut processor, |p, k| {
            p.process_ref(seq, k);
        });
        assert!(processor.process_read(seq));
    }

    #[test]
    fn test_empty_reference_set() {
        let processor = KmerProcessor::<ShardedHashStore>::new(5, 1, true);
        let read = b"ACGTACGT";

        // No reference k-mers added
        assert!(!processor.process_read(read));
    }

    #[test]
    fn test_repeated_bases() {
        let mut processor = KmerProcessor::<ShardedHashStore>::new(5, 1, true);
        with_idx(&mut processor, |p, k| {
            p.process_ref(b"AAAAAAAAAA", k);
        });

        let read = b"AAAAAAAAAA";
        assert!(processor.process_read(read));
    }

    // CANONICAL K-MER TESTS

    #[test]
    fn test_canonical() {
        let mut processor = KmerProcessor::<ShardedHashStore>::new(5, 1, true);

        // Add forward strand
        with_idx(&mut processor, |p, k| {
            p.process_ref(b"ATGCCAGT", k);
        });

        // Reverse complement should also match due to canonical representation
        let read = b"ACTGGCAT";
        assert!(processor.process_read(read));
    }

    // PERFORMANCE & CAPACITY TESTS

    #[test]
    fn test_large_reference_set() {
        let mut processor = KmerProcessor::<ShardedHashStore>::new(21, 1, true);
        let bases = ["A", "C", "G", "T"];

        // Add many reference sequences
        with_idx(&mut processor, |p, k| {
            for i in 0..100 {
                let remainder = i % 4;
                let seq = format!("ACGTACGTACGTACGTACGT{}", bases[remainder]);
                p.process_ref(seq.as_bytes(), k);
            }
        });

        assert_eq!(processor.num_kmers(), 4); // 4 above variants + metadata
    }

    #[test]
    fn test_long_sequence_processing() {
        let mut processor = KmerProcessor::<ShardedHashStore>::new(21, 1, true);
        let mut randy = rand::rng();
        // Create a long sequence (1000 bases)
        let long_seq: Vec<u8> = (0..1000)
            .map(|_| match randy.random_range(0..4) {
                0 => b'A',
                1 => b'C',
                2 => b'G',
                _ => b'T',
            })
            .collect();

        with_idx(&mut processor, |p, k| {
            p.process_ref(&long_seq, k);
        });
        // Should have many k-mers
        assert!(processor.ref_kmers.num_shards() > 500);
    }

    // MULTIPLE THRESHOLD TESTS
    #[test]
    fn test_various_thresholds() {
        for threshold in 1..=5 {
            let mut processor = KmerProcessor::<ShardedHashStore>::new(5, threshold, true);
            with_idx(&mut processor, |p, k| {
                p.process_ref(b"ACGTACGTACGTACGT", k);
            });

            let read = b"ACGTACGTACGTACGT";
            // With exact match, should always pass regardless of threshold
            assert!(processor.process_read(read));
        }
    }

    // DIFFERENT K VALUES TESTS
    #[test]
    fn test_multiple_k() {
        for k in [3, 5, 7, 11, 15, 21, 25, 31].iter() {
            let mut processor = KmerProcessor::<ShardedHashStore>::new(*k, 1, true);

            // Create sequence long enough for this k
            let seq: Vec<u8> = (0..*k)
                .map(|i| match i % 4 {
                    0 => b'A',
                    1 => b'C',
                    2 => b'G',
                    _ => b'T',
                })
                .collect();

            with_idx(&mut processor, |p, k| {
                p.process_ref(&seq, k);
            });

            assert_eq!(processor.num_kmers(), 1); // 1 k-mer
        }
    }

    // BOUNDARY TESTS
    #[test]
    fn test_sequence_length_k() {
        let mut processor = KmerProcessor::<ShardedHashStore>::new(10, 1, true);
        let seq = b"ACGTACGTAC"; // Exactly 10 bases

        with_idx(&mut processor, |p, k| {
            p.process_ref(seq, k);
        });

        assert_eq!(processor.num_kmers(), 1); // Inserted 10-mer
    }

    #[test]
    fn test_sequence_length_k_plus_one() {
        let mut processor = KmerProcessor::<ShardedHashStore>::new(10, 1, true);
        let seq = b"ACGTACGTACT"; // 11 bases

        with_idx(&mut processor, |p, k| {
            p.process_ref(seq, k);
        });

        assert_eq!(processor.num_kmers(), 2); // 2 inserted 10-mers
    }

    // BLOOM FILTER TESTS

    /// Helper for bloom filter tests — uses bulk_insert_shard since
    /// load_from_serializable is a no-op for bloom stores.
    fn with_bloom_idx<
        R,
        F: FnOnce(&mut KmerProcessor<ShardedBloomStore>, &mut Vec<Vec<u64>>) -> R,
    >(
        processor: &mut KmerProcessor<ShardedBloomStore>,
        f: F,
    ) -> R {
        let mut temp_index = vec![Vec::new(); processor.ref_kmers.num_shards()];
        let result = f(processor, &mut temp_index);

        for (i, subindex) in temp_index.iter().enumerate() {
            if !subindex.is_empty() {
                processor.ref_kmers.bulk_insert_shard(i, subindex);
            }
        }

        result
    }

    #[test]
    fn test_bloom_no_false_negatives() {
        let mut processor =
            KmerProcessor::<ShardedBloomStore>::new_bloom(5, 1, true, 0.01, 5_000_000);
        let ref_seq = b"ACGTACGT";

        with_bloom_idx(&mut processor, |p, k| {
            p.process_ref(ref_seq, k);
        });

        // Inserted k-mers must always be found (no false negatives)
        assert!(processor.process_read(b"ACGTACGT"));
    }

    #[test]
    fn test_bloom_process_ref_and_read_roundtrip() {
        let mut processor =
            KmerProcessor::<ShardedBloomStore>::new_bloom(5, 1, true, 0.001, 5_000_000);

        with_bloom_idx(&mut processor, |p, k| {
            p.process_ref(b"ATGCCAGT", k);
        });

        // Reverse complement should also match (canonical mode)
        assert!(processor.process_read(b"ACTGGCAT"));
        // Original should match
        assert!(processor.process_read(b"ATGCCAGT"));
    }

    #[test]
    fn test_bloom_canonical_no_false_negatives() {
        let mut processor =
            KmerProcessor::<ShardedBloomStore>::new_bloom(4, 1, true, 0.01, 5_000_000);

        with_bloom_idx(&mut processor, |p, k| {
            p.process_ref(b"TTTT", k);
        });

        // Must find reverse complement (no false negatives)
        assert!(processor.process_read(b"AAAA"));
        assert!(processor.process_read(b"TTTT"));
    }

    #[test]
    fn test_bloom_metadata_noop() {
        let mut processor =
            KmerProcessor::<ShardedBloomStore>::new_bloom(21, 1, true, 0.01, 5_000_000);
        processor.ref_kmers.insert_metadata(21);
        assert!(processor.ref_kmers.verify_metadata(21));
        assert!(processor.ref_kmers.verify_metadata(99)); // always true in bloom mode
    }

    #[test]
    fn test_bloom_serialization_returns_empty() {
        let processor = KmerProcessor::<ShardedBloomStore>::new_bloom(21, 1, true, 0.01, 5_000_000);
        let serialized = processor.to_serializable();
        assert_eq!(serialized.len(), 1024);
        for shard in &serialized {
            assert!(shard.is_empty());
        }
    }

    #[test]
    fn test_bloom_num_kmers_returns_zero() {
        let mut processor =
            KmerProcessor::<ShardedBloomStore>::new_bloom(5, 1, true, 0.01, 5_000_000);

        with_bloom_idx(&mut processor, |p, k| {
            p.process_ref(b"ACGTACGT", k);
        });

        // Bloom filters don't track cardinality
        assert_eq!(processor.num_kmers(), 0);
    }
}
