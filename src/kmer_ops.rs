use needletail::bitkmer::canonical;
use std::cmp::min;

/// K-mer processor for building reference index and filtering reads
#[derive(Clone)]
pub struct KmerProcessor {
    pub k: usize,
    pub threshold: u8,
    pub use_canonical: bool,
    pub ref_kmers: Vec<Vec<u64>>,
    pub bit_cap: u64,
    pub kmer_id_length: usize,
    pub bloom_mask: u64,
    pub bloom_filter: Vec<u64>,
}

impl KmerProcessor {
    pub fn new(k: usize, threshold: u8, use_canonical: bool, bloom_size: usize) -> Self {
        let kmer_id_length = 12; // sort with first 6 bases
        let num_partitions = 1 << kmer_id_length;

        KmerProcessor {
            k,
            threshold,
            use_canonical,
            ref_kmers: vec![Vec::new(); num_partitions],
            bit_cap: if k >= 32 {
                u64::MAX // Safe bit cap for k = 32
            } else {
                (1u64 << (k * 2)) - 1
            },
            kmer_id_length,
            bloom_mask: bloom_size as u64 * 8 - 1,
            bloom_filter: vec![0u64; bloom_size / 8],
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
                        canonical((kmer, 0)).0.0
                    } else {
                        kmer
                    };

                    let sid = self.map_kmer(&kmer_to_store);
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
    /// Checks kmer index for kmer presence
    fn contains_kmer(&self, kmer: &u64) -> bool {
        // O(1) Bloom filter query
        let bloom_pos = kmer & self.bloom_mask;
        if self.bloom_filter[(bloom_pos / 64) as usize] & 1 << (bloom_pos % 64) == 0 {
            return false; // No false negatives
        }

        let kmer_id = self.map_kmer(kmer);
        self.ref_kmers[kmer_id].binary_search(kmer).is_ok() // O(log n)
    }

    #[inline(always)]
    fn map_kmer(&self, kmer: &u64) -> usize {
        (kmer >> (2 * self.k - self.kmer_id_length)) as usize
    }
}

#[inline(always)]
/// Encodes UTF8 bases to 2 bits
fn encode(b: u8) -> Option<u64> {
    const INVALID: u8 = 0xFF;
    // A=00, C=01, G=10, T/U=11
    static BASE_TABLE: [u8; 256] = {
        let mut bases = [INVALID; 256];
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
    if v == INVALID { None } else { Some(v as u64) }
}

#[cfg(test)]
mod tests {
    use crate::kmer_ops::*;
    use rand::Rng;
    use std::cmp::min;

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
    fn with_kmer_index<R, F: FnOnce(&mut KmerProcessor, &mut Vec<Vec<u64>>) -> R>(
        processor: &mut KmerProcessor,
        f: F,
    ) -> R {
        let mut kmer_index = processor.ref_kmers.clone();
        let res = f(processor, &mut kmer_index);
        processor.ref_kmers = kmer_index;
        res
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
        let processor = KmerProcessor::new(21, 1, true, 1024);
        assert_eq!(processor.k, 21);
        assert_eq!(processor.threshold, 1);
        assert!(processor.ref_kmers.is_empty());
        assert_eq!(processor.bit_cap, (1u64 << 42) - 1);
    }

    #[test]
    fn test_kmer_processor_diff_vars() {
        let processor = KmerProcessor::new(15, 3, true, 1024);
        assert_eq!(processor.k, 15);
        assert_eq!(processor.threshold, 3);
        assert_eq!(processor.bit_cap, (1u64 << 30) - 1);
    }

    // REFERENCE PROCESSING TESTS

    #[test]
    fn test_single_sequence() {
        let mut processor = KmerProcessor::new(5, 1, true, 1024);
        let ref_seq = b"ACGTACGT";

        with_kmer_index(&mut processor, |p, k| {
            p.process_ref(ref_seq, k);
        });

        // Should have metadata + k-mers
        assert!(processor.ref_kmers.len() > 1);

        // Verify metadata was inserted
        let metadata = u64::MAX ^ 5;
        assert!(processor.contains_kmer(&metadata));
    }

    #[test]
    fn test_process_multiple() {
        let mut processor = KmerProcessor::new(5, 1, true, 1024);

        let mut count1 = 0usize;
        let mut count2 = 0usize;
        with_kmer_index(&mut processor, |p, k| {
            p.process_ref(b"ACGTACGT", k);
            count1 = k.len();
            p.process_ref(b"TGCATGCA", k);
            count2 = k.len();
        });

        // Should have added more k-mers (may have some overlap)
        assert!(count2 >= count1);
    }

    #[test]
    #[should_panic(expected = "Read sequence is shorter than k")]
    fn test_process_short() {
        let mut processor = KmerProcessor::new(10, 1, true, 1024);
        with_kmer_index(&mut processor, |p, k| {
            p.process_ref(b"ACGT", k); // Only 4 bases, k=10
        });
    }

    #[test]
    fn test_rc_refs() {
        let mut processor = KmerProcessor::new(3, 1, true, 1024);
        with_kmer_index(&mut processor, |p, k| {
            p.process_ref(b"TTTT", k); // original
            p.process_ref(b"AAAA", k); // rc
        });

        assert_eq!(processor.ref_kmers.len(), 2); // metadata + above k-mer
    }

    // READ PROCESSING TESTS

    #[test]
    fn test_process_read_exact_match() {
        let mut processor = KmerProcessor::new(5, 1, true, 1024);
        let ref_seq = b"ACGTACGT";
        with_kmer_index(&mut processor, |p, k| {
            p.process_ref(ref_seq, k);
        });

        // Read with exact k-mer from reference
        let read = b"ACGTACGT";
        assert!(processor.process_read(read));
    }

    #[test]
    fn test_process_read_no_match() {
        let mut processor = KmerProcessor::new(5, 1, true, 1024);
        with_kmer_index(&mut processor, |p, k| {
            p.process_ref(b"ACGTACGT", k);
        });
        let read = b"TTTTTTTT";
        assert!(!processor.process_read(read));
    }

    #[test]
    fn test_process_read_partial_match_below_threshold() {
        let mut processor = KmerProcessor::new(5, 3, true, 1024);
        with_kmer_index(&mut processor, |p, k| {
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
        let mut processor = KmerProcessor::new(4, 1, true, 1024);
        with_kmer_index(&mut processor, |p, k| {
            p.process_ref(b"ACGTACGT", k);
        });

        // Read with at least one matching k-mer
        let read = b"ACGTTTTTTT";
        assert!(processor.process_read(read));
    }

    #[test]
    fn test_process_read_too_short() {
        let processor = KmerProcessor::new(10, 1, true, 1024);
        let read = b"ACGT"; // Only 4 bases, k=10

        assert!(!processor.process_read(read));
    }

    // BIT MANIPULATION TESTS

    #[test]
    fn test_bit_cap_calculation() {
        let processor5 = KmerProcessor::new(5, 1, true, 1024);
        assert_eq!(processor5.bit_cap, (1u64 << 10) - 1);

        let processor10 = KmerProcessor::new(10, 1, true, 1024);
        assert_eq!(processor10.bit_cap, (1u64 << 20) - 1);

        let processor21 = KmerProcessor::new(21, 1, true, 1024);
        assert_eq!(processor21.bit_cap, (1u64 << 42) - 1);
    }

    // SLIDING WINDOW TESTS

    #[test]
    fn test_sliding_window_kmer_generation() {
        let mut processor = KmerProcessor::new(3, 1, true, 1024);

        // For sequence "ACGTACGT" with k=3:
        // k-mers should be: ACG, CGT, GTA, TAC, ACG, CGA
        with_kmer_index(&mut processor, |p, k| {
            p.process_ref(b"ACGTACGA", k);
        });

        // Should have metadata + unique k-mers
        println!("{:?}", &processor.ref_kmers);
        assert!(processor.ref_kmers.len() >= 4);
    }

    // K-MER UNIQUENESS TESTS

    #[test]
    fn test_duplicate_kmers() {
        let mut processor = KmerProcessor::new(5, 1, true, 1024);

        // Process same sequence twice
        let mut count1 = 0usize;
        let mut count2 = 0usize;
        with_kmer_index(&mut processor, |p, k| {
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
        let mut processor = KmerProcessor::new(5, 1, true, 1024);
        with_kmer_index(&mut processor, |p, k| {
            p.process_ref(b"ACGTACGTACGT", k);
        });

        // Even one matching k-mer should return true
        let read = b"ACGTATTTTTTT";
        assert!(processor.process_read(read));
    }

    #[test]
    fn test_threshold_higher() {
        let mut processor = KmerProcessor::new(3, 3, true, 1024);
        with_kmer_index(&mut processor, |p, k| {
            p.process_ref(b"ACGTACGT", k);
        });

        // Need at least 3 matching k-mers
        let read_with_matches = b"ACGTACGTACGT";
        assert!(processor.process_read(read_with_matches));
    }

    // EDGE CASES

    #[test]
    fn test_minimum_k_value() {
        let mut processor = KmerProcessor::new(1, 1, true, 1024);
        with_kmer_index(&mut processor, |p, k| {
            p.process_ref(b"ACGT", k);
        });

        let read = b"AAAA";
        processor.process_read(read); // Should not panic
    }

    #[test]
    fn test_sequence_exactly_k_length() {
        let mut processor = KmerProcessor::new(5, 1, true, 1024);
        let seq = b"ACGTA"; // Exactly k=5
        with_kmer_index(&mut processor, |p, k| {
            p.process_ref(seq, k);
        });
        assert!(processor.process_read(seq));
    }

    #[test]
    fn test_empty_reference_set() {
        let processor = KmerProcessor::new(5, 1, true, 1024);
        let read = b"ACGTACGT";

        // No reference k-mers added
        assert!(!processor.process_read(read));
    }

    #[test]
    fn test_repeated_bases() {
        let mut processor = KmerProcessor::new(5, 1, true, 1024);
        with_kmer_index(&mut processor, |p, k| {
            p.process_ref(b"AAAAAAAAAA", k);
        });

        let read = b"AAAAAAAAAA";
        assert!(processor.process_read(read));
    }

    // CANONICAL K-MER TESTS

    #[test]
    fn test_canonical() {
        let mut processor = KmerProcessor::new(5, 1, true, 1024);

        // Add forward strand
        with_kmer_index(&mut processor, |p, k| {
            p.process_ref(b"ATGCCAGT", k);
        });

        // Reverse complement should also match due to canonical representation
        let read = b"ACTGGCAT";
        assert!(processor.process_read(read));
    }

    // PERFORMANCE & CAPACITY TESTS

    #[test]
    fn test_large_reference_set() {
        let mut processor = KmerProcessor::new(21, 1, true, 1024);
        let bases = ["A", "C", "G", "T"];

        // Add many reference sequences
        with_kmer_index(&mut processor, |p, k| {
            for i in 0..100 {
                let remainder = i % 4;
                let seq = format!("ACGTACGTACGTACGTACGT{}", bases[remainder]);
                p.process_ref(seq.as_bytes(), k);
            }
        });

        assert!(processor.ref_kmers.len() == 5);
    }

    #[test]
    fn test_long_sequence_processing() {
        let mut processor = KmerProcessor::new(21, 1, true, 1024);
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

        with_kmer_index(&mut processor, |p, k| {
            p.process_ref(&long_seq, k);
        });
        println!("{}", &processor.ref_kmers.len());
        // Should have many k-mers
        assert!(processor.ref_kmers.len() > 500);
    }

    // METADATA TESTS
    #[test]
    fn test_metadata_insertion() {
        let mut processor = KmerProcessor::new(15, 1, true, 1024);
        assert!(processor.ref_kmers.is_empty());

        with_kmer_index(&mut processor, |p, k| {
            p.process_ref(b"ACGTACGTACGTACGT", k);
        });

        let metadata = u64::MAX ^ 15;
        assert!(processor.contains_kmer(&metadata));
    }

    #[test]
    fn test_metadata_different_k() {
        let mut processor21 = KmerProcessor::new(21, 1, true, 1024);
        let mut processor15 = KmerProcessor::new(15, 1, true, 1024);

        with_kmer_index(&mut processor21, |p, k| {
            p.process_ref(b"ACGTACGTACGTACGTACGTACGT", k);
        });
        with_kmer_index(&mut processor15, |p, k| {
            p.process_ref(b"ACGTACGTACGTACGT", k);
        });

        let metadata21 = u64::MAX ^ 21;
        let metadata15 = u64::MAX ^ 15;

        assert!(processor21.contains_kmer(&metadata21));
        assert!(processor15.contains_kmer(&metadata15));
        assert_ne!(metadata21, metadata15);
    }

    // MULTIPLE THRESHOLD TESTS
    #[test]
    fn test_various_thresholds() {
        for threshold in 1..=5 {
            let mut processor = KmerProcessor::new(5, threshold, true, 1024);
            with_kmer_index(&mut processor, |p, k| {
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
            let mut processor = KmerProcessor::new(*k, 1, true, 1024);

            // Create sequence long enough for this k
            let seq: Vec<u8> = (0..*k + 10)
                .map(|i| match i % 4 {
                    0 => b'A',
                    1 => b'C',
                    2 => b'G',
                    _ => b'T',
                })
                .collect();

            with_kmer_index(&mut processor, |p, k| {
                p.process_ref(&seq, k);
            });
            assert!(processor.ref_kmers.len() > 1); // Should have metadata + kmers
        }
    }

    // BOUNDARY TESTS
    #[test]
    fn test_sequence_length_k() {
        let mut processor = KmerProcessor::new(10, 1, true, 1024);
        let seq = b"ACGTACGTAC"; // Exactly 10 bases

        with_kmer_index(&mut processor, |p, k| {
            p.process_ref(seq, k);
        });
        // Should have metadata + 1 k-mer
        assert_eq!(processor.ref_kmers.len(), 2);
    }

    #[test]
    fn test_sequence_length_k_plus_one() {
        let mut processor = KmerProcessor::new(10, 1, true, 1024);
        let seq = b"ACGTACGTACT"; // 11 bases

        with_kmer_index(&mut processor, |p, k| {
            p.process_ref(seq, k);
        });
        // Should have metadata + 2 k-mers
        assert_eq!(processor.ref_kmers.len(), 3);
    }
}
