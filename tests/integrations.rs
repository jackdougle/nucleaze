use assert_cmd::Command;
use bincode::{config, decode_from_std_read};
use predicates::prelude::*;
use std::fs::{self, File};
use std::io::{BufReader, Write};
use std::path::Path;
use tempfile::TempDir;

fn create_fastq(path: &Path, records: &[(&str, &str, &str)]) -> std::io::Result<()> {
    let mut file = File::create(path)?;
    for (id, seq, qual) in records {
        writeln!(file, "@{}", id)?;
        writeln!(file, "{}", seq)?;
        writeln!(file, "+")?;
        writeln!(file, "{}", qual)?;
    }
    Ok(())
}

fn create_fasta(path: &Path, records: &[(&str, &str)]) -> std::io::Result<()> {
    let mut file = File::create(path)?;
    for (id, seq) in records {
        writeln!(file, ">{}", id)?;
        writeln!(file, "{}", seq)?;
    }
    Ok(())
}

fn nucleaze_cmd() -> Command {
    Command::from_std(std::process::Command::new(assert_cmd::cargo::cargo_bin!(
        "nucleaze"
    )))
}

#[test]
fn test_basic_filtering_unpaired() {
    let temp = TempDir::new().unwrap();
    let ref_path = temp.path().join("ref.fa");
    let reads_path = temp.path().join("reads.fq");
    let matched_path = temp.path().join("matched.fq");
    let unmatched_path = temp.path().join("unmatched.fq");

    // Reference: simple ACGT pattern
    create_fasta(&ref_path, &[("ref1", "ACGTACGTACGTACGTACGTA")]).unwrap();

    // Reads: 1 match, 1 total mismatch, 1 partial mismatch (suffix differs)
    create_fastq(
        &reads_path,
        &[
            ("read1", "ACGTACGTACGTACGTACGTA", "IIIIIIIIIIIIIIIIIIIII"), // Match
            ("read2", "TTTTTTTTTTTTTTTTTTTTT", "IIIIIIIIIIIIIIIIIIIII"), // Mismatch
            ("read3", "ACGTACGTACGTTTTTTTTTT", "IIIIIIIIIIIIIIIIIIIII"), // Mismatch (partial)
        ],
    )
    .unwrap();

    nucleaze_cmd()
        .arg("--in")
        .arg(&reads_path)
        .arg("--ref")
        .arg(&ref_path)
        .arg("--outm")
        .arg(&matched_path)
        .arg("--outu")
        .arg(&unmatched_path)
        .arg("--k")
        .arg("21")
        .assert()
        .success();

    let matched = fs::read_to_string(&matched_path).unwrap();
    assert!(matched.contains("read1"));
    assert!(!matched.contains("read2"));
    assert!(!matched.contains("read3"));

    let unmatched = fs::read_to_string(&unmatched_path).unwrap();
    assert!(!unmatched.contains("read1"));
    assert!(unmatched.contains("read2"));
    assert!(unmatched.contains("read3"));
}

#[test]
fn test_filtering_ambiguous_sequences() {
    let temp = TempDir::new().unwrap();
    let ref_path = temp.path().join("ref.fa");
    let reads_path = temp.path().join("reads.fq");
    let matched_path = temp.path().join("matched.fq");
    let unmatched_path = temp.path().join("unmatched.fq");

    create_fasta(&ref_path, &[("ref1", "ACGTAAAAAAAAAAAAAAAAA")]).unwrap();

    // 'N' bases should break the k-mer chain
    create_fastq(
        &reads_path,
        &[("read1", "ACGTNNNNNNNNNNNNNNNNN", "IIIIIIIIIIIIIIIIIIIII")],
    )
    .unwrap();

    nucleaze_cmd()
        .arg("--in")
        .arg(&reads_path)
        .arg("--ref")
        .arg(&ref_path)
        .arg("--outm")
        .arg(&matched_path)
        .arg("--outu")
        .arg(&unmatched_path)
        .arg("--k")
        .arg("10")
        .assert()
        .success();

    let matched = fs::read_to_string(&matched_path).unwrap();
    assert!(!matched.contains("read1"));

    let unmatched = fs::read_to_string(&unmatched_path).unwrap();
    assert!(unmatched.contains("read1"));
}

#[test]
fn test_interleaved_input_separate_output() {
    let temp = TempDir::new().unwrap();
    let ref_path = temp.path().join("ref.fa");
    let reads_path = temp.path().join("interleaved.fq");
    let matched1_path = temp.path().join("m1.fq");
    let matched2_path = temp.path().join("m2.fq");
    let unmatched1_path = temp.path().join("u1.fq");
    let unmatched2_path = temp.path().join("u2.fq");

    create_fasta(&ref_path, &[("ref1", "AAAAA")]).unwrap();

    // R1 matches, R2 mismatches -> Both should go to matched because they are paired
    // R3 mismatches, R4 mismatches -> Both to unmatched
    create_fastq(
        &reads_path,
        &[
            ("read1/1", "AAAAA", "IIIII"), // Match
            ("read1/2", "TTTTT", "IIIII"), // Mismatch, but pair matches
            ("read2/1", "CCCCC", "IIIII"), // Mismatch
            ("read2/2", "GGGGG", "IIIII"), // Mismatch
        ],
    )
    .unwrap();

    nucleaze_cmd()
        .arg("--in")
        .arg(&reads_path)
        .arg("--ref")
        .arg(&ref_path)
        .arg("--outm")
        .arg(&matched1_path)
        .arg("--outm2")
        .arg(&matched2_path)
        .arg("--outu")
        .arg(&unmatched1_path)
        .arg("--outu2")
        .arg(&unmatched2_path)
        .arg("-i") // Critical flag
        .arg("-k")
        .arg("5")
        .assert()
        .success();

    // Check Pair 1 (Matched)
    let m1 = fs::read_to_string(&matched1_path).unwrap();
    let m2 = fs::read_to_string(&matched2_path).unwrap();
    assert!(m1.contains("read1/1"));
    assert!(m2.contains("read1/2"));

    // Check Pair 2 (Unmatched)
    let u1 = fs::read_to_string(&unmatched1_path).unwrap();
    let u2 = fs::read_to_string(&unmatched2_path).unwrap();
    assert!(u1.contains("read2/1"));
    assert!(u2.contains("read2/2"));
}

#[test]
fn test_paired_input_files() {
    let temp = TempDir::new().unwrap();
    let ref_path = temp.path().join("ref.fa");
    let r1_path = temp.path().join("r1.fq");
    let r2_path = temp.path().join("r2.fq");
    let m1_path = temp.path().join("m1.fq");
    let m2_path = temp.path().join("m2.fq");
    let u1_path = temp.path().join("u1.fq");
    let u2_path = temp.path().join("u2.fq");

    create_fasta(&ref_path, &[("ref", "AAAAA")]).unwrap();

    // Pair 1: R1 matches
    create_fastq(
        &r1_path,
        &[("p1_1", "AAAAA", "IIIII"), ("p2_1", "TTTTT", "IIIII")],
    )
    .unwrap();
    // Pair 2: Neither matches
    create_fastq(
        &r2_path,
        &[("p1_2", "CCCCC", "IIIII"), ("p2_2", "GGGGG", "IIIII")],
    )
    .unwrap();

    nucleaze_cmd()
        .arg("--in")
        .arg(&r1_path)
        .arg("--in2")
        .arg(&r2_path)
        .arg("--ref")
        .arg(&ref_path)
        .arg("--outm")
        .arg(&m1_path)
        .arg("--outm2")
        .arg(&m2_path)
        .arg("--outu")
        .arg(&u1_path)
        .arg("--outu2")
        .arg(&u2_path)
        .arg("--k")
        .arg("5")
        .assert()
        .success();

    assert!(fs::read_to_string(&m1_path).unwrap().contains("p1_1"));
    assert!(fs::read_to_string(&m2_path).unwrap().contains("p1_2"));
    assert!(fs::read_to_string(&u1_path).unwrap().contains("p2_1"));
    assert!(fs::read_to_string(&u2_path).unwrap().contains("p2_2"));
}

#[test]
fn test_missing_interleaved_flag_failure() {
    // This tests the logic that detects missing --interinput flag
    // based on arg combinations
    let temp = TempDir::new().unwrap();
    let ref_path = temp.path().join("ref.fa");
    let reads_path = temp.path().join("reads.fq");

    create_fasta(&ref_path, &[("ref", "AAAAA")]).unwrap();
    create_fastq(&reads_path, &[("r1", "AAAAA", "IIIII")]).unwrap();

    nucleaze_cmd()
        .arg("--in")
        .arg(&reads_path)
        .arg("--ref")
        .arg(&ref_path)
        // Providing paired OUTPUTs for single INPUT implies interleaved,
        // but without --interinput flag, this should fail/panic.
        .arg("--outm")
        .arg("m1.fq")
        .arg("--outm2")
        .arg("m2.fq")
        .arg("--outu")
        .arg("u1.fq")
        .arg("--outu2")
        .arg("u2.fq")
        .assert()
        .failure();
}

// --- New Architecture Tests (Bloom Filter & Serialization) ---

#[test]
fn test_serialized_reference_persistence() {
    // This test ensures the new Vec+Bloom system can serialize and deserialize correctly.
    // NOTE: This test verifies that the binary format works as expected.
    let temp = TempDir::new().unwrap();
    let ref_path = temp.path().join("ref.fa");
    let bin_path = temp.path().join("ref.bin");
    let reads_path = temp.path().join("reads.fq");
    let matched_path = temp.path().join("matched.fq");
    let unmatched_path = temp.path().join("unmatched.fq");

    create_fasta(&ref_path, &[("ref", "ACGTACGT")]).unwrap();
    create_fastq(&reads_path, &[("r1", "ACGTACGT", "IIIIIIII")]).unwrap();

    // 1. Run to generate the binary reference
    nucleaze_cmd()
        .arg("--in")
        .arg(&reads_path)
        .arg("--ref")
        .arg(&ref_path)
        .arg("--saveref")
        .arg(&bin_path) // Save the index
        .arg("--outm")
        .arg(&matched_path)
        .arg("--outu")
        .arg(&unmatched_path)
        .arg("--k")
        .arg("8")
        .assert()
        .success();

    assert!(bin_path.exists(), "Binary reference file was not created");

    // Clear outputs
    fs::remove_file(&matched_path).unwrap();

    // 2. Run using the binary reference
    // This stresses the load_serialized_kmers function.
    // IF the Bloom filter is not rebuilt on load, this step will find 0 matches!
    nucleaze_cmd()
        .arg("--in")
        .arg(&reads_path)
        .arg("--binref")
        .arg(&bin_path) // Load the index
        // No --ref argument here
        .arg("--outm")
        .arg(&matched_path)
        .arg("--outu")
        .arg(&unmatched_path)
        .arg("--k")
        .arg("8")
        .assert()
        .success();

    let matched = fs::read_to_string(&matched_path).unwrap();
    assert!(
        matched.contains("r1"),
        "Failed to match read using loaded binary reference"
    );
}

// --- Parameter Validation Tests ---

#[test]
fn test_minhits_threshold() {
    let temp = TempDir::new().unwrap();
    let ref_path = temp.path().join("ref.fa");
    let reads_path = temp.path().join("reads.fq");
    let matched_path = temp.path().join("matched.fq");
    let unmatched_path = temp.path().join("unmatched.fq");

    create_fasta(&ref_path, &[("ref", "AAAAA")]).unwrap();
    // Read has 5 bases.
    // k=4.
    // kmers: AAAA, AAAA (2 kmers).
    create_fastq(&reads_path, &[("r1", "AAAAA", "IIIII")]).unwrap();

    // With minhits=2, it should match (both kmers match)
    nucleaze_cmd()
        .arg("--in")
        .arg(&reads_path)
        .arg("--ref")
        .arg(&ref_path)
        .arg("--outm")
        .arg(&matched_path)
        .arg("--outu")
        .arg(&unmatched_path)
        .arg("--k")
        .arg("4")
        .arg("--minhits")
        .arg("2")
        .assert()
        .success();

    assert!(fs::read_to_string(&matched_path).unwrap().contains("r1"));

    fs::remove_file(&matched_path).unwrap();

    // With minhits=3, it should fail (only 2 kmers total)
    nucleaze_cmd()
        .arg("--in")
        .arg(&reads_path)
        .arg("--ref")
        .arg(&ref_path)
        .arg("--outm")
        .arg(&matched_path)
        .arg("--outu")
        .arg(&unmatched_path)
        .arg("--k")
        .arg("4")
        .arg("--minhits")
        .arg("3")
        .assert()
        .success(); // Command succeeds, but read is unmatched

    assert!(!fs::read_to_string(&matched_path).unwrap().contains("r1"));
    assert!(fs::read_to_string(&unmatched_path).unwrap().contains("r1"));
}

#[test]
#[cfg(target_os = "linux")]
fn test_memory_allocation_limits() {
    let temp = TempDir::new().unwrap();
    let ref_path = temp.path().join("ref.fa");
    let reads_path = temp.path().join("reads.fq");

    create_fasta(&ref_path, &[("ref", "ACGT")]).unwrap();
    create_fastq(&reads_path, &[("r1", "ACGT", "IIII")]).unwrap();

    // Impossible memory limit
    nucleaze_cmd()
        .arg("--in")
        .arg(&reads_path)
        .arg("--ref")
        .arg(&ref_path)
        .arg("--k")
        .arg("4")
        .arg("--maxmem")
        .arg("1B")
        .assert()
        .failure();

    // Reasonable memory limit
    nucleaze_cmd()
        .arg("--in")
        .arg(&reads_path)
        .arg("--ref")
        .arg(&ref_path)
        .arg("--k")
        .arg("4")
        .arg("--maxmem")
        .arg("1G")
        .assert()
        .success();
}

#[test]
fn test_canonical_kmers() {
    let temp = TempDir::new().unwrap();
    let ref_path = temp.path().join("ref.fa");
    let reads_path = temp.path().join("reads.fq");
    let matched_path = temp.path().join("matched.fq");
    let unmatched_path = temp.path().join("unmatched.fq");

    // Ref: All T's
    create_fasta(&ref_path, &[("ref", "TTTTT")]).unwrap();
    // Read: All A's (Reverse complement of TTTTT)
    create_fastq(&reads_path, &[("r1", "AAAAA", "IIIII")]).unwrap();

    // 1. Without canonical, AAAAA != TTTTT
    nucleaze_cmd()
        .arg("--in")
        .arg(&reads_path)
        .arg("--ref")
        .arg(&ref_path)
        .arg("--outm")
        .arg(&matched_path)
        .arg("--outu")
        .arg(&unmatched_path)
        .arg("--k")
        .arg("5")
        .assert()
        .success();

    assert!(fs::read_to_string(&unmatched_path).unwrap().contains("r1"));

    // 2. With canonical, AAAAA == TTTTT (rc)
    nucleaze_cmd()
        .arg("--in")
        .arg(&reads_path)
        .arg("--ref")
        .arg(&ref_path)
        .arg("--outm")
        .arg(&matched_path)
        .arg("--outu")
        .arg(&unmatched_path)
        .arg("--k")
        .arg("5")
        .arg("--canonical")
        .assert()
        .success();

    assert!(fs::read_to_string(&matched_path).unwrap().contains("r1"));
}

#[test]
fn test_duplicate_file_args_error() {
    let temp = TempDir::new().unwrap();
    let p = temp.path().join("f.fq");
    create_fastq(&p, &[("r", "A", "I")]).unwrap();

    // Same file for input 1 and 2
    nucleaze_cmd()
        .arg("--in")
        .arg(&p)
        .arg("--in2")
        .arg(&p)
        .arg("--ref")
        .arg(&p)
        .assert()
        .failure();
}

#[test]
fn test_metadata_kmer_functional() {
    let temp = TempDir::new().unwrap();
    let ref_path = temp.path().join("ref.fa");
    let reads_path = temp.path().join("reads.fq");
    let matched_path = temp.path().join("matched.fq");
    let unmatched_path = temp.path().join("unmatched.fq");
    let saveref_path = temp.path().join("index.bin");

    create_fasta(&ref_path, &[("ref1", "ACGTACGTACGT")]).unwrap();
    create_fastq(&reads_path, &[("r1", "ACGTACGTACGT", "IIIIIIIIIIII")]).unwrap();

    let k: i32 = 5;

    nucleaze_cmd()
        .arg("--in")
        .arg(&reads_path)
        .arg("--ref")
        .arg(&ref_path)
        .arg("--outm")
        .arg(&matched_path)
        .arg("--outu")
        .arg(&unmatched_path)
        .arg("--k")
        .arg(&k.to_string())
        .arg("--saveref")
        .arg(&saveref_path)
        .assert()
        .success();

    let f = File::open(&saveref_path).unwrap();
    let mut reader = BufReader::new(f);

    let ref_kmers: Vec<Vec<u64>> = decode_from_std_read(&mut reader, config::standard().with_fixed_int_encoding()).unwrap();

    let size_metadata = u64::MAX ^ k as u64;
    assert!(ref_kmers[0].contains(&size_metadata));

    nucleaze_cmd()
        .arg("--in")
        .arg(&reads_path)
        .arg("--ref")
        .arg(&ref_path)
        .arg("--outm")
        .arg(&matched_path)
        .arg("--outu")
        .arg(&unmatched_path)
        .arg("--k")
        .arg("6".to_string()) // metadata is set for k = 5
        .arg("--binref")
        .arg(&saveref_path)
        .assert()
        .failure(); 
}

#[test]
fn test_ref_with_short_seqs_errors() {
    let temp = TempDir::new().unwrap();
    let ref_path = temp.path().join("ref.fa");
    let reads_path = temp.path().join("reads.fq");

    // Reference sequences shorter than k (k=10)
    create_fasta(&ref_path, &[("short1", "ACGT")]).unwrap();
    create_fastq(&reads_path, &[("r1", "ACGT", "IIII")]).unwrap();

    nucleaze_cmd()
        .arg("--in")
        .arg(&reads_path)
        .arg("--ref")
        .arg(&ref_path)
        .arg("--k")
        .arg("10")
        .assert()
        .failure();
}

#[test]
fn test_interleaved() {
    let temp = TempDir::new().unwrap();
    let ref_path = temp.path().join("ref.fa");
    let reads_path = temp.path().join("reads1.fq");
    let matched_path = temp.path().join("matched1.fq");
    let unmatched_path = temp.path().join("unmatched1.fq");

    create_fasta(&ref_path, &[("ref1", "ACGTACGTACGTACGTACGTA")]).unwrap();

    create_fastq(
        &reads_path,
        &[
            ("read1", "ACGTACGTACGTACGTACGTA", "IIIIIIIIIIIIIIIIIIIII"),
            ("read2", "TTTTTTTTTTTTTTTTTTTTT", "IIIIIIIIIIIIIIIIIIIII"),
            ("read3", "TTTTTTTTTTTTTTTTTTTTT", "IIIIIIIIIIIIIIIIIIIII"),
            ("read4", "TTTTTTTTTTTTTTTTTTTTT", "IIIIIIIIIIIIIIIIIIIII"),
        ],
    )
    .unwrap();

    Command::from_std(std::process::Command::new(assert_cmd::cargo::cargo_bin!(
        "nucleaze"
    )))
    .arg("--in")
    .arg(&reads_path)
    .arg("--ref")
    .arg(&ref_path)
    .arg("--outm")
    .arg(&matched_path)
    .arg("--outu")
    .arg(&unmatched_path)
    .arg("--interinput")
    .assert()
    .success();

    assert!(matched_path.exists());

    let matched_content = fs::read_to_string(&matched_path).unwrap();
    assert!(matched_content.contains("read1"));
    assert!(matched_content.contains("ACGTACGTACGTACGTACGTA"));
    assert!(matched_content.contains("read2"));
    assert!(matched_content.contains("TTTTTTTTTTTTTTTTTTTTT"));

    let unmatched_content = fs::read_to_string(&unmatched_path).unwrap();
    assert!(unmatched_content.contains("read3"));
    assert!(unmatched_content.contains("TTTTTTTTTTTTTTTTTTTTT"));
    assert!(unmatched_content.contains("read4"));
    assert!(unmatched_content.contains("TTTTTTTTTTTTTTTTTTTTT"));
}

#[test]
fn test_inter_in_paired_out() {
    let temp = TempDir::new().unwrap();
    let ref_path = temp.path().join("ref.fa");
    let reads_path = temp.path().join("reads1.fq");
    let matched1_path = temp.path().join("matched1.fq");
    let matched2_path = temp.path().join("matched2.fq");
    let unmatched1_path = temp.path().join("unmatched1.fq");
    let unmatched2_path = temp.path().join("unmatched2.fq");

    create_fasta(&ref_path, &[("ref1", "ACGTACGTACGTACGTACGTA")]).unwrap();

    create_fastq(
        &reads_path,
        &[
            ("read1", "ACGTACGTACGTACGTACGTA", "IIIIIIIIIIIIIIIIIIIII"),
            ("read2", "TTTTTTTTTTTTTTTTTTTTT", "IIIIIIIIIIIIIIIIIIIII"),
            ("read3", "TTTTTTTTTTTTTTTTTTTTT", "IIIIIIIIIIIIIIIIIIIII"),
            ("read4", "TTTTTTTTTTTTTTTTTTTTT", "IIIIIIIIIIIIIIIIIIIII"),
        ],
    )
    .unwrap();

    Command::from_std(std::process::Command::new(assert_cmd::cargo::cargo_bin!(
        "nucleaze"
    )))
    .arg("--in")
    .arg(&reads_path)
    .arg("--ref")
    .arg(&ref_path)
    .arg("--outm")
    .arg(&matched1_path)
    .arg("--outm2")
    .arg(&matched2_path)
    .arg("--outu")
    .arg(&unmatched1_path)
    .arg("--outu2")
    .arg(&unmatched2_path)
    .arg("--interinput")
    .assert()
    .success();

    assert!(matched1_path.exists());
    assert!(matched2_path.exists());

    let matched1_content = fs::read_to_string(&matched1_path).unwrap();
    assert!(matched1_content.contains("read1"));
    assert!(matched1_content.contains("ACGTACGTACGTACGTACGTA"));

    let matched2_content = fs::read_to_string(&matched2_path).unwrap();
    assert!(matched2_content.contains("read2"));
    assert!(matched2_content.contains("TTTTTTTTTTTTTTTTTTTTT"));

    let unmatched1_content = fs::read_to_string(&unmatched1_path).unwrap();
    assert!(unmatched1_content.contains("read3"));
    assert!(unmatched1_content.contains("TTTTTTTTTTTTTTTTTTTTT"));

    let unmatched2_content = fs::read_to_string(&unmatched2_path).unwrap();
    assert!(unmatched2_content.contains("read4"));
    assert!(unmatched2_content.contains("TTTTTTTTTTTTTTTTTTTTT"));
}

#[test]
fn test_inter_in_paired_out_no_flag() {
    let temp = TempDir::new().unwrap();
    let ref_path = temp.path().join("ref.fa");
    let reads_path = temp.path().join("reads1.fq");
    let matched1_path = temp.path().join("matched1.fq");
    let matched2_path = temp.path().join("matched2.fq");
    let unmatched1_path = temp.path().join("unmatched1.fq");
    let unmatched2_path = temp.path().join("unmatched2.fq");

    create_fasta(&ref_path, &[("ref1", "ACGTACGTACGTACGTACGTA")]).unwrap();

    create_fastq(
        &reads_path,
        &[
            ("read1", "ACGTACGTACGTACGTACGTA", "IIIIIIIIIIIIIIIIIIIII"),
            ("read2", "TTTTTTTTTTTTTTTTTTTTT", "IIIIIIIIIIIIIIIIIIIII"),
            ("read3", "TTTTTTTTTTTTTTTTTTTTT", "IIIIIIIIIIIIIIIIIIIII"),
            ("read4", "TTTTTTTTTTTTTTTTTTTTT", "IIIIIIIIIIIIIIIIIIIII"),
        ],
    )
    .unwrap();

    Command::from_std(std::process::Command::new(assert_cmd::cargo::cargo_bin!(
        "nucleaze"
    )))
    .arg("--in")
    .arg(&reads_path)
    .arg("--ref")
    .arg(&ref_path)
    .arg("--outm")
    .arg(&matched1_path)
    .arg("--outm2")
    .arg(&matched2_path)
    .arg("--outu")
    .arg(&unmatched1_path)
    .arg("--outu2")
    .arg(&unmatched2_path)
    .assert()
    .failure();
}

#[test]
fn test_paired_in_inter_out() {
    let temp = TempDir::new().unwrap();
    let ref_path = temp.path().join("ref.fa");
    let reads1_path = temp.path().join("reads1.fq");
    let reads2_path = temp.path().join("reads2.fq");
    let matched_path = temp.path().join("matched1.fq");
    let unmatched_path = temp.path().join("unmatched1.fq");

    create_fasta(&ref_path, &[("ref1", "ACGTACGTACGTACGTACGTA")]).unwrap();

    create_fastq(
        &reads1_path,
        &[
            ("read1/1", "ACGTACGTACGTACGTACGTA", "IIIIIIIIIIIIIIIIIIIII"),
            ("read2/1", "TTTTTTTTTTTTTTTTTTTTT", "IIIIIIIIIIIIIIIIIIIII"),
            ("read3/1", "ACGTACGTACGTTTTTTTTTT", "IIIIIIIIIIIIIIIIIIIII"),
        ],
    )
    .unwrap();

    create_fastq(
        &reads2_path,
        &[
            ("read1/2", "CCCCCCCCCCCCCCCCCCCCC", "IIIIIIIIIIIIIIIIIIIII"),
            ("read2/2", "AAAAAAAAAAAAAAAAAAAAA", "IIIIIIIIIIIIIIIIIIIII"),
            ("read3/2", "TTTTTTTTTTTTACGTACGTA", "IIIIIIIIIIIIIIIIIIIII"),
        ],
    )
    .unwrap();

    Command::from_std(std::process::Command::new(assert_cmd::cargo::cargo_bin!(
        "nucleaze"
    )))
    .arg("--in")
    .arg(&reads1_path)
    .arg("--in2")
    .arg(&reads2_path)
    .arg("--ref")
    .arg(&ref_path)
    .arg("--outm")
    .arg(&matched_path)
    .arg("--outu")
    .arg(&unmatched_path)
    .assert()
    .success();

    assert!(matched_path.exists());
    assert!(matched_path.exists());

    let matched_content = fs::read_to_string(&matched_path).unwrap();
    assert!(matched_content.contains("read1/1"));
    assert!(matched_content.contains("ACGTACGTACGTACGTACGTA"));
    assert!(matched_content.contains("read1/2"));
    assert!(matched_content.contains("CCCCCCCCCCCCCCCCCCCCC"));

    let unmatched_content = fs::read_to_string(&unmatched_path).unwrap();
    assert!(unmatched_content.contains("read2/1"));
    assert!(unmatched_content.contains("TTTTTTTTTTTTTTTTTTTTT"));
    assert!(unmatched_content.contains("read2/2"));
    assert!(unmatched_content.contains("AAAAAAAAAAAAAAAAAAAAA"));
    assert!(unmatched_content.contains("read3/1"));
    assert!(unmatched_content.contains("ACGTACGTACGTTTTTTTTTT"));
    assert!(unmatched_content.contains("read3/2"));
    assert!(unmatched_content.contains("TTTTTTTTTTTTACGTACGTA"));
}

#[test]
fn test_serialized_reference() {
    let temp = TempDir::new().unwrap();
    let ref_path = temp.path().join("ref.fa");
    let bin_path = temp.path().join("ref.bin");
    let reads_path = temp.path().join("reads.fq");
    let matched_path = temp.path().join("matched.fq");
    let unmatched_path = temp.path().join("unmatched.fq");

    create_fasta(&ref_path, &[("ref1", "ACGTACGTACGTACGTACGTG")]).unwrap();

    create_fastq(
        &reads_path,
        &[("read1", "ACGTACGTACGTACGTACGTG", "IIIIIIIIIIIIIIIIIIIII")],
    )
    .unwrap();

    // First run: create serialized reference
    Command::from_std(std::process::Command::new(assert_cmd::cargo::cargo_bin!(
        "nucleaze"
    )))
    .arg("--in")
    .arg(&reads_path)
    .arg("--ref")
    .arg(&ref_path)
    .arg("--saveref")
    .arg(&bin_path)
    .arg("--outm")
    .arg(&matched_path)
    .arg("--outu")
    .arg(&unmatched_path)
    .assert()
    .success();

    // Verify binary file was created
    assert!(bin_path.exists());

    // Second run: use serialized reference
    fs::remove_file(&matched_path).unwrap();
    fs::remove_file(&unmatched_path).unwrap();

    Command::from_std(std::process::Command::new(assert_cmd::cargo::cargo_bin!(
        "nucleaze"
    )))
    .arg("--in")
    .arg(&reads_path)
    .arg("--ref")
    .arg(&ref_path)
    .arg("--binref")
    .arg(&bin_path)
    .arg("--outm")
    .arg(&matched_path)
    .arg("--outu")
    .arg(&unmatched_path)
    .assert()
    .success()
    .stdout(predicate::str::contains("Loaded"));
}

#[test]
fn test_different_k_values() {
    let temp = TempDir::new().unwrap();
    let ref_path = temp.path().join("ref.fa");
    let reads_path = temp.path().join("reads.fq");

    create_fasta(&ref_path, &[("ref1", "ACGTACGTACGTACGTACGTC")]).unwrap();

    create_fastq(
        &reads_path,
        &[("read1", "ACGTACGTACGTACGTACGTC", "IIIIIIIIIIIIIIIIIIIII")],
    )
    .unwrap();

    for k in [5, 11, 15, 21].iter() {
        let matched = temp.path().join(format!("matched_k{}.fq", k));
        let unmatched = temp.path().join(format!("unmatched_k{}.fq", k));

        Command::from_std(std::process::Command::new(assert_cmd::cargo::cargo_bin!(
            "nucleaze"
        )))
        .arg("--in")
        .arg(&reads_path)
        .arg("--ref")
        .arg(&ref_path)
        .arg("--outm")
        .arg(&matched)
        .arg("--outu")
        .arg(&unmatched)
        .arg("--k")
        .arg(k.to_string())
        .assert()
        .success();
    }
}

#[test]
fn test_duplicate_read_files() {
    let temp = TempDir::new().unwrap();
    let ref_path = temp.path().join("ref.fa");
    let reads_path = temp.path().join("reads.fq");
    let reads_path_dup = temp.path().join("reads.fq");
    let matched_path = temp.path().join("matched.fq");
    let unmatched_path = temp.path().join("unmatched.fq");

    create_fasta(&ref_path, &[("ref1", "ACGTACGTACGTACGTACGTC")]).unwrap();

    create_fastq(
        &reads_path,
        &[("read1", "ACGTACGTACGTACGTACGTC", "IIIIIIIIIIIIIIIIIIIII")],
    )
    .unwrap();

    create_fastq(
        &reads_path_dup,
        &[("read1", "ACGTACGTACGTACGTACGTC", "IIIIIIIIIIIIIIIIIIIII")],
    )
    .unwrap();

    Command::from_std(std::process::Command::new(assert_cmd::cargo::cargo_bin!(
        "nucleaze"
    )))
    .arg("--in")
    .arg(&reads_path)
    .arg("--in2")
    .arg(&reads_path_dup)
    .arg("--ref")
    .arg(&ref_path)
    .arg("--outm")
    .arg(&matched_path)
    .arg("--outu")
    .arg(&unmatched_path)
    .assert()
    .failure();
}

#[test]
fn test_fasta_output_format() {
    let temp = TempDir::new().unwrap();
    let ref_path = temp.path().join("ref.fa");
    let reads_path = temp.path().join("reads.fq");
    let matched_fasta_path = temp.path().join("matched.fa");
    let matched_fastq_path = temp.path().join("matched.fq");
    let unmatched_fasta_path = temp.path().join("unmatched.fa");
    let unmatched_fastq_path = temp.path().join("unmatched.fq");

    create_fasta(&ref_path, &[("ref1", "ACGTACGTACGTACGTACGTA")]).unwrap();

    create_fastq(
        &reads_path,
        &[("read1", "ACGTACGTACGTACGTACGTA", "IIIIIIIIIIIIIIIIIIIII")],
    )
    .unwrap();

    Command::from_std(std::process::Command::new(assert_cmd::cargo::cargo_bin!(
        "nucleaze"
    )))
    .arg("--in")
    .arg(&reads_path)
    .arg("--ref")
    .arg(&ref_path)
    .arg("--outm")
    .arg(&matched_fasta_path)
    .arg("--outu")
    .arg(&unmatched_fasta_path)
    .assert()
    .success();

    let matched = fs::read_to_string(&matched_fasta_path).unwrap();
    assert!(matched.starts_with(">"));
    assert!(!matched.contains("@"));
    assert!(!matched.contains("+"));
    assert!(!matched.contains("IIIIIIIIIIIIIIIIIIIII"));

    Command::from_std(std::process::Command::new(assert_cmd::cargo::cargo_bin!(
        "nucleaze"
    )))
    .arg("--in")
    .arg(&reads_path)
    .arg("--ref")
    .arg(&ref_path)
    .arg("--outm")
    .arg(&matched_fastq_path)
    .arg("--outu")
    .arg(&unmatched_fastq_path)
    .assert()
    .success();

    let matched = fs::read_to_string(&matched_fastq_path).unwrap();
    assert!(matched.starts_with("@"));
    assert!(matched.contains("IIIIIIIIIIIIIIIIIIIII"));
    assert!(!matched.contains(">"));
}

#[test]
fn test_missing_required_args() {
    Command::from_std(std::process::Command::new(assert_cmd::cargo::cargo_bin!(
        "nucleaze"
    )))
    .assert()
    .failure();
}

#[test]
fn test_threads_argument() {
    let temp = TempDir::new().unwrap();
    let ref_path = temp.path().join("ref.fa");
    let reads_path = temp.path().join("reads.fq");
    let matched_path = temp.path().join("matched.fq");
    let unmatched_path = temp.path().join("unmatched.fq");

    create_fasta(&ref_path, &[("ref1", "ACGTACGTACGTACGTACGTG")]).unwrap();

    create_fastq(
        &reads_path,
        &[("read1", "ACGTACGTACGTACGTACGTG", "IIIIIIIIIIIIIIIIIIIII")],
    )
    .unwrap();

    Command::from_std(std::process::Command::new(assert_cmd::cargo::cargo_bin!(
        "nucleaze"
    )))
    .arg("--in")
    .arg(&reads_path)
    .arg("--ref")
    .arg(&ref_path)
    .arg("--outm")
    .arg(&matched_path)
    .arg("--outu")
    .arg(&unmatched_path)
    .arg("--threads")
    .arg("2")
    .assert()
    .success();
}
