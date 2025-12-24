use assert_cmd::Command;
use predicates::prelude::*;
use std::fs::{self, File};
use std::io::Write;
use tempfile::TempDir;

// Helper function to create test FASTQ files
fn create_fastq(path: &str, records: &[(&str, &str, &str)]) -> std::io::Result<()> {
    let mut file = File::create(path)?;
    for (id, seq, qual) in records {
        writeln!(file, "@{}", id)?;
        writeln!(file, "{}", seq)?;
        writeln!(file, "+")?;
        writeln!(file, "{}", qual)?;
    }
    Ok(())
}

// Helper for FASTA files
fn create_fasta(path: &str, records: &[(&str, &str)]) -> std::io::Result<()> {
    let mut file = File::create(path)?;
    for (id, seq) in records {
        writeln!(file, ">{}", id)?;
        writeln!(file, "{}", seq)?;
    }
    Ok(())
}

#[test]
fn test_basic_filtering_unpaired() {
    let temp = TempDir::new().unwrap();
    let ref_path = temp.path().join("ref.fa");
    let reads_path = temp.path().join("reads.fq");
    let matched_path = temp.path().join("matched.fq");
    let unmatched_path = temp.path().join("unmatched.fq");

    // Create reference with known sequences
    create_fasta(
        ref_path.to_str().unwrap(),
        &[("ref1", "ACGTACGTACGTACGTACGTA")],
    )
    .unwrap();

    // Create reads: some match, some don't
    create_fastq(
        reads_path.to_str().unwrap(),
        &[
            ("read1", "ACGTACGTACGTACGTACGTA", "IIIIIIIIIIIIIIIIIIIII"), // match
            ("read2", "TTTTTTTTTTTTTTTTTTTTT", "IIIIIIIIIIIIIIIIIIIII"), // total non-match
            ("read3", "ACGTACGTACGTTTTTTTTTT", "IIIIIIIIIIIIIIIIIIIII"), // partial non-match
        ],
    )
    .unwrap();

    assert!(ref_path.exists());
    assert!(reads_path.exists());

    // Run nucleaze
    Command::from_std(std::process::Command::new(assert_cmd::cargo::cargo_bin!(
        "nucleaze"
    )))
    .arg("--in")
    .arg(reads_path.to_str().unwrap())
    .arg("--ref")
    .arg(ref_path.to_str().unwrap())
    .arg("--outm")
    .arg(matched_path.to_str().unwrap())
    .arg("--outu")
    .arg(unmatched_path.to_str().unwrap())
    .arg("--k")
    .arg("21")
    .assert()
    .success();

    // Verify output files exist
    assert!(matched_path.exists());
    assert!(unmatched_path.exists());

    // Verify content
    let matched_content = fs::read_to_string(&matched_path).unwrap();
    assert!(matched_content.contains("read1"));
    assert!(matched_content.contains("ACGTACGTACGTACGTACGTA"));

    let unmatched_content = fs::read_to_string(&unmatched_path).unwrap();
    assert!(unmatched_content.contains("read2"));
    assert!(unmatched_content.contains("read3"));
}

#[test]
fn test_filtering_ambiguous_sequences() {
    let temp = TempDir::new().unwrap();
    let ref_path = temp.path().join("ref.fa");
    let reads_path = temp.path().join("reads.fq");
    let matched_path = temp.path().join("matched.fq");
    let unmatched_path = temp.path().join("unmatched.fq");

    create_fasta(
        ref_path.to_str().unwrap(),
        &[("ref1", "ACGTAAAAAAAAAAAAAAAAA")],
    )
    .unwrap();

    create_fastq(
        reads_path.to_str().unwrap(),
        &[("read1", "ACGTNNNNNNNNNNNNNNNNN", "IIIIIIIIIIIIIIIIIIIII")],
    )
    .unwrap();

    Command::from_std(std::process::Command::new(assert_cmd::cargo::cargo_bin!(
        "nucleaze"
    )))
    .arg("--in")
    .arg(reads_path.to_str().unwrap())
    .arg("--ref")
    .arg(ref_path.to_str().unwrap())
    .arg("--outm")
    .arg(matched_path.to_str().unwrap())
    .arg("--outu")
    .arg(unmatched_path.to_str().unwrap())
    .arg("--k")
    .arg("10")
    .assert()
    .success();

    assert!(matched_path.exists());
    assert!(unmatched_path.exists());

    let matched_content = fs::read_to_string(&matched_path).unwrap();
    assert!(!matched_content.contains("read1"));

    let unmatched_content = fs::read_to_string(&unmatched_path).unwrap();
    assert!(unmatched_content.contains("read1"));
}

#[test]
fn test_interleaved() {
    let temp = TempDir::new().unwrap();
    let ref_path = temp.path().join("ref.fa");
    let reads_path = temp.path().join("reads1.fq");
    let matched_path = temp.path().join("matched1.fq");
    let unmatched_path = temp.path().join("unmatched1.fq");

    create_fasta(
        ref_path.to_str().unwrap(),
        &[("ref1", "ACGTACGTACGTACGTACGTA")],
    )
    .unwrap();

    create_fastq(
        reads_path.to_str().unwrap(),
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
    .arg(reads_path.to_str().unwrap())
    .arg("--ref")
    .arg(ref_path.to_str().unwrap())
    .arg("--outm")
    .arg(matched_path.to_str().unwrap())
    .arg("--outu")
    .arg(unmatched_path.to_str().unwrap())
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

    create_fasta(
        ref_path.to_str().unwrap(),
        &[("ref1", "ACGTACGTACGTACGTACGTA")],
    )
    .unwrap();

    create_fastq(
        reads_path.to_str().unwrap(),
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
    .arg(reads_path.to_str().unwrap())
    .arg("--ref")
    .arg(ref_path.to_str().unwrap())
    .arg("--outm")
    .arg(matched1_path.to_str().unwrap())
    .arg("--outm2")
    .arg(matched2_path.to_str().unwrap())
    .arg("--outu")
    .arg(unmatched1_path.to_str().unwrap())
    .arg("--outu2")
    .arg(unmatched2_path.to_str().unwrap())
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

    create_fasta(
        ref_path.to_str().unwrap(),
        &[("ref1", "ACGTACGTACGTACGTACGTA")],
    )
    .unwrap();

    create_fastq(
        reads_path.to_str().unwrap(),
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
    .arg(reads_path.to_str().unwrap())
    .arg("--ref")
    .arg(ref_path.to_str().unwrap())
    .arg("--outm")
    .arg(matched1_path.to_str().unwrap())
    .arg("--outm2")
    .arg(matched2_path.to_str().unwrap())
    .arg("--outu")
    .arg(unmatched1_path.to_str().unwrap())
    .arg("--outu2")
    .arg(unmatched2_path.to_str().unwrap())
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

    create_fasta(
        ref_path.to_str().unwrap(),
        &[("ref1", "ACGTACGTACGTACGTACGTA")],
    )
    .unwrap();

    create_fastq(
        reads1_path.to_str().unwrap(),
        &[
            ("read1", "ACGTACGTACGTACGTACGTA", "IIIIIIIIIIIIIIIIIIIII"),
            ("read3", "TTTTTTTTTTTTTTTTTTTTT", "IIIIIIIIIIIIIIIIIIIII"),
        ],
    )
    .unwrap();

    create_fastq(
        reads2_path.to_str().unwrap(),
        &[
            ("read2", "CCCCCCCCCCCCCCCCCCCCC", "IIIIIIIIIIIIIIIIIIIII"),
            ("read4", "GGGGGGGGGGGGGGGGGGGGG", "IIIIIIIIIIIIIIIIIIIII"),
        ],
    )
    .unwrap();

    Command::from_std(std::process::Command::new(assert_cmd::cargo::cargo_bin!(
        "nucleaze"
    )))
    .arg("--in")
    .arg(reads1_path.to_str().unwrap())
    .arg("--in2")
    .arg(reads2_path.to_str().unwrap())
    .arg("--ref")
    .arg(ref_path.to_str().unwrap())
    .arg("--outm")
    .arg(matched_path.to_str().unwrap())
    .arg("--outu")
    .arg(unmatched_path.to_str().unwrap())
    .assert()
    .success();

    assert!(matched_path.exists());
    assert!(matched_path.exists());

    let matched_content = fs::read_to_string(&matched_path).unwrap();
    assert!(matched_content.contains("read1"));
    assert!(matched_content.contains("ACGTACGTACGTACGTACGTA"));
    assert!(matched_content.contains("read2"));
    assert!(matched_content.contains("CCCCCCCCCCCCCCCCCCCCC"));

    let unmatched_content = fs::read_to_string(&unmatched_path).unwrap();
    assert!(unmatched_content.contains("read3"));
    assert!(unmatched_content.contains("TTTTTTTTTTTTTTTTTTTTT"));
    assert!(unmatched_content.contains("read4"));
    assert!(unmatched_content.contains("GGGGGGGGGGGGGGGGGGGGG"));
}

#[test]
fn test_paired_reads() {
    let temp = TempDir::new().unwrap();
    let ref_path = temp.path().join("ref.fa");
    let reads1_path = temp.path().join("reads1.fq");
    let reads2_path = temp.path().join("reads2.fq");
    let matched1_path = temp.path().join("matched1.fq");
    let matched2_path = temp.path().join("matched2.fq");
    let unmatched1_path = temp.path().join("unmatched1.fq");
    let unmatched2_path = temp.path().join("unmatched2.fq");

    create_fasta(
        ref_path.to_str().unwrap(),
        &[("ref1", "ACGTACGTACGTACGTACGTA")],
    )
    .unwrap();

    create_fastq(
        reads1_path.to_str().unwrap(),
        &[
            ("read1", "ACGTACGTACGTACGTACGTA", "IIIIIIIIIIIIIIIIIIIII"),
            ("read3", "GGGGGGGGGGGGGGGGGGGGG", "IIIIIIIIIIIIIIIIIIIII"),
        ],
    )
    .unwrap();

    create_fastq(
        reads2_path.to_str().unwrap(),
        &[
            ("read2", "GGGGGGGGGGGGGGGGGGGGG", "IIIIIIIIIIIIIIIIIIIII"),
            ("read4", "GGGGGGGGGGGGGGGGGGGGG", "IIIIIIIIIIIIIIIIIIIII"),
        ],
    )
    .unwrap();

    Command::from_std(std::process::Command::new(assert_cmd::cargo::cargo_bin!(
        "nucleaze"
    )))
    .arg("--in")
    .arg(reads1_path.to_str().unwrap())
    .arg("--in2")
    .arg(reads2_path.to_str().unwrap())
    .arg("--ref")
    .arg(ref_path.to_str().unwrap())
    .arg("--outm")
    .arg(matched1_path.to_str().unwrap())
    .arg("--outm2")
    .arg(matched2_path.to_str().unwrap())
    .arg("--outu")
    .arg(unmatched1_path.to_str().unwrap())
    .arg("--outu2")
    .arg(unmatched2_path.to_str().unwrap())
    .assert()
    .success();

    assert!(matched1_path.exists());
    assert!(matched2_path.exists());

    let matched1_content = fs::read_to_string(&matched1_path).unwrap();
    assert!(matched1_content.contains("read1"));
    assert!(matched1_content.contains("ACGTACGTACGTACGTACGTA"));

    let matched2_content = fs::read_to_string(&matched2_path).unwrap();
    assert!(matched2_content.contains("read2"));
    assert!(matched2_content.contains("GGGGGGGGGGGGGGGGGGGGG"));

    let unmatched1_content = fs::read_to_string(&unmatched1_path).unwrap();
    assert!(unmatched1_content.contains("read3"));
    assert!(unmatched1_content.contains("GGGGGGGGGGGGGGGGGGGGG"));

    let unmatched2_content = fs::read_to_string(&unmatched2_path).unwrap();
    assert!(unmatched2_content.contains("read4"));
    assert!(unmatched2_content.contains("GGGGGGGGGGGGGGGGGGGGG"));
}

#[test]
fn test_serialized_reference() {
    let temp = TempDir::new().unwrap();
    let ref_path = temp.path().join("ref.fa");
    let bin_path = temp.path().join("ref.bin");
    let reads_path = temp.path().join("reads.fq");
    let matched_path = temp.path().join("matched.fq");
    let unmatched_path = temp.path().join("unmatched.fq");

    create_fasta(
        ref_path.to_str().unwrap(),
        &[("ref1", "ACGTACGTACGTACGTACGTG")],
    )
    .unwrap();

    create_fastq(
        reads_path.to_str().unwrap(),
        &[("read1", "ACGTACGTACGTACGTACGTG", "IIIIIIIIIIIIIIIIIIIII")],
    )
    .unwrap();

    // First run: create serialized reference
    Command::from_std(std::process::Command::new(assert_cmd::cargo::cargo_bin!(
        "nucleaze"
    )))
    .arg("--in")
    .arg(reads_path.to_str().unwrap())
    .arg("--ref")
    .arg(ref_path.to_str().unwrap())
    .arg("--saveref")
    .arg(bin_path.to_str().unwrap())
    .arg("--outm")
    .arg(matched_path.to_str().unwrap())
    .arg("--outu")
    .arg(unmatched_path.to_str().unwrap())
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
    .arg(reads_path.to_str().unwrap())
    .arg("--ref")
    .arg(ref_path.to_str().unwrap())
    .arg("--binref")
    .arg(bin_path.to_str().unwrap())
    .arg("--outm")
    .arg(matched_path.to_str().unwrap())
    .arg("--outu")
    .arg(unmatched_path.to_str().unwrap())
    .assert()
    .success()
    .stdout(predicate::str::contains("Loaded"));
}

#[test]
fn test_different_k_values() {
    let temp = TempDir::new().unwrap();
    let ref_path = temp.path().join("ref.fa");
    let reads_path = temp.path().join("reads.fq");

    create_fasta(
        ref_path.to_str().unwrap(),
        &[("ref1", "ACGTACGTACGTACGTACGTC")],
    )
    .unwrap();

    create_fastq(
        reads_path.to_str().unwrap(),
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
        .arg(reads_path.to_str().unwrap())
        .arg("--ref")
        .arg(ref_path.to_str().unwrap())
        .arg("--outm")
        .arg(matched.to_str().unwrap())
        .arg("--outu")
        .arg(unmatched.to_str().unwrap())
        .arg("--k")
        .arg(k.to_string())
        .assert()
        .success();
    }
}

// #[test]
// fn test_memory_allocation() {
//     let temp = TempDir::new().unwrap();
//     let ref_path = temp.path().join("ref.fa");
//     let reads_path = temp.path().join("reads.fq");

//     // Create reference with known sequences
//     create_fasta(
//         ref_path.to_str().unwrap(),
//         &[("ref1", "ACGTACGTACGTACGTACGTA")],
//     )
//     .unwrap();

//     // Create reads: some match, some don't
//     create_fastq(
//         reads_path.to_str().unwrap(),
//         &[
//             ("read1", "ACGTACGTACGTACGTACGTA", "IIIIIIIIIIIIIIIIIIIII"), // match
//             ("read2", "TTTTTTTTTTTTTTTTTTTTT", "IIIIIIIIIIIIIIIIIIIII"), // total non-match
//             ("read3", "ACGTACGTACGTTTTTTTTTT", "IIIIIIIIIIIIIIIIIIIII"), // partial non-match
//         ],
//     )
//     .unwrap();

//     assert!(ref_path.exists());
//     assert!(reads_path.exists());

//     // Run nucleaze with no memory allocated
//     Command::cargo_bin("nucleaze")
//         .unwrap()
//         .arg("--in")
//         .arg(reads_path.to_str().unwrap())
//         .arg("--ref")
//         .arg(ref_path.to_str().unwrap())
//         .arg("--k")
//         .arg("21")
//         .arg("--maxmem")
//         .arg("1B")
//         .assert()
//         .failure();

//     // Run nucleaze with enough memory allocated
//     Command::cargo_bin("nucleaze")
//         .unwrap()
//         .arg("--in")
//         .arg(reads_path.to_str().unwrap())
//         .arg("--ref")
//         .arg(ref_path.to_str().unwrap())
//         .arg("--k")
//         .arg("21")
//         .arg("--maxmem")
//         .arg("1G")
//         .assert()
//         .success();
// }
#[test]
fn test_duplicate_read_files() {
    let temp = TempDir::new().unwrap();
    let ref_path = temp.path().join("ref.fa");
    let reads_path = temp.path().join("reads.fq");
    let reads_path_dup = temp.path().join("reads.fq");
    let matched_path = temp.path().join("matched.fq");
    let unmatched_path = temp.path().join("unmatched.fq");

    create_fasta(
        ref_path.to_str().unwrap(),
        &[("ref1", "ACGTACGTACGTACGTACGTC")],
    )
    .unwrap();

    create_fastq(
        reads_path.to_str().unwrap(),
        &[("read1", "ACGTACGTACGTACGTACGTC", "IIIIIIIIIIIIIIIIIIIII")],
    )
    .unwrap();

    create_fastq(
        reads_path_dup.to_str().unwrap(),
        &[("read1", "ACGTACGTACGTACGTACGTC", "IIIIIIIIIIIIIIIIIIIII")],
    )
    .unwrap();

    Command::from_std(std::process::Command::new(assert_cmd::cargo::cargo_bin!(
        "nucleaze"
    )))
    .arg("--in")
    .arg(reads_path.to_str().unwrap())
    .arg("--in2")
    .arg(reads_path_dup.to_str().unwrap())
    .arg("--ref")
    .arg(ref_path.to_str().unwrap())
    .arg("--outm")
    .arg(matched_path.to_str().unwrap())
    .arg("--outu")
    .arg(unmatched_path.to_str().unwrap())
    .assert()
    .failure();
}

#[test]
fn test_minhits_threshold() {
    let temp = TempDir::new().unwrap();
    let ref_path = temp.path().join("ref.fa");
    let reads_path = temp.path().join("reads.fq");
    let matched_path = temp.path().join("matched.fq");
    let unmatched_path = temp.path().join("unmatched.fq");

    create_fasta(
        ref_path.to_str().unwrap(),
        &[("ref1", "ACGTACGTACGTACGTACGTA")],
    )
    .unwrap();

    create_fastq(
        reads_path.to_str().unwrap(),
        &[("read1", "ACGTATTTTTTTTTTTTTTTT", "IIIIIIIIIIIIIIIIIIIII")],
    )
    .unwrap();

    // With minhits=1, should match
    Command::from_std(std::process::Command::new(assert_cmd::cargo::cargo_bin!(
        "nucleaze"
    )))
    .arg("--in")
    .arg(reads_path.to_str().unwrap())
    .arg("--ref")
    .arg(ref_path.to_str().unwrap())
    .arg("--outm")
    .arg(matched_path.to_str().unwrap())
    .arg("--outu")
    .arg(unmatched_path.to_str().unwrap())
    .arg("--k")
    .arg("4")
    .arg("--minhits")
    .arg("2")
    .assert()
    .success();

    let matched = fs::read_to_string(&matched_path).unwrap();
    assert!(matched.contains("read1"));
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

    create_fasta(
        ref_path.to_str().unwrap(),
        &[("ref1", "ACGTACGTACGTACGTACGTA")],
    )
    .unwrap();

    create_fastq(
        reads_path.to_str().unwrap(),
        &[("read1", "ACGTACGTACGTACGTACGTA", "IIIIIIIIIIIIIIIIIIIII")],
    )
    .unwrap();

    Command::from_std(std::process::Command::new(assert_cmd::cargo::cargo_bin!(
        "nucleaze"
    )))
    .arg("--in")
    .arg(reads_path.to_str().unwrap())
    .arg("--ref")
    .arg(ref_path.to_str().unwrap())
    .arg("--outm")
    .arg(matched_fasta_path.to_str().unwrap())
    .arg("--outu")
    .arg(unmatched_fasta_path.to_str().unwrap())
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
    .arg(reads_path.to_str().unwrap())
    .arg("--ref")
    .arg(ref_path.to_str().unwrap())
    .arg("--outm")
    .arg(matched_fastq_path.to_str().unwrap())
    .arg("--outu")
    .arg(unmatched_fastq_path.to_str().unwrap())
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

    create_fasta(
        ref_path.to_str().unwrap(),
        &[("ref1", "ACGTACGTACGTACGTACGTG")],
    )
    .unwrap();

    create_fastq(
        reads_path.to_str().unwrap(),
        &[("read1", "ACGTACGTACGTACGTACGTG", "IIIIIIIIIIIIIIIIIIIII")],
    )
    .unwrap();

    Command::from_std(std::process::Command::new(assert_cmd::cargo::cargo_bin!(
        "nucleaze"
    )))
    .arg("--in")
    .arg(reads_path.to_str().unwrap())
    .arg("--ref")
    .arg(ref_path.to_str().unwrap())
    .arg("--outm")
    .arg(matched_path.to_str().unwrap())
    .arg("--outu")
    .arg(unmatched_path.to_str().unwrap())
    .arg("--threads")
    .arg("2")
    .assert()
    .success();
}
