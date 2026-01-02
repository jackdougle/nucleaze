mod core;
mod kmer_ops;

use clap::Parser;
use rlimit::{Resource, setrlimit};
use std::io;
use std::time::Instant;

const ABOUT: &str = "Nucleaze 1.3.0
Written by Jack Douglass
Last modified January 2nd, 2026

Nucleaze compares DNA sequences from input file to DNA sequences from reference
 file using k-mer analysis. Splits up reference file sequences into k-mers of
 specified length to build k-mer index, then compares each input read sequence
 for matching k-mers. If a read sequence has >= minhits matching k-mers, it
 will be printed as a match. Very memory-efficient and performant. Processes
 paired reads in two files or as a single interleaved file.

USAGE: nucleaze --in <reads file> --ref <ref file> ...

INPUT PARAMETERS
    --in <file>         Input FASTA/FASTQ file containing reads to be filtered.
                        Use 'stdin.fq' or 'stdin' to pipe from stdin.
    --in2 <file>        Second input file for 2nd pair of reads. 
                        Must be same length as main input file.
    --ref <file>        (-r) Reference FASTA/FASTQ file containing sequences to
                        build reference k-mer index. Program will serialize
                        reference k-mers and build to --binref path for future
                        use. Not necessary if '--binref <file>' is provided.
    --saveref <file>    (-s) Path at which to store serialized k-mer index if
                        no valid serialized k-mer index is provided.
    --binref <file>     (-b) Binary file containing serialized k-mer index,
                        increases performance. Nucleaze makes this
                        automatically based on '--ref' file if '--saveref'
                        <file> is provided.

OUTPUT PARAMETERS: use 'stdout.fa / stdout.fq to pipe to stdout'
    --outm <file>       Output file for reads that have >= minhits k-mer 
                        matches to reference k-mers. FASTA format if .fa, .fna,
                        or .fasta, FASTQ format otherwise.
    --outu <file>       Output file for reads that have < minhits k-mer matches
                        to reference k-mers. FASTA format if .fa, .fna, or
                        .fasta, FASTQ format otherwise.
    --outm2 <file>      Output file for 2nd pair of matched reads.
    --outu2 <file>      Output file for 2nd pair of unmatched reads.

MEMORY & PERFORMANCE PARAMETERS
    --k 21              (-k) K-mer size (number of bases per k-mer). Ranges
                        from 1-31, larger k-mers will have less matches.
    --minhits 1         Minimum number of k-mer matches a read sequence must 
                        have to be considered a match.
    --threads auto      (-t) Number of threads to use for parallel processing.
                        Program will use all available threads by default.
    --maxmem auto       (-m) Maximum memory to use, in human readable format. 
                        '--maxmem 5G' will specify 5 gigabytes, '--maxmem 200M'
                        will specify 200 megabytes. No memory limit by default.
    --interinput        (-i) Enable flag to input as interleaved paired-end
                        reads, omit flag for unpaired reads.
    --order             (-o) Enable flag to get read outputs ordered by
                        sequence ID.
    --canonical         (-c) K-mers are stored and compared in canonical form 
                        (lowest of forward and reverse-complement).
    --bloomsize 16M     (-f) Memory size of Bloom filter in human-readable
                        format. '-l 64M' delegates 64 MB of memory to the Bloom
                        filter. 16 Mb by default.

Function and usage documentation at /README.md.
Contact jack.gdouglass@gmail.com for any questions or issues encountered.
";

#[derive(Parser)]
#[command(version, override_help = ABOUT)]
struct Args {
    /// Amount of bases in a k-mer
    #[arg(short, long)]
    k: Option<usize>,

    /// Min number of k-mer hits to match a read
    #[arg(long)]
    minhits: Option<u8>,

    /// Max amount of threads to use
    #[arg(short, long)]
    threads: Option<usize>,

    /// Memory cap in human-readable format
    #[arg(short, long)]
    maxmem: Option<String>,

    /// FASTA/FASTQ path for reference sequences
    #[arg(short, long)]
    r#ref: Option<String>,

    /// FASTA/FASTQ path for read sequences
    #[arg(long)]
    r#in: String,

    /// FASTA/FASTQ path for 2nd pair of reads
    #[arg(long)]
    in2: Option<String>,

    /// Path to store serialized reference index
    #[arg(short, long)]
    saveref: Option<String>,

    /// Binary file containing serialized ref k-mers
    #[arg(short, long)]
    binref: Option<String>,

    /// Output file of matched reads
    #[arg(long)]
    outm: Option<String>,

    /// Output file of unmatched reads
    #[arg(long)]
    outu: Option<String>,

    /// Output file for second pair of matched reads
    #[arg(long)]
    outm2: Option<String>,

    /// Output file for second pair of unmatched reads
    #[arg(long)]
    outu2: Option<String>,

    /// Enabling flag signals interleaved input
    #[arg(short, long)]
    interinput: bool,

    /// Enabling flag causes ordered output
    #[arg(short, long)]
    order: bool,

    /// Enabling flag uses canonical k-mers for comparison
    #[arg(short, long)]
    canonical: bool,

    /// Size of Bloom filter in human-readable format
    #[arg(short = 'f', long)]
    bloomsize: Option<String>,
}

fn main() -> io::Result<()> {
    let start_time = Instant::now();

    let args = Args::parse();

    let version = env!("CARGO_PKG_VERSION");
    let user_args: Vec<String> = std::env::args().skip(1).collect();
    println!("Nucleaze {} [{}]", version, user_args.join(" "));

    if let Some(maxmem_str) = &args.maxmem {
        match parse_memory_size(maxmem_str) {
            Ok(size_bytes) => {
                let soft_limit = size_bytes;
                let hard_limit = size_bytes;
                println!("Setting memory limit to {} bytes", size_bytes);
                if let Err(e) = setrlimit(Resource::AS, soft_limit, hard_limit) {
                    eprintln!("Failed to set memory limit: {}", e);
                    std::process::exit(1);
                }
            }
            Err(e) => {
                eprintln!("Could not parse --maxmem value: {}, exiting.", e);
                std::process::exit(1);
            }
        }
    }

    validate_args(&args)?;

    core::run(args, start_time)?;

    Ok(())
}

fn parse_memory_size(input: &str) -> Result<u64, String> {
    let s = input.trim().to_uppercase();

    // Split into (number, suffix)
    let (num_str, multiplier) = match s.strip_suffix('G') {
        Some(n) => (n, 1024_u64.pow(3)),
        None => match s.strip_suffix('M') {
            Some(n) => (n, 1024_u64.pow(2)),
            None => match s.strip_suffix('K') {
                Some(n) => (n, 1024_u64),
                None => match s.strip_suffix('B') {
                    Some(n) => (n, 1),
                    None => (s.as_str(), 1), // no unit
                },
            },
        },
    };

    num_str
        .parse::<u64>()
        .map(|n| n * multiplier)
        .map_err(|_| format!("Invalid size format: '{}'", input))
}

/// Validate command-line arguments to catch common errors early
fn validate_args(args: &Args) -> io::Result<()> {
    // Prevent using same file for both inputs
    if let Some(ref in2) = args.in2 {
        if args.r#in == *in2 {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                "both inputs (--in and --in2) cannot be the same file",
            ));
        }
    }

    // Prevent using same file for matched and unmatched outputs
    if let (Some(outm), Some(outu)) = (&args.outm, &args.outu) {
        if outm == outu {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                "matches (--outm) and non-matches (--outu) cannot have the same output path",
            ));
        }
    }

    // Prevent using same file for paired matched and unmatched outputs
    if let (Some(outm2), Some(outu2)) = (&args.outm2, &args.outu2) {
        if outm2 == outu2 {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                "Matches (--outm2) and non-matches (--outu2) cannot have the same output path",
            ));
        }
    }

    Ok(())
}
