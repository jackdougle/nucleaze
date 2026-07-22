//! I/O for reference indexing and read processing operations using k-mers
use nucleaze::kmer_ops::{HashShards, KmerProcessor, KmerStore, MiniBloom, encode, per_kmer_fpr};

use std::cmp::Ordering;
use std::collections::BinaryHeap;
use std::error::Error;
use std::fmt::{Display, Formatter, Result as FormatResult};
use std::fs::{File, metadata, remove_file};
use std::io::{BufReader, BufWriter, Result as IOResult, Write, stdin};
use std::mem;
use std::process::exit;
use std::str::from_utf8_unchecked;
use std::sync::{
    Arc, Mutex,
    atomic::{AtomicU32, AtomicU64, Ordering as AtomicOrdering},
    mpsc::{Receiver, SyncSender, sync_channel},
};
use std::thread;
use std::time::Instant;

use bincode::{config, decode_from_std_read, encode_into_std_write};
use crossbeam::channel::{Sender, bounded};
use needletail::{parse_fastx_file, parse_fastx_reader};
use rayon::prelude::*;
use rustc_hash::FxHashSet;

/// Source of input reads - either a file path or stdin
pub enum InputSource {
    File(String),
    Stdin,
}

impl Display for InputSource {
    fn fmt(&self, f: &mut Formatter<'_>) -> FormatResult {
        match self {
            InputSource::File(path) => write!(f, "{}", path),
            InputSource::Stdin => write!(f, "stdin"),
        }
    }
}

/// Limits how many read arenas may exist before their results are written.
/// `acquire` waits when no arena slot is available.
struct ArenaSlotPool {
    slot_sender: Sender<()>,
    slot_receiver: crossbeam::channel::Receiver<()>,
}

impl ArenaSlotPool {
    fn new(capacity: usize) -> Self {
        let (slot_sender, slot_receiver) = bounded::<()>(capacity);
        for _ in 0..capacity {
            slot_sender.send(()).expect("initialize arena slot pool"); // every slot starts free
        }
        ArenaSlotPool {
            slot_sender,
            slot_receiver,
        }
    }

    fn acquire(&self) -> ArenaSlotGuard {
        self.slot_receiver
            .recv()
            .expect("arena slot pool disconnected");
        ArenaSlotGuard {
            slot_sender: self.slot_sender.clone(),
        }
    }
}

/// Stored in each chunk to keep one arena slot unavailable.
/// Its `Drop` implementation returns the slot when the chunk is dropped.
struct ArenaSlotGuard {
    slot_sender: Sender<()>,
}

impl Drop for ArenaSlotGuard {
    fn drop(&mut self) {
        // `acquire` removed one entry, so the bounded channel has room to return it.
        // Ignore errors after the slot receiver is dropped during shutdown.
        let _ = self.slot_sender.send(());
    }
}

/// A read arena, its offsets, and its k-mer match results.
struct SequenceChunk {
    id: u32,
    data_arena: Vec<u8>, // contiguous bytes for every read in this chunk
    offsets: Vec<(u32, u32, u32, u32, u32, u32)>, // (id_start, id_len, seq_start, seq_len, qual_start, qual_len)
    matches: Vec<bool>,                           // k-mer match results
    _arena_slot: ArenaSlotGuard,                  // returns the slot when dropped
}

// Compare only IDs; the arena slot does not affect output order.
impl PartialEq for SequenceChunk {
    fn eq(&self, other: &Self) -> bool {
        self.id == other.id
    }
}

impl Eq for SequenceChunk {}

// Reverse ID comparison so `BinaryHeap::pop` returns the lowest ID first.
impl Ord for SequenceChunk {
    fn cmp(&self, other: &Self) -> Ordering {
        other.id.cmp(&self.id) // min-heap for sequential processing
    }
}

impl PartialOrd for SequenceChunk {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

/// Processing reads against a reference k-mer index
pub fn run(args: crate::Args, start_time: Instant) -> IOResult<()> {
    let memory_limited = args.maxmem.is_some();
    let available_threads = num_cpus::get();
    let num_threads = args
        .threads
        .unwrap_or(available_threads)
        .min(available_threads);

    rayon::ThreadPoolBuilder::new()
        .num_threads(num_threads)
        .stack_size(12 * 1024 * 1024)
        .build_global()
        .expect("Could not build Rayon Pool with specified thread amount");

    let k = args.k.unwrap_or(21);
    let min_hits = args.minhits.unwrap_or(1);

    let use_canonical = args.canonical;
    let fpr = args.fpr.unwrap_or_default();

    let ref_path = args.r#ref.unwrap_or_default();
    let bin_kmers_path = args.binref.unwrap_or_default();
    let save_kmers_path = args.saveref.unwrap_or_default();

    // Copy settings needed after input and output paths are moved.
    let ordered_output = args.order;
    let interleaved_input = args.interinput;
    let input = match args.r#in.as_deref() {
        None | Some("-") => InputSource::Stdin,
        Some(path) => InputSource::File(path.to_string()),
    };
    let in2 = args.in2;
    let outm = args.outm;
    let outu = args.outu;
    let outm2 = args.outm2;
    let outu2 = args.outu2;

    let ref_size = if !ref_path.is_empty() {
        let meta = metadata(&ref_path)?;
        meta.len()
    } else if !bin_kmers_path.is_empty() {
        let meta = metadata(&bin_kmers_path)?;
        meta.len() / 8 // serialized files are ~8x larger
    } else {
        0
    };

    if fpr > 0.0 {
        // Load a saved Bloom index, or build a new one if loading fails.
        let kmer_fpr = per_kmer_fpr(fpr as f64, k);

        if !bin_kmers_path.is_empty() {
            let nkb_path = format!("{}.nkb", bin_kmers_path);
            match MiniBloom::load(&nkb_path) {
                Ok((bloom, loaded_k, kmer_count)) => {
                    if loaded_k as usize != k {
                        eprintln!(
                            "Error: bloom file has k={} but --k {} was specified",
                            loaded_k, k
                        );
                        exit(1);
                    }
                    println!(
                        "Loaded ~{} k-mer(s) from {} (bloom: {} bits, {} hashes)",
                        kmer_count, nkb_path, bloom.num_bits, bloom.num_hashes,
                    );
                    let processor = Arc::new(KmerProcessor::new(k, min_hits, use_canonical, bloom));
                    let read_fn: Arc<dyn Fn(&[u8]) -> bool + Send + Sync> = {
                        let p = processor.clone();
                        Arc::new(move |seq: &[u8]| p.process_read_bloom(seq))
                    };

                    return run_inner(
                        read_fn,
                        input,
                        in2,
                        outm,
                        outu,
                        outm2,
                        outu2,
                        ordered_output,
                        interleaved_input,
                        start_time,
                        num_threads,
                    );
                }
                Err(e) => {
                    eprintln!("Could not load bloom index from {}: {}", nkb_path, e);
                    eprintln!("Building from scratch.");
                }
            }
        }

        let store = MiniBloom::new(ref_size as usize, kmer_fpr);
        let processor = KmerProcessor::new(k, min_hits, use_canonical, store);

        match get_reference_kmers_bloom(&ref_path, &processor, num_threads) {
            Ok(kmer_count) => {
                let (avg, min_f, max_f) = processor.ref_kmers.fill_stats();
                eprintln!(
                    "Bloom fill: avg={:.1}%, min={:.1}%, max={:.1}%",
                    avg * 100.0,
                    min_f * 100.0,
                    max_f * 100.0,
                );
                println!("Added {} k-mer(s) from {}", kmer_count, ref_path);
            }
            Err(e) => {
                eprintln!("Error loading reference sequences: {}", e);
                exit(1);
            }
        }

        if !save_kmers_path.is_empty() {
            let nkb_path = format!("{}.nkb", save_kmers_path);
            match processor
                .ref_kmers
                .save(&nkb_path, k as u64, processor.num_kmers() as u64)
            {
                Ok(()) => println!("Saved bloom filter to {}", nkb_path),
                Err(e) => eprintln!("Could not save bloom filter: {}", e),
            }
        }

        let processor = Arc::new(processor);
        let read_fn: Arc<dyn Fn(&[u8]) -> bool + Send + Sync> = {
            let p = processor.clone();
            Arc::new(move |seq: &[u8]| p.process_read_bloom(seq))
        };

        run_inner(
            read_fn,
            input,
            in2,
            outm,
            outu,
            outm2,
            outu2,
            ordered_output,
            interleaved_input,
            start_time,
            num_threads,
        )
    } else {
        // Exact hash mode supports loading and saving an index.
        let mut store = HashShards::new(0, 0.0);

        // Try loading a saved hash index before building a new one.
        let loaded = try_deserialize_hash(&bin_kmers_path, &mut store, k);
        let mut processor = KmerProcessor::new(k, min_hits, use_canonical, store);

        if loaded {
            println!(
                "Loaded {} k-mer(s) from {}",
                processor.num_kmers() - 1,
                bin_kmers_path
            );
        } else {
            if !bin_kmers_path.is_empty() {
                eprintln!("Invalid serialized reference file, building from scratch.");
            }
            match get_reference_kmers(&ref_path, &mut processor, num_threads, !memory_limited) {
                Ok(()) => println!(
                    "Added {} k-mer(s) from {}",
                    processor.num_kmers() - 1,
                    ref_path
                ),
                Err(e) => {
                    eprintln!("Error loading reference sequences: {}", e);
                    exit(1);
                }
            }
            if !save_kmers_path.is_empty() {
                match serialize_hash_kmers(&save_kmers_path, &processor) {
                    Ok(()) => println!("Saved serialized k-mers to {}", save_kmers_path),
                    Err(e) => eprintln!("Could not serialize reference k-mers: {}", e),
                }
            } else {
                println!("K-mer index not serialized");
            }
        }

        let processor = Arc::new(processor);
        let read_fn: Arc<dyn Fn(&[u8]) -> bool + Send + Sync> = {
            let p = processor.clone();
            Arc::new(move |seq: &[u8]| p.process_read(seq))
        };

        run_inner(
            read_fn,
            input,
            in2,
            outm,
            outu,
            outm2,
            outu2,
            ordered_output,
            interleaved_input,
            start_time,
            num_threads,
        )
    }
}

fn run_inner(
    read_fn: Arc<dyn Fn(&[u8]) -> bool + Send + Sync>,
    input: InputSource,
    in2: Option<String>,
    outm: Option<String>,
    outu: Option<String>,
    outm2: Option<String>,
    outu2: Option<String>,
    ordered_output: bool,
    interleaved_input: bool,
    start_time: Instant,
    num_threads: usize,
) -> IOResult<()> {
    let indexing_time = start_time.elapsed().as_secs_f32();
    println!("Indexing time:\t\t{:.3} seconds\n", indexing_time);

    let process_mode = detect_mode(&in2, &outm2, &outu2, interleaved_input);

    let in2_path = in2.unwrap_or_default();

    let outm_path = outm.unwrap_or(String::from("/dev/null"));
    let outu_path = outu.unwrap_or(String::from("/dev/null"));

    let outm2_path = outm2.unwrap_or(String::from("/dev/null"));
    let outu2_path = outu2.unwrap_or(String::from("/dev/null"));

    println!(
        "Using {} threads to process reads from {}",
        num_threads, input
    );

    match process_reads(
        input,
        in2_path,
        read_fn,
        &outm_path,
        &outu_path,
        &outm2_path,
        &outu2_path,
        process_mode,
        ordered_output,
    ) {
        Ok((mseq_count, mbase_count, useq_count, ubase_count)) => {
            let read_count = mseq_count + useq_count;
            let matched_percent = if read_count > 0 {
                (mseq_count as f32 / read_count as f32) * 100.0
            } else {
                0.0
            };
            let unmatched_percent = if read_count > 0 {
                (useq_count as f32 / read_count as f32) * 100.0
            } else {
                0.0
            };

            let base_count = mbase_count + ubase_count;
            let mbase_percent = if base_count > 0 {
                (mbase_count as f32 / base_count as f32) * 100.0
            } else {
                0.0
            };
            let ubase_percent = if base_count > 0 {
                (ubase_count as f32 / base_count as f32) * 100.0
            } else {
                0.0
            };

            let end_time = start_time.elapsed().as_secs_f32();
            println!("Processing time:\t{:.3} seconds", end_time - indexing_time);

            println!(
                "\nInput:\t\t\t{} reads         \t{} bases",
                read_count, base_count
            );
            println!(
                "Matches:\t\t{} reads ({:.2}%) \t{} bases ({:.2}%)",
                mseq_count, matched_percent, mbase_count, mbase_percent
            );
            println!(
                "Nonmatches:\t\t{} reads ({:.2}%) \t{} bases ({:.2}%)\n",
                useq_count, unmatched_percent, ubase_count, ubase_percent
            );
        }
        Err(e) => {
            eprintln!("\nError processing read sequences: {}", e);
        }
    }

    Ok(())
}

/// Attempt to load a pre-built hash index. Returns true on success.
fn try_deserialize_hash(path: &str, store: &mut HashShards, k: usize) -> bool {
    if path.is_empty() {
        return false;
    }
    let Ok(file) = File::open(path) else {
        return false;
    };
    let mut reader = BufReader::new(file);
    let Ok(data): Result<Vec<Vec<u64>>, _> =
        decode_from_std_read(&mut reader, config::standard().with_fixed_int_encoding())
    else {
        return false;
    };

    store.from_serialized(data);

    // Check that the saved index contains the metadata k-mer for this k value.
    let sentinel = u64::MAX ^ k as u64;
    let sentinel_present = store.contains(&sentinel);
    if !sentinel_present {
        store.clear();
    }
    sentinel_present
}

/// Build k-mer index from reference FASTA/FASTQ file
fn get_reference_kmers<S: KmerStore>(
    ref_path: &str,
    processor: &mut KmerProcessor<S>,
    num_threads: usize,
    coarse_batches: bool,
) -> Result<(), Box<dyn Error>> {
    let ref_meta = metadata(&ref_path)?;
    if ref_meta.is_file() && ref_meta.len() == 0 {
        return Err("reference file is empty".into());
    }

    const SUBINDEX_COUNT: usize = 1024;
    let merged_idx: Arc<Vec<Mutex<FxHashSet<u64>>>> = Arc::new(
        (0..SUBINDEX_COUNT)
            .map(|_| Mutex::new(FxHashSet::default()))
            .collect(),
    );

    let (sender, receiver) = bounded::<Vec<u8>>(16);

    spawn_reference_reader(
        ref_path,
        sender,
        processor.k,
        num_threads,
        true,
        coarse_batches,
    )
    .expect("k-mer extraction failed");

    (0..num_threads).into_par_iter().for_each(|_| {
        let mut local_idx = vec![Vec::with_capacity(64); SUBINDEX_COUNT];

        while let Ok(seq) = receiver.recv() {
            processor.process_ref(&seq, &mut local_idx);

            for (i, local_subidx) in local_idx.iter_mut().enumerate() {
                if !local_subidx.is_empty() {
                    let mut merged_subidx = merged_idx[i].lock().unwrap();
                    const LARGE_BATCH: usize = 4096;

                    if local_subidx.len() <= LARGE_BATCH {
                        // Reserve small batches while holding the shard lock to avoid
                        // repeated hash-table growth.
                        merged_subidx.reserve(local_subidx.len());
                        for &kmer in local_subidx.iter() {
                            merged_subidx.insert(kmer);
                        }
                    } else {
                        // Sample large batches to estimate unique k-mers before
                        // reserving space, reducing excess allocation for repeats.
                        const SAMPLE_SIZE: usize = 64;
                        let sample_len = local_subidx.len().min(SAMPLE_SIZE);
                        let first_kmer = local_subidx[0];

                        if local_subidx[..sample_len]
                            .iter()
                            .all(|&kmer| kmer == first_kmer)
                        {
                            for &kmer in local_subidx.iter() {
                                merged_subidx.insert(kmer);
                            }
                            local_subidx.clear();
                            continue;
                        }

                        let len_before_sample = merged_subidx.len();

                        for &kmer in &local_subidx[..sample_len] {
                            merged_subidx.insert(kmer);
                        }

                        let sample_uniques = merged_subidx.len() - len_before_sample;
                        let remaining = local_subidx.len() - sample_len;
                        if sample_uniques * 2 >= sample_len {
                            let estimated_uniques = remaining
                                .saturating_mul(sample_uniques)
                                .div_ceil(sample_len);
                            merged_subidx.reserve(estimated_uniques);
                        }

                        for &kmer in &local_subidx[sample_len..] {
                            merged_subidx.insert(kmer);
                        }
                    }

                    local_subidx.clear();
                }
            }
        }
    });

    let shards = Arc::try_unwrap(merged_idx)
        .expect("merged_idx still has multiple owners")
        .into_iter()
        .map(|m| m.into_inner().unwrap())
        .collect();
    processor.ref_kmers.absorb_shards(shards);

    if processor.num_kmers() == 0 {
        return Err(format!("reference file(s) contained no usable k-mers").into());
    }

    processor.insert_kmer(&(u64::MAX ^ processor.k as u64)); // insert metadata

    Ok(())
}

/// Build a Bloom index with parallel atomic inserts into the shared filter.
fn get_reference_kmers_bloom(
    ref_path: &str,
    processor: &KmerProcessor<MiniBloom>,
    num_threads: usize,
) -> Result<u64, Box<dyn Error>> {
    let ref_meta = metadata(ref_path)?;
    if ref_meta.is_file() && ref_meta.len() == 0 {
        return Err("reference file is empty".into());
    }

    let kmer_counter = Arc::new(AtomicU64::new(0));
    let (sender, receiver) = bounded::<Vec<u8>>(16);
    spawn_reference_reader(ref_path, sender, processor.k, num_threads, false, false)?;

    // Rayon workers update shared atomic Bloom blocks without a merge step.
    let counter_ref = kmer_counter.clone();
    (0..num_threads).into_par_iter().for_each(|_| {
        let mut count = 0u64;
        while let Ok(seq) = receiver.recv() {
            count += processor.insert_kmers_bloom(&seq);
        }
        counter_ref.fetch_add(count, AtomicOrdering::Relaxed);
    });

    let total_kmers = kmer_counter.load(AtomicOrdering::Relaxed);
    if total_kmers == 0 {
        return Err("reference file(s) contained no usable k-mers".into());
    }

    // The Bloom file header includes k, so no metadata k-mer is needed.

    Ok(total_kmers)
}

// A dedicated parser thread divides long records into chunks for worker threads.
// Chunks overlap by k-1 bases so k-mers spanning a boundary are included.
fn spawn_reference_reader(
    path: &str,
    sender: Sender<Vec<u8>>,
    k: usize,
    num_threads: usize,
    exact_mode: bool,
    coarse_exact_batches: bool,
) -> Result<(), Box<dyn Error>> {
    let path = path.to_string();

    thread::spawn(move || {
        let mut reader = parse_fastx_file(&path).expect("FASTA open failed");
        let overlap = k.saturating_sub(1);

        while let Some(record) = reader.next() {
            let record = match record {
                Ok(r) => r,
                Err(_) => return,
            };
            let seq = record.seq();
            let mut start = 0;
            let low_complexity = exact_mode && reference_sample_is_low_complexity(&seq, k);

            if low_complexity {
                if let Some(period) = reference_repeat_period(&seq) {
                    // For period p, the first p + k - 1 bases contain every distinct
                    // k-mer. Additional periods produce only duplicates.
                    let compact_len = (period + k.saturating_sub(1)).min(seq.len());
                    if sender.send(seq[..compact_len].to_vec()).is_err() {
                        return;
                    }
                    continue;
                }
            }
            let chunk_bases = if num_threads == 1 {
                seq.len().max(1)
            } else if coarse_exact_batches && !low_complexity {
                // Records up to 8 MiB use one chunk. Longer records use about one
                // chunk per worker, or 16 MiB chunks when that split exceeds 8 MiB.
                if seq.len() <= 8 << 20 {
                    seq.len().max(1)
                } else {
                    let per_worker = seq.len().div_ceil(num_threads);
                    if per_worker > 8 << 20 {
                        16 << 20
                    } else {
                        per_worker
                    }
                }
            } else {
                // Atomic Bloom inserts need no hash-set merge. Smaller chunks give
                // Rayon more tasks to distribute among worker threads.
                (seq.len() / num_threads.saturating_mul(2).max(1)).clamp(1 << 18, 1 << 20)
            };

            loop {
                let end = (start + chunk_bases).min(seq.len());
                if sender.send(seq[start..end].to_vec()).is_err() || end == seq.len() {
                    break;
                }
                start = end - overlap;
            }
        }
    });

    Ok(())
}

/// Sample k-mers across a record to detect low-complexity repeats.
/// The result selects a chunk size without a second full scan.
fn reference_sample_is_low_complexity(seq: &[u8], k: usize) -> bool {
    const MAX_SAMPLES: usize = 1024;
    const MIN_VALID_SAMPLES: usize = 64;

    if k == 0 || seq.len() < k {
        return false;
    }

    let prefix_len = seq.len().min(4096);
    if seq[..prefix_len].iter().all(|&base| base == seq[0]) {
        return true;
    }

    let windows = seq.len() - k + 1;
    let sample_count = windows.min(MAX_SAMPLES);
    let mut valid_samples = 0usize;
    let mut distinct = FxHashSet::default();

    for sample in 0..sample_count {
        let start = if sample_count == 1 {
            0
        } else {
            sample * (windows - 1) / (sample_count - 1)
        };
        let mut kmer = 0u64;
        let mut valid = true;

        for &base in &seq[start..start + k] {
            let Some(bits) = encode(base) else {
                valid = false;
                break;
            };
            kmer = (kmer << 2) | bits;
        }

        if valid {
            valid_samples += 1;
            distinct.insert(kmer);
        }
    }

    valid_samples >= MIN_VALID_SAMPLES && distinct.len() * 4 < valid_samples
}

/// Find an exact repeat period of at most 64 bases, test 409 bases first.
fn reference_repeat_period(seq: &[u8]) -> Option<usize> {
    const MAX_PERIOD: usize = 64;
    const PREFIX_BASES: usize = 4096;

    if seq.len() < 2 {
        return None;
    }

    let prefix_len = seq.len().min(PREFIX_BASES);
    let max_period = MAX_PERIOD.min(prefix_len / 2);

    'periods: for period in 1..=max_period {
        for i in 0..prefix_len {
            let Some(base) = encode(seq[i]) else {
                continue 'periods;
            };
            let Some(expected) = encode(seq[i % period]) else {
                continue 'periods;
            };
            if base != expected {
                continue 'periods;
            }
        }

        for i in prefix_len..seq.len() {
            let Some(base) = encode(seq[i]) else {
                continue 'periods;
            };
            let Some(expected) = encode(seq[i % period]) else {
                continue 'periods;
            };
            if base != expected {
                continue 'periods;
            }
        }

        return Some(period);
    }

    None
}

/// Save hash index to binary file.
fn serialize_hash_kmers(
    path: &str,
    processor: &KmerProcessor<HashShards>,
) -> Result<(), Box<dyn Error>> {
    let file = File::create(path)?;
    let mut writer = BufWriter::new(file);
    encode_into_std_write(
        &processor.serialize_kmers(),
        &mut writer,
        config::standard().with_fixed_int_encoding(),
    )?;
    Ok(())
}

#[derive(PartialEq, Clone, Copy, Default)]
enum ProcessMode {
    #[default]
    Unpaired, // Single-end reads
    Paired,           // Paired-end in two files, output to two files
    PairedInInterOut, // Paired-end in two files, interleaved output
    InterInPairedOut, // Interleaved input, paired-end output
    Interleaved,      // Interleaved input and output
}

/// Determine input/output mode based on file arguments
fn detect_mode(
    reads2_path: &Option<String>,
    matched2_path: &Option<String>,
    unmatched2_path: &Option<String>,
    interleaved_input: bool,
) -> ProcessMode {
    if reads2_path.is_some() {
        assert!(
            !interleaved_input,
            "Please disable the --interinput flag if providing 2 input files"
        );
        if matched2_path.is_none() && unmatched2_path.is_none() {
            println!(
                "Forcing interleaved output because paired input was specified for single output files"
            );
            ProcessMode::PairedInInterOut
        } else {
            assert!(
                matched2_path.is_some(),
                "Please add a second matched output path using: --outm2 <file>"
            );
            assert!(
                unmatched2_path.is_some(),
                "Please add a second unmatched output path using: --outu2 <file>"
            );
            println!("Input and output is processed as paired");
            ProcessMode::Paired
        }
    } else if interleaved_input {
        if matched2_path.is_none() && unmatched2_path.is_none() {
            println!("Input and output is processed as interleaved");
            ProcessMode::Interleaved
        } else {
            println!("Processing interleaved input and paired output");
            ProcessMode::InterInPairedOut
        }
    } else if matched2_path.is_some() || unmatched2_path.is_some() {
        panic!("Please enable the --interinput flag for 1 input file with paired output files");
    } else {
        println!("Input and output are processed as unpaired");
        ProcessMode::Unpaired
    }
}

/// Process reads from input file(s), filter by k-mer matches, and write to output file(s)
fn process_reads(
    input: InputSource,
    reads2_path: String,
    read_fn: Arc<dyn Fn(&[u8]) -> bool + Send + Sync>,
    matched_path: &str,
    unmatched_path: &str,
    matched2_path: &str,
    unmatched2_path: &str,
    process_mode: ProcessMode,
    ordered_output: bool,
) -> Result<(u64, u64, u64, u64), Box<dyn Error + Send + Sync>> {
    if let InputSource::File(ref path) = input {
        let meta = metadata(path)?;
        if meta.is_file() && meta.len() == 0 {
            return Err("reads file is empty".into());
        }
    }

    // Use a bounded sync channel to send processed chunks to the writing thread.
    let (chunk_sender, chunk_receiver): (SyncSender<SequenceChunk>, Receiver<SequenceChunk>) =
        sync_channel(20);

    // Limit arenas retained by workers, the output queue, and ordered output.
    // Capacity is one arena per Rayon worker plus 24 additional arenas.
    let arena_slot_pool = ArenaSlotPool::new(rayon::current_num_threads() + 24);

    let chunk_pos = Arc::new(AtomicU32::new(0));

    // Get each output format from its file name.
    let matched_filetype = matched_path.rsplit('.').next().unwrap_or("").to_string();
    let unmatched_filetype = unmatched_path.rsplit('.').next().unwrap_or("").to_string();
    let matched2_filetype = matched2_path.rsplit('.').next().unwrap_or("").to_string();
    let unmatched2_filetype = unmatched2_path.rsplit('.').next().unwrap_or("").to_string();

    // Check which output paths mean standard output.
    let matched_stdout = matched_path == "stdout" || matched_path.starts_with("stdout.");
    let unmatched_stdout = unmatched_path == "stdout" || unmatched_path.starts_with("stdout.");
    let matched2_stdout = matched2_path == "stdout" || matched2_path.starts_with("stdout.");
    let unmatched2_stdout = unmatched2_path == "stdout" || unmatched2_path.starts_with("stdout.");

    let parallel_sender = chunk_sender.clone();
    let parellel_chunk_pos = chunk_pos.clone();

    // Parser thread fills arenas and sends them for k-mer matching.
    let worker_thread = thread::spawn(move || -> Result<(), Box<dyn Error + Send + Sync>> {
        let mut reader = match input {
            InputSource::Stdin => parse_fastx_reader(stdin()),
            InputSource::File(ref path) => parse_fastx_file(path),
        }?;

        const CHUNK_SIZE: usize = 10_000;
        const ARENA_CAPACITY: usize = CHUNK_SIZE * 500;

        // Build each chunk in one byte arena and store offsets into that arena.
        let mut arena: Vec<u8> = Vec::with_capacity(ARENA_CAPACITY);
        let mut offsets: Vec<(u32, u32, u32, u32, u32, u32)> = Vec::with_capacity(CHUNK_SIZE);

        // Submit a full arena for processing. Store its slot guard in the result until the
        // chunk is written and dropped, including time in the ordered-output buffer.
        let process_arena = |local_arena: Vec<u8>,
                             local_offsets: Vec<(u32, u32, u32, u32, u32, u32)>,
                             arena_slot: ArenaSlotGuard| {
            let read_fn = read_fn.clone();
            let sender = parallel_sender.clone();
            let current_chunk_pos = parellel_chunk_pos.fetch_add(1, AtomicOrdering::SeqCst);

            rayon::spawn(move || {
                // Each Rayon task matches one arena sequentially, while Rayon runs
                // multiple arena tasks in parallel. This avoids nested parallelism.
                let matches: Vec<bool> = local_offsets
                    .iter()
                    .map(|(_, _, seq_start, seq_len, _, _)| {
                        let seq =
                            &local_arena[*seq_start as usize..(*seq_start + *seq_len) as usize];
                        read_fn(seq)
                    })
                    .collect();

                let chunk = SequenceChunk {
                    id: current_chunk_pos,
                    data_arena: local_arena,
                    offsets: local_offsets,
                    matches,
                    _arena_slot: arena_slot,
                };

                // `SyncSender::send` blocks when the 20-chunk output queue is full.
                let _ = sender.send(chunk);
            });
        };

        // Append a read's ID, sequence, and quality to the arena and record offsets.
        let push_record =
            |id: &[u8],
             seq: &[u8],
             qual: &[u8],
             arena: &mut Vec<u8>,
             offsets: &mut Vec<(u32, u32, u32, u32, u32, u32)>| {
                let id_start = arena.len() as u32;
                arena.extend_from_slice(id);
                let id_len = (arena.len() as u32) - id_start;

                let seq_start = arena.len() as u32;
                arena.extend_from_slice(seq);
                let seq_len = (arena.len() as u32) - seq_start;

                let qual_start = arena.len() as u32;
                arena.extend_from_slice(qual);
                let qual_len = (arena.len() as u32) - qual_start;

                offsets.push((id_start, id_len, seq_start, seq_len, qual_start, qual_len));
            };

        // Call `acquire` before filling an arena. It waits while the maximum number
        // of arenas are being processed, queued, or written.
        let mut arena_slot = arena_slot_pool.acquire();

        // Read sequences in the selected single or paired mode.
        if process_mode == ProcessMode::Unpaired
            || process_mode == ProcessMode::Interleaved
            || process_mode == ProcessMode::InterInPairedOut
        {
            while let Some(record) = reader.next() {
                let record = match record {
                    Ok(r) => r,
                    Err(_) => return Err("no usable reads found".into()),
                };
                push_record(
                    record.id(),
                    &record.seq(),
                    record.qual().unwrap_or(b""),
                    &mut arena,
                    &mut offsets,
                );

                if offsets.len() == CHUNK_SIZE {
                    // Move the full arena to Rayon and replace it with an empty arena
                    // so the parser can continue reading.
                    let local_arena = mem::replace(&mut arena, Vec::with_capacity(ARENA_CAPACITY));
                    let local_offsets = mem::replace(&mut offsets, Vec::with_capacity(CHUNK_SIZE));
                    process_arena(local_arena, local_offsets, arena_slot);
                    arena_slot = arena_slot_pool.acquire(); // reserve the next arena slot
                }
            }
        } else {
            // Read both files together in paired modes.
            let mut reader2 = parse_fastx_file(&reads2_path)?;

            loop {
                let rec1 = reader.next();
                let rec2 = reader2.next();

                match (rec1, rec2) {
                    (Some(r1), Some(r2)) => {
                        let record1 = r1?;
                        let record2 = r2?;

                        push_record(
                            record1.id(),
                            &record1.seq(),
                            record1.qual().unwrap_or(b""),
                            &mut arena,
                            &mut offsets,
                        );
                        push_record(
                            record2.id(),
                            &record2.seq(),
                            record2.qual().unwrap_or(b""),
                            &mut arena,
                            &mut offsets,
                        );

                        if offsets.len() == CHUNK_SIZE {
                            let local_arena =
                                mem::replace(&mut arena, Vec::with_capacity(ARENA_CAPACITY));
                            let local_offsets =
                                mem::replace(&mut offsets, Vec::with_capacity(CHUNK_SIZE));
                            process_arena(local_arena, local_offsets, arena_slot);
                            arena_slot = arena_slot_pool.acquire(); // reserve the next arena slot
                        }
                    }
                    (Some(_), None) => {
                        eprintln!(
                            "Warning: --in has more reads than --in2. \
                             Processing stopped at the end of the shorter file."
                        );
                        break;
                    }
                    (None, Some(_)) => {
                        eprintln!(
                            "Warning: --in2 has more reads than --in. \
                             Processing stopped at the end of the shorter file."
                        );
                        break;
                    }
                    (None, None) => break,
                }
            }
        }

        // Send the final partially filled arena. If it is empty, dropping
        // `arena_slot` returns the unused slot to the pool.
        if !offsets.is_empty() {
            process_arena(arena, offsets, arena_slot);
        }

        Ok(())
    });

    // Drop the original sender; the receiver closes after all Rayon clones finish.
    drop(chunk_sender);

    let mut matched_writer: BufWriter<File> = BufWriter::new(File::create(matched_path)?);
    let mut unmatched_writer: BufWriter<File> =
        BufWriter::with_capacity(4_000_000, File::create(unmatched_path)?);

    // Open second output files when writing paired reads separately.
    let mut m2_writer =
        if process_mode == ProcessMode::InterInPairedOut || process_mode == ProcessMode::Paired {
            Some(BufWriter::new(File::create(matched2_path)?))
        } else {
            None
        };
    let mut u2_writer =
        if process_mode == ProcessMode::InterInPairedOut || process_mode == ProcessMode::Paired {
            Some(BufWriter::new(File::create(unmatched2_path)?))
        } else {
            None
        };

    let mseq_count = Arc::new(AtomicU64::new(0));
    let mbase_count = Arc::new(AtomicU64::new(0));
    let useq_count = Arc::new(AtomicU64::new(0));
    let ubase_count = Arc::new(AtomicU64::new(0));

    // Create ID, sequence, and quality slices from the contiguous arena without
    // copying their bytes.
    let mut chunk_output = |chunk: &SequenceChunk| -> Result<(), Box<dyn Send + Sync + Error>> {
        let arena = &chunk.data_arena;
        let offsets = &chunk.offsets;
        let matches = &chunk.matches;

        // Use stored offsets to borrow this read's ID, sequence, and quality slices.
        let get_read = |pos: usize| {
            let (id_s, id_l, seq_s, seq_l, qual_s, qual_l) = offsets[pos];
            (
                &arena[id_s as usize..(id_s + id_l) as usize],
                &arena[seq_s as usize..(seq_s + seq_l) as usize],
                &arena[qual_s as usize..(qual_s + qual_l) as usize],
            )
        };

        if process_mode == ProcessMode::Unpaired {
            for i in 0..offsets.len() {
                let (id, seq, qual) = get_read(i);
                if matches[i] {
                    write_read(
                        &mut matched_writer,
                        id,
                        seq,
                        qual,
                        &matched_filetype,
                        matched_stdout,
                    )?;
                    mseq_count.fetch_add(1, AtomicOrdering::Relaxed);
                    mbase_count.fetch_add(seq.len() as u64, AtomicOrdering::Relaxed);
                } else {
                    write_read(
                        &mut unmatched_writer,
                        id,
                        seq,
                        qual,
                        &unmatched_filetype,
                        unmatched_stdout,
                    )?;
                    useq_count.fetch_add(1, AtomicOrdering::Relaxed);
                    ubase_count.fetch_add(seq.len() as u64, AtomicOrdering::Relaxed);
                }
            }
        } else {
            // Handle two reads at a time. If either matches, write both to the
            // matched output.
            let num_reads = offsets.len();
            let reads_to_process = if num_reads % 2 != 0 {
                // Skip a final read that has no partner.
                eprintln!(
                    "Warning: Odd number of reads ({}) in interleaved pairs mode. \
                     The last unpaired read will be skipped.",
                    num_reads
                );
                num_reads - 1
            } else {
                num_reads
            };
            for i in (0..reads_to_process).step_by(2) {
                let has_match = matches[i] || matches[i + 1];
                let (id1, seq1, qual1) = get_read(i);
                let (id2, seq2, qual2) = get_read(i + 1);

                let (w1, w2, count, bases) = if has_match {
                    (
                        &mut matched_writer,
                        m2_writer.as_mut(),
                        &mseq_count,
                        &mbase_count,
                    )
                } else {
                    (
                        &mut unmatched_writer,
                        u2_writer.as_mut(),
                        &useq_count,
                        &ubase_count,
                    )
                };

                write_read(
                    w1,
                    id1,
                    seq1,
                    qual1,
                    if has_match {
                        &matched_filetype
                    } else {
                        &unmatched_filetype
                    },
                    if has_match {
                        matched_stdout
                    } else {
                        unmatched_stdout
                    },
                )?;

                if let Some(w2_real) = w2 {
                    write_read(
                        w2_real,
                        id2,
                        seq2,
                        qual2,
                        if has_match {
                            &matched2_filetype
                        } else {
                            &unmatched2_filetype
                        },
                        if has_match {
                            matched2_stdout
                        } else {
                            unmatched2_stdout
                        },
                    )?;
                } else {
                    write_read(
                        w1,
                        id2,
                        seq2,
                        qual2,
                        if has_match {
                            &matched_filetype
                        } else {
                            &unmatched_filetype
                        },
                        if has_match {
                            matched_stdout
                        } else {
                            unmatched_stdout
                        },
                    )?;
                }

                count.fetch_add(2, AtomicOrdering::Relaxed);
                bases.fetch_add((seq1.len() + seq2.len()) as u64, AtomicOrdering::Relaxed);
            }
        }
        Ok(())
    };

    if ordered_output {
        // Store out-of-order chunks in a `BinaryHeap`. Reversed `Ord` makes `pop`
        // return the chunk with the lowest ID.
        let mut out_of_order_buffer: BinaryHeap<SequenceChunk> = BinaryHeap::new();
        let mut next_chunk_id = 0;
        const MAX_BUFFERED_CHUNKS: usize = 1000;

        for chunk in chunk_receiver {
            if chunk.id == next_chunk_id {
                // This is the next chunk to write.
                chunk_output(&chunk)?;
                next_chunk_id += 1;

                // Write any following chunks that are already waiting.
                while let Some(buffered) = out_of_order_buffer.peek() {
                    if buffered.id == next_chunk_id {
                        let buffered = out_of_order_buffer.pop().unwrap();
                        chunk_output(&buffered)?;
                        next_chunk_id += 1;
                    } else {
                        break;
                    }
                }
            } else {
                // Save this chunk until earlier chunks arrive.
                out_of_order_buffer.push(chunk);
                if out_of_order_buffer.len() > MAX_BUFFERED_CHUNKS {
                    return Err(Box::from("Too many out-of-order chunks buffered."));
                }
            }
        }
        // Write any chunks that are still waiting.
        while let Some(buffered) = out_of_order_buffer.pop() {
            if buffered.id == next_chunk_id {
                chunk_output(&buffered)?;
                next_chunk_id += 1;
            } else {
                return Err(Box::from("Missing chunk in ordered output stream."));
            }
        }
    } else {
        // Write each chunk as soon as it is ready.
        for chunk in chunk_receiver {
            chunk_output(&chunk)?;
        }
    }

    matched_writer.flush()?;
    unmatched_writer.flush()?;
    if let Some(mut w) = m2_writer {
        w.flush()?;
    }
    if let Some(mut w) = u2_writer {
        w.flush()?;
    }

    worker_thread.join().unwrap()?;

    // Clean up temporary stdout files
    for path in [
        &matched_path,
        unmatched_path,
        matched2_path,
        unmatched2_path,
    ] {
        if path.starts_with("stdout") {
            let _ = remove_file(path);
        }
    }

    Ok((
        mseq_count.load(AtomicOrdering::Relaxed),
        mbase_count.load(AtomicOrdering::Relaxed),
        useq_count.load(AtomicOrdering::Relaxed),
        ubase_count.load(AtomicOrdering::Relaxed),
    ))
}

/// Write a single read to file in FASTA or FASTQ format
fn write_read(
    writer: &mut BufWriter<File>,
    id: &[u8],
    sequence: &[u8],
    quality: &[u8],
    format: &str,
    stdout: bool,
) -> Result<(), Box<dyn Send + Sync + Error>> {
    if stdout {
        unsafe {
            let id = from_utf8_unchecked(id);
            let seq = from_utf8_unchecked(sequence);

            if format == "fa" || format == "fna" || format == "fasta" {
                println!(">{}\n{}", id, seq);
            } else {
                let qual = from_utf8_unchecked(quality);
                println!("@{}\n{}\n+\n{}", id, seq, qual);
            }
        }
    } else {
        if ["fa", "fna", "fasta"].contains(&format) {
            writer.write_all(b">")?;
            writer.write_all(id)?;
            writer.write_all(b"\n")?;
            writer.write_all(sequence)?;
            writer.write_all(b"\n")?;
        } else {
            writer.write_all(b"@")?;
            writer.write_all(id)?;
            writer.write_all(b"\n")?;
            writer.write_all(sequence)?;
            writer.write_all(b"\n+\n")?;
            writer.write_all(quality)?;
            writer.write_all(b"\n")?;
        }
    }

    Ok(())
}
