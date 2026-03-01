//! I/O for reference indexing and read processing operations using k-mers
use crate::kmer_ops::{HashShards, KmerProcessor, KmerStore, MiniBloom};

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

/// A chunk of sequences with their match results
#[derive(Eq, PartialEq)]
struct SequenceChunk {
    id: u32,
    data_arena: Vec<u8>,                          // raw bytes for all sequences
    offsets: Vec<(u32, u32, u32, u32, u32, u32)>, // (id_start, id_len, seq_start, seq_len, qual_start, qual_len)
    matches: Vec<bool>,                           // k-mer match results
}

// Implement Ord for BinaryHeap to support ordered output
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

    // Extract remaining args fields before dispatch (avoids partial move)
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
        // Bloom mode — no serialization support
        if !bin_kmers_path.is_empty() {
            eprintln!("Warning: --binref is not supported in bloom filter mode. Ignoring.");
        }
        if !save_kmers_path.is_empty() {
            eprintln!("Warning: --saveref is not supported in bloom filter mode. Ignoring.");
        }
        let store = MiniBloom::new(ref_size as usize, fpr as f64);
        let mut processor = KmerProcessor::new(k, min_hits, use_canonical, store);

        match get_reference_kmers(&ref_path, &mut processor, num_threads) {
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

        run_inner(
            processor,
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
        // Hash mode — full serialization support
        let mut store = HashShards::new(0, 0.0);

        // Attempt to load pre-built index directly onto the store before going generic
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
            match get_reference_kmers(&ref_path, &mut processor, num_threads) {
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

        run_inner(
            processor,
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

fn run_inner<S: KmerStore + 'static>(
    processor: KmerProcessor<S>,
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

    let kmer_processor = Arc::new(processor);
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
        kmer_processor,
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

    // Verify metadata sentinel
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

    spawn_reader(ref_path, sender).expect("k-mer extraction failed");

    (0..num_threads).into_par_iter().for_each(|_| {
        let mut local_idx = vec![Vec::with_capacity(64); SUBINDEX_COUNT];

        while let Ok(seq) = receiver.recv() {
            processor.process_ref(&seq, &mut local_idx);

            for (i, local_subidx) in local_idx.iter_mut().enumerate() {
                if !local_subidx.is_empty() {
                    let mut merged_subidx = merged_idx[i].lock().unwrap();

                    for &kmer in local_subidx.iter() {
                        merged_subidx.insert(kmer);
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

// Delegate single thread to parse ref file and send read sequences
fn spawn_reader(path: &str, sender: Sender<Vec<u8>>) -> Result<(), Box<dyn Error>> {
    let path = path.to_string();

    thread::spawn(move || {
        let mut reader = parse_fastx_file(&path).expect("FASTA open failed");

        while let Some(record) = reader.next() {
            let record = match record {
                Ok(r) => r,
                Err(_) => return,
            };

            sender.send(record.seq().to_vec()).unwrap();
        }
    });

    Ok(())
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
fn process_reads<S: KmerStore + 'static>(
    input: InputSource,
    reads2_path: String,
    processor: Arc<KmerProcessor<S>>,
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

    // Channel for passing chunks from reader thread to writer
    let (chunk_sender, chunk_receiver): (SyncSender<SequenceChunk>, Receiver<SequenceChunk>) =
        sync_channel(20);

    let chunk_pos = Arc::new(AtomicU32::new(0));

    // Extract file extensions
    let matched_filetype = matched_path.rsplit('.').next().unwrap_or("").to_string();
    let unmatched_filetype = unmatched_path.rsplit('.').next().unwrap_or("").to_string();
    let matched2_filetype = matched2_path.rsplit('.').next().unwrap_or("").to_string();
    let unmatched2_filetype = unmatched2_path.rsplit('.').next().unwrap_or("").to_string();

    // Check if output is to stdout
    let matched_stdout = matched_path == "stdout" || matched_path.starts_with("stdout.");
    let unmatched_stdout = unmatched_path == "stdout" || unmatched_path.starts_with("stdout.");
    let matched2_stdout = matched2_path == "stdout" || matched2_path.starts_with("stdout.");
    let unmatched2_stdout = unmatched2_path == "stdout" || unmatched2_path.starts_with("stdout.");

    let parallel_sender = chunk_sender.clone();
    let parellel_chunk_pos = chunk_pos.clone();

    // Worker thread: reads sequences and dispatches chunks to Rayon for parallel k-mer processing
    let worker_thread = thread::spawn(move || -> Result<(), Box<dyn Error + Send + Sync>> {
        let mut reader = match input {
            InputSource::Stdin => parse_fastx_reader(stdin()),
            InputSource::File(ref path) => parse_fastx_file(path),
        }?;

        const CHUNK_SIZE: usize = 10_000;
        const ARENA_CAPACITY: usize = CHUNK_SIZE * 500;

        // Current chunk being built
        let mut arena: Vec<u8> = Vec::with_capacity(ARENA_CAPACITY);
        let mut offsets: Vec<(u32, u32, u32, u32, u32, u32)> = Vec::with_capacity(CHUNK_SIZE);

        // Dispatch filled arena to Rayon for parallel k-mer processing
        let process_arena =
            |local_arena: Vec<u8>, local_offsets: Vec<(u32, u32, u32, u32, u32, u32)>| {
                let processor = processor.clone();
                let sender = parallel_sender.clone();
                let current_chunk_pos = parellel_chunk_pos.fetch_add(1, AtomicOrdering::SeqCst);

                rayon::spawn(move || {
                    // Parallel k-mer matching: map over offsets to create slices into arena
                    // This avoids copying sequence data
                    let matches: Vec<bool> = local_offsets
                        .par_iter()
                        .map(|(_, _, seq_start, seq_len, _, _)| {
                            let seq =
                                &local_arena[*seq_start as usize..(*seq_start + *seq_len) as usize];
                            processor.process_read(seq)
                        })
                        .collect();

                    let chunk = SequenceChunk {
                        id: current_chunk_pos,
                        data_arena: local_arena,
                        offsets: local_offsets,
                        matches,
                    };

                    // Blocks if channel is full (backpressure)
                    let _ = sender.send(chunk);
                });
            };

        // Add one record to the arena
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

        // Read and chunk sequences based on mode
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
                    // Double buffering: swap buffers and dispatch filled one
                    let local_arena = mem::replace(&mut arena, Vec::with_capacity(ARENA_CAPACITY));
                    let local_offsets = mem::replace(&mut offsets, Vec::with_capacity(CHUNK_SIZE));
                    process_arena(local_arena, local_offsets);
                }
            }
        } else {
            // Paired modes: read from two files simultaneously
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
                            process_arena(local_arena, local_offsets);
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

        // Process remaining sequences
        if !offsets.is_empty() {
            process_arena(arena, offsets);
        }

        Ok(())
    });

    // Close channel from sender side so receiver knows when to stop
    drop(chunk_sender);

    let mut matched_writer: BufWriter<File> = BufWriter::new(File::create(matched_path)?);
    let mut unmatched_writer: BufWriter<File> =
        BufWriter::with_capacity(4_000_000, File::create(unmatched_path)?);

    // Optional writers for paired output splitting
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

    // Write a SequenceChunk to disk, reconstructing reads from the arena
    let mut chunk_output = |chunk: &SequenceChunk| -> Result<(), Box<dyn Send + Sync + Error>> {
        let arena = &chunk.data_arena;
        let offsets = &chunk.offsets;
        let matches = &chunk.matches;

        // Slice arena to reconstruct individual read
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
            // Paired-end: process reads in pairs (stride 2)
            // If either read matches, both are written to matched output
            let num_reads = offsets.len();
            let reads_to_process = if num_reads % 2 != 0 {
                // Odd number of reads in paired/interleaved mode - skip the last unpaired read
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
        // Ordered mode: use min-heap to reorder chunks as they arrive
        let mut out_of_order_buffer: BinaryHeap<SequenceChunk> = BinaryHeap::new();
        let mut next_chunk_id = 0;
        const MAX_BUFFERED_CHUNKS: usize = 1000;

        for chunk in chunk_receiver {
            if chunk.id == next_chunk_id {
                // Chunk arrived in order
                chunk_output(&chunk)?;
                next_chunk_id += 1;

                // Check if subsequent chunks are already buffered
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
                // Chunk arrived early, buffer it
                out_of_order_buffer.push(chunk);
                if out_of_order_buffer.len() > MAX_BUFFERED_CHUNKS {
                    return Err(Box::from("Too many out-of-order chunks buffered."));
                }
            }
        }
        // Drain remaining buffered chunks
        while let Some(buffered) = out_of_order_buffer.pop() {
            if buffered.id == next_chunk_id {
                chunk_output(&buffered)?;
                next_chunk_id += 1;
            } else {
                return Err(Box::from("Missing chunk in ordered output stream."));
            }
        }
    } else {
        // Unordered mode: write chunks as they arrive
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
