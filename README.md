![Topic](https://img.shields.io/badge/bioinformatics-sequence_decontamination-blue)
![GitHub Workflow Status](https://img.shields.io/github/actions/workflow/status/jackdougle/nucleaze/rust.yml)
[![License: MIT](https://img.shields.io/badge/license-MIT-yellow.svg)](LICENSE) 

# **Nucleaze 🧬**
A high-performance Rust tool for filtering DNA/RNA reads based on a set of reference k-mers.
Inspired by [BBDuk](https://archive.jgi.doe.gov/data-and-tools/software-tools/bbtools/bb-tools-user-guide/bbduk-guide/) by Brian Bushnell. Provides performant and memory-efficient read processing with support for both paired and unpaired FASTA/FASTQ files, with multiple files or in interleaved format.  

---

## **Features and Default Behavior**

**K-mer based read filtering**:  
- Reads are compared to reference sequences by matching k-mers.
- If a read sequence has at least x k-mers also found in reference dataset, it is a match
  - x is 1 by default, changed with `--minhits <int>`

**Piping**:  
- Use `--in stdin` to pipe from stdin
- use `--outm`/`outu`/`outm2`/`outu2` `stdout.fa`/`stdout.fq` to pipe results to stdout

**Paired reads support**:  
- Paired inputs and outputs can be specified by adding more input/output files
- Interleaved inputs or outputs, signify interleaved input with `--interinput`
- Automatic detection of input/output mode

**Multithreading with Rayon**:  
- Adjustable thread count via `--threads` argument  
- Defaults to all available CPU cores

**Memory Limit**:  
- Specify maximum memory usage with `--maxmem <String>` (e.g., `5G` for 5 gigabytes, `500M` for 500 megabytes)  

**Automatic Reference Indexing**:  
- Builds a serialized reference k-mer index using Bincode if `--binref <file>` is provided from references provided with `--ref <file>`
- Uses saved index on subsequent runs if `--binref <file>` is included

---

## **Installation**

### **1. Install Rust**
If using UNIX, run this command and follow the ensuing instructions:

```bash
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
```

If using Windows, download the correct installer from [Rustup](https://rustup.rs/#).

### **2. Download the release executable**

```bash
brew install nucleaze
```

### **3. Run program with necessary arguments**
Nucleaze requires `--in` and at least one of `--ref` or `--binref` to be provided.

See more parameter documentation at **[./src/main.rs](/src/main.rs)**

---

## **Example Usage**
```bash
./nucleaze --in reads.fq --ref refs.fa --outm matched.fq --outu unmatched.fq --k 21
```

This command:
1. Builds 21-mer index from `refs.fa` sequences
2. Reads input reads from `reads.fq` into chunks of size 10,000
3. Processes each read into 21-mers and checks against reference index
4. Outputs matched reads to `matched.fq` and unmatched reads to `unmatched.fq`

---

## **License**

This project is licensed under the MIT License, see [LICENSE](LICENSE) for details. There is lots of room for improvement here so new additions or suggestions are welcome!

---

## **Crates Used**

- [Needletail](https://github.com/onecodex/needletail) — FASTA/FASTQ file parsing and bitkmer operations
- [Bincode](https://sr.ht/~stygianentity/bincode/) — K-mer hashset serialization/deserialization
- [Rayon](https://github.com/rayon-rs/rayon) — Multithreading
- [Clap](https://github.com/clap-rs/clap) — CLI
- [Num-Cpus](https://github.com/seanmonstar/num_cpus) — detection of available threads
- [Sysinfo](https://github.com/GuillaumeGomez/sysinfo) — system memory information
- [Crossbeam](https://github.com/crossbeam-rs/crossbeam) — asynchronous channels

---

#### Please email jack.gdouglass@gmail.com with any questions or feature requests.
