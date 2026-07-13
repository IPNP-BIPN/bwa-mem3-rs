//! `bwa-mem3 mem` subcommand.
//!
//! Phase 0: load the reference metadata for the SAM header, stream reads in fixed-size (`-K`)
//! batches, and emit every read as unmapped. Real alignment lands in phases 3+.

use std::io::{BufWriter, Write};
use std::path::PathBuf;

use clap::Args;

use bwa_core::MemOpt;
use bwa_index::BntSeq;
use bwa_io::{sam, FastqReader, SqRecord};

#[derive(Args)]
pub struct MemArgs {
    /// Number of threads (phase 0 is single-threaded; accepted but ignored).
    #[arg(short = 't', default_value_t = 1)]
    pub threads: i32,
    /// Process INT input bases per batch (fixes batch boundaries for reproducibility).
    #[arg(short = 'K')]
    pub k_batch: Option<i64>,
    /// Index prefix: the FASTA path that was indexed.
    pub index_prefix: PathBuf,
    /// Reads in FASTQ (phase 0: single-end / R1 only).
    pub reads: PathBuf,
}

pub fn run(args: MemArgs, argv: &[String]) -> anyhow::Result<()> {
    let opt = MemOpt::default();
    let k_batch = args
        .k_batch
        .unwrap_or(opt.chunk_size * i64::from(args.threads))
        .max(1) as usize;

    let bns = BntSeq::load(&args.index_prefix)?;
    let sqs: Vec<SqRecord> = bns
        .contigs
        .iter()
        .map(|c| SqRecord {
            name: c.name.clone(),
            len: i64::from(c.len),
        })
        .collect();

    let stdout = std::io::stdout();
    let mut out = BufWriter::new(stdout.lock());

    let cl = argv.join(" ");
    sam::write_header(
        &mut out,
        &sqs,
        "bwa-mem3",
        "bwa-mem3",
        env!("CARGO_PKG_VERSION"),
        &cl,
    )?;

    // Phase 0: read in fixed-size batches, emit every read as unmapped (FLAG 4).
    let mut reader = FastqReader::from_path(&args.reads)?;
    loop {
        let batch = reader.next_batch(k_batch)?;
        if batch.is_empty() {
            break;
        }
        for rec in &batch {
            sam::write_unmapped(&mut out, &rec.name, &rec.seq, rec.qual.as_deref())?;
        }
    }
    out.flush()?;
    Ok(())
}
