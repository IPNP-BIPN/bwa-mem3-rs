//! `bwa-mem3 mem` subcommand.
//!
//! Phase 6: seed -> chain -> extend -> best region -> `reg2aln` (exact CIGAR + NM/MD). MAPQ and
//! secondary/XA handling follow.

use std::io::{BufWriter, Write};
use std::path::PathBuf;

use clap::Args;

use bwa_core::{dna, MemOpt};
use bwa_index::{BntSeq, FmIndex};
use bwa_io::{sam, FastqReader, SqRecord};
use bwa_mem::{align_read_se, cigar_string, mem_approx_mapq_se, reg2aln};

#[derive(Args)]
pub struct MemArgs {
    /// Number of threads (phase 6 is single-threaded; accepted but ignored).
    #[arg(short = 't', default_value_t = 1)]
    pub threads: i32,
    /// Process INT input bases per batch (fixes batch boundaries for reproducibility).
    #[arg(short = 'K')]
    pub k_batch: Option<i64>,
    /// Index prefix: the FASTA path that was indexed.
    pub index_prefix: PathBuf,
    /// Reads in FASTQ (phase 6: single-end / R1 only).
    pub reads: PathBuf,
}

pub fn run(args: MemArgs, argv: &[String]) -> anyhow::Result<()> {
    let opt = MemOpt::default();
    let k_batch = args
        .k_batch
        .unwrap_or(opt.chunk_size * i64::from(args.threads))
        .max(1) as usize;

    let fm = FmIndex::load(&args.index_prefix)?;
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

    let mut reader = FastqReader::from_path(&args.reads)?;
    let mut read_id = 0u64;
    loop {
        let batch = reader.next_batch(k_batch)?;
        if batch.is_empty() {
            break;
        }
        for rec in &batch {
            let codes: Vec<u8> = rec.seq.iter().map(|&b| dna::nt4(b)).collect();
            let regs = align_read_se(&fm, &bns, &opt, &codes, read_id);
            read_id += 1;
            // After marking, regs[0] is the highest-scoring primary region.
            let best = regs.first().filter(|r| r.score >= opt.t);

            match best {
                Some(r) => {
                    let aln = reg2aln(&fm, &bns, &opt, codes.len() as i32, &codes, r);
                    let mapq = mem_approx_mapq_se(&opt, r);
                    let rname = &bns.contigs[aln.rid as usize].name;
                    let flag = if aln.is_rev { 16 } else { 0 };
                    let cigar = cigar_string(&aln.cigar);
                    let tags = format!(
                        "NM:i:{}\tMD:Z:{}\tAS:i:{}\tXS:i:{}",
                        aln.nm, aln.md, aln.score, aln.sub
                    );
                    let pos = aln.pos + 1;
                    if aln.is_rev {
                        let seq = dna::revcomp_ascii(&rec.seq);
                        let qual = rec.qual.as_ref().map(|q| {
                            let mut v = q.clone();
                            v.reverse();
                            v
                        });
                        sam::write_mapped_se(
                            &mut out,
                            &rec.name,
                            flag,
                            rname,
                            pos,
                            mapq,
                            &cigar,
                            &seq,
                            qual.as_deref(),
                            &tags,
                        )?;
                    } else {
                        sam::write_mapped_se(
                            &mut out,
                            &rec.name,
                            flag,
                            rname,
                            pos,
                            mapq,
                            &cigar,
                            &rec.seq,
                            rec.qual.as_deref(),
                            &tags,
                        )?;
                    }
                }
                None => sam::write_unmapped(&mut out, &rec.name, &rec.seq, rec.qual.as_deref())?,
            }
        }
    }
    out.flush()?;
    Ok(())
}
