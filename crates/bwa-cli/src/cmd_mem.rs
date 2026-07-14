//! `bwa-mem3 mem` subcommand.
//!
//! Phase 6: seed -> chain -> extend -> best region -> `reg2aln` (exact CIGAR + NM/MD). MAPQ and
//! secondary/XA handling follow.

use std::io::{BufWriter, Write};
use std::path::PathBuf;

use clap::Args;

use bwa_core::{dna, MemOpt};
use bwa_index::{BntSeq, FmIndex};
use bwa_io::{sam, FastqReader, PairedFastqReader, SqRecord};
use bwa_mem::{
    align_read_dedup, align_read_se, cigar_string, mem_approx_mapq_se, mem_pestat, mem_sam_pe,
    reg2aln, MemAlnReg,
};

#[derive(Args)]
pub struct MemArgs {
    /// Number of threads (single-threaded for now; accepted but ignored).
    #[arg(short = 't', default_value_t = 1)]
    pub threads: i32,
    /// Process INT input bases per batch (fixes batch boundaries for reproducibility).
    #[arg(short = 'K')]
    pub k_batch: Option<i64>,
    /// Index prefix: the FASTA path that was indexed.
    pub index_prefix: PathBuf,
    /// Reads in FASTQ (R1, or the only file for single-end).
    pub reads: PathBuf,
    /// Optional mate reads (R2): triggers paired-end mode.
    pub reads2: Option<PathBuf>,
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

    if let Some(reads2) = args.reads2.clone() {
        run_pe(&fm, &bns, &opt, &args.reads, &reads2, k_batch, &mut out)?;
        out.flush()?;
        return Ok(());
    }

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

/// Paired-end driver: per batch, align+dedup both ends of every pair, estimate insert sizes
/// (`mem_pestat`), then emit paired SAM (`mem_sam_pe`). The pair index is global across batches (for
/// the `hash` tie-break), matching bwa-mem2's `(n_processed>>1)+i`.
#[allow(clippy::too_many_arguments)]
fn run_pe<W: std::io::Write>(
    fm: &FmIndex,
    bns: &BntSeq,
    opt: &MemOpt,
    reads1: &std::path::Path,
    reads2: &std::path::Path,
    k_batch: usize,
    out: &mut W,
) -> anyhow::Result<()> {
    let mut reader = PairedFastqReader::from_paths(reads1, reads2)?;
    let mut pair_id = 0u64;
    loop {
        let batch = reader.next_batch(k_batch)?;
        if batch.is_empty() {
            break;
        }
        // nt4-encode and align+dedup both ends of every pair (interleaved: 2i=R1, 2i+1=R2).
        let mut codes: Vec<(Vec<u8>, Vec<u8>)> = Vec::with_capacity(batch.len());
        let mut regs: Vec<Vec<MemAlnReg>> = Vec::with_capacity(batch.len() * 2);
        for (r1, r2) in &batch {
            let c1: Vec<u8> = r1.seq.iter().map(|&b| dna::nt4(b)).collect();
            let c2: Vec<u8> = r2.seq.iter().map(|&b| dna::nt4(b)).collect();
            regs.push(align_read_dedup(fm, bns, opt, &c1));
            regs.push(align_read_dedup(fm, bns, opt, &c2));
            codes.push((c1, c2));
        }
        let pes = mem_pestat(opt, bns.l_pac, &regs);
        for (i, (r1, r2)) in batch.iter().enumerate() {
            let mut a1 = std::mem::take(&mut regs[2 * i]);
            let mut a2 = std::mem::take(&mut regs[2 * i + 1]);
            let names = [r1.name.clone(), r2.name.clone()];
            let (c1, c2) = &codes[i];
            let seqs = [c1.as_slice(), c2.as_slice()];
            let quals = [r1.qual.as_deref(), r2.qual.as_deref()];
            mem_sam_pe(
                fm, bns, opt, &pes, pair_id, &names, &seqs, &quals, &mut a1, &mut a2, out,
            )?;
            pair_id += 1;
        }
    }
    Ok(())
}
