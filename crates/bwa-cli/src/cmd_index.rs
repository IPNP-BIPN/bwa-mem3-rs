//! `bwa-mem3 index` subcommand: build the FMD index (byte-identical to `bwa-mem2 index`).

use std::path::PathBuf;
use std::time::Instant;

use clap::Args;

#[derive(Args)]
pub struct IndexArgs {
    /// FASTA reference to index.
    pub fasta: PathBuf,
}

pub fn run(args: IndexArgs) -> anyhow::Result<()> {
    let t0 = Instant::now();
    bwa_index::build_index(&args.fasta)?;
    eprintln!(
        "[bwa-mem3 index] built index for {} in {:.3}s",
        args.fasta.display(),
        t0.elapsed().as_secs_f64()
    );
    Ok(())
}
