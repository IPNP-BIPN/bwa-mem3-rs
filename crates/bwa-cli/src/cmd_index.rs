//! `bwa-mem3 index` subcommand (phase 1; currently a stub).

use clap::Args;

#[derive(Args)]
pub struct IndexArgs {
    /// FASTA reference to index.
    pub fasta: String,
}

pub fn run(_args: IndexArgs) -> anyhow::Result<()> {
    anyhow::bail!(
        "`bwa-mem3 index` is not implemented yet (phase 1); use `bwa-mem2 index` meanwhile"
    )
}
