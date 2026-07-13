//! `bwa-mem3` command-line entry point.

use clap::{Parser, Subcommand};

mod cmd_index;
mod cmd_mem;

#[derive(Parser)]
#[command(
    name = "bwa-mem3",
    version,
    about = "Native Rust reimplementation of bwa-mem2"
)]
struct Cli {
    #[command(subcommand)]
    cmd: Cmd,
}

#[derive(Subcommand)]
enum Cmd {
    /// Build the FMD index from a FASTA reference (phase 1; not yet implemented).
    Index(cmd_index::IndexArgs),
    /// Align reads to an indexed reference.
    Mem(cmd_mem::MemArgs),
}

fn main() -> anyhow::Result<()> {
    // Capture the raw command line for the @PG CL tag before clap consumes it.
    let argv: Vec<String> = std::env::args().collect();
    match Cli::parse().cmd {
        Cmd::Index(args) => cmd_index::run(args),
        Cmd::Mem(args) => cmd_mem::run(args, &argv),
    }
}
