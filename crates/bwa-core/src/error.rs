//! Crate-wide error type.

use thiserror::Error;

/// Errors produced across bwa-mem3-rs crates.
#[derive(Debug, Error)]
pub enum Error {
    /// Underlying I/O failure.
    #[error("I/O error: {0}")]
    Io(#[from] std::io::Error),
    /// A malformed or unexpected index file.
    #[error("index format error: {0}")]
    IndexFormat(String),
    /// A malformed FASTQ/FASTA input.
    #[error("sequence input error: {0}")]
    Fastq(String),
    /// Anything else, with context.
    #[error("{0}")]
    Other(String),
}

/// Convenience alias.
pub type Result<T> = std::result::Result<T, Error>;
