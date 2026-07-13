//! Sequence I/O: FASTQ input (needletail) and hand-formatted SAM output.
//!
//! SAM is written by hand (not via a library) so we control the exact bytes, which is required for
//! the bit-identity goal against bwa-mem2.

pub mod fastq;
pub mod sam;

pub use fastq::{FastqReader, Record};
pub use sam::{write_header, write_unmapped, SqRecord};
