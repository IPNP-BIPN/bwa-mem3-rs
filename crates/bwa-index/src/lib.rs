//! FMD index construction and loading.
//!
//! Phase 0 implements only reference-metadata parsing (`.ann`/`.amb`) needed for the SAM header.
//! Index construction (`build`) and the FM traversal (`fmindex`) arrive in phases 1-2.

pub mod bntseq;

pub use bntseq::{Amb, BntSeq, Contig};
