//! Seed extension via banded Smith-Waterman.
//!
//! Phase 0 defines only the backend trait; the scalar implementation (the bit-identity source of
//! truth, mirroring the patched `bandedSWA.cpp`) lands in phase 5, with NEON/Metal backends later.

/// A batched banded Smith-Waterman backend.
///
/// The scalar backend is authoritative; SIMD (NEON) and GPU (Metal) backends must produce
/// byte-identical integer results to it.
pub trait SwBackend {
    /// Short backend name, e.g. `"scalar"`, `"neon"`, `"metal"`.
    fn name(&self) -> &'static str;
}
