//! glibc `srand48`/`lrand48` reproduction (48-bit linear congruential generator).
//!
//! bwa-mem2 randomizes ambiguous (N) bases in `.pac` with `lrand48() & 3` after `srand48(11)`.
//! Reproducing glibc's generator exactly is required for byte-identical `.pac` on references that
//! contain N runs. (N-free references never call it.)

/// The 48-bit LCG state, matching glibc's `drand48` family.
pub struct Rand48 {
    state: u64, // only low 48 bits used
}

const MULT: u64 = 0x5DEECE66D;
const ADD: u64 = 0xB;
const MASK48: u64 = 0xFFFF_FFFF_FFFF;

impl Rand48 {
    /// Equivalent to glibc `srand48(seed)`: state = `(seed << 16) | 0x330E`.
    pub fn srand48(seed: i64) -> Self {
        let state = (((seed as u64) & 0xFFFF_FFFF) << 16) | 0x330E;
        Rand48 {
            state: state & MASK48,
        }
    }

    fn step(&mut self) -> u64 {
        self.state = self.state.wrapping_mul(MULT).wrapping_add(ADD) & MASK48;
        self.state
    }

    /// Equivalent to glibc `lrand48()`: the high 31 bits of the next state, in `[0, 2^31)`.
    pub fn lrand48(&mut self) -> u32 {
        (self.step() >> 17) as u32
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn matches_libc_lrand48_seed11() {
        // Reference values from the system libc: srand48(11); lrand48() x8 (the oracle ran here).
        let mut r = Rand48::srand48(11);
        let got: Vec<u32> = (0..8).map(|_| r.lrand48()).collect();
        assert_eq!(
            got,
            vec![
                1609868485, 1074594562, 470884846, 2128573038, 960673312, 346697164, 303961605,
                444770020
            ]
        );
    }
}
