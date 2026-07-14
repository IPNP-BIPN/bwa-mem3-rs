//! Nucleotide encoding, matching bwa-mem2's `nst_nt4_table`.

/// Maps an ASCII byte to a 2-bit base code (A=0, C=1, G=2, T=3), or 4 for anything else (N).
pub const NT4_TABLE: [u8; 256] = build_nt4_table();

const fn build_nt4_table() -> [u8; 256] {
    let mut t = [4u8; 256];
    t[b'A' as usize] = 0;
    t[b'a' as usize] = 0;
    t[b'C' as usize] = 1;
    t[b'c' as usize] = 1;
    t[b'G' as usize] = 2;
    t[b'g' as usize] = 2;
    t[b'T' as usize] = 3;
    t[b't' as usize] = 3;
    t
}

/// Encode a base to its 2-bit code (4 = N / other).
#[inline]
pub fn nt4(b: u8) -> u8 {
    NT4_TABLE[b as usize]
}

/// Complement of a 2-bit base code (0<->3, 1<->2); 4 (N) stays 4.
#[inline]
pub fn comp2(c: u8) -> u8 {
    if c < 4 {
        3 - c
    } else {
        4
    }
}

/// Reverse-complement an ASCII nucleotide sequence (A<->T, C<->G, case-normalized to upper;
/// non-ACGT bytes pass through unchanged). Used for reverse-strand SAM SEQ output.
pub fn revcomp_ascii(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .rev()
        .map(|&b| match b {
            b'A' | b'a' => b'T',
            b'C' | b'c' => b'G',
            b'G' | b'g' => b'C',
            b'T' | b't' => b'A',
            other => other,
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn encodes_bases() {
        assert_eq!(nt4(b'A'), 0);
        assert_eq!(nt4(b'c'), 1);
        assert_eq!(nt4(b'G'), 2);
        assert_eq!(nt4(b't'), 3);
        assert_eq!(nt4(b'N'), 4);
        assert_eq!(nt4(b'-'), 4);
    }

    #[test]
    fn complements() {
        assert_eq!(comp2(nt4(b'A')), nt4(b'T'));
        assert_eq!(comp2(nt4(b'C')), nt4(b'G'));
        assert_eq!(comp2(4), 4);
    }
}
