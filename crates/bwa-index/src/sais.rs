//! Suffix array construction by induced sorting (SA-IS), pure Rust.
//!
//! [`suffix_array_with_sentinel`] returns the suffix array of a byte string with an implicit `$`
//! sentinel that sorts before every real symbol. This matches bwa-mem2's convention exactly: the
//! returned vector has length `n+1`, element `0` is `n` (the sentinel suffix), and `[1..]` is the
//! suffix array of the real string (equivalent to `saisxx(s, sa+1, n)` followed by `sa[0]=n`).

const EMPTY: usize = usize::MAX;

/// Suffix array of `s` including the sentinel row first.
///
/// Result length is `s.len()+1`; `result[0] == s.len()` and `result[1..]` is the suffix array of
/// `s` (shorter suffixes sort first on ties, i.e. the standard suffix array).
pub fn suffix_array_with_sentinel(s: &[u8]) -> Vec<i64> {
    let n = s.len();
    // Map bytes to 1..=256 and append the 0 sentinel, so 0 is uniquely smallest.
    let mut t: Vec<usize> = Vec::with_capacity(n + 1);
    for &b in s {
        t.push(b as usize + 1);
    }
    t.push(0);
    sa_is(&t, 257).into_iter().map(|x| x as i64).collect()
}

#[inline]
fn is_lms(is_s: &[bool], i: usize) -> bool {
    i > 0 && is_s[i] && !is_s[i - 1]
}

/// Bucket start offsets (`out[c]` = first index of symbol `c`'s bucket).
fn bucket_starts(s: &[usize], k: usize) -> Vec<usize> {
    let mut c = vec![0usize; k];
    for &x in s {
        c[x] += 1;
    }
    let mut sum = 0;
    for slot in c.iter_mut() {
        let cnt = *slot;
        *slot = sum;
        sum += cnt;
    }
    c
}

/// Bucket end offsets (`out[c]` = one past the last index of symbol `c`'s bucket).
fn bucket_ends(s: &[usize], k: usize) -> Vec<usize> {
    let mut c = vec![0usize; k];
    for &x in s {
        c[x] += 1;
    }
    let mut sum = 0;
    for slot in c.iter_mut() {
        sum += *slot;
        *slot = sum;
    }
    c
}

/// Induce L-type then S-type suffixes from the LMS suffixes already placed in `sa`.
fn induce(sa: &mut [usize], s: &[usize], is_s: &[bool], k: usize) {
    let n = s.len();
    let mut starts = bucket_starts(s, k);
    for i in 0..n {
        let j = sa[i];
        if j != EMPTY && j != 0 && !is_s[j - 1] {
            let c = s[j - 1];
            sa[starts[c]] = j - 1;
            starts[c] += 1;
        }
    }
    let mut ends = bucket_ends(s, k);
    for i in (0..n).rev() {
        let j = sa[i];
        if j != EMPTY && j != 0 && is_s[j - 1] {
            let c = s[j - 1];
            ends[c] -= 1;
            sa[ends[c]] = j - 1;
        }
    }
}

fn lms_substr_equal(s: &[usize], is_s: &[bool], a: usize, b: usize) -> bool {
    if a == b {
        return true;
    }
    let n = s.len();
    let mut i = 0usize;
    loop {
        let pa = a + i;
        let pb = b + i;
        if pa >= n || pb >= n {
            return pa >= n && pb >= n;
        }
        if s[pa] != s[pb] || is_s[pa] != is_s[pb] {
            return false;
        }
        if i > 0 {
            let al = is_lms(is_s, pa);
            let bl = is_lms(is_s, pb);
            if al || bl {
                return al && bl;
            }
        }
        i += 1;
    }
}

fn sa_is(s: &[usize], k: usize) -> Vec<usize> {
    let n = s.len();
    let mut sa = vec![EMPTY; n];
    if n == 1 {
        sa[0] = 0;
        return sa;
    }

    // Classify S-type (true) / L-type positions.
    let mut is_s = vec![false; n];
    is_s[n - 1] = true;
    for i in (0..n - 1).rev() {
        is_s[i] = s[i] < s[i + 1] || (s[i] == s[i + 1] && is_s[i + 1]);
    }

    // Step 1: place LMS suffixes at bucket ends, then induce.
    let mut ends = bucket_ends(s, k);
    for i in (0..n).rev() {
        if is_lms(&is_s, i) {
            let c = s[i];
            ends[c] -= 1;
            sa[ends[c]] = i;
        }
    }
    induce(&mut sa, s, &is_s, k);

    // Step 2: name LMS substrings in sorted order.
    let mut lms_sorted = Vec::new();
    for &p in sa.iter() {
        if p != EMPTY && is_lms(&is_s, p) {
            lms_sorted.push(p);
        }
    }
    let mut names = vec![EMPTY; n];
    let mut cur = 0usize;
    names[lms_sorted[0]] = 0;
    let mut prev = lms_sorted[0];
    for &p in lms_sorted.iter().skip(1) {
        if !lms_substr_equal(s, &is_s, prev, p) {
            cur += 1;
        }
        names[p] = cur;
        prev = p;
    }
    let num_names = cur + 1;

    // Reduced string in position order.
    let mut lms_positions = Vec::new();
    let mut reduced = Vec::new();
    for (i, &name) in names.iter().enumerate() {
        if is_lms(&is_s, i) {
            lms_positions.push(i);
            reduced.push(name);
        }
    }

    // Step 3: suffix array of the reduced string.
    let reduced_sa = if num_names == reduced.len() {
        let mut rsa = vec![0usize; reduced.len()];
        for (i, &name) in reduced.iter().enumerate() {
            rsa[name] = i;
        }
        rsa
    } else {
        sa_is(&reduced, num_names)
    };

    // Step 4: re-place LMS suffixes in sorted order, then induce the final SA.
    for x in sa.iter_mut() {
        *x = EMPTY;
    }
    let mut ends = bucket_ends(s, k);
    for &idx in reduced_sa.iter().rev() {
        let p = lms_positions[idx];
        let c = s[p];
        ends[c] -= 1;
        sa[ends[c]] = p;
    }
    induce(&mut sa, s, &is_s, k);

    sa
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Naive suffix array (with sentinel-first), for validation.
    fn naive(s: &[u8]) -> Vec<i64> {
        let n = s.len();
        let mut idx: Vec<usize> = (0..n).collect();
        idx.sort_by(|&a, &b| s[a..].cmp(&s[b..]));
        let mut out = vec![n as i64];
        out.extend(idx.into_iter().map(|x| x as i64));
        out
    }

    fn check(s: &[u8]) {
        assert_eq!(
            suffix_array_with_sentinel(s),
            naive(s),
            "mismatch for {s:?}"
        );
    }

    #[test]
    fn known_strings() {
        check(b"");
        check(b"a");
        check(b"aa");
        check(b"banana");
        check(b"mississippi");
        check(b"abracadabra");
        check(&[0, 0, 0, 1, 0, 1, 1]);
        check(&[3, 2, 1, 0, 3, 2, 1, 0]);
    }

    #[test]
    fn random_small_alphabet() {
        // Deterministic pseudo-random strings over {0,1,2,3} (the DNA alphabet), various lengths.
        let mut state: u64 = 0x1234_5678_9abc_def0;
        let mut next = || {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            state
        };
        for _ in 0..300 {
            let len = (next() % 200) as usize;
            let s: Vec<u8> = (0..len).map(|_| (next() % 4) as u8).collect();
            check(&s);
        }
    }

    #[test]
    fn random_byte_alphabet() {
        let mut state: u64 = 0xdead_beef_cafe_babe;
        let mut next = || {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            state
        };
        for _ in 0..100 {
            let len = (next() % 150) as usize;
            let s: Vec<u8> = (0..len).map(|_| (next() % 256) as u8).collect();
            check(&s);
        }
    }
}
