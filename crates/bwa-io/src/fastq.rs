//! FASTQ reading with bwa-compatible QNAME derivation and fixed-size (`-K`) batching.

use bwa_core::{Error, Result};
use needletail::{parse_fastx_file, FastxReader};

/// One read: SAM QNAME, sequence, and (optional) quality string.
pub struct Record {
    pub name: String,
    pub seq: Vec<u8>,
    pub qual: Option<Vec<u8>>,
}

/// Streaming FASTQ reader.
pub struct FastqReader {
    inner: Box<dyn FastxReader>,
}

impl FastqReader {
    /// Open a FASTQ (optionally gzipped) file.
    pub fn from_path<P: AsRef<std::path::Path>>(path: P) -> Result<Self> {
        let inner = parse_fastx_file(path).map_err(|e| Error::Fastq(e.to_string()))?;
        Ok(Self { inner })
    }

    /// Pull the next read, or `None` at EOF.
    pub fn next_record(&mut self) -> Result<Option<Record>> {
        match self.inner.next() {
            None => Ok(None),
            Some(rec) => {
                let rec = rec.map_err(|e| Error::Fastq(e.to_string()))?;
                let name = qname_from_id(rec.id());
                let seq = rec.seq().into_owned();
                let qual = rec.qual().map(<[u8]>::to_vec);
                Ok(Some(Record { name, seq, qual }))
            }
        }
    }

    /// Read a batch whose cumulative sequence length reaches at least `k_batch` bytes (or EOF).
    /// Fixing this boundary (via `-K`) is what makes downstream per-batch statistics reproducible.
    pub fn next_batch(&mut self, k_batch: usize) -> Result<Vec<Record>> {
        let mut batch = Vec::new();
        let mut bases = 0usize;
        while bases < k_batch {
            match self.next_record()? {
                None => break,
                Some(rec) => {
                    bases += rec.seq.len();
                    batch.push(rec);
                }
            }
        }
        Ok(batch)
    }
}

/// Derive the SAM QNAME from a FASTQ id line, mirroring bwa: take the field up to the first
/// whitespace, then trim a trailing `/<digit>` (bwa's `trim_readno`).
fn qname_from_id(id: &[u8]) -> String {
    let end = id
        .iter()
        .position(u8::is_ascii_whitespace)
        .unwrap_or(id.len());
    let mut s = &id[..end];
    if s.len() > 2 && s[s.len() - 2] == b'/' && s[s.len() - 1].is_ascii_digit() {
        s = &s[..s.len() - 2];
    }
    String::from_utf8_lossy(s).into_owned()
}

#[cfg(test)]
mod tests {
    use super::qname_from_id;

    #[test]
    fn qname_strips_comment_and_readno() {
        assert_eq!(qname_from_id(b"read1 some comment"), "read1");
        assert_eq!(qname_from_id(b"read1/1"), "read1");
        assert_eq!(qname_from_id(b"read1/2 desc"), "read1");
        assert_eq!(
            qname_from_id(b"20:2000000-2200000_50861_51313_0:0:0_0:1:0_0"),
            "20:2000000-2200000_50861_51313_0:0:0_0:1:0_0"
        );
    }
}
