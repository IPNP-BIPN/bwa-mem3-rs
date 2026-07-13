//! SAM concordance: compare two SAM files field-by-field and summarise agreement.
//!
//! This is the inner-loop oracle-diff used to gate every phase. It compares primary alignments
//! (skipping secondary/supplementary) keyed by QNAME + read-in-pair bit.

use std::collections::HashMap;
use std::path::Path;

use serde::Serialize;

/// A minimal primary-alignment record (the fields we gate on).
#[derive(Debug, Clone, PartialEq)]
pub struct Rec {
    pub flag: u32,
    pub rname: String,
    pub pos: i64,
    pub mapq: u32,
    pub cigar: String,
}

/// A single field-level divergence example.
#[derive(Debug, Serialize)]
pub struct Divergence {
    pub key: String,
    pub field: String,
    pub oracle: String,
    pub ours: String,
}

/// Concordance summary between an oracle SAM and ours.
#[derive(Debug, Default, Serialize)]
pub struct Report {
    pub oracle_records: usize,
    pub our_records: usize,
    pub compared: usize,
    pub flag_match: usize,
    pub rname_pos_match: usize,
    pub mapq_match: usize,
    pub cigar_match: usize,
    pub all_fields_match: usize,
    pub only_in_oracle: usize,
    pub only_in_ours: usize,
    pub examples: Vec<Divergence>,
}

/// Parse primary alignments (skip `@` headers, secondary `0x100`, supplementary `0x800`), keyed
/// by QNAME plus the read-in-pair bit so PE mates do not collide.
pub fn parse_primary(path: &Path) -> std::io::Result<HashMap<String, Rec>> {
    let text = std::fs::read_to_string(path)?;
    let mut map = HashMap::new();
    for line in text.lines() {
        if line.starts_with('@') || line.is_empty() {
            continue;
        }
        let mut f = line.split('\t');
        let qname = f.next().unwrap_or("");
        let flag: u32 = f.next().unwrap_or("0").parse().unwrap_or(0);
        if flag & 0x100 != 0 || flag & 0x800 != 0 {
            continue;
        }
        let rname = f.next().unwrap_or("*").to_string();
        let pos: i64 = f.next().unwrap_or("0").parse().unwrap_or(0);
        let mapq: u32 = f.next().unwrap_or("0").parse().unwrap_or(0);
        let cigar = f.next().unwrap_or("*").to_string();
        let mate = if flag & 0x80 != 0 { 2 } else { 1 };
        map.insert(
            format!("{qname}/{mate}"),
            Rec {
                flag,
                rname,
                pos,
                mapq,
                cigar,
            },
        );
    }
    Ok(map)
}

/// Compare two parsed SAMs, producing a concordance report (up to `max_examples` divergences).
pub fn compare(
    oracle: &HashMap<String, Rec>,
    ours: &HashMap<String, Rec>,
    max_examples: usize,
) -> Report {
    let mut r = Report {
        oracle_records: oracle.len(),
        our_records: ours.len(),
        ..Default::default()
    };
    for (key, o) in oracle {
        let Some(u) = ours.get(key) else {
            r.only_in_oracle += 1;
            continue;
        };
        r.compared += 1;
        let flag_ok = o.flag == u.flag;
        let pos_ok = o.rname == u.rname && o.pos == u.pos;
        let mapq_ok = o.mapq == u.mapq;
        let cigar_ok = o.cigar == u.cigar;
        r.flag_match += usize::from(flag_ok);
        r.rname_pos_match += usize::from(pos_ok);
        r.mapq_match += usize::from(mapq_ok);
        r.cigar_match += usize::from(cigar_ok);
        if flag_ok && pos_ok && mapq_ok && cigar_ok {
            r.all_fields_match += 1;
        } else if r.examples.len() < max_examples {
            let (field, oracle_v, ours_v) = if !pos_ok {
                (
                    "RNAME/POS",
                    format!("{}:{}", o.rname, o.pos),
                    format!("{}:{}", u.rname, u.pos),
                )
            } else if !cigar_ok {
                ("CIGAR", o.cigar.clone(), u.cigar.clone())
            } else if !mapq_ok {
                ("MAPQ", o.mapq.to_string(), u.mapq.to_string())
            } else {
                ("FLAG", o.flag.to_string(), u.flag.to_string())
            };
            r.examples.push(Divergence {
                key: key.clone(),
                field: field.to_string(),
                oracle: oracle_v,
                ours: ours_v,
            });
        }
    }
    r.only_in_ours = ours.keys().filter(|k| !oracle.contains_key(*k)).count();
    r
}
