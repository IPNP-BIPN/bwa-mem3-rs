//! Byte-identity gate for the indexer: building the tiny fixture must reproduce the five index
//! files produced by `bwa-mem2 index` exactly.

use std::path::{Path, PathBuf};

fn tiny_dir() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../../testdata/tiny")
}

#[test]
fn tiny_index_is_byte_identical() {
    let src = tiny_dir().join("tiny.fa");
    let dir = std::env::temp_dir().join(format!("bwamem3_idx_{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let dst = dir.join("tiny.fa");
    std::fs::copy(&src, &dst).unwrap();

    bwa_index::build_index(&dst).expect("build_index");

    for ext in ["pac", "ann", "amb", "0123", "bwt.2bit.64"] {
        let ours = std::fs::read(dir.join(format!("tiny.fa.{ext}"))).unwrap();
        let oracle = std::fs::read(tiny_dir().join(format!("tiny.fa.{ext}"))).unwrap();
        if ours != oracle {
            let pos = ours.iter().zip(&oracle).position(|(a, b)| a != b);
            std::fs::remove_dir_all(&dir).ok();
            panic!(
                ".{ext}: differs at byte {:?} (len ours={}, oracle={})",
                pos,
                ours.len(),
                oracle.len()
            );
        }
    }
    std::fs::remove_dir_all(&dir).ok();
}
