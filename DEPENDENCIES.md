# Dependances

Versions exactes gelees dans `Cargo.lock` (pin de record). Justifications :

| Crate | Usage | Licence |
|---|---|---|
| `clap` | parsing CLI (`bwa-mem3 index`/`mem`) | MIT/Apache-2.0 |
| `thiserror` | type d'erreur de `bwa-core` | MIT/Apache-2.0 |
| `anyhow` | erreurs applicatives (CLI) | MIT/Apache-2.0 |
| `needletail` | lecture FASTQ/FASTA (hot path, analogue `kseq.h`) | MIT |
| `memmap2` | mmap de l'index (phase 2+) | MIT/Apache-2.0 |
| `rayon` | parallelisme (phase 8, differe) | MIT/Apache-2.0 |
| `serde`, `serde_json` | rapport de concordance JSON (`sam-diff`) | MIT/Apache-2.0 |
| `sha2` | shasum de determinisme dans le harnais | MIT/Apache-2.0 |

Dev-only (jamais dans le binaire livre) : `noodles-*` (validation SAM/BAM cote test), `rust-bio` /
`block-aligner` (bootstrap de seeding/SW).
