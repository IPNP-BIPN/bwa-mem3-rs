# bwa-mem3-rs

Réimplémentation native Rust, complète et bit-identique, de l'aligneur de reads courts
[bwa-mem2](https://github.com/bwa-mem2/bwa-mem2) (indexeur inclus), avec accélération GPU (Metal)
prévue en phase finale.

Projet dans la lignée de STAR-rs / piscem-rs. Objectif d'acceptation : sortie **octet-identique**
(index et SAM) au binaire `bwa-mem2` de référence.

## Oracle de référence

La parité vise le binaire `bwa-mem2` Homebrew installé (`bwa-mem2 version` => `2.3`), qui est un build
**patché** : upstream tag `v2.3` (rév `7aa5ff6c…`, source 2.2.1), plus `fastmap.patch` et
`bandedSWA.cpp.patch`, SIMD via sse2neon (voie SSE 128-bit traduite NEON). Voir le plan et
`scripts/setup_reference.sh`.

## Structure

Workspace Cargo de crates `bwa-*` (une par étage du pipeline). Binaire : `bwa-mem3` (`index`, `mem`).

| Crate | Rôle |
|---|---|
| `bwa-core` | Types, constantes, options d'alignement |
| `bwa-io` | FASTA/FASTQ in, sortie SAM (formatée à la main) |
| `bwa-index` | Construction + chargement de l'index FMD |
| `bwa-seed` | Seeding SMEM (FM-index) |
| `bwa-chain` | Chaînage des graines |
| `bwa-extend` | Extension Smith-Waterman bandée (scalaire ; trait `SwBackend`) |
| `bwa-sam` | Marquage primaire, MAPQ, CIGAR, tags |
| `bwa-mem` | Glue du pipeline (SE + PE) |
| `bwa-cli` | Binaire `bwa-mem3` |
| `bwa-diff` | Concordance SAM (`sam-diff`) |
| `bwa-gpu` | Backend GPU Metal du SW (phase finale) |

## Développement

```sh
cargo build --release
bash scripts/check.sh        # fmt + clippy + tests
bash scripts/oracle_diff.sh  # diff SAM vs bwa-mem2 natif
```

## Statut

En cours, phase 0 (squelette). Voir `ROADMAP.md`.

## Licence

MIT (comme bwa-mem2).
