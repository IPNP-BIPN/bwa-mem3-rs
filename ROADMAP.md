# Roadmap

Une phase = une branche. Commits frequents ; merge vers `main` quand le gate de la phase passe.
Cible d'acceptation : index et SAM **octet-identiques** au binaire `bwa-mem2` 2.3 patche (oracle).

| Phase | Branche | But | Gate |
|---|---|---|---|
| 0 | `phase0-skeleton` | CLI `mem`/`index` (index stub), FASTQ -> SAM non-mappe, en-tete correct, harnais | `@SQ` octet-identique + JSON de concordance produit |
| 1 | `phase1-indexer` | `bwa-mem3 index` : `.pac/.ann/.amb/.bwt.2bit.64/.0123` | `cmp` octet vs `bwa-mem2 index` |
| 2 | `phase2-index-load` | chargeur + `get_occ`/`get_sa`/`backward_ext` + validation | assertions + cross-check occ/SMEM |
| 3 | `phase3-seeding` | SMEM + reseed + filtrage occ | seeds == oracle |
| 4 | `phase4-chaining` | `mem_chain` + `mem_chain_flt` | chaines == attendu |
| 5 | `phase5-extension` | SW bande scalaire (`SwBackend`) + `mem_chain2aln` | `AS` + coords == oracle |
| 6 | `phase6-se-sam` | primaire/MAPQ/CIGAR/tags | **SAM SE octet-identique** |
| 7 | `phase7-pe` | `mem_pestat`/`mem_matesw`/`mem_pair` | **SAM PE octet-identique** |
| 8 | `phase8-scale` | GRCh38 complet, rayon, (NEON optionnel) | ~100% concordance reads reels |
| 9 | `phase9-gpu` | backend Metal du SW (entier -> bit-identique) | identique au scalaire + speedup |
| 10+ | | gate GIAB `hap.py`/`vcfeval` ; packaging | |

Statut : **phase 1 terminee** (index octet-identique tiny + 2Mb via notre SA-IS). Prochaine : phase 2 (chargement + traversee FM-index).
