# Divergences connues vs oracle bwa-mem2 2.3

Suivi de la traine de parite. Chaque entree : champ concerne, cause, statut, plan.

## Acceptees (par conception)

- **`@PG`** : notre sortie emet `ID:bwa-mem3 PN:bwa-mem3 VN:<ver> CL:<notre argv>`, l'oracle emet
  `bwa-mem2`. Exclu du gate d'octet-identite (on compare `@SQ` + lignes d'alignement). Decision finale
  (spoof eventuel de `@PG`) reportee en fin de projet.

## En cours

- **Parité des régions sous-optimales (`XS:i`, et indirectement `csub`/`sub_n` donc parfois MAPQ),
  ~6-7% des reads, SE et PE.** Sur 5000 reads/paires wgsim (chr20 2 Mb) :
  - **SE** : byte-identique sur **4658/5000** lignes ; FLAG/RNAME/POS/CIGAR 100%, MAPQ 99.5%.
  - **PE** : byte-identique sur **9329/10000** enregistrements ; le cœur (FLAG/POS/CIGAR/RNEXT/
    PNEXT/TLEN/MC) est quasi 100%, `mem_pestat` est **identique bit-à-bit** à l'oracle (mêmes
    percentiles/bornes), `mem_pair` + MAPQ combinée vérifiés.

  **Cause racine** (diagnostiquée) : notre cascade de seeds (SMEM 3 rondes + reseed) ne produit
  pas toujours **exactement** le même ensemble de régions sous-optimales que l'oracle. Exemple SE
  `_514496` : l'oracle trouve une région sous-optimale à score 44 (via reseed), nous n'obtenons
  qu'un seed partiel de 26 bp non étendu → `XS` 44 vs 26. Les écarts vont **dans les deux sens**
  (on trouve parfois une sous-optimale plus haute), donc ce n'est pas un simple biais de bande SW ;
  c'est la **complétude/parité exacte du SMEM**. `mem_patch_reg` est déjà porté (branche de fusion
  de `mem_sort_dedup_patch`) mais ne se déclenche pas sur ces reads courts.
  **Impact** : nul sur position/CIGAR des cas testés ; `XS` (et rarement MAPQ via `csub`) diffèrent.
  **Plan** : audit dédié de la parité SMEM (comparer nos seeds à `bwa-mem2 fastmap`), puis
  tie-breaks des tris de score égal.

- **`XA:Z` (hits alternatifs, `mem_gen_alt`) non porté.** ~37 enregistrements PE (et l'équivalent
  SE) où l'oracle émet un tag `XA:Z` listant des positions alternatives. Sous-système `mem_gen_alt`
  non encore implémenté (commun SE/PE). **Plan** : porter `mem_gen_alt` (phase 8).

- **Mate rescue (`mem_matesw` / `ksw_align2`) non porté.** Sur données concordantes il ne s'exécute
  pas (toutes orientations « skip »), donc l'impact ici est ~11 enregistrements (paires
  discordantes / bouts manquants). Requis pour la parité sur données réelles. **Plan** : porter
  `ksw_align2` puis `mem_matesw`.
