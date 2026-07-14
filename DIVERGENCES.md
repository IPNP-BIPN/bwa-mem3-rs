# Divergences connues vs oracle bwa-mem2 2.3

Suivi de la traine de parite. Chaque entree : champ concerne, cause, statut, plan.

## Acceptees (par conception)

- **`@PG`** : notre sortie emet `ID:bwa-mem3 PN:bwa-mem3 VN:<ver> CL:<notre argv>`, l'oracle emet
  `bwa-mem2`. Exclu du gate d'octet-identite (on compare `@SQ` + lignes d'alignement). Decision finale
  (spoof eventuel de `@PG`) reportee en fin de projet.

## En cours

- **`XS:i` (score sous-optimal `sub`) sur ~7% des reads SE.** Sur 5000 reads wgsim (chr20 2 Mb),
  la sortie est **byte-identique sur 4658/5000 lignes** ; les FLAG/RNAME/POS/CIGAR sont à 100%,
  MAPQ à 99.5%. Les divergences restantes sont **uniquement sur `XS`** (dans les deux sens), et
  viennent de deux sources non encore portées :
  - la branche de fusion **`mem_patch_reg`** de `mem_sort_dedup_patch` (fusion de régions séparées
    par un gap), qui modifie l'ensemble des régions candidates et donc `sub`/`csub` ;
  - les tie-breaks exacts des tris (`ks_introsort` non stable vs tri stable Rust) sur régions/SMEM
    de score égal.
  **Impact** : nul sur position/CIGAR/MAPQ des cas testés ; seul le tag `XS` diffère. **Plan** :
  porter `mem_patch_reg` puis auditer les tie-breaks.
