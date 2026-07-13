# Divergences connues vs oracle bwa-mem2 2.3

Suivi de la traine de parite. Chaque entree : champ concerne, cause, statut, plan.

## Acceptees (par conception)

- **`@PG`** : notre sortie emet `ID:bwa-mem3 PN:bwa-mem3 VN:<ver> CL:<notre argv>`, l'oracle emet
  `bwa-mem2`. Exclu du gate d'octet-identite (on compare `@SQ` + lignes d'alignement). Decision finale
  (spoof eventuel de `@PG`) reportee en fin de projet.

## En cours

_(rien pour l'instant : phase 0)_
