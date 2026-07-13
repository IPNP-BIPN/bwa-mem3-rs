#!/usr/bin/env bash
# Build the end-to-end test data under work/ : a small reference slice, its bwa-mem2 index, and
# reproducible simulated reads (wgsim, fixed seed, truth positions encoded in read names).
#
# Env overrides: REF (bgzipped FASTA), REGION (samtools region), NREADS.
set -euo pipefail
cd "$(dirname "$0")/.."

REF="${REF:-/Users/benjamin/Database/VEP_database/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz}"
REGION="${REGION:-20:2000000-4000000}"   # mid-contig, no telomeric N runs
NREADS="${NREADS:-5000}"
W=work
mkdir -p "$W"

echo "[extract] $REGION from $REF"
samtools faidx "$REF" "$REGION" > "$W/region.fa"

echo "[index] bwa-mem2 index"
bwa-mem2 index "$W/region.fa" 2>/dev/null

echo "[wgsim] $NREADS pairs (seed 11)"
wgsim -S 11 -N "$NREADS" -1 150 -2 150 -e 0.005 -r 0.001 \
  "$W/region.fa" "$W/r1.fq" "$W/r2.fq" >/dev/null 2>&1

echo "testdata ready in $W/ (region $REGION, $NREADS read pairs)"
