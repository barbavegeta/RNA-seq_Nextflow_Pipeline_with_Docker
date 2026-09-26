#!/usr/bin/env bash
# Ensembl GRCh38 annotation, genome (STAR route) and transcript sequences (Salmon route).
# STAR needs the genome and GTF uncompressed; Salmon reads the gzipped transcripts directly.
# Safe to re-run: interrupted downloads resume where they stopped.
set -euo pipefail
REL=${ENSEMBL_RELEASE:-112}
OUT=data/airway/ref
mkdir -p "$OUT"

fetch() {  # fetch URL FILE : resume + retry until the file is complete and passes gzip -t
  local url=$1 f=$2 n=0
  until curl -fL --retry 20 --retry-all-errors --retry-delay 10 --connect-timeout 30 \
             -C - -o "$f" "$url" && gzip -t "$f" 2>/dev/null; do
    n=$((n+1)); [ $n -ge 30 ] && { echo "giving up on $url" >&2; return 1; }
    if ! gzip -t "$f" 2>/dev/null && [ -s "$f" ]; then echo "  incomplete, resuming (attempt $n)"; fi
    sleep 15
  done
}

FA=Homo_sapiens.GRCh38.dna.primary_assembly.fa
GTF=Homo_sapiens.GRCh38.${REL}.gtf
BASE=https://ftp.ensembl.org/pub/release-${REL}
[ -s "$OUT/$FA" ]  || { fetch "$BASE/fasta/homo_sapiens/dna/$FA.gz" "$OUT/$FA.gz" && gunzip -f "$OUT/$FA.gz"; }
[ -s "$OUT/$GTF" ] || { fetch "$BASE/gtf/homo_sapiens/$GTF.gz" "$OUT/$GTF.gz" && gunzip -f "$OUT/$GTF.gz"; }

# Transcript sequences for the low-memory Salmon route (protein-coding cDNA + non-coding RNA, ~100 MB)
TX=Homo_sapiens.GRCh38.transcripts.fa.gz
if [ ! -s "$OUT/$TX" ]; then
  fetch "$BASE/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz" "$OUT/cdna.fa.gz"
  fetch "$BASE/fasta/homo_sapiens/ncrna/Homo_sapiens.GRCh38.ncrna.fa.gz" "$OUT/ncrna.fa.gz"
  cat "$OUT/cdna.fa.gz" "$OUT/ncrna.fa.gz" > "$OUT/$TX" && rm "$OUT/cdna.fa.gz" "$OUT/ncrna.fa.gz"
fi
ls -lh "$OUT"
