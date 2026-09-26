#!/usr/bin/env bash
# Airway RNA-seq (GEO GSE52778, Himes et al. 2014) from ENA.
#
# Default: stream only the first SUBSET_READS read pairs of each sample (3 million),
# so ~2-3 GB is downloaded instead of ~20 GB. The dexamethasone response is strong
# enough to detect at this depth.   FULL=1 bash airway/download_fastq.sh  -> whole files.
# Safe to re-run: finished files are skipped, broken ones are retried.
set -euo pipefail
OUT=data/airway/fastq
SUBSET_READS=${SUBSET_READS:-3000000}
mkdir -p "$OUT"

urls_for() {
  curl -fsS --retry 10 --retry-all-errors \
    "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=$1&result=read_run&fields=fastq_ftp&format=tsv" \
    | tail -n 1 | cut -f2 | tr ';' '\n' | grep -E '_[12]\.fastq\.gz$'   # paired mates only
}

stream_subset() {  # URL FILE : first N read pairs; a shorter file is kept whole
  local url=$1 f=$2 want=$((4 * SUBSET_READS)) n=0
  while :; do
    set +o pipefail   # head closes the pipe early on purpose
    curl -fsSL --connect-timeout 30 "$url" 2>/dev/null | gzip -dc 2>/dev/null \
      | head -n "$want" | gzip -1 > "$f.part"
    local st=("${PIPESTATUS[@]}")
    set -o pipefail
    got=$(gzip -dc "$f.part" | wc -l)
    # enough reads, or curl and gunzip both finished cleanly (the file is simply shorter)
    if [ "$got" -eq "$want" ] || { [ "${st[0]}" -eq 0 ] && [ "${st[1]}" -eq 0 ] && [ "$got" -gt 0 ]; }; then
      mv "$f.part" "$f"; [ "$got" -lt "$want" ] && echo "  file has only $((got / 4)) reads; kept all"
      return 0
    fi
    n=$((n+1)); echo "  got $got/$want lines, retrying ($n)"; [ $n -ge 10 ] && return 1; sleep 15
  done
}

fetch_full() {
  local url=$1 f=$2 n=0
  until curl -fL --retry 20 --retry-all-errors --retry-delay 10 -C - -o "$f" "$url" && gzip -t "$f"; do
    n=$((n+1)); [ $n -ge 30 ] && return 1; echo "  resuming $f ($n)"; sleep 15
  done
}

for run in $(tail -n +2 airway/samples_airway.tsv | cut -f1); do
  for u in $(urls_for "$run"); do
    f="$OUT/$(basename "$u")"
    [ -s "$f" ] && gzip -t "$f" 2>/dev/null && { echo "have $f"; continue; }
    echo "$run -> $f"
    if [ "${FULL:-0}" = 1 ]; then fetch_full "https://$u" "$f"; else stream_subset "https://$u" "$f"; fi
  done
done
ls -lh "$OUT"
