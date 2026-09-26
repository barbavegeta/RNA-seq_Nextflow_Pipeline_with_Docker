#!/usr/bin/env bash
# Checks the things that commonly stop the airway run. Run from the repository root.
ok=1
say() { printf '%-28s %s\n' "$1" "$2"; }

jv=$(java -version 2>&1 | grep -m1 ' version ' | sed -E 's/.*version "([0-9]+).*/\1/')
if [ -n "$jv" ] && [ "$jv" -ge 17 ] 2>/dev/null; then say "Java" "OK ($jv)"; else
  say "Java" "PROBLEM: need 17+, found '${jv:-none}' ($(command -v java))"; ok=0; fi
command -v nextflow >/dev/null && say "Nextflow" "OK" || { say "Nextflow" "PROBLEM: not on PATH"; ok=0; }
if docker info >/dev/null 2>&1; then
  for img in rnaseq-tools:0.1.0 rnaseq-r:0.1.0 $( [ "${ALIGNER:-salmon}" = salmon ] && echo rnaseq-salmon:0.1.0 ); do
    docker image inspect "$img" >/dev/null 2>&1 && say "image $img" "OK" || { say "image $img" "PROBLEM: not built"; ok=0; }
  done
else say "Docker" "PROBLEM: not running (enable WSL integration in Docker Desktop)"; ok=0; fi
[ -f main.nf ] && [ -d containers ] && say "location" "OK (repo root)" || { say "location" "PROBLEM: run from the repository root"; ok=0; }
case "$PWD" in /mnt/*) say "filesystem" "WARNING: on the Windows drive ($PWD); move the repo under ~ for speed";; *) say "filesystem" "OK";; esac
mem=$(awk '/MemTotal/ {printf "%d", $2/1024/1024}' /proc/meminfo)
if [ "${ALIGNER:-salmon}" = salmon ]; then
  [ "$mem" -ge 6 ] && say "RAM" "OK (${mem} GB; Salmon needs ~4 GB)" || { say "RAM" "PROBLEM: ${mem} GB visible to WSL"; ok=0; }
elif [ "$mem" -ge 30 ]; then say "RAM" "OK (${mem} GB, full STAR index possible)";
elif [ "$mem" -ge 12 ]; then say "RAM" "OK-ish (${mem} GB: low-memory STAR index, relies on swap; Salmon is safer)";
else say "RAM" "PROBLEM: ${mem} GB is too little for a STAR index; use ALIGNER=salmon"; ok=0; fi
if [ "${ALIGNER:-salmon}" = salmon ]; then REFFILE=data/airway/ref/Homo_sapiens.GRCh38.transcripts.fa.gz
else REFFILE=data/airway/ref/Homo_sapiens.GRCh38.dna.primary_assembly.fa; fi
[ -s "$REFFILE" ] && say "reference" "OK ($(basename "$REFFILE"))" || { say "reference" "missing: bash airway/download_reference.sh"; ok=0; }
n=$(ls data/airway/fastq/*_1.fastq.gz 2>/dev/null | wc -l)
[ "$n" -eq 8 ] && say "FASTQ" "OK (8 samples)" || { say "FASTQ" "found $n/8: bash airway/download_fastq.sh"; ok=0; }
[ $ok = 1 ] && echo "All checks passed." || echo "Fix the PROBLEM lines, then run: bash airway/run_airway.sh"
