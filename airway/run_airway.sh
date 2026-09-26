#!/usr/bin/env bash
# Full airway run.
#   ALIGNER=salmon (default): transcript quantification, ~4 GB RAM, laptop-friendly
#   ALIGNER=star:             genome alignment, needs ~16 GB free RAM for the index (32 GB machine recommended)
set -euo pipefail
ALIGNER=${ALIGNER:-salmon}
bash airway/preflight.sh || true
REL=${ENSEMBL_RELEASE:-112}
REF=data/airway/ref
if [ "$ALIGNER" = star ]; then
  ALIGN_ARGS=(--aligner star --genome_fasta $REF/Homo_sapiens.GRCh38.dna.primary_assembly.fa
              --star_index_args '--genomeSAsparseD 3 --limitGenomeGenerateRAM 12000000000')
else
  ALIGN_ARGS=(--aligner salmon --transcript_fasta $REF/Homo_sapiens.GRCh38.transcripts.fa.gz)
fi
nextflow run main.nf -profile docker \
  --samplesheet airway/samplesheet_airway.csv \
  --design airway/design_airway.tsv \
  --deseq2_formula '~ cell + condition' \
  --deseq2_reference control \
  --genome_gtf $REF/Homo_sapiens.GRCh38.${REL}.gtf \
  "${ALIGN_ARGS[@]}" \
  --threads 4 \
  --outdir results_airway \
  -with-report results_airway/pipeline_report.html \
  -with-timeline results_airway/timeline.html \
  -resume
python airway/check_known_genes.py \
  --results results_airway/deseq2_results.tsv \
  --gtf $REF/Homo_sapiens.GRCh38.${REL}.gtf \
  --out results_airway
