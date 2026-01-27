# RNA-seq Nextflow Pipeline (Docker)

Minimal end-to-end RNA-seq pipeline:
FASTQ -> FastQC -> cutadapt -> STAR -> sorted BAM -> featureCounts -> counts.tsv

## Quickstart (test data)

```bash
bash tests/get_test_data.sh

docker build -t rnaseq-tools:0.1.0 containers/rnaseq-tools
docker build -t rnaseq-r:0.1.0 containers/rnaseq-r

nextflow run main.nf -profile docker \
  --samplesheet tests/samplesheet.test.csv \
  --genome_fasta tests/ref/genome.fasta \
  --genome_gtf tests/ref/genes.gtf \
  --outdir results_test
```
