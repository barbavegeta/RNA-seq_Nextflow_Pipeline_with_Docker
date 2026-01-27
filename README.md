# RNA-seq Nextflow Pipeline (Docker)

Minimal end-to-end RNA-seq pipeline:
**FASTQ -> FastQC -> cutadapt -> STAR -> sorted BAM -> featureCounts -> MultiQC → DESeq2**

This repository includes a small test setup so you can run an end-to-end check locally.

---

## Requirements

- **Nextflow** (tested with `25.10.3`)
- **Docker**
- **Git**

---

## Quick start (test run)

From the repo root:

```bash
# (Optional) build the local containers (only needed if you use local Docker builds)
docker build -t rnaseq-tools:0.1.0 containers/rnaseq-tools
docker build -t rnaseq-r:0.1.0 containers/rnaseq-r

# run the test
nextflow run main.nf \
  --samplesheet tests/samplesheet.test.csv \
  --genome_fasta tests/ref/genome.fasta \
  --genome_gtf tests/ref/genes.gtf \
  --design tests/design.test.tsv \
  --outdir results_test \
  -with-trace  "trace-$(date +%Y%m%d-%H%M%S).txt" \
  -with-report "report-$(date +%Y%m%d-%H%M%S).html" \
  -with-timeline "timeline-$(date +%Y%m%d-%H%M%S).html"
```

### Note on report/timeline filenames
If you reuse the same `report.html` / `timeline.html`, Nextflow will refuse to overwrite them.  
Using timestamped names (as above) avoids that.

---

## Inputs

### Sample sheet (`--samplesheet`)
CSV with columns:

- `sample` — unique sample ID
- `fastq_1` — path to R1 FASTQ.gz
- `fastq_2` — path to R2 FASTQ.gz

Example: `tests/samplesheet.test.csv`

### Design matrix (`--design`)
TSV with columns:

- `sample`
- `condition` (e.g. `WT` / `KO`)

Example: `tests/design.test.tsv`

### Reference
- `--genome_fasta` — reference genome FASTA
- `--genome_gtf` — annotation GTF

Test reference in: `tests/ref/`

---

## Outputs

Output folder: `--outdir` (example: `results_test/`)

Key outputs:
- `counts.tsv` — gene counts (featureCounts)
- `qc/multiqc_report.html` — aggregated QC report
- `deseq2_results.tsv` — DESeq2 differential expression results
- `deseq2_ma_plot.png` — MA plot
- `deseq2_pca_plot.png` — PCA plot

---

## Repository structure

- `main.nf` — Nextflow pipeline (DSL2)
- `nextflow.config` — pipeline configuration
- `bin/deseq2.R` — DESeq2 script
- `containers/` — container build context(s)
- `tests/` — test samplesheet/design/reference and helper scripts

---

## Notes / limitations

- The test reference and test FASTQs are intentionally small; they are for validating the workflow wiring, not biological interpretation.
- For real datasets, update references, resources (cpus/memory), and output locations accordingly.
