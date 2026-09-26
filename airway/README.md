# Real-data run: dexamethasone response in airway smooth muscle cells

This runs the pipeline on the classic **airway** RNA-seq dataset (GEO [GSE52778](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE52778), Himes et al. 2014, *PLoS One*): primary human airway smooth muscle cells from **4 donors**, each **untreated or treated with 1 µM dexamethasone** for 18 h. That gives 8 paired-end samples. The biology is well understood, so the results can be checked against known glucocorticoid-response genes.

## Design

The four donor cell lines (N61311, N052611, N080611, N061011) differ a lot at baseline, so the model accounts for donor:

```
~ cell + condition        (condition: control = reference, dex)
```

With a plain `~ condition` model the donor differences end up in the residual variance and many real dex-response genes are lost.

## Setup (once)

1. **Work inside WSL's own filesystem, not on the Windows drive.** Nextflow is very slow on `/mnt/c/...` and file locking can fail there:
```bash
   mkdir -p ~/projects && cd ~/projects
   git clone https://github.com/barbavegeta/RNA-seq_Nextflow_Pipeline_with_Docker.git rnaseq && cd rnaseq
```
2. **Java 17+ for Nextflow.** Conda environments often carry an older Java that Nextflow picks up. Use a small dedicated environment:
```bash
   conda create -n nextflow -c conda-forge -c bioconda nextflow openjdk=17
   conda activate nextflow
   java -version        # must say 17 or later
```
3. **Docker images** (Docker Desktop with WSL integration switched on):
```bash
   docker build -t rnaseq-tools:0.1.0 containers/rnaseq-tools
   docker build -t rnaseq-r:0.1.0 containers/rnaseq-r
   docker build -t rnaseq-salmon:0.1.0 containers/rnaseq-salmon
```
4. **Memory for WSL.** Only needed for the STAR route. By default WSL sees only half of the laptop's RAM. The STAR index needs about 16 GB with the low-memory setting. Create `C:\Users\<you>\.wslconfig` with
```
   [wsl2]
   memory=13GB
   swap=20GB
```
   (adjust to about 75% of your RAM), then run `wsl --shutdown` in PowerShell and reopen Ubuntu.

`bash airway/preflight.sh` checks all of this and tells you what is still missing.

## Run

```bash
bash airway/download_reference.sh   # ~1 GB; resumes automatically if the connection drops
bash airway/download_fastq.sh       # first 3M read pairs per sample (~2-3 GB); FULL=1 for whole files (~20 GB)
bash airway/run_airway.sh           # preflight, pipeline, known-gene check
```

Both download scripts can simply be re-run after an interruption: finished files are skipped, partial ones resume. The subset mode takes the *first* reads of each file rather than a random sample, which is fine here because read order in these FASTQs is not related to which gene a read comes from.

**STAR or Salmon.** The pipeline can quantify expression in two ways, chosen with `ALIGNER`:

| | `ALIGNER=salmon` (default) | `ALIGNER=star` |
|---|---|---|
| Method | transcript quantification (selective alignment) | genome alignment + featureCounts |
| Memory | ~4 GB | ~16 GB free for the low-memory index; 32 GB machine recommended |
| Time (3M read pairs) | index ~20-30 min, ~2 min per sample | index 1-2 h+, ~10 min per sample |
| Output | gene counts summed from transcript estimates (as tximport) | gene counts from reads over exons |

Both give gene-level counts for the same DESeq2 model, and the dexamethasone genes should come out either way. On a 16 GB laptop use Salmon; to try STAR on a bigger machine run `ALIGNER=star bash airway/run_airway.sh`.

## Results (3M read pairs per sample, Salmon route)

`airway/check_known_genes.py` runs at the end of `run_airway.sh`. It adds gene symbols and checks eight well-known dexamethasone-induced genes. All eight were recovered (log2FC > 1, padj < 0.05):

| Gene | log2FC | padj |
|---|---|---|
| FKBP5 | 3.79 | 9.6e-21 |
| TSC22D3 (GILZ) | 3.18 | 2.7e-10 |
| ZBTB16 | 7.36 | 9.6e-19 |
| PER1 | 3.17 | 1.2e-28 |
| KLF15 | 4.44 | 6.5e-24 |
| DUSP1 | 2.99 | 3.4e-81 |
| CRISPLD2 | 2.69 | 3.0e-39 |
| SPARCL1 | 4.64 | 2.7e-38 |

Salmon mapped 93.9–95.3% of fragments in every sample. 615 genes pass padj < 0.05 and |log2FC| > 1 (350 up, 265 down). The whole run took 9.5 minutes on a 16 GB laptop, with FastQC and Cutadapt cached from an earlier attempt. The "65,083 genes tested" in the log counts every Ensembl gene, including the roughly 35,000 with no reads, which DESeq2 excludes from the adjusted p-values.

Kept in `airway/results/`: `known_gene_check.tsv`, `top20_genes.tsv` and `salmon_mapping.tsv`. The full outputs go to `results_airway/`: counts, DESeq2 results with symbols, volcano, PCA and MA plots, and MultiQC. That folder is not committed.

In the PCA, samples separate by treatment on PC1 (42%), and PC2 mainly separates donor N080611 from the others.

## Pipeline changes made for real data

- `--deseq2_formula` and `--deseq2_reference`: paired/blocked designs, and an explicit reference level. Previously the reference was whichever level came first alphabetically, so the sign of log2FC could silently flip.
- `--star_index_args` / `--sjdb_overhang`: low-memory human index.
- `--count_read_pairs` (default on): featureCounts now counts fragments, not individual reads, for paired-end data, which is what DESeq2 expects.
