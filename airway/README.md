# Real-data run: dexamethasone response in airway smooth muscle cells

This runs the pipeline on the classic **airway** RNA-seq dataset (GEO [GSE52778](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE52778), Himes et al. 2014, *PLoS One*): primary human airway smooth muscle cells from **4 donors**, each **untreated or treated with 1 µM dexamethasone** for 18 h. That gives 8 paired-end samples. The biology is well understood, so the results can be checked against known glucocorticoid-response genes.

## Design

The four donor cell lines (N61311, N052611, N080611, N061011) differ a lot at baseline, so the model accounts for donor:

```
~ cell + condition        (condition: control = reference, dex)
```

With a plain `~ condition` model the donor differences end up in the residual variance and many real dex-response genes are lost.

## Setup (once)

1. **Work inside WSL's own filesystem, not on the Windows drive.** Nextflow and STAR are very slow on `/mnt/c/...` and file locking can fail there:
   ```bash
   mkdir -p ~/projects && cp -r "/mnt/c/Users/<you>/Downloads/RNA-seq_Nextflow_Pipeline_with_Docker-main" ~/projects/rnaseq
   cd ~/projects/rnaseq
   ```
2. **Put the airway files into the repository root** (not a subfolder), overwriting `main.nf`, `nextflow.config` and `bin/deseq2.R`:
   ```bash
   cp -r rnaseq-airway-update/airway rnaseq-airway-update/main.nf rnaseq-airway-update/nextflow.config .
   cp rnaseq-airway-update/bin/deseq2.R bin/
   ```
3. **Java 17+ for Nextflow.** Conda environments often carry an older Java that Nextflow picks up. Use a small dedicated environment:
   ```bash
   conda create -n nextflow -c conda-forge -c bioconda nextflow openjdk=17
   conda activate nextflow
   java -version        # must say 17 or later
   ```
4. **Docker images** (Docker Desktop with WSL integration switched on):
   ```bash
   docker build -t rnaseq-tools:0.1.0 containers/rnaseq-tools
   docker build -t rnaseq-r:0.1.0 containers/rnaseq-r
   ```
5. **Memory for WSL.** By default WSL sees only half of the laptop's RAM. The STAR index needs about 16 GB with the low-memory setting. Create `C:\Users\<you>\.wslconfig` with
   ```
   [wsl2]
   memory=14GB
   swap=16GB
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

Both give gene-level counts for the same DESeq2 model, and the dexamethasone genes should come out either way. On a 16 GB laptop use Salmon; to try STAR on a bigger machine run `ALIGNER=star bash airway/run_airway.sh`. Salmon needs its own small image: `docker build -t rnaseq-salmon:0.1.0 containers/rnaseq-salmon`.

## What to expect

`airway/check_known_genes.py` (run at the end of `run_airway.sh`) adds gene symbols and checks eight well-known dexamethasone-induced genes: **FKBP5, TSC22D3 (GILZ), ZBTB16, PER1, KLF15, DUSP1, CRISPLD2, SPARCL1**. All should be strongly up-regulated (log2FC > 1, padj < 0.05). CRISPLD2 was the headline finding of the original paper. The script also writes:

- `results_airway/volcano_airway.png`: volcano plot with the known genes labelled (put this at the top of the README)
- `results_airway/top20_genes.tsv`: most significant genes with symbols
- `results_airway/known_gene_check.tsv`: the pass/fail table
- the pipeline's own MultiQC report, PCA and MA plots

In the PCA you should see samples separate by treatment on one axis and by donor on the other.

## Pipeline changes made for real data

- `--deseq2_formula` and `--deseq2_reference`: paired/blocked designs, and an explicit reference level. Previously the reference was whichever level came first alphabetically, so the sign of log2FC could silently flip.
- `--star_index_args` / `--sjdb_overhang`: low-memory human index.
- `--count_read_pairs` (default on): featureCounts now counts fragments, not individual reads, for paired-end data, which is what DESeq2 expects.
