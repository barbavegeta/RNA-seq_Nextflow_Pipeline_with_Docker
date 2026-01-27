#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process FASTQC {
  tag "${sample}"
  cpus 1
  input:
    tuple val(sample), path(r1), path(r2), val(condition)
  output:
    tuple val(sample), path("*_fastqc.zip"), path("*_fastqc.html")
  script:
    """
    fastqc -t ${task.cpus} ${r1} ${r2}
    """
}

process CUTADAPT {
  tag "${sample}"
  cpus params.threads
  input:
    tuple val(sample), path(r1), path(r2), val(condition)
  output:
    tuple val(sample),
          path("${sample}_trim_R1.fastq.gz"),
          path("${sample}_trim_R2.fastq.gz"),
          val(condition),
          path("${sample}.cutadapt.log")
  script:
    """
    cutadapt -j ${task.cpus} \
      -q 20 -m 20 \
      -o ${sample}_trim_R1.fastq.gz -p ${sample}_trim_R2.fastq.gz \
      ${r1} ${r2} > ${sample}.cutadapt.log
    """
}

process STAR_INDEX {
  tag "star_index"
  cpus params.threads
  input:
    path fasta
    path gtf
  output:
    path "star_index"
  script:
    """
    mkdir -p star_index
    STAR --runThreadN ${task.cpus} \
      --runMode genomeGenerate \
      --genomeDir star_index \
      --genomeFastaFiles ${fasta} \
      --sjdbGTFfile ${gtf} \
      --sjdbOverhang 99
    """
}

process STAR_ALIGN {
  tag "${sample}"
  cpus params.threads
  input:
    path index_dir
    tuple val(sample), path(r1), path(r2), val(condition), path(cutlog)
  output:
    tuple val(sample),
          path("${sample}.Aligned.out.bam"),
          val(condition),
          path("${sample}.star.log"),
          path("${sample}.Log.final.out")
  script:
    """
    STAR --runThreadN ${task.cpus} \
      --genomeDir ${index_dir} \
      --readFilesIn ${r1} ${r2} \
      --readFilesCommand zcat \
      --outSAMtype BAM Unsorted \
      --outFileNamePrefix ${sample}. \
      > ${sample}.star.log 2>&1
    """
}

process SORT_INDEX_BAM {
  tag "${sample}"
  cpus 2
  input:
    tuple val(sample), path(bam), val(condition), path(starlog), path(finalout)
  output:
    tuple val(sample), path("${sample}.sorted.bam"), path("${sample}.sorted.bam.bai"), val(condition)
  script:
    """
    samtools sort -@ ${task.cpus} -o ${sample}.sorted.bam ${bam}
    samtools index ${sample}.sorted.bam
    """
}

process FEATURECOUNTS {
  publishDir params.outdir, mode: 'copy'

  cpus 4

  input:
    path genome_gtf
    path bam_map
    path bam_files

  output:
    path "counts.tsv", emit: counts
    path "featureCounts.log", emit: log

  shell:
  '''
  set -euo pipefail

  bam_map="!{bam_map}"
  gtf="!{genome_gtf}"

  mapfile -t bam_files < <(cut -f2 "$bam_map")

  featureCounts \
    -T !{task.cpus} \
    -p -B -C \
    -a "$gtf" \
    -o counts.raw.tsv \
    -t exon \
    -g gene_id \
    "${bam_files[@]}" \
    > featureCounts.log 2>&1

  grep -v '^#' counts.raw.tsv > counts.nohash.tsv

  hdr="gene_id"
  while IFS=$'\t' read -r sid bam; do
    hdr="${hdr}"$'\t'"${sid}"
  done < "$bam_map"
  printf "%s\n" "$hdr" > counts.tsv

  tail -n +2 counts.nohash.tsv | cut -f1,7- >> counts.tsv
  '''
}



process MULTIQC {
  tag "multiqc"
  cpus 1
  publishDir "${params.outdir}/qc", mode: 'copy'
  input:
    path(qc_files)
  output:
    path "multiqc_report.html"
    path "multiqc_data"
  script:
    """
    multiqc -o . ${qc_files.join(' ')}
    """
}

process DESEQ2 {
  publishDir params.outdir, mode: 'copy'

  cpus 1

  input:
    path counts
    path design
    path deseq_script

  output:
    path "deseq2_results.tsv", emit: results
    path "deseq2_ma_plot.png", emit: ma_plot
    path "deseq2_pca_plot.png", emit: pca_plot

  shell:
  '''
  Rscript "!{deseq_script}" \
    --counts "!{counts}" \
    --design "!{design}" \
    --out_results deseq2_results.tsv \
    --out_ma deseq2_ma_plot.png \
    --out_pca deseq2_pca_plot.png
  '''

}



workflow {

  if( !params.samplesheet )  error "Provide --samplesheet"
  if( !params.genome_fasta ) error "Provide --genome_fasta"
  if( !params.genome_gtf )   error "Provide --genome_gtf"

  def samplesheet = file(params.samplesheet)
  if( !samplesheet.exists() ) error "Samplesheet not found: ${samplesheet}"

  def fasta = file(params.genome_fasta)
  def gtf   = file(params.genome_gtf)

  def samples_ch = Channel
    .from(samplesheet)
    .splitCsv(header:true)
    .map { row ->
      def sample = row.sample.toString()
      def r1 = file(row.fastq_1.toString())
      def r2 = file(row.fastq_2.toString())
      def condition = row.containsKey('condition') ? row.condition.toString() : 'NA'
      tuple(sample, r1, r2, condition)
    }

  def fq      = FASTQC(samples_ch)
  def trimmed = CUTADAPT(samples_ch)
  def idx     = STAR_INDEX(fasta, gtf)
  def aligned = STAR_ALIGN(idx, trimmed)
  def bams    = SORT_INDEX_BAM(aligned)

  // keep only sample + bam
  def bam_pairs = bams.map { sample, bam, bai, condition -> tuple(sample, bam) }

  // bam_map.tsv will contain filenames that will exist in the FEATURECOUNTS task dir
  def bam_map = bam_pairs
    .map { sid, bam -> "${sid}\t${bam.getName()}" }
    .collectFile(name: 'bam_map.tsv', newLine: true)

  // stage all BAMs into the FEATURECOUNTS task (portable across containers/executors)
  def bam_files = bam_pairs.map { sid, bam -> bam }.collect()

def fc = FEATURECOUNTS(gtf, bam_map, bam_files)

  def fastqc_files  = fq.flatMap { sample, zips, htmls -> (zips + htmls) }
  def cutadapt_logs = trimmed.map { sample, r1, r2, condition, cutlog -> cutlog }
  def star_final    = aligned.map { sample, bam, condition, starlog, finalout -> finalout }

  def qc_files_ch = fastqc_files
    .mix(cutadapt_logs)
    .mix(star_final)
    .mix(fc.log)
    .collect()

  MULTIQC(qc_files_ch)

  if( params.design && !params.skip_deseq2 ) {
    def design = file(params.design)
    if( !design.exists() ) error "Design file not found: ${design}"

    def deseq_script = file("${projectDir}/bin/deseq2.R")
    if( !deseq_script.exists() ) error "DESeq2 script not found: ${deseq_script}"

    DESEQ2(fc.counts, design, deseq_script)
  }
}