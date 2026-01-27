#!/bin/bash -ue
cutadapt -j 4       -q 20 -m 20       -o KO_REP1_trim_R1.fastq.gz -p KO_REP1_trim_R2.fastq.gz       SRR6357072_1.fastq.gz SRR6357072_2.fastq.gz > KO_REP1.cutadapt.log
