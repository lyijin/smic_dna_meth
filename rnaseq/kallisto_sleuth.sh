#!/bin/bash

# > kallisto_sleuth.sh <
#
# Bash script documents the commands used to run kallisto (index/quant) and
# the sleuth R script.

~/tools/kallisto_linux-v0.44.0/kallisto index -i smic_cds Smic.genome.annotation.CDS.longest.sorted.fa
mkdir results
for a in ../trimmed_reads/*R1.trim.fastq; do b=`echo $a | grep -Po '\d\d\w'` && ~/tools/kallisto_linux-v0.44.0/kallisto quant -i smic_cds -o results/${b} --bias --rf-stranded -t 8 -b 100 ${a} ${a/R1/R2} 2> ${b}.kallisto.log; done; wait

Rscript sleuth_analysis.R
