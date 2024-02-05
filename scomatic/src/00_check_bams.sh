#!/bin/bash

for bam in /lustre/scratch126/casm/team268im/mp34/scRNAseq_data/RA_ZhangEtal2023/cellranger_output/*/possorted_genome_bam.bam ; do
  echo $bam
  java -jar ~/bin/picard.jar ValidateSamFile \
    I=$bam \
    MODE=SUMMARY
  echo
done > out/Zhang2023/check_bams.txt