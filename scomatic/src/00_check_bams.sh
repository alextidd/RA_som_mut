#!/bin/bash
# cd /lustre/scratch125/casm/team268im/at31/RA_som_mut/scomatic/ ; . ~/.bashrc ; mamba activate jupy ; /software/team205/bin/jsub lsf -q week -n check_bams -c 1 -m 12g -l log "bash src/00_check_bams.sh" | bsub

wd=/lustre/scratch125/casm/team268im/at31/RA_som_mut/scomatic/
cd $wd

# env
. ~/.bashrc
mamba activate jupy

for bam in /lustre/scratch126/casm/team268im/mp34/scRNAseq_data/RA_ZhangEtal2023/cellranger_output/*/possorted_genome_bam.bam ; do
  echo $bam
  java -jar ~/bin/picard.jar ValidateSamFile I=$bam MODE=SUMMARY 2>&1 
  echo
done > out/Zhang2023/check_bams.txt