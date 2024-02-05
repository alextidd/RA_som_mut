#!/bin/bash

java -jar ~/bin/picard.jar ValidateSamFile \
      I=/lustre/scratch126/casm/team268im/mp34/scRNAseq_data/RA_ZhangEtal2023/cellranger_output/BRI-419/possorted_genome_bam.bam \
      MODE=SUMMARY