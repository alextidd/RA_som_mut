#!/bin/bash
# cd /lustre/scratch125/casm/team268im/at31/RA_som_mut/scomatic/ ; . ~/.bashrc ; mamba activate jupy ; /software/team205/bin/jsub lsf -q week -n targeted_mut_calling -c 1 -m 8g -l log "bash src/05_run_targeted_mutation_calling.sh" | bsub

# wd
wd=/lustre/scratch125/casm/team268im/at31/RA_som_mut/scomatic/
out_dir=${wd}out/Zhang2023/targeted_mutation_calling/
mkdir -p $out_dir
cd $wd

# conda activate an environment which has docopt installed so jsub works
. ~/.bashrc
mamba activate jupy

# get mappings of celltype-specific bams only
cat ${wd}out/Zhang2023/coverage/mappings.csv | 
awk -F',' -v OFS=',' '$3 != "NA" {print}' \
> ${out_dir}/mappings.csv

# # run target_mutation_calling.nf on celltype bams only
# /software/team205/nextflow-23.04.1-all run nextflow/targeted_mutation_calling.nf \
#   --mappings ${out_dir}/mappings.csv \
#   --window 5 \
#   --out_dir ${out_dir} \
#   -c config/LSF.config \
#   -w work/ \
#   -resume

# run on all impactful TNFRSF14 mutations
/software/team205/nextflow-23.04.1-all run nextflow/targeted_mutation_calling.nf \
  --mappings ${out_dir}/mappings.csv \
  --recursites $wd/data/driver_genes/TNFRSF14_withingene_sites_GRCh38.tsv \
  --window 1 \
  --out_dir ${out_dir}/TNFRSF14/ \
  -c config/LSF.config \
  -w work/ 
  
  