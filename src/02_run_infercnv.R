#!/bin/bash
# cd /lustre/scratch125/casm/team268im/at31/RA_som_mut/ ; . ~/.bashrc ; mamba activate jupy ; /software/team205/bin/jsub lsf -q week -n infercnv -c 1 -m 2g -l log "bash src/02_run_infercnv.sh" | bsub

set -eo pipefail

# wd
wd=/lustre/scratch125/casm/team268im/at31/RA_som_mut/
cd $wd
mkdir -p work/ out/infercnv/ log/

# conda activate an environment which has docopt installed so jsub works
. ~/.bashrc
mamba activate jupy

# run infercnv
/software/team205/nextflow-23.04.1-all run nextflow/infercnv.nf \
    --output_dir out/infercnv \
    --mappings data/infercnv/mappings.csv \
    --gene_order_file data/infercnv/hg38_gencode_v27.txt \
    -c /nfs/team205/kp9/nextflow/scomatic/LSF.config  \
    -w work/ \
    -resume
    
