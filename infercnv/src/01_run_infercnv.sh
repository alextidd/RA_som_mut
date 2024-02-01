#!/bin/bash
# cd /lustre/scratch125/casm/team268im/at31/RA_som_mut/infercnv/ ; . ~/.bashrc ; mamba activate jupy ; /software/team205/bin/jsub lsf -q week -n infercnv -c 1 -m 2g -l log "bash src/01_run_infercnv.sh" | bsub

set -eo pipefail

# wd
wd=/lustre/scratch125/casm/team268im/at31/RA_som_mut/infercnv/
cd $wd
mkdir -p work/ out/ log/

# conda activate an environment which has docopt installed so jsub works
. ~/.bashrc
mamba activate jupy

# # run infercnv
# /software/team205/nextflow-23.04.1-all run nextflow/infercnv.nf \
#     -c config/Zhang2023.config \
#     -c config/infercnv.config \
#     -c /nfs/team205/kp9/nextflow/scomatic/LSF.config  \
#     -w work/ \
#     -resume
    
# run infercnv wo immune
/software/team205/nextflow-23.04.1-all run nextflow/infercnv.nf \
    -c config/Zhang2023_wo_immune.config \
    -c config/infercnv.config \
    -c /nfs/team205/kp9/nextflow/scomatic/LSF.config  \
    -w work/ \
    -resume