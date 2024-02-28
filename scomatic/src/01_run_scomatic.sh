#!/bin/bash
# cd /lustre/scratch125/casm/team268im/at31/RA_som_mut/scomatic/ ; . ~/.bashrc ; mamba activate jupy ; /software/team205/bin/jsub lsf -q week -n scomatic -c 1 -m 2g -l log "bash src/01_run_scomatic.sh" | bsub

set -eo pipefail

# wd
wd=/lustre/scratch125/casm/team268im/at31/RA_som_mut/scomatic/
cd $wd
config_dir=/nfs/team205/kp9/nextflow/scomatic/
mkdir -p work/ log/

# conda activate an environment which has docopt installed so jsub works
. ~/.bashrc
mamba activate jupy

# use scomatic to get celltype bams for driver coverage
# set all thresholds such that no cells are removed and we get full unfiltered
# celltype-specific bams
/software/team205/nextflow-23.04.1-all run nextflow/scomatic.nf \
    -entry STEP1 \
    --publish_celltype_bams true \
    --output_dir out/Zhang2023/coverage/ \
    --max_nM 1000000 \
    --max_NH 1000000 \
    --min_MQ 0 \
    --n_trim 0 \
    --min_dp 0 \
    --min_cc 0 \
    --projectDir $PWD \
    -params-file $config_dir/GEX.json \
    -c config/Zhang2023.config \
    -c $config_dir/LSF.config \
    -w work/ \
    --location local

# # run scomatic
# /software/team205/nextflow-23.04.1-all run nextflow/scomatic.nf \
#     --projectDir $PWD \
#     -params-file $config_dir/GEX.json \
#     -c config/Zhang2023.config \
#     -c $config_dir/LSF.config \
#     -w work/ \
#     --output_dir out/Zhang2023/default_thresholds/ \
#     --location local \
#     -resume
# 
# # get single cell genotypes
# /software/team205/nextflow-23.04.1-all run nextflow/scomatic.nf \
#     -entry genotypes \
#     --projectDir $PWD \
#     -params-file $config_dir/GEX.json \
#     -c config/Zhang2023.config \
#     -c $config_dir/LSF.config \
#     -w work/ \
#     --output_dir out/Zhang2023/default_thresholds/ \
#     --location local \
#     -resume
