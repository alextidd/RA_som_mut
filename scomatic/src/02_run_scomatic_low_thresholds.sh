#!/bin/bash
# cd /lustre/scratch125/casm/team268im/at31/RA_som_mut/scomatic/ ; . ~/.bashrc ; mamba activate jupy ; /software/team205/bin/jsub lsf -q week -n scomatic_lo -c 1 -m 2g -l log "bash src/02_run_scomatic_low_thresholds.sh" | bsub

set -eo pipefail

# wd
wd=/lustre/scratch125/casm/team268im/at31/RA_som_mut/scomatic/
cd $wd
config_dir=/nfs/team205/kp9/nextflow/scomatic/
mkdir -p work/ out/Zhang2023/default_thresholds/ log/

# conda activate an environment which has docopt installed so jsub works
. ~/.bashrc
mamba activate jupy

# run scomatic
/software/team205/nextflow-23.04.1-all run nextflow/scomatic.nf \
    --projectDir $PWD \
    --min_dp 1 \
    --output_dir '/lustre/scratch125/casm/team268im/at31/RA_som_mut/scomatic/out/Zhang2023/low_thresholds/' \
    -params-file $config_dir/GEX.json \
    -c config/Zhang2023.config \
    -c $config_dir/LSF.config \
    -w work/ \
    --location local \
    -resume

# get single cell genotypes
/software/team205/nextflow-23.04.1-all run nextflow/scomatic.nf \
    -entry genotypes \
    --projectDir $PWD \
    --min_dp 1 \
    --output_dir '/lustre/scratch125/casm/team268im/at31/RA_som_mut/scomatic/out/Zhang2023/low_thresholds/' \
    -params-file $config_dir/GEX.json \
    -c config/Zhang2023.config \
    -c $config_dir/LSF.config \
    -w work/ \
    --location local \
    -resume
