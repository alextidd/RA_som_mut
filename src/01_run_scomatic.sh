#!/bin/bash
# cd /lustre/scratch125/casm/team268im/at31/RA_som_mut/ ; . ~/.bashrc ; mamba activate jupy ; /software/team205/bin/jsub lsf -q week -n scomatic -c 1 -m 2g -l log "bash src/01_run_scomatic.sh" | bsub

set -eo pipefail

# wd
wd=/lustre/scratch125/casm/team268im/at31/RA_som_mut/
cd $wd
config_dir=/nfs/team205/kp9/nextflow/scomatic/
mkdir -p work/ out/scomatic/ log/

# conda activate an environment which has docopt installed so jsub works
. ~/.bashrc
mamba activate jupy

(
    cd out/scomatic/

    # run scomatic
    /software/team205/nextflow-23.04.1-all run \
        /lustre/scratch125/casm/team268im/at31/nextflow/scomatic.nf \
        --projectDir $PWD \
        -params-file $config_dir/GEX.json \
        -c ../../config/run.config \
        -c $config_dir/LSF.config \
        -w $wd/work/ \
        --location local \
        -resume
    
    # get single cell genotypes
    /software/team205/nextflow-23.04.1-all run \
        /lustre/scratch125/casm/team268im/at31/nextflow/scomatic.nf \
        -entry genotypes \
        --projectDir $PWD \
        -params-file $config_dir/GEX.json \
        -c ../../config/run.config \
        -c $config_dir/LSF.config \
        -w $wd/work/ \
        --location local \
        -resume
)