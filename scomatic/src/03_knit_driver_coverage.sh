#!/bin/bash
# cd /lustre/scratch125/casm/team268im/at31/RA_som_mut/scomatic/ ; . ~/.bashrc ; mamba activate jupy ; /software/team205/bin/jsub lsf -q week -n driver_coverage -c 1 -m 20g -l log "bash src/03_knit_driver_coverage.sh" | bsub

# wd
wd=/lustre/scratch125/casm/team268im/at31/RA_som_mut/scomatic/
cd $wd

# conda activate an environment which has docopt installed so jsub works
. ~/.bashrc
mamba activate jupy

# knit
Rscript -e "rmarkdown::render('reports/driver_coverage.Rmd', output_file = 'driver_coverage.html', output_dir = 'reports/', params = list(rerun = T)"