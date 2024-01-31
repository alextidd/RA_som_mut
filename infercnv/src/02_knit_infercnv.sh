#!/bin/bash
# cd /lustre/scratch125/casm/team268im/at31/RA_som_mut/infercnv/ ; . ~/.bashrc ; mamba activate jupy ; /software/team205/bin/jsub lsf -q week -n knit_infercnv -c 1 -m 20g -l log "bash src/02_knit_infercnv.sh" | bsub

# wd
wd=/lustre/scratch125/casm/team268im/at31/RA_som_mut/infercnv/
cd $wd

# conda activate an environment which has docopt installed so jsub works
. ~/.bashrc
mamba activate jupy

# knit
Rscript -e "rmarkdown::render('reports/infercnv.Rmd', output_file = 'infercnv.html', output_dir = 'reports/', params = list(rerun = F))"
