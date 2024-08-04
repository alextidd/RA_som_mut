#!/bin/bash
# cd /lustre/scratch125/casm/team268im/at31/RA_som_mut/infercnv/ ; . ~/.bashrc ; mamba activate jupy ; /software/team205/bin/jsub lsf -q week -n knit_infercnv -c 1 -m 20g -l log "bash src/02_knit_infercnv.sh" | bsub

# wd
wd=/lustre/scratch125/casm/team268im/at31/RA_som_mut/infercnv/
cd $wd

# conda activate an environment which has docopt installed so jsub works
. ~/.bashrc
mamba activate jupy

# # knit by_celltype
# mkdir out/Zhang2023/by_celltype/summary/
# Rscript -e "rmarkdown::render('reports/infercnv.Rmd', params = list(infercnv_dir = 'out/Zhang2023/by_celltype/', rerun = F), output_file = 'infercnv.html', output_dir = 'out/Zhang2023/by_celltype/summary/')"
# 
# # knit by_cluster
# mkdir out/Zhang2023/by_cluster/summary/
# Rscript -e "rmarkdown::render('reports/infercnv.Rmd', params = list(infercnv_dir = 'out/Zhang2023/by_cluster/', rerun = F), output_file = 'infercnv.html', output_dir = 'out/Zhang2023/by_cluster/summary/')"

# knit by_stromal_cluster
mkdir out/Zhang2023/by_stromal_cluster/summary/
Rscript -e "rmarkdown::render('reports/infercnv.Rmd', params = list(infercnv_dir = 'out/Zhang2023/by_stromal_cluster/', rerun = F), output_file = 'infercnv.html', output_dir = 'out/Zhang2023/by_stromal_cluster/summary/')"
