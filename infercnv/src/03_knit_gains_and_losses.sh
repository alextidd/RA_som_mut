#!/bin/bash
# cd /lustre/scratch125/casm/team268im/at31/RA_som_mut/infercnv/ ; . ~/.bashrc ; mamba activate jupy ; /software/team205/bin/jsub lsf -q week -n knit_gains_and_losses -c 1 -m 70g -l log "bash src/03_knit_gains_and_losses.sh" | bsub

# wd
wd=/lustre/scratch125/casm/team268im/at31/RA_som_mut/infercnv/
cd $wd

# conda activate an environment which has docopt installed so jsub works
. ~/.bashrc
mamba activate jupy

# knit by_celltype
mkdir -p out/Zhang2023/by_celltype/summary/
Rscript -e "rmarkdown::render('reports/gains_and_losses.Rmd', params = list(sce = 'data/Zhang2023/sce/sce.rds', cell_group_name = 'cell_type', infercnv_dir = 'out/Zhang2023/by_celltype/', cache_dir = 'out/Zhang2023/by_celltype/summary/gains_and_losses_cache/', rerun = F), output_file = 'gains_and_losses.html', output_dir = 'out/Zhang2023/by_celltype/summary/')"

# knit by_stromal_cluster
mkdir -p out/Zhang2023/by_stromal_cluster/summary/
Rscript -e "rmarkdown::render('reports/gains_and_losses.Rmd', params = list(sce = 'data/Zhang2023/sce/Stromal_cell_sce.rds', cell_group_name = 'cluster_name', infercnv_dir = 'out/Zhang2023/by_stromal_cluster/', cache_dir = 'out/Zhang2023/by_stromal_cluster/summary/gains_and_losses_cache/', rerun = F), output_file = 'gains_and_losses.html', output_dir = 'out/Zhang2023/by_stromal_cluster/summary/')"

