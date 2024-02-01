#!/bin/bash
# cd /lustre/scratch125/casm/team268im/at31/RA_som_mut/scomatic/ ; . ~/.bashrc ; mamba activate jupy ; /software/team205/bin/jsub lsf -q week -n knit_driver_coverage -c 1 -m 30g -l log "bash src/03_knit_driver_coverage.sh" | bsub

# wd
wd=/lustre/scratch125/casm/team268im/at31/RA_som_mut/scomatic/
cd $wd

# conda activate an environment which has docopt installed so jsub works
. ~/.bashrc
mamba activate jupy

# knit
Rscript -e "rmarkdown::render('reports/driver_coverage.Rmd', output_file = 'driver_coverage.html', output_dir = 'reports/', params = list(rerun = F))"


# run scomatic
(
  echo 'id,bam,out_dir' ;
  cat data/Zhang2023/mappings.csv |
  sed 1d | head -3 | cut -d, -f 1,2 |
  awk -F'/' '{print $0",/lustre/scratch125/casm/team268im/at31/RA_som_mut/scomatic/out/Zhang2023/coverage/"$10"/"}' ;
) | cat > data/test_mappings.csv

/software/team205/nextflow-23.04.1-all run nextflow/driver_coverage.nf \
    --mappings data/test_mappings.csv \
    --drivers /lustre/scratch125/casm/team268im/at31/RA_som_mut/scomatic/data/driver_genes/driver_gene_coords_for_coverage.tsv \
    --location local \
    -w work/ \
    -resume
