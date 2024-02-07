#!/bin/bash
# cd /lustre/scratch125/casm/team268im/at31/RA_som_mut/scomatic/ ; . ~/.bashrc ; mamba activate jupy ; /software/team205/bin/jsub lsf -q week -n knit_driver_coverage -c 1 -m 8g -l log "bash src/03_knit_driver_coverage.sh" | bsub

# wd
wd=/lustre/scratch125/casm/team268im/at31/RA_som_mut/scomatic/
cd $wd

# conda activate an environment which has docopt installed so jsub works
. ~/.bashrc
mamba activate jupy

# # knit
# Rscript -e "rmarkdown::render('reports/driver_coverage.Rmd', output_file = 'driver_coverage.html', output_dir = 'reports/', params = list(rerun = F))"

# # create mappings 
# bams_dir=/lustre/scratch126/casm/team268im/mp34/scRNAseq_data/RA_ZhangEtal2023/cellranger_output/
# ct_bams_dir=/lustre/scratch125/casm/team268im/at31/RA_som_mut/scomatic/out/Zhang2023/default_thresholds/
# out_dir=/lustre/scratch125/casm/team268im/at31/RA_som_mut/scomatic/out/Zhang2023/coverage/
# (
#   echo 'id,bam,celltype' ;
#   while read -r id ; do
#     # check if ct bams directory exists
#     if [ -d $ct_bams_dir/${id/-/_}/ ] ; then 
#       ct_bams=($(ls $ct_bams_dir/${id/-/_}/celltype_bams/*.bam))
#       # check if scomatic successfully ran
#       if [ ${#ct_bams[*]} -gt 0 ] ; then
#         echo "$id,$bams_dir/$id/possorted_genome_bam.bam,NA" ;
#         for file in ${ct_bams[@]} ; do
#           ct=$(basename $file | cut -d. -f2)
#           echo "$id,$file,$ct" ;
#         done
#       fi
#     fi
#   done < \
#   <(sed 1d /lustre/scratch125/casm/team268im/at31/RA_som_mut/infercnv/data/Zhang2023/mappings.csv | cut -d, -f1 | sort -u)
# ) | cat > data/Zhang2023/coverage_mappings.csv

# # test mappings
# cat data/Zhang2023/coverage_mappings.csv | head -10 > data/Zhang2023/coverage_mappings_test.csv

# run
/software/team205/nextflow-23.04.1-all run nextflow/driver_coverage.nf \
    --mappings data/Zhang2023/coverage_mappings.csv \
    --drivers /lustre/scratch125/casm/team268im/at31/RA_som_mut/scomatic/data/driver_genes/driver_gene_coords_for_coverage.tsv \
    --ref_cds /lustre/scratch125/casm/team268im/at31/RA_som_mut/scomatic/data/dndscv/ref_cds.rds \
    --out_dir /lustre/scratch125/casm/team268im/at31/RA_som_mut/scomatic/out/Zhang2023/coverage/ \
    -c /nfs/team205/kp9/nextflow/scomatic/LSF.config  \
    --location local \
    -w work/ \
    -resume
