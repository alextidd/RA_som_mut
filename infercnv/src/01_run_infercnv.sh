#!/bin/bash
# cd /lustre/scratch125/casm/team268im/at31/RA_som_mut/infercnv/ ; bsub -q week -M 2000 -R 'select[mem>=2000] span[hosts=1] rusage[mem=2000]' -J infercnv -o log "bash src/01_run_infercnv.sh"

# get singularity module
module load singularity

# wd
wd=/lustre/scratch125/casm/team268im/at31/RA_som_mut/infercnv/
cd $wd
mkdir -p work/ out/ log/

# run infercnv on celltypes, running all individuals together
# TODO: decide on cluster_by_groups 
# arg definition from inferCNV - If observations are defined according to groups 
# (ie. patients), each group of cells will be clustered separately. 
# (default = FALSE, instead will use k_obs_groups setting)
nextflow run nextflow/infercnv.nf \
  --out_dir out/Zhang2023/by_celltype_run_all/ \
  --mappings data/Zhang2023/stromal_mappings_run_all.csv \
  --annotations data/Zhang2023/annotations.tsv \
  --annotation_col "celltype" \
  --analysis_mode "samples" \
  --cluster_by_groups "FALSE" \
  -c config/infercnv.config \
  -c /lustre/scratch125/casm/team268im/at31/RA_som_mut/scomatic/config/LSF.config \
  -w work/ \
  -resume

# # run infercnv on celltypes
# /software/team205//nextflow-23.04.1-all run nextflow/infercnv.nf \
#   --out_dir $wd/out/Zhang2023/by_celltype/ \
#   --mappings $wd/data/Zhang2023/mappings.csv \
#   --annotations $wd/data/Zhang2023/annotations.tsv \
#   --annotation_col 'celltype' \
#   -c config/infercnv.config \
#   -c /nfs/team205/kp9/nextflow/scomatic/LSF.config  \
#   -w work/ \
#   -resume
    
# # run infercnv on clusters
# nextflow run nextflow/infercnv.nf \
#   --out_dir $wd/out/Zhang2023/by_cluster/ \
#   --mappings $wd/data/Zhang2023/mappings.csv \
#   --annotations $wd/data/Zhang2023/annotations.tsv \
#   --annotation_col 'cluster' \
#   -c config/infercnv.config \
#   -c /nfs/team205/kp9/nextflow/scomatic/LSF.config  \
#   -w work/ \
#   -resume

# # run infercnv on stromal clusters
# cat $wd/data/Zhang2023/annotations.tsv |
# awk -F'\t' '{if ((NR == 1) || ($2 == "Stromal cell")) {print} }' \
# > $wd/data/Zhang2023/stromal_annotations.tsv
# nextflow run nextflow/infercnv.nf \
#   --out_dir $wd/out/Zhang2023/by_stromal_cluster/ \
#   --mappings $wd/data/Zhang2023/mappings.csv \
#   --annotations $wd/data/Zhang2023/stromal_annotations.tsv \
#   --annotation_col 'cluster' \
#   -c config/infercnv.config \
#   -c /nfs/team205/kp9/nextflow/scomatic/LSF.config  \
#   -w work/ \
#   -resume
  
# # run infercnv on stromal clusters, running all individuals together
# cat $wd/data/Zhang2023/stromal_annotations.tsv |
# awk -F'\t' 'BEGIN{OFS="\t";} NR == 1 {print} ; NR > 1 {print $1,$2,$3,"all"}' \
# > $wd/data/Zhang2023/stromal_annotations_run_all.tsv

# cat $wd/data/Zhang2023/mappings.csv |
# awk -F',' 'BEGIN{OFS=",";} NR == 1 {print} ; NR == 2 {print "all",$2}' \
# > $wd/data/Zhang2023/stromal_mappings_run_all.csv

# nextflow run nextflow/infercnv.nf \
#   --out_dir $wd/out/Zhang2023/by_stromal_cluster_run_all/ \
#   --mappings $wd/data/Zhang2023/stromal_mappings_run_all.csv \
#   --annotations $wd/data/Zhang2023/stromal_annotations_run_all.tsv \
#   --annotation_col 'cluster' \
#   --analysis_mode 'samples' \
#   --cluster_by_groups 'FALSE' \
#   -c config/infercnv.config \
#   -c /nfs/team205/kp9/nextflow/scomatic/LSF.config  \
#   -w work/ \
#   -resume

# # run test
# nextflow run nextflow/infercnv.nf \
#   --out_dir $wd/out/test/ \
#   --mappings <((head -1 $wd/data/Zhang2023/mappings.csv; grep 'BRI-527\|BRI-581' $wd/data/Zhang2023/mappings.csv)| cat) \
#   --annotations $wd/data/Zhang2023/annotations.tsv \
#   --annotation_col 'celltype' \
#   -c config/infercnv.config \
#   -c /nfs/team205/kp9/nextflow/scomatic/LSF.config  \
#   -w work/test/
