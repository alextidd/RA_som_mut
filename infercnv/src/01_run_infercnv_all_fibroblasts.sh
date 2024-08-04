#!/bin/bash
# cd /lustre/scratch125/casm/team268im/at31/RA_som_mut/infercnv/ ; bsub -q week -M 2000 -R 'select[mem>=2000] span[hosts=1] rusage[mem=2000]' -J infercnv_all -o log "bash src/01_run_infercnv_all_fibroblasts.sh"

# run infercnv on celltypes, running all individuals together
# TODO: decide on cluster_by_groups value
# arg definition from inferCNV - If observations are defined according to groups 
# (ie. patients), each group of cells will be clustered separately. 
# (default = FALSE, instead will use k_obs_groups setting)
# run infercnv on fibroblast clusters, running all individuals together

# get singularity module
module load singularity

# run infercnv on all fibroblasts, running all individuals together
# must partition work subdirectories so concomitant jobs don't cause lock issues
mkdir -p out/Zhang2023/all_fibroblasts/ work/all_fibroblasts/

# generate annotations, filter to only fibroblasts
cat data/Zhang2023/annotations.tsv |
awk -F'\t' -v OFS="\t" \
  '{if (NR == 1) {print $1,$2,$3,"id"} ; if(($2 == "Stromal cell") && $3 ~ /^F-/) {print $1,$2,$3,"all"} }' \
> out/Zhang2023/all_fibroblasts/annotations.tsv

# generate mappings for all
echo -e "id,raw_counts_matrix\nall,out/Zhang2023/sce/counts.rds" \
> out/Zhang2023/all_fibroblasts/mappings.csv

# run infercnv 
# (parameters recommended by mp34@sanger.ac.uk)
# (/lustre/scratch126/casm/team268im/mp34/analysis/synovium_scRNA/infercnv_analysis.R)
nextflow run nf-infercnv \
  --out_dir out/Zhang2023/all_fibroblasts/ \
  -w work/all_fibroblasts/ \
  --mappings out/Zhang2023/all_fibroblasts/mappings.csv \
  --annotations out/Zhang2023/all_fibroblasts/annotations.tsv \
  --annotation_col cluster \
  --analysis_mode subclusters \
  --cluster_by_groups TRUE \
  --cutoff 0.1 \
  --window_length 151 \
  --noise_logistic FALSE \
  --denoise TRUE \
  --HMM TRUE \
  --HMM_transition_prob 0.000001 \
  --HMM_report_by subcluster \
  --sd_amplifier 0.65 \
  --useRaster FALSE \
  --output_format pdf \
  -c config/infercnv.config \
  -resume \
  -N at31@sanger.ac.uk

# # run infercnv on fibroblast clusters, running all individuals together, 
# # clustered by fibroblast subtype
# cat data/Zhang2023/annotations.tsv |
# awk -F'\t' -v OFS="\t" \
#   '{if (NR == 1) {print $1,$2,$3,"id"} ; if(($2 == "Stromal cell") && $3 ~ /^F-/) {print $1,$2,$3,$3} }' \
# > out/Zhang2023/all_fibroblasts/annotations.tsv

# cat data/Zhang2023/mappings.csv |
# awk -F',' 'BEGIN{OFS=",";} NR == 1 {print} ; NR == 2 {print "all",$2}' \
# > data/Zhang2023/stromal_mappings_run_all.csv

# nextflow run nextflow/infercnv.nf \
#   --out_dir out/Zhang2023/by_stromal_cluster_run_all/ \
#   --mappings data/Zhang2023/stromal_mappings_run_all.csv \
#   --annotations data/Zhang2023/stromal_annotations_run_all.tsv \
#   --annotation_col 'cluster' \
#   --analysis_mode 'samples' \
#   --cluster_by_groups 'FALSE' \
#   -c config/infercnv.config \
#   -c /nfs/team205/kp9/nextflow/scomatic/LSF.config  \
#   -w work/ \
#   -resume

# # run infercnv on celltypes
# /software/team205//nextflow-23.04.1-all run nextflow/infercnv.nf \
#   --out_dir out/Zhang2023/by_celltype/ \
#   --mappings data/Zhang2023/mappings.csv \
#   --annotations data/Zhang2023/annotations.tsv \
#   --annotation_col 'celltype' \
#   -c config/infercnv.config \
#   -c /nfs/team205/kp9/nextflow/scomatic/LSF.config  \
#   -w work/ \
#   -resume
    
# # run infercnv on clusters
# nextflow run nextflow/infercnv.nf \
#   --out_dir out/Zhang2023/by_cluster/ \
#   --mappings data/Zhang2023/mappings.csv \
#   --annotations data/Zhang2023/annotations.tsv \
#   --annotation_col 'cluster' \
#   -c config/infercnv.config \
#   -c /nfs/team205/kp9/nextflow/scomatic/LSF.config  \
#   -w work/ \
#   -resume

# # run infercnv on stromal clusters
# cat data/Zhang2023/annotations.tsv |
# awk -F'\t' '{if ((NR == 1) || ($2 == "Stromal cell")) {print} }' \
# > data/Zhang2023/stromal_annotations.tsv
# nextflow run nextflow/infercnv.nf \
#   --out_dir out/Zhang2023/by_stromal_cluster/ \
#   --mappings data/Zhang2023/mappings.csv \
#   --annotations data/Zhang2023/stromal_annotations.tsv \
#   --annotation_col 'cluster' \
#   -c config/infercnv.config \
#   -c /nfs/team205/kp9/nextflow/scomatic/LSF.config  \
#   -w work/ \
#   -resume
  
# # run infercnv on stromal clusters, running all individuals together
# cat data/Zhang2023/stromal_annotations.tsv |
# awk -F'\t' 'BEGIN{OFS="\t";} NR == 1 {print} ; NR > 1 {print $1,$2,$3,"all"}' \
# > data/Zhang2023/stromal_annotations_run_all.tsv

# cat data/Zhang2023/mappings.csv |
# awk -F',' 'BEGIN{OFS=",";} NR == 1 {print} ; NR == 2 {print "all",$2}' \
# > data/Zhang2023/stromal_mappings_run_all.csv

# nextflow run nextflow/infercnv.nf \
#   --out_dir out/Zhang2023/by_stromal_cluster_run_all/ \
#   --mappings data/Zhang2023/stromal_mappings_run_all.csv \
#   --annotations data/Zhang2023/stromal_annotations_run_all.tsv \
#   --annotation_col 'cluster' \
#   --analysis_mode 'samples' \
#   --cluster_by_groups 'FALSE' \
#   -c config/infercnv.config \
#   -c /nfs/team205/kp9/nextflow/scomatic/LSF.config  \
#   -w work/ \
#   -resume

# # run test
# nextflow run nextflow/infercnv.nf \
#   --out_dir out/test/ \
#   --mappings <((head -1 data/Zhang2023/mappings.csv; grep 'BRI-527\|BRI-581' data/Zhang2023/mappings.csv)| cat) \
#   --annotations data/Zhang2023/annotations.tsv \
#   --annotation_col 'celltype' \
#   -c config/infercnv.config \
#   -c /nfs/team205/kp9/nextflow/scomatic/LSF.config  \
#   -w work/test/
