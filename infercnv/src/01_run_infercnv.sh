#!/bin/bash
# cd /lustre/scratch125/casm/team268im/at31/RA_som_mut/infercnv/ ; . ~/.bashrc ; mamba activate jupy ; ~/bin/jsub lsf -q week -n infercnv -c 1 -m 2g -l log "bash src/01_run_infercnv.sh" | bsub

set -eo pipefail

# wd
wd=/lustre/scratch125/casm/team268im/at31/RA_som_mut/infercnv/
cd $wd
mkdir -p work/ out/ log/

# conda activate an environment which has docopt installed so jsub works
. ~/.bashrc
mamba activate jupy

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
# /software/team205/nextflow-23.04.1-all run nextflow/infercnv.nf \
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
# /software/team205/nextflow-23.04.1-all run nextflow/infercnv.nf \
#   --out_dir $wd/out/Zhang2023/by_stromal_cluster/ \
#   --mappings $wd/data/Zhang2023/mappings.csv \
#   --annotations $wd/data/Zhang2023/stromal_annotations.tsv \
#   --annotation_col 'cluster' \
#   -c config/infercnv.config \
#   -c /nfs/team205/kp9/nextflow/scomatic/LSF.config  \
#   -w work/ \
#   -resume
  
# run infercnv on stromal clusters, running all individuals together
cat $wd/data/Zhang2023/stromal_annotations.tsv |
awk -F'\t' 'BEGIN{OFS="\t";} NR == 1 {print} ; NR > 1 {print $1,$2"_"$4,$2"_"$4,"all"}' \
> $wd/data/Zhang2023/stromal_annotations_run_all.tsv
/software/team205/nextflow-23.04.1-all run nextflow/infercnv.nf \
  --out_dir $wd/out/Zhang2023/by_stromal_cluster_run_all/ \
  --mappings $wd/data/Zhang2023/mappings.csv \
  --annotations $wd/data/Zhang2023/stromal_annotations_run_all.tsv \
  --annotation_col 'cluster' \
  -c config/infercnv.config \
  -c /nfs/team205/kp9/nextflow/scomatic/LSF.config  \
  -w work/ \
  -resume

# # run test
# /software/team205/nextflow-23.04.1-all run nextflow/infercnv.nf \
#   --out_dir $wd/out/test/ \
#   --mappings <((head -1 $wd/data/Zhang2023/mappings.csv; grep 'BRI-527\|BRI-581' $wd/data/Zhang2023/mappings.csv)| cat) \
#   --annotations $wd/data/Zhang2023/annotations.tsv \
#   --annotation_col 'celltype' \
#   -c config/infercnv.config \
#   -c /nfs/team205/kp9/nextflow/scomatic/LSF.config  \
#   -w work/test/
