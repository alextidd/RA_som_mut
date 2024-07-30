#!/bin/bash
# cd /lustre/scratch125/casm/team268im/at31/RA_som_mut/infercnv/ ; bsub -q week -M 2000 -R 'select[mem>=2000] span[hosts=1] rusage[mem=2000]' -J infercnv_by_celltype -o log "bash src/01_run_infercnv.sh"

# get singularity module
module load singularity

# run infercnv by celltype
mkdir -p out/Zhang2023/by_celltype_test/

# generate annotations, differentiate mural cells and fibroblasts within stromal compartment
cat data/Zhang2023/annotations.tsv |
awk -F'\t' '{
    if (NR == 1) {
      print $1, $2, "id"
    } else if ($2 == "Stromal cell") {
        if ($3 ~ /^F-/) {
            print $1, "Fibroblast", $4
        } else if ($3 ~ /^M-/) {
            print $1, "Mural cell", $4
        } 
    } else {
        print $1, $2, $4
    }
}' OFS='\t' \
> out/Zhang2023/by_celltype_test/annotations.tsv

# generate mappings for all
cat data/Zhang2023/mappings.csv |
awk -F, '{ if (NR==1) { print } ; if (NR>1) {print $1",out/Zhang2023/sce/counts.rds"}}' \
> out/Zhang2023/by_celltype_test/mappings.csv

# run infercnv
nextflow run nf-infercnv \
  --out_dir out/Zhang2023/by_celltype_test/ \
  --mappings out/Zhang2023/by_celltype_test/mappings.csv \
  --annotations out/Zhang2023/by_celltype_test/annotations.tsv \
  --annotation_col celltype \
  --analysis_mode subclusters \
  --cluster_by_groups FALSE \
  --cutoff 0.1 \
  --window_length 151 \
  --noise_logistic FALSE \
  --denoise TRUE \
  --HMM TRUE \
  --HMM_transition_prob 0.000001 \
  --HMM_report_by subcluster \
  --sd_amplifier 0.65 \
  --useRaster FALSE \
  -c config/infercnv.config \
  -resume \
  -N at31@sanger.ac.uk
