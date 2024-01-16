# Zhange et al. 2019 dataset  (AMP Phase 1)
/lustre/scratch126/casm/team268im/mp34/analysis/synovium_scRNA/SeuratCelSeq_final_seurat_obj.rds

# Zhang et al. 2023 dataset
/lustre/scratch126/casm/team268im/mp34/scRNAseq_data/RA_ZhangEtal2023/
  |- processed_output/  # annotated reference objects for each cell type (as rds files)
  |- cellranger_output/ # cellranger output files per patient sample

# Moritz's inferCNV analysis
/lustre/scratch126/casm/team268im/mp34/analysis/synovium_scRNA/
  |- infercnv_analysis.R      # script with inferCNV command
  |- infercnv_prepare_input.R # script preparing matrix and seurat object to be able to run inferCNV
  |- bsub_inferCNV.sh         # submission script for running inferCNV on farm5

# Moritz's coverage analysis
/lustre/scratch126/casm/team268im/mp34/scRNAseq_data/RA_ZhangEtal2023/
  |- src/    # post-processing
  |- plots/  # visualisation
