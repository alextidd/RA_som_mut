# Investigating somatic mutations and copy number alterations in scRNA-seq data from rheumatoid arthritis

## Main aims

1. Evaluate the cell type-specific coverage across each of the genes of interest.
2. Call known somatic mutations from tNanoSeq in the scRNA-seq data.
3. Perform de novo mutation calling in the different cell types (SComatic / cellSNP).
4. Perform copy number calling in fibroblasts and different fine-grained cell states in 2023 dataset to validate chr7 finding in 2019 dataset.
5. Are there differences in the transcriptomes of cells with and without specific mutations?

## Plan

* Check coverage of genes of interest in the dataset
* Run inferCNV at patient level on the dataset, using the average expression profile across all cells as the 'normal'.
* Run SComatic on the dataset and look for validation of tNanoSeq hotspots in the dataset.

## Resources

### Zhang et al. 2019 dataset  (AMP Phase 1)
```
/lustre/scratch126/casm/team268im/mp34/analysis/synovium_scRNA/SeuratCelSeq_final_seurat_obj.rds
```

### Zhang et al. 2023 dataset
```
/lustre/scratch126/casm/team268im/mp34/scRNAseq_data/RA_ZhangEtal2023/
  |- processed_output/  # annotated reference objects for each cell type (as rds files)
  |- cellranger_output/ # cellranger output files per patient sample
```

### Moritz's inferCNV analysis
```
/lustre/scratch126/casm/team268im/mp34/analysis/synovium_scRNA/
  |- infercnv_analysis.R      # script with inferCNV command
  |- infercnv_prepare_input.R # script preparing matrix and seurat object to be able to run inferCNV
  |- bsub_inferCNV.sh         # submission script for running inferCNV on farm5
```

### Moritz's coverage analysis
```
/lustre/scratch126/casm/team268im/mp34/scRNAseq_data/RA_ZhangEtal2023/
  |- src/    # post-processing
  |- plots/  # visualisation
```
