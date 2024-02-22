# Zhang et al. 2023 performed CITE-seq on rheumatoid arthritis samples.
# The processed scRNA-seq output from this was saved to various files.
# Here I will construct these outputs into SingleCellExperiment objects.

# libraries
library(magrittr)
library(SummarizedExperiment)
library(SingleCellExperiment)

# dirs
wd <- 
  ifelse(Sys.info()['nodename'] == "mib119104s", '~/Volumes/', '') |>
  paste0('/lustre/scratch125/casm/team268im/at31/RA_som_mut/infercnv/') 
setwd(wd)
data_dir <- '/lustre/scratch126/casm/team268im/mp34/scRNAseq_data/RA_ZhangEtal2023/processed_output/'
out_dir <- 'data/Zhang2023/sce/'
dir.create(out_dir)

# function: convert Zhang 2023 object to sce
get_sce <- function(obj, assays) {
  
  # get cell order from metadata
  cell_order <- obj$meta_data$cell
  
  # filter and reorder to cells that pass qc from counts
  assays_pass <-
    assays %>%
    purrr::map(~ .x[rownames(assays$logcounts), cell_order])
  
  # populate sce
  sce <- 
    SingleCellExperiment::SingleCellExperiment(
      assays = assays_pass,
      colData = obj$meta_data
    )
  
  # name, reorder, and add umap embeddings
  reducedDim(sce, 'UMAP_uwot') <- obj$umap$embedding
  
  # return
  sce
}

# all cells ----

# load object
obj <- 
  readRDS(file.path(data_dir, 'all_cells_reference.rds'))

# load rna data
assays <-
  list(
    counts = 
      readRDS(file.path(data_dir, 'raw_mRNA_count_matrix.rds')),
    logcounts = 
      readRDS(file.path(data_dir, 'qc_mRNA_314011cells_log_normalized_matrix.rds')))

# generate sce
sce <-
  get_sce(obj, assays)

# save
saveRDS(sce, 'data/Zhang2023/sce/sce.rds')

# stromal cells ----

# load obj
obj <-
  readRDS(file.path(data_dir, 'fibroblast_reference.rds'))

# fix cell type names (Fibroblast -> Stromal cell), because it also includes 
# mural cells
obj$meta_data$cell_type <- 'Stromal cell'

# generate sce
sce <-
  get_sce(obj, assays)

# save
saveRDS(sce, 'data/Zhang2023/sce/Stromal_cell_sce.rds')

# fibroblasts ----

# subset to fibroblasts only
sce <-
  sce[, sce$cluster_name != 'Mu-0: Mural']

# save
saveRDS(sce, 'data/Zhang2023/sce/Fibroblast_sce.rds')
