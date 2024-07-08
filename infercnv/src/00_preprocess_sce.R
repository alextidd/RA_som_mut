# Zhang et al. 2023 performed CITE-seq on rheumatoid arthritis samples.
# The processed scRNA-seq output from this was saved to various files.
# Here I will construct these outputs into SingleCellExperiment objects.

# libraries
library(magrittr)
library(SummarizedExperiment)
library(SingleCellExperiment)

# dirs
data_dir <- "/lustre/scratch126/casm/team268im/mp34/scRNAseq_data/RA_ZhangEtal2023/processed_output/"
out_dir <- "out/Zhang2023/sce/"

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
  reducedDim(sce, "UMAP_uwot") <- obj$umap$embedding

  # return
  sce
}

# all cells ----

# load object
obj <- readRDS(file.path(data_dir, "all_cells_reference.rds"))

# load rna data (raw + normalised)
counts_rds <- file.path(data_dir, "raw_mRNA_count_matrix.rds")
logcounts_rds <- file.path(data_dir, "qc_mRNA_314011cells_log_normalized_matrix.rds")
assays <- list(counts = readRDS(counts_rds), logcounts = readRDS(logcounts_rds))

# generate sce
sce <- get_sce(obj, assays)

# save sce
saveRDS(sce, file.path(out_dir, "sce.rds"))

# save QC-passing raw counts
saveRDS(sce@assays@data$counts, file.path(out_dir, "counts.rds"))

# stromal cells ----

# load obj
obj <- readRDS(file.path(data_dir, "fibroblast_reference.rds"))

# fix cell type names (Fibroblast -> Stromal cell), because it also includes
# mural cells
obj$meta_data$cell_type <- "Stromal cell"

# generate sce
sce <-
  get_sce(obj, assays)

# save
saveRDS(sce, "data/Zhang2023/sce/Stromal_cell_sce.rds")

# fibroblasts ----

# subset to fibroblasts only
sce <- sce[, sce$cluster_name != "Mu-0: Mural"]

# save
saveRDS(sce, "data/Zhang2023/sce/Fibroblast_sce.rds")