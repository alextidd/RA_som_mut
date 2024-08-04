# libraries
library(magrittr)
library(Seurat)

# dir
wd <- 
  ifelse(Sys.info()['nodename'] == "mib119104s", '~/Volumes/', '') |> 
  paste0('/lustre/scratch125/casm/team268im/at31/RA_som_mut/infercnv/') 
setwd(wd)

# read in sce, convert to seurat
sce <-
  readRDS('data/Zhang2023/sce/sce.rds')
seu <- as.Seurat(sce)
Idents(seu) <- seu$cell_type
sce <- NULL

# find top markers for each celltype
markers <-
  FindAllMarkers(seu)

# save the markers
markers %>%
  tibble::as_tibble() %>%
  readr::write_tsv('out/Zhang2023/by_celltype/summary/markers/celltype_markers.tsv')