# libraries
library(magrittr)

# dirs
setwd('/lustre/scratch125/casm/team268im/at31/RA_som_mut')
dat_dir <- '/lustre/scratch126/casm/team268im/mp34/scRNAseq_data/RA_ZhangEtal2023/'
cellranger_dir <- paste0(dat_dir, 'cellranger_output/')
processed_dir <- paste0(dat_dir, 'processed_output/')

# mappings CSV (sample_id,bam_file,id)
list.files(
        cellranger_dir,
        pattern = '^BRI\\-',
        include.dirs = T) %>%
    tibble::enframe(value = 'path_id') %>%
    dplyr::transmute(
        sample_id = gsub('-', '_', path_id),
        bam_file = paste0(cellranger_dir, path_id, '/possorted_genome_bam.bam'),
        id = gsub('-', '_', path_id)) %>%
    readr::write_csv('data/mappings.csv')

# celltypes TSV (Index Cell_type) ----
dat <-
    readRDS(paste0(processed_dir, 'all_cells_reference.rds'))$meta_data
aliases <-
    dat %>%
    dplyr::transmute(
        sample_id = sample %>% gsub('-', '_', .),
        # fix indices - must be formatted SAMPLE_BARCODE-1
        cell = cell %>% gsub('-', '_', .) %>% paste0('-1'),
        # fix cell types - only letters, numbers and underscores
        celltype = cell_type,
        celltype_alias = cell_type %>% gsub(' |/', '_', .))
aliases %>%
    dplyr::select(Index = cell, Cell_type = celltype_alias) %>%
    readr::write_tsv('data/celltypes.tsv')