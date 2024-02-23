# libraries
library(magrittr)

# dirs
wd <- 
  ifelse(Sys.info()['nodename'] == "mib119104s", '~/Volumes/', '') %>% 
  paste0('/lustre/scratch125/casm/team268im/at31/RA_som_mut/scomatic/') 
setwd(wd)
dat_dir <- '/lustre/scratch126/casm/team268im/mp34/scRNAseq_data/RA_ZhangEtal2023/'
cellranger_dir <- paste0(dat_dir, 'cellranger_output/')
processed_dir <- paste0(dat_dir, 'processed_output/')
dir.create('data/')

# Zhang 2023 ----
# mappings CSV (sample_id,bam_file,id) -----
dir.create('data/Zhang2023/')
mappings <-
  list.files(
      cellranger_dir,
      pattern = 'possorted_genome_bam.bam$',
      recursive = T, full.names = T) %>%
  tibble::tibble(bam_file = .) %>%
  # get parent directory (= id)
  dplyr::mutate(
    bam_file = normalizePath(bam_file),
    sample_id = bam_file %>% dirname() %>% basename() %>% gsub('-', '_', .),
    id = sample_id
  ) 

# check for truncation
truncated_bams <- 
  readr::read_tsv('data/Zhang2023/truncated_bams.tsv')

# filter out
mappings <-
  mappings %>%
  dplyr::filter(!bam_file %in% truncated_bams$bam_file)

mappings %>%
    readr::write_csv('data/Zhang2023/mappings.csv')

# celltypes TSV (Index Cell_type) -----
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
    readr::write_tsv('data/Zhang2023/celltypes.tsv')

# celltype counts
dat %>%
    dplyr::select(sample_id = sample,
                  id = sample, 
                  celltype = cell_type) %>%
    dplyr::group_by(sample_id, id, celltype) %>%
    dplyr::count(name = 'n_cells_per_id_per_celltype') %>%
    readr::write_tsv('data/Zhang2023/celltype_counts.tsv')
