library(magrittr)

# dirs
wd <- '/lustre/scratch125/casm/team268im/at31/RA_som_mut'
setwd(wd)
dat_dir <- '/lustre/scratch126/casm/team268im/mp34/scRNAseq_data/RA_ZhangEtal2023/'
cellranger_dir <- paste0(dat_dir, 'cellranger_output/')
processed_dir <- paste0(dat_dir, 'processed_output/')
dir.create('data/infercnv/annotations', recursive = T)

# mappings CSV (sample_id,raw_counts_matrix,annotations,id)
list.files(
        cellranger_dir,
        pattern = '^BRI\\-',
        include.dirs = T) %>%
    tibble::enframe(value = 'sample_id') %>%
    dplyr::transmute(
        sample_id,
        raw_counts_matrix = paste0(processed_dir, 'raw_mRNA_count_matrix.rds'),
	annotations = paste0(wd, '/data/infercnv/annotations/', sample_id, '.tsv'),
        id = sample_id) %>%
    readr::write_csv('data/infercnv/mappings.csv')

# annotations TSVs (cell celltype) ----
dat <-
    readRDS(paste0(processed_dir, 'all_cells_reference.rds'))$meta_data
dat %>%
    dplyr::select(cell, cell_type, sample) %>%
    dplyr::group_by(sample) %>%
    dplyr::group_split() %>%
    purrr::walk(function(df) {
        sample <- unique(df$sample)
	df %>%
            dplyr::select(cell, cell_type) %>%
	    readr::write_tsv(
		paste0('data/infercnv/annotations/', sample, '.tsv'),
		col_names = F
	    )
        })

