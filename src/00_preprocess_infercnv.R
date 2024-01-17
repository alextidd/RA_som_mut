library(magrittr)
setwd('/lustre/scratch125/casm/team268im/at31/RA_som_mut')
dir.create('data/infercnv/annotations', recursive = T)

readRDS('/lustre/scratch126/casm/team268im/mp34/scRNAseq_data/RA_ZhangEtal2023/processed_output/all_cells_reference.rds')$meta_data %>%
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

