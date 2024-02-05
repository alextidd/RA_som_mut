# libraries
library(magrittr)

# dirs
wd <- 
  ifelse(Sys.info()['nodename'] == "mib119104s", '~/Volumes/', '') |>
  paste0('/lustre/scratch125/casm/team268im/at31/RA_som_mut/infercnv/') 
setwd(wd)

# Zhang 2023 ---
dat_dir <- '/lustre/scratch126/casm/team268im/mp34/scRNAseq_data/RA_ZhangEtal2023/'
cellranger_dir <- paste0(dat_dir, 'cellranger_output/')
processed_dir <- paste0(dat_dir, 'processed_output/')

# annotations TSVs (cell celltype) ----

# load annotations
annotations <-
  readRDS(paste0(processed_dir, 'all_cells_reference.rds'))$meta_data %>%
  dplyr::select(cell, 
                celltype = cell_type, 
                cluster = cluster_name, 
                sample_id = sample) 

# all celltypes
dir.create('data/Zhang2023/annotations/celltypes/', recursive = T)
annotations %>%
  {split(., .$sample_id)} %>%
  purrr::walk2(.x = names(.), .y = ., function(sample_id, df) {
      df %>%
          dplyr::select(cell, celltype) %>%
          readr::write_tsv(
            paste0('data/Zhang2023/annotations/celltypes/', sample_id, '.tsv'),
            col_names = F
          )
  })

# wo immune celltypes
dir.create('data/Zhang2023/annotations/celltypes_wo_immune/', recursive = T)
immune_cts <- c('NK', 'B cell/plasma cell', 'Myeloid cell', 'T cell')
annotations %>%
  dplyr::filter(!(celltype %in% immune_cts)) %>%
  {split(., .$sample_id)} %>%
  purrr::walk2(.x = names(.), .y = ., function(sample_id, df) {
    df %>%
      dplyr::select(cell, celltype) %>%
      readr::write_tsv(
        paste0('data/Zhang2023/annotations/celltypes_wo_immune/', sample_id, '.tsv'),
        col_names = F
      )
  })

# clusters
dir.create('data/Zhang2023/annotations/clusters')
annotations %>%
  {split(., .$sample_id)} %>%
  purrr::walk2(.x = names(.), .y = ., function(sample_id, df) {
    df %>%
      dplyr::select(cell, cluster) %>%
      readr::write_tsv(
        paste0('data/Zhang2023/annotations/clusters/', sample_id, '.tsv'),
        col_names = F
      )
  })

# mappings CSV (sample_id,raw_counts_matrix,annotations,id) ----

# all celltypes
dir.create('data/Zhang2023/mappings')
list('celltypes', 'celltypes_wo_immune', 'clusters') %>% 
  purrr::walk(function(i) {
    
    # get files
    files <-
      list.files(
        cellranger_dir,
        pattern = '^BRI\\-',
        include.dirs = T) %>%
        tibble::tibble(sample_id = .)  %>%
      # write file paths
      dplyr::transmute(
        sample_id,
        raw_counts_matrix = paste0(processed_dir, 'raw_mRNA_count_matrix.rds'),
        annotations = paste0(wd, '/data/Zhang2023/annotations/', i, '/', sample_id, '.tsv'),
        id = sample_id)
    
    # check files exist
    missing <- 
      files %>%
      tidyr::pivot_longer(c('raw_counts_matrix', 'annotations')) %>%
      dplyr::filter(!file.exists(value))
    if(nrow(missing) > 0) {
      message('Files missing!\n', paste(missing$value, collapse = '\n'))
    }
    
    # write files
    files %>%
      dplyr::filter(
        file.exists(annotations),
        file.exists(raw_counts_matrix)
      ) %>%
      readr::write_csv(
        paste0('data/Zhang2023/mappings/', i, '.csv')
      )
    
  })

# Zhang 2019 ---
# mp34 previously ran inferCNV on RA samples from the Zhang 2019 paper. 
# Here we will run this smaller dataset again as a test.
dir.create('data/Zhang2019/annotations/', recursive = T)

dat_dir_2019 <- '/lustre/scratch126/casm/team268im/mp34/analysis/synovium_scRNA/Celseq_wo_normal/'
annots_2019 <-
    file.path(dat_dir_2019, 'cellAnnotations_Celseq_wo_normal.txt')
matrix_2019 <-
    file.path(dat_dir_2019, 'sc.10x.counts_Celseq_wo_normal.RData')
load(matrix_2019)

# mappings CSV ----
mtx_file <- paste0(wd, 'data/Zhang2019/raw_counts_matrix.tsv')
saveRDS(mat.matrix.sparse, mtx_file)
tibble::tibble(sample_id = 'all',
               raw_counts_matrix = mtx_file,
               annotations = annots_2019,
               id = 'all') %>%
    readr::write_csv('data/Zhang2019/mappings.csv')
