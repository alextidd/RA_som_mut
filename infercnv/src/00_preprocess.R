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

# all celltypes
dir.create('data/Zhang2023/annotations/', recursive = T)
annotations <-
  readRDS(paste0(processed_dir, 'all_cells_reference.rds'))$meta_data %>%
  dplyr::select(cell, celltype = cell_type, sample_id = sample) 
annotations %>%
  {split(., .$sample_id)} %>%
  purrr::walk2(.x = names(.), .y = ., function(sample_id, df) {
      df %>%
          dplyr::select(cell, celltype) %>%
          readr::write_tsv(
            paste0('data/Zhang2023/annotations/', sample_id, '.tsv'),
            col_names = F
          )
  })

# wo immune celltypes
dir.create('data/Zhang2023/annotations_wo_immune/', recursive = T)
immune_cts <- c('NK', 'B cell/plasma cell', 'Myeloid cell', 'T cell')
annotations %>%
  dplyr::filter(!(celltype %in% immune_cts)) %>%
  {split(., .$sample_id)} %>%
  purrr::walk2(.x = names(.), .y = ., function(sample_id, df) {
    df %>%
      dplyr::select(cell, celltype) %>%
      readr::write_tsv(
        paste0('data/Zhang2023/annotations_wo_immune/', sample_id, '.tsv'),
        col_names = F
      )
  })

# mappings CSV (sample_id,raw_counts_matrix,annotations,id) ----

# all celltypes
mappings <-
  list.files(
    cellranger_dir,
    pattern = '^BRI\\-',
    include.dirs = T) %>%
  tibble::tibble(sample_id = .)  %>%
  # write file paths
  dplyr::transmute(
    sample_id,
    raw_counts_matrix = paste0(processed_dir, 'raw_mRNA_count_matrix.rds'),
    annotations = paste0(wd, '/data/Zhang2023/annotations/', sample_id, '.tsv'),
    id = sample_id)

# wo immune celltypes
mappings_wo_immune <- 
  mappings %>%
  dplyr::mutate(
    annotations = paste0(wd, '/data/Zhang2023/annotations_wo_immune/', sample_id, '.tsv')
  ) 

# check annotations files exist
checks <-
  tibble::tibble(
    file = c(mappings$raw_counts_matrix, 
             mappings$annotations,
             mappings_wo_immune$raw_counts_matrix,
             mappings_wo_immune$annotations)) %>%
  dplyr::mutate(exists = file.exists(file)) %>%
  dplyr::filter(!exists)
if(nrow(checks) > 0) {
  message('Files missing!\n', paste(checks$file, collapse = '\n'))
}

# write files
filter_files <- function(df) {
  df %>% 
    dplyr::filter(file.exists(annotations), 
                  file.exists(raw_counts_matrix))
}
mappings %>%
  filter_files() %>%
  readr::write_csv('data/Zhang2023/mappings.csv')
mappings_wo_immune %>%
  filter_files() %>%
  readr::write_csv('data/Zhang2023/mappings_wo_immune.csv')

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
