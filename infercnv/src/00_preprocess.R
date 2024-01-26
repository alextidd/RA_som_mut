# libraries
library(magrittr)

# dirs
wd <- 
  ifelse(Sys.info()['nodename'] == "mib119104s", '~/Volumes/', '') %>%
  paste0('/lustre/scratch125/casm/team268im/at31/RA_som_mut/infercnv/') 
setwd(wd)

# Zhang 2023 ---
dat_dir <- '/lustre/scratch126/casm/team268im/mp34/scRNAseq_data/RA_ZhangEtal2023/'
cellranger_dir <- paste0(dat_dir, 'cellranger_output/')
processed_dir <- paste0(dat_dir, 'processed_output/')
dir.create('data/Zhang2023/annotations/', recursive = T)

# annotations TSVs (cell celltype) ----
annotations <-
  readRDS(paste0(processed_dir, 'all_cells_reference.rds'))$meta_data %>%
  dplyr::select(cell, cell_type, sample_id = sample) 
annotations %>%
  dplyr::group_by(sample_id) %>%
  dplyr::group_split() %>%
  purrr::walk(function(df) {
      sample_id <- unique(df$sample_id)
      df %>%
          dplyr::select(cell, cell_type) %>%
          readr::write_tsv(
            paste0('data/Zhang2023/annotations/', sample_id, '.tsv'),
            col_names = F
          )
  })

# mappings CSV (sample_id,raw_counts_matrix,annotations,id) ----
mappings <-
  list.files(
    cellranger_dir,
    pattern = '^BRI\\-',
    include.dirs = T) %>%
  tibble::enframe(value = 'sample_id') %>%
  # get only sample_ids present in the annotations
  dplyr::filter(sample_id %in% annotations$sample_id) %>%
  # write file paths
  dplyr::transmute(
    sample_id,
    raw_counts_matrix = paste0(processed_dir, 'raw_mRNA_count_matrix.rds'),
    annotations = paste0(wd, '/data/Zhang2023/annotations/', sample_id, '.tsv'),
    id = sample_id)
mappings %>%
  readr::write_csv('data/Zhang2023/mappings.csv')

# check all files exist
checks <-
  tibble::tibble(
    file = c(mappings$raw_counts_matrix, mappings$annotations)) %>%
  dplyr::mutate(exists = file.exists(file)) %>%
  dplyr::filter(!exists)
if(nrow(checks) > 0) {
  message('Files missing!\n', paste(checks$file, collapse = '\n'))
  }

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
