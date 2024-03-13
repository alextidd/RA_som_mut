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

# clinical metadata
clinical <- 
  readr::read_csv('data/Zhang2023/syn26710600/AMP-RA.SLE_clinical.csv') %>%
  dplyr::rename(id2 = individualID) %>%
  dplyr::mutate(age = as.numeric(ageEnrollment))
readr::read_tsv('data/Zhang2023/syn26710600/CTAP_donor_mapping.tsv') %>%
  dplyr::rename(id = donor, id2 = subject_id) %>%
  dplyr::left_join(clinical) %>%
  readr::write_tsv('data/Zhang2023/clinical_metadata.tsv')

# annotations TSVs (cell celltype) ----
readRDS(paste0(processed_dir, 'all_cells_reference.rds'))$meta_data %>%
  dplyr::select(cell, 
                celltype = cell_type, 
                cluster = cluster_name, 
                id = sample) %>%
  readr::write_tsv('data/Zhang2023/annotations.tsv')

# mappings CSV (id,raw_counts_matrix) ----
files <-
  list.files(
    cellranger_dir,
    pattern = '^BRI\\-',
    include.dirs = T) %>%
  tibble::tibble(id = .)  %>%
  # write file paths
  dplyr::transmute(
    id,
    raw_counts_matrix = paste0(processed_dir, 'raw_mRNA_count_matrix.rds')) 

# check for missing files
missing <- 
  files %>%
  dplyr::filter(!file.exists(raw_counts_matrix))
if (nrow(missing) > 0) {
  message('Files missing!\n', paste(missing$raw_counts_matrix, collapse = '\n'))
}

# write mappings
files %>%
  # matrix exists
  dplyr::filter(file.exists(raw_counts_matrix)) %>%
  # id in annotations file
  dplyr::filter(id %in% annotations$id) %>%
  readr::write_csv(
    paste0('data/Zhang2023/mappings.csv')
  )

######################################
# # Zhang 2019 ---
# # mp34 previously ran inferCNV on RA samples from the Zhang 2019 paper. 
# # Here we will run this smaller dataset again as a test.
# dir.create('data/Zhang2019/annotations/', recursive = T)
# 
# dat_dir_2019 <- '/lustre/scratch126/casm/team268im/mp34/analysis/synovium_scRNA/Celseq_wo_normal/'
# annots_2019 <-
#     file.path(dat_dir_2019, 'cellAnnotations_Celseq_wo_normal.txt')
# matrix_2019 <-
#     file.path(dat_dir_2019, 'sc.10x.counts_Celseq_wo_normal.RData')
# load(matrix_2019)
# 
# # mappings CSV ----
# mtx_file <- paste0(wd, 'data/Zhang2019/raw_counts_matrix.tsv')
# saveRDS(mat.matrix.sparse, mtx_file)
# tibble::tibble(sample_id = 'all',
#                raw_counts_matrix = mtx_file,
#                annotations = annots_2019,
#                id = 'all') %>%
#     readr::write_csv('data/Zhang2019/mappings.csv')
