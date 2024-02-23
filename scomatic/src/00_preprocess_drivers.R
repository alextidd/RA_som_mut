# libraries
library(magrittr)
library(ggplot2)
library(biomaRt)
library(org.Hs.eg.db)

# dirs
wd <- 
  ifelse(Sys.info()['nodename'] == "mib119104s", '~/Volumes/', '') %>% 
  paste0('/lustre/scratch125/casm/team268im/at31/RA_som_mut/scomatic/') 
setwd(wd)

# Genes of interest in rheumatoid arthritis (hereafter referred to as 'drivers') 
# are defined as genes that are found to undergo significant selection (q-value 
# < 0.2) according to the whole exome and targeted NanoSeq data. The coverage of 
# these sites in the dataset, across samples and celltypes, will be checked.
# First, we read in the dNdS results and filter by q-value.

# read in dnds results
muts <-
  'data/dnds_nanoseq/' %>%
  list.files(pattern = 'nanoseq') %>%
  purrr::map(function(assay) {
    list.files(paste0('data/dnds_nanoseq/', assay), 
               pattern = '_dnds_sel_cv.tsv', 
               full.names = T) %>%
      purrr::map(function(file) {
        file %>%
          readr::read_tsv(show_col_types = F) %>%
          dplyr::mutate(
            assay = assay,
            celltype = file %>% 
              basename() %>% 
              gsub('_dnds_sel_cv.tsv', '', .))
      }) %>%
      dplyr::bind_rows()
  }) %>%
  dplyr::bind_rows() %>%
  dplyr::select(assay, celltype, everything())

# get genes with strong evidence (q-value)
drivers <-
  muts %>%
  dplyr::filter(qglobal_cv < 0.2) %>%
  dplyr::rename(gene = gene_name)

# manually write aliases for genes that are not matched
# HIST1H3B -> H3C2
# HIST1H1E -> H1-4
# CRIPAK   -> This record has been withdrawn by NCBI, after discussions with 
#             CCDS collaborators. It was decided that this locus is not an 
#             independent gene.
drivers <-
  drivers %>%
  dplyr::mutate(
    gene_fixed = dplyr::case_when(
      gene == 'HIST1H3B' ~ 'H3C2',
      gene == 'HIST1H1E' ~ 'H1-4',
      TRUE ~ gene
    )
  )

# plot drivers
drivers %>%
  dplyr::distinct(celltype, gene, assay, qglobal_cv) %>%
  dplyr::mutate(neglog10q = -log10(qglobal_cv)) %>%
  ggplot(aes(x = tidytext::reorder_within(gene, -neglog10q, celltype), 
             y = neglog10q, fill = celltype)) +
  geom_col(position = 'dodge') +
  theme_classic() +
  theme(axis.text.x = element_text(angle = -90, hjust = 0)) +
  facet_grid(assay ~ celltype, scales = 'free_x', space = 'free_x') +
  tidytext::scale_x_reordered() +
  labs(x = 'gene', y = '-log10(p)')

# In order to calculate their coverage in the BAMs, we must get their positions.

# load ref
library(biomaRt)
ensembl <- useEnsembl(biomart = 'ensembl', dataset = 'hsapiens_gene_ensembl')

# get positions
drivers_pos <-
  getBM(attributes = c('hgnc_symbol', 'chromosome_name', 'start_position', 'end_position'),
        filters = 'hgnc_symbol', values = unique(drivers$gene),
        mart = ensembl) %>%
  tibble::as_tibble() %>% 
  # fix colnames
  dplyr::transmute(
    chr = chromosome_name, 
    start = start_position,
    end = end_position,
    gene = hgnc_symbol) %>%
  # remove weird chromosomes
  dplyr::filter(chr %in% c(as.character(1:22, 'X', 'Y')))

# save bed
drivers_pos %>%
  dplyr::transmute(
    paste0('chr', chr), start, end, gene
  ) %>%
  readr::write_tsv('data/driver_genes/driver_gene_coords_for_coverage.bed',
                   col_names = F)

# save positions
drivers_pos %>% 
  dplyr::group_by(gene) %>%
  dplyr::summarise(coords = 
                     paste0(unique(chr), ':', 
                            min(start), '-', 
                            max(end))) %>%
  readr::write_tsv('data/driver_genes/driver_gene_coords_for_coverage.tsv')

# get reference object ----

# load rda 
refcds_38 <- '/lustre/scratch125/casm/team268im/fa8/119/DNDSCV_COVARIATES/refcds_GRCh38-GencodeV18+Appris.rda'
load(refcds_38)

# set gene names
ref_cds <-
  RefCDS %>%
  setNames(RefCDS %>% purrr::map(~ .x$gene_name) %>% unlist())

# save 
ref_cds %>%
  saveRDS('data/dndscv/ref_cds.rds')
