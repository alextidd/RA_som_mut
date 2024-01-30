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

# load ref
ensembl <- useEnsembl(biomart = 'ensembl',
                      dataset = 'hsapiens_gene_ensembl')

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

# get genes with strong evidence
drivers <-
  muts %>%
  dplyr::filter(qglobal_cv < 0.2) %>%
  dplyr::rename(gene = gene_name)

# get positions
drivers_pos <-
  getBM(attributes = c('hgnc_symbol', 'chromosome_name', 'start_position', 'end_position'),
        filters = 'hgnc_symbol', values = unique(drivers$gene),
        mart = ensembl)

# manually get aliases for genes that are not matched
unique(drivers$gene[!drivers$gene %in% drivers_pos$hgnc_symbol])
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

# get positions again
drivers_pos <-
  getBM(attributes = c('hgnc_symbol', 'ensembl_transcript_id',
                       'chromosome_name', 'start_position', 'end_position'),
        filters = 'hgnc_symbol', values = unique(drivers$gene_fixed),
        mart = ensembl)

# check if it was fixed
unique(drivers$gene_fixed[!drivers$gene_fixed %in% drivers_pos$hgnc_symbol])

# get exon positions
drivers_exons <-
  getBM(attributes = c('ensembl_transcript_id', 'exon_chrom_start', 'exon_chrom_end'),
        filters = 'hgnc_symbol', values = unique(drivers$gene_fixed),
        mart = ensembl)

# combine spans + exons
drivers_info <-
  dplyr::left_join(
    drivers_pos,
    drivers_exons
  ) %>%
  tibble::as_tibble() %>%
  dplyr::transmute(
    gene = hgnc_symbol,
    ensembl_transcript_id,
    chr = chromosome_name,
    start = start_position,
    end = end_position,
    exon_start = exon_chrom_start,
    exon_end = exon_chrom_end
  ) %>%
  # filter weird chromosomes
  dplyr::filter(chr %in% c(as.character(1:22, 'X', 'Y'))) %>%
  # number exons
  dplyr::group_by(gene, ensembl_transcript_id) %>%
  dplyr::arrange(exon_start, .by_group = T) %>%
  dplyr::mutate(exon = paste0('exon', dplyr::row_number()))

# plot genes with strong evidence
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

# save coords for coverage
drivers_info %>% 
  dplyr::group_by(gene) %>%
  dplyr::summarise(coords = paste0(unique(chr), ':', min(start), '-', max(end))) %>%
  readr::write_tsv('data/driver_genes/driver_gene_coords_for_coverage.tsv')
