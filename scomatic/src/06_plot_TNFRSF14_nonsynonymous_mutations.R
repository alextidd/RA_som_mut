#!/usr/bin/env Rscript

# wd
wd <- 
  ifelse(Sys.info()['nodename'] == "mib119104s", '~/Volumes/', '') |> 
  paste0('/lustre/scratch125/casm/team268im/at31/RA_som_mut/scomatic/')
setwd(wd)

#libraries
library(magrittr)
library(GenomicFeatures)

# get recursites in TNFRSF14
recursites <-
  readr::read_tsv('data/driver_genes/TNFRSF14_withingene_sites_GRCh38.tsv')

# get all reads
reads <-
  readRDS('out/Zhang2023/targeted_mutation_calling/TNFRSF14/summary/reads.rds')

# get exon coords
exons <- 
  readRDS('data/dndscv/ref_cds.rds')$TNFRSF14 %>%
  {tibble::tibble(
    chr = as.numeric(.$chr),
    start = .$intervals_cds[, 1],
    end = .$intervals_cds[, 2],
    gene = 'TNFRSF14') %>%
  dplyr::mutate(feature = paste0('exon', dplyr::row_number()))}
introns <-
  gene_spans %>%
  dplyr::group_by(chr, gene) %>%
  dplyr::transmute(
    new_end = start - 1, start = dplyr::lag(end) + 1) %>%
  dplyr::filter(!is.na(start)) %>%
  dplyr::transmute(
    start, end = new_end,
    feature = paste0('intron', dplyr::row_number()))

# get all features
features <-
  dplyr::bind_rows(exons, introns) %>%
  dplyr::arrange(start) %>%
  dplyr::mutate(feature = forcats::fct_reorder(feature, dplyr::row_number())) 
features_all_positions <-
  features %>%
  dplyr::mutate(pos = purrr::map2(start, end, ~ .x:.y)) %>%
  tidyr::unnest(cols = 'pos') %>%
  dplyr::select(-start, -end)

# alts - colours and order
nt_cols <-
  c('INS' = '#e6ab02',
    'DEL' = '#1c2d78',
    'A' = '#1b9e77', 
    'C' = '#d95f02', 
    'G' = '#7570b3', 
    'T' = '#e7298a',
    'ref' = 'grey')

# impacts - colours and order
imp_cols <-
  c('Start_loss' = '#e6ab02',
    'Missense' = '#e7298a',
    'Nonsense' = '#1b9e77', 
    'Essential_Splice' = '#d95f02', 
    'Stop_loss' = '#7570b3')

# plot mutations
p_dat <-
  reads %>%
  dplyr::bind_rows() %>%
  dplyr::filter(pos == mut_pos) %>%
  dplyr::left_join(recursites, relationship = 'many-to-many') %>%
  dplyr::filter(name == mut) %>%
  dplyr::group_by(chr, pos, ref, status, aachange, celltype, gene, impact) %>%
  dplyr::summarise(value = as.integer(sum(value))) %>%
  dplyr::mutate(
    aa_change_label = dplyr::case_when(value > 0 ~ aachange,
                                 TRUE ~ NA)) %>%
  dplyr::group_by(celltype) %>%
  dplyr::filter(sum(value) > 0) %>%
  dplyr::left_join(features_all_positions %>% dplyr::select(-gene))

pdf('test.pdf', height = 15, width = 20)

# exons and introns
p_dat %>%
  ggplot(aes(x = pos, y = value, colour = impact)) + 
  # plot line at y = 0
  geom_segment(aes(x = -Inf, xend = Inf, y = 0, yend = 0),
               colour = 'black') +
  # plot all detected nonsense mutations
  geom_segment(aes(x = pos, xend = pos, y = 0, yend = value)) +
  geom_point() +
  theme_classic() +
  theme(panel.grid.major.y = element_line(),
        panel.border = element_blank(),
        axis.line.x = element_blank()) +
  scale_colour_manual(values = imp_cols) +
  scale_x_continuous(guide = guide_axis(angle = -45), n.breaks = 2) +
  ggh4x::facet_grid2(
    celltype ~ feature, scales = 'free', space = 'free') +
  labs(x = 'TNFRSF14', 
       title = 'All nonsynonymous mutations detected in TNFRSF14') 

# exons only  
p_dat %>%
  dplyr::filter(grepl('exon', feature), value > 0) %>%
  ggplot(aes(x = pos, y = value, colour = impact)) + 
  # plot line at y = 0
  geom_segment(aes(x = -Inf, xend = Inf, y = 0, yend = 0),
               colour = 'black') +
  # plot all detected nonsense mutations
  geom_segment(aes(x = pos, xend = pos, y = 0, yend = value)) +
  geom_point() +
  ggrepel::geom_label_repel(
    data = p_dat %>% dplyr::group_by(celltype) %>% dplyr::slice_max(value, n = 3),
    aes(label = aa_change_label),
    min.segment.length = 0, nudge_y = 0.2,
    show.legend = F) +
  theme_classic() +
  theme(panel.grid.major.y = element_line(),
        panel.border = element_blank(),
        axis.line.x = element_blank()) +
  scale_colour_manual(values = imp_cols) +
  scale_x_continuous(guide = guide_axis(angle = -45), n.breaks = 2) +
  ggh4x::facet_grid2(
    celltype ~ feature, scales = 'free', space = 'free') +
  labs(x = 'TNFRSF14', 
       title = 'All nonsynonymous mutations detected in exons of TNFRSF14') 

dev.off()

# file for cbioportal mutation mapper
reads %>%
  dplyr::bind_rows() %>%
  dplyr::filter(pos == mut_pos) %>%
  dplyr::left_join(recursites, relationship = 'many-to-many') %>%
  dplyr::filter(name == mut, value > 0) %>%
  dplyr::transmute(
    Sample_ID = paste0(id, '_', celltype),
    Chromosome = chr,
    Start_Position = pos,
    End_Position = pos,
    Reference_Allele = ref,
    Variant_Allele = status,
    n = value) %>% 
  tidyr::uncount(n) %>%
  readr::write_tsv('out/Zhang2023/targeted_mutation_calling/TNFRSF14/summary/muts_for_cbioportal.tsv')