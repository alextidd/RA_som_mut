library(magrittr)
library(ggplot2)

obj <- readRDS('/lustre/scratch126/casm/team268im/mp34/scRNAseq_data/RA_ZhangEtal2023/processed_output/all_cells_reference.rds')

cnvs <- 
  list.files(
    'out/Zhang2023/by_celltype/',
    pattern = 'map_metadata_from_infercnv.txt',
    recursive = T
  ) %>%
  purrr::set_names(., dirname(.)) %>%
  purrr::map(function(file) {
    paste0('out/Zhang2023/by_celltype/', file) %>% 
      read.table(sep = '\t') %>% 
      tibble::as_tibble(rownames = 'cell')
    }) %>%
  dplyr::bind_rows(.id = 'id')

metadata <-
  obj$meta_data %>%
  tibble::as_tibble() %>%
  dplyr::mutate(id = sample)

umap <-
  obj$umap$embedding %>%
  tibble::as_tibble() %>%
  dplyr::bind_cols(metadata) %>%
  dplyr::left_join(cnvs, by = c('id', 'cell'))

long_umap <-
  umap %>%
  dplyr::filter(!is.na(top_loss_1)) %>%
  tidyr::pivot_longer(
    cols = tidyselect::starts_with(c('has_', 'proportion_', 'top_'))
  ) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(
    chr = name %>%
      gsub('.*_', '', .) %>%
      {forcats::fct_reorder(., gsub('chr', '', .) %>% as.numeric())},
    event = sub(paste0(chr, '_'), '', name),
    measure = sub(paste0(chr, '_', event), '', name))
  
# plot gains and losses
plot_gains_and_losses <- 
  function(long_umap, annot) {
    long_umap %>%
      dplyr::rename(annotation = tidyselect::all_of(annot)) %>%
      dplyr::filter(measure == 'has', event != 'cnv') %>%
      dplyr::group_by(id, chr, event, annotation) %>%
      dplyr::summarise(n_cells = dplyr::n(),
                       n_events = sum(value),
                       prop_of_cells = n_events / n_cells) %>%
      ggplot(aes(x = id, y = chr, alpha = prop_of_cells, fill = event)) +
      geom_tile() +
      theme_classic() +
      theme(panel.spacing = unit(0, "lines")) +
      scale_fill_manual(values = c('loss' = scales::muted('blue'),
                                   'dupli' = scales::muted('red'))) +
      scale_alpha(range = c(0, 1)) +
      facet_grid(event ~ annotation)
  }

pdf('test.pdf', width = 30)
umap %>%
  dplyr::filter(!is.na(has_dupli_chr7)) %>%
  ggplot(aes(x = UMAP1, y = UMAP2, colour = cell_type, alpha = has_dupli_chr7)) +
  geom_point(size = 0.1) +
  theme_classic() +
  coord_fixed()
plot_gains_and_losses(long_umap, 'cell_type')
plot_gains_and_losses(long_umap, 'cluster_name')

long_umap %>%
  dplyr::filter(measure == 'proportion', event != 'cnv') %>%
  dplyr::group_by(id, chr, event, cell_type) %>%
  dplyr::summarise(avg_prop_of_cells = mean(value)) %>%
  ggplot(aes(x = id, y = chr, alpha = avg_prop_of_cells, fill = event)) +
  geom_tile() +
  theme_classic() +
  theme(panel.spacing = unit(0, "lines")) +
  scale_fill_manual(values = c('loss' = scales::muted('blue'),
                               'dupli' = scales::muted('red'))) +
  scale_alpha(range = c(0, 1)) +
  facet_grid(event ~ cell_type)
dev.off()

# plot proportions


