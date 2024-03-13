#!/usr/bin/env Rscript

#libraries
library(magrittr)
library(GenomicRanges)

# directory
wd <- 
  ifelse(Sys.info()['nodename'] == "mib119104s", '~/Volumes/', '') |> 
  paste0('/lustre/scratch125/casm/team268im/at31/RA_som_mut/scomatic/')
setwd(wd)

# reference files
fasta <- '/lustre/scratch125/casm/team268im/al28/bed_ref/GRCh38_full_analysis_set_plus_decoy_hla_genome.fa'
hg19_to_hg38 <- '/lustre/scratch125/casm/team268im/at31/RA_som_mut/scomatic/data/liftOver/hg19ToHg38.over.chain'

# recursites from tnanoseq ----

# files
tnanoseq <- '/lustre/scratch125/casm/team268im/lh22/synovium_nanoseq/ex_t_dndsout/dnds_nano_muts.tsv'
tgenes <- 'data/driver_genes/Sanger_TERT-v4_TE-95148282_hg19_highstringencyfilter_buccal_gene_list.tsv'

# read targeted nanoseq 
tnanoseq <- readr::read_tsv(tnanoseq)

# read genes of interest
target_genes <- readr::read_tsv(tgenes)

# run dndscv
dnds <-
  dndscv::dndscv(
    unique(tnanoseq[, 1:5]),
    max_muts_per_gene_per_sample = Inf,
    max_coding_muts_per_sample = Inf,
    gene_list = target_genes$gene,
    outmats = T)

# run sitednds
sitednds <-
  dndscv::sitednds(dnds)

# get the significant genes
signif_genes <-
  dnds$sel_cv %>%
  tibble::as_tibble() %>%
  dplyr::filter(qglobal_cv < 0.1)

# get the sites with positively selected mutations
recursites <-
  sitednds$recursites %>%
  tibble::as_tibble() %>%
  # fix 1d table columns -> numeric columns
  dplyr::mutate(dplyr::across(where(is.table), as.numeric))

# convert to granges (1-based!)
gr <- GRanges(
  seqnames = paste0('chr', recursites$chr),
  ranges = IRanges(start = recursites$pos, end = recursites$pos),
  gene = recursites$gene,
  pos_GRCh37 = recursites$pos)

# get chain file
chain <-
  rtracklayer::import.chain(hg19_to_hg38)

# lift over
lifted <-
  rtracklayer::liftOver(gr, chain) %>%
  unlist() %>%
  tibble::as_tibble() %>%
  dplyr::transmute(chr = gsub('chr', '', seqnames), 
                   pos = start, gene, pos_GRCh37)

# update positions
recursites <-
  recursites %>%
  dplyr::rename(pos_GRCh37 = pos) %>%
  dplyr::left_join(lifted, by = c('chr', 'gene', 'pos_GRCh37')) %>%
  dplyr::select(chr, pos, gene, ref, mut, pos_GRCh37, everything())

# save
recursites %>%
  readr::write_tsv('out/Zhang2023/targeted_mutation_calling/summary/recursites.tsv')

# nonsense sites in TNFRSF14
library(dndscv)
source('src/withingenednds.R')

# read in tnfrsf14
tnfrsf14 <-
  readr::read_tsv('data/driver_genes/TNFRSF14_withingene_sites.tsv') %>%
  # get all possible nonsense mutations
  dplyr::filter(!is.na(impact), impact != 'Synonymous')

# lift over GRCh37 -> GRCh38
# convert to granges (1-based!)
gr <- 
  tnfrsf14 %>%
  dplyr::distinct(chr, pos, gene) %>%
  {GRanges(
    seqnames = paste0('chr', .$chr),
    ranges = IRanges(
      start = .$pos, end = .$pos),
    gene = .$gene,
    pos_GRCh37 = .$pos) 
  }
  
# get chain file
chain <-
  rtracklayer::import.chain(hg19_to_hg38)

# lift over
lifted <-
  rtracklayer::liftOver(gr, chain) %>%
  unlist() %>%
  tibble::as_tibble() %>%
  dplyr::transmute(chr = as.numeric(gsub('chr', '', seqnames)), 
                   pos = start, gene, pos_GRCh37)

# update positions
tnfrsf14_lifted <-
  tnfrsf14 %>%
  dplyr::rename(pos_GRCh37 = pos) %>%
  dplyr::left_join(lifted, by = c('chr', 'gene', 'pos_GRCh37')) %>%
  dplyr::select(chr, pos, gene, ref, mut, pos_GRCh37, everything())

# save
tnfrsf14_lifted %>%
  dplyr::mutate(gene = paste0(gene, '_', pos)) %>%
  readr::write_tsv('data/driver_genes/TNFRSF14_withingene_sites_GRCh38.tsv')


