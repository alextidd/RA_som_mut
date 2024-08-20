# Load the libraries
library(httr)
library(jsonlite)
library(ggplot2)
library(cowplot)
library(magrittr)

# Function to fetch exon and intron data from Ensembl (GRCh37) for all transcripts of a gene
get_all_transcripts_data_grch37 <- function(gene_name) {
  ensembl_url <-
    paste0("https://grch37.rest.ensembl.org/lookup/symbol/homo_sapiens/",
           gene_name, "?expand=1")
  response <- GET(ensembl_url, content_type("application/json"))

  if (status_code(response) != 200) {
    stop("Failed to retrieve data from Ensembl API")
  }

  gene_data <- fromJSON(content(response, as = "text", encoding = "UTF-8"))

  all_transcripts <- list()

  names(gene_data$Transcript$Exon) <- gene_data$Transcript$id
  for (transcript in names(gene_data$Transcript$Exon)) {
    exons <- gene_data$Transcript$Exon[[transcript]]

    if (nrow(exons) > 1) {
      introns <- data.frame(
        start = exons[-1, "end"] + 1,
        end = exons[-length(exons$end), "start"] - 1,
        Type = "Intron"
      )
    }

    exons$Type <- "Exon"
    exons <- exons[, c("start", "end", "Type")]
    introns <- introns[, c("start", "end", "Type")]

    transcript_data <- rbind(exons, introns)
    transcript_data$transcript_id <- transcript

    all_transcripts[[transcript]] <- transcript_data
  }

  out <- do.call(rbind, all_transcripts)

  # Differentiate canonical transcript
  canonical_transcript_id <-
    gsub("\\..*", "", gene_data$canonical_transcript)
  out$canonical <- 
    out$transcript_id == canonical_transcript_id
  canonical <- out %>%
    dplyr::filter(canonical == TRUE) %>%
    dplyr::mutate(Type = "Canonical")
  out <- dplyr::bind_rows(out, canonical) %>%
    dplyr::mutate(
      Type = factor(Type, levels = c("Canonical", "Intron", "Exon")),
      size = dplyr::case_match(Type, "Canonical" ~ 1, "Exon" ~ 0.9,
                               "Intron" ~ 0.5))

  # Add biotype information
  out <-
    gene_data$Transcript %>%
    dplyr::select(transcript_id = id, biotype) %>%
    {dplyr::left_join(out, .)} %>%
    dplyr::mutate(biotype = dplyr::case_when(Type == "Canonical" ~ "Canonical",
                                             TRUE ~ biotype))

  return(out)
}

# Function to Fetch Protein Domain Data from InterPro API
get_protein_domains <- function(gene_name) {
  # Convert the gene name to a UniProt ID using the Ensembl API
  query <- list("gene_exact" = gene_name, "reviewed" = "true",
                "organism_id" = "9606")
  uniprot_id <- queryup::query_uniprot(query)$Entry

  # Fetch protein domains from InterPro using the UniProt ID
  interpro_url <- paste0("https://www.ebi.ac.uk/proteins/api/proteins/",
                         uniprot_id)
  response <- GET(interpro_url, content_type("application/json"))
  protein_data <- fromJSON(content(response, as = "text", encoding = "UTF-8"))

  domains <-
    protein_data$features %>%
    dplyr::filter(type == "DOMAIN") %>%
    dplyr::transmute(start = begin, end, Type = description,
                     length = protein_data$sequence$length) %>%
    readr::type_convert()

  return(domains)
}

# Function to plot the exon/intron structure for all transcripts and protein domains
plot_gene_tracks <- function(gene_name, annotmuts) {
  # Get exon/intron data (GRCh37) for all transcripts
  all_transcripts_data <- get_all_transcripts_data_grch37(gene_name)

  # Get protein domain data
  protein_domains <- get_protein_domains(gene_name)

  # Get canonical exon boundaries to align with CDS
  exons <-
    all_transcripts_data %>%
    tibble::as_tibble() %>%
    dplyr::filter(Type == "Exon", canonical == TRUE) %>%
    dplyr::arrange(start) %>%
    dplyr::mutate(interval = end - start,
                  cumulative = cumsum(interval),
                  cds_start = cumulative - interval + 1,
                  cds_end = cumulative,
                  midpoint = (cds_start + cds_end) / 2,
                  exon_number = dplyr::row_number(),
                  alternating = as.character(exon_number %% 2))

  # TODO: Plot mutations from dNdScv output
  # try using ggrepel to place aachange labels
  mut_plot <-
    annotmuts %>%
    dplyr::filter(gene == gene_name) %>%
    dplyr::distinct() %>%
    dplyr::group_by(chr, pos, ref, mut) %>%
    dplyr::mutate(n_donors = dplyr::n_distinct(donor)) %>%
    ggplot(aes(x = pos, xend = pos, y = n_donors, yend = 0,
               label = aachange, colour = impact)) +
    geom_segment() +
    ggrepel::geom_label_repel()

  # Plot for Exon/Intron Structure of All Transcripts
  exon_intron_plot <- 
    ggplot(dplyr::arrange(all_transcripts_data, Type),
           aes(x = start, xend = end, y = transcript_id, yend = transcript_id,
           size = Type, colour = biotype)) +
    geom_segment() +
    labs(title = paste("Exon/Intron Structure of All Transcripts for",
                       gene_name),
         x = "Genomic Position", y = "Transcript ID") +
    scale_size_manual(values = c("Exon" = 3, "Intron" = 1, "Canonical" = 3)) +
    theme_minimal() +
    theme(legend.position = "bottom") +
    guides(size = FALSE)

  # Plot for Exon Boundaries in the Canonical CDS
  exon_boundaries_plot <-
    ggplot(exons, aes(xmin = cds_start, xmax = cds_end, ymin = -1, ymax = 1,
           fill = alternating)) +
    geom_rect() +
    scale_fill_manual(values = c("1" = "lightblue", "0" = "darkblue")) +
    theme_minimal() +
    labs(title = paste("Exon Boundaries in the Canonical", gene_name, "CDS"),
         x = "CDS Position", y = "") +
    theme(axis.text.x = element_text(), legend.position = "none")

  # Plot for Protein Domains
  protein_domain_plot <-
  ggplot(protein_domains, aes(x = start, xend = end, y = Type, yend = Type)) +
    geom_segment(size = 5, color = "blue") +
    labs(title = paste("Protein Domains for", gene_name),
         x = "Amino Acid Position", y = "") +
    theme_minimal() +
    theme(axis.ticks.y = element_blank(),
          axis.title.y = element_blank()) +
    xlim(0, unique(protein_domains$length))

  # Combine the two plots using cowplot
  combined_plot <- plot_grid(exon_intron_plot, exon_boundaries_plot,
                             protein_domain_plot, ncol = 1)

  print(combined_plot)
}

# load the annotated mutations from nanoseq output (`nano_muts`)
lh22_dir <- "/lustre/scratch125/casm/team268im/lh22/Rheumatoid_Arthritis_Paper/"
load(file.path(lh22_dir, "Figures/Figure_2/LCM_NanoSeq_Overlap/GOI",
               "pileup_muts.rda"))
nano_muts <- nano_muts %>% tibble::as_tibble()

# run dNdScv
dndscv_in <- nano_muts %>% dplyr::select(sampleID, chr, pos, ref, mut)
dndscv_out <- dndscv::dndscv(dndscv_in)
annotmuts <-
  nano_muts %>%
  dplyr::select(sampleID, donor) %>%
  {dplyr::left_join(dndscv_out$annotmuts, .)}

# Example usage for any gene
plot_gene_tracks("NF1", dndscv_out_w_donor)
