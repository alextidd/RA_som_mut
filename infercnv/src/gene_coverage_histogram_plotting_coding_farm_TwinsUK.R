library("tidyverse")
library("dndscv")
library("scales")
library("patchwork")

args <- commandArgs(TRUE)

gene_num <- as.numeric(args[1])
datasets <- c("plate_014_SQPP-1106-C_targeted_buccal_v3","plate_017_SQPP-1640-T_targeted_buccal_v2","plate_018_SQPP-1728-D_targeted_buccal_v2","plate_019_SQPP-1734-B_targeted_buccal_v2",
              "plate_020_SQPP-1737-E_targeted_buccal_v2","plate_021_SQPP-2297-F_targeted_buccal","plate_022_SQPP-2906-H_targeted_buccal","plate_024_SQPP_4022-L_targeted_buccal",
              "plate_029_SQPP-13328-P_targeted_buccal","plate_035_SQPP-18868-M_targeted_buccal","plate_036_SQPP-18870-G_targeted_buccal","plate_038_SQPP-19909-I_targeted_buccal",
              "plate_040_SQPP-20086-J_targeted_buccal","plate_042_SQPP-22755-I_targeted_buccal")

twins_uk_sample_list <- read.table("/lustre/scratch126/casm/team268im/al28/TwinsUK_Manuscript/TwinsUKSampleList.csv", sep = ",", stringsAsFactors = F, header = F)
twins_uk_sample_list <- twins_uk_sample_list$V1

exclude_sample_list <- read.table("/lustre/scratch126/casm/team268im/al28/TwinsUK_Manuscript/SAMPLES2EXCLUDE.2023-11-15.tsv", sep = "\t", stringsAsFactors=F, header = T)
twins_uk_sample_list <- twins_uk_sample_list[which(!(twins_uk_sample_list) %in% rownames(exclude_sample_list))]
twins_uk_sample_list <- twins_uk_sample_list[which(!(grepl("Unknown",twins_uk_sample_list)))]

buccal_genes <- read.table("/lustre/scratch126/casm/team268im/al28/TwinsUK_Manuscript/Sanger_TERT-v4_TE-95148282_hg19_highstringencyfilter_buccal_gene_list.tsv", sep = "\t", stringsAsFactors = F, header = T)
gene <- buccal_genes$gene[gene_num]

muts <- read.table("/lustre/scratch126/casm/team268im/al28/TwinsUK_Manuscript/2023-12-21_TwinsUK_CompleteTNanoSeqCalls_ExcludeSamplesQposArtefact_IndelsFiltered.v3.tsv", sep = "\t", stringsAsFactors=F, header = T)
muts <- muts[which(muts$sampleID %in% twins_uk_sample_list),]

dndsout <- dndscv(muts[,1:5], gene_list = buccal_genes$gene, max_muts_per_gene_per_sample = Inf, max_coding_muts_per_sample = Inf, outmats = T)

data(list = "refcds_hg19", package = "dndscv")

ref_cds_gene_list <- NULL

for(i in 1:length(RefCDS)){
  ref_cds_gene_list <- c(ref_cds_gene_list,RefCDS[[i]]$gene_name)
}

for(i in 1:length(datasets)){
  if(i == 1){
    gene_coverage <- read.table(paste0("/lustre/scratch126/casm/team268im/al28/targeted_nanoseq/",datasets[i],"/Analysis/per_gene_cov/",gene,".tsv"),sep = "\t", stringsAsFactors=F, header=T)
    gene_coverage <- gene_coverage[,c(1:(ncol(gene_coverage) - 1))]
  }else{
    temp_coverage <- read.table(paste0("/lustre/scratch126/casm/team268im/al28/targeted_nanoseq/",datasets[i],"/Analysis/per_gene_cov/",gene,".tsv"),sep = "\t", stringsAsFactors=F, header=T)
    temp_coverage <- temp_coverage[,c(5:(ncol(temp_coverage) - 1))]
    gene_coverage <- cbind(gene_coverage,temp_coverage)
  }
}

gene_coverage <- gene_coverage[,c(1:4,which(colnames(gene_coverage) %in% twins_uk_sample_list)),]

gene_coverage <- gene_coverage[which(gene_coverage$coding == 1),]

if(buccal_genes$strand[which(buccal_genes$gene == gene)] == 1){
  gene_coverage$idx <- 1:nrow(gene_coverage)
}else if(buccal_genes$strand[which(buccal_genes$gene == gene)] == -1){
  gene_coverage$idx <- nrow(gene_coverage):1
}

gene_coverage <- gene_coverage[order(gene_coverage$idx, decreasing = F),]

gene_coverage$total <- NA
gene_coverage$total <- rowSums(gene_coverage[,which(substr(colnames(gene_coverage),1,2) == "PD")])

gene_coverage$pos_transitions <- NA
for(i in 1:nrow(gene_coverage)){
  if(i == 1){
    gene_coverage$pos_transitions[i] <- 1
  }else if(abs(gene_coverage$pos[i] - gene_coverage$pos[i - 1]) == 1){
    gene_coverage$pos_transitions[i] <- gene_coverage$pos_transitions[i - 1]
  }else{
    gene_coverage$pos_transitions[i] <- gene_coverage$pos_transitions[i - 1] + 1
  }
}

gene_pos_transition_table <- as.data.frame(array(data = NA, dim = c(max(gene_coverage$pos_transitions),4)))
colnames(gene_pos_transition_table) <- c("start_pos","end_pos","start_idx","end_idx")

for(i in 1:nrow(gene_pos_transition_table)){
  gene_pos_transition_table$start_pos[i] <- gene_coverage$pos[min(which(gene_coverage$pos_transitions == i))]
  gene_pos_transition_table$end_pos[i] <- gene_coverage$pos[max(which(gene_coverage$pos_transitions == i))]
  gene_pos_transition_table$start_idx[i] <- gene_coverage$idx[min(which(gene_coverage$pos_transitions == i))]
  gene_pos_transition_table$end_idx[i] <- gene_coverage$idx[max(which(gene_coverage$pos_transitions == i))]
}

gene_coverage_smooth <- as.data.frame(array(data = NA, dim = c(nrow(gene_coverage),9)))
colnames(gene_coverage_smooth) <- c("chr","pos","idx","total","pos_transitions","idx_first","idx_last","bin_size","bin_average")
gene_coverage_smooth$chr <- gene_coverage$chr
gene_coverage_smooth$pos <- gene_coverage$pos
gene_coverage_smooth$idx <- gene_coverage$idx
gene_coverage_smooth$total <- gene_coverage$total
gene_coverage_smooth$pos_transitions <- gene_coverage$pos_transitions

gene_cov_bin_size <- 50

for(i in 1:nrow(gene_coverage_smooth)){
  gene_coverage_smooth$idx_first[i] <- max(min(gene_coverage_smooth$idx[which(gene_coverage_smooth$pos_transitions == gene_coverage_smooth$pos_transitions[i])]),gene_coverage_smooth$idx[i] - (gene_cov_bin_size/2))
  gene_coverage_smooth$idx_last[i] <- min(max(gene_coverage_smooth$idx[which(gene_coverage_smooth$pos_transitions == gene_coverage_smooth$pos_transitions[i])]),gene_coverage_smooth$idx[i] + (gene_cov_bin_size/2))
  gene_coverage_smooth$gene_cov_bin_size[i] <- (gene_coverage_smooth$idx_last[i] - gene_coverage_smooth$idx_first[i]) + 1
  gene_coverage_smooth$bin_average[i] <- sum(gene_coverage_smooth$total[gene_coverage_smooth$idx_first[i]:gene_coverage_smooth$idx_last[i]]) / gene_coverage_smooth$gene_cov_bin_size[i]
}

gene_muts <- dndsout$annotmuts[which(dndsout$annotmuts$gene == gene),]

if(nrow(gene_muts) > 0){
  
  for(i in 1:nrow(gene_muts)){
    if(!(gene_muts$pos[i] %in% gene_coverage$pos)){
      diff_in_pos <- abs(gene_coverage$pos - gene_muts$pos[i])
      gene_muts$pos[i] <- gene_coverage$pos[which(diff_in_pos == min(diff_in_pos))][1]
    }
  }
  
  gene_start_loss <- which(gene_muts$impact != "Synonymous" & substr(gene_muts$aachange,2,nchar(gene_muts$aachange)-1) == "1")
  gene_muts$impact[gene_start_loss] <- "Start_loss"
  
  mut_types <- names(table(gene_muts$impact))
  
  fixed_bin_size <- 3
  variable_bin_size <- max(round(buccal_genes$CDS_length[gene_num] / 1500),1) * 3
  
  cosmic_muts <- read.table(paste0("/lustre/scratch126/casm/team268im/al28/TwinsUK_Manuscript/COSMIC_v99_WholeGenomeExome_FullAnnotated_chr",unique(gene_muts$chr),".tsv"), sep = "\t", stringsAsFactors=F, header = T)
  cosmic_gene_muts <- cosmic_muts[which(cosmic_muts$gene == gene),]
  
  if(nrow(cosmic_gene_muts) > 0){
    cosmic_start_loss <- which(cosmic_gene_muts$impact != "Synonymous" & substr(cosmic_gene_muts$aachange,2,nchar(cosmic_gene_muts$aachange)-1) == "1")
    cosmic_gene_muts$impact[cosmic_start_loss] <- "Start_loss"
    
    mut_types <- unique(c(mut_types,names(table(cosmic_gene_muts$impact))))
    
    for(i in 1:nrow(cosmic_gene_muts)){
      if(!(cosmic_gene_muts$pos[i] %in% gene_coverage$pos)){
        diff_in_pos <- abs(gene_coverage$pos - cosmic_gene_muts$pos[i])
        cosmic_gene_muts$pos[i] <- gene_coverage$pos[which(diff_in_pos == min(diff_in_pos))][1]
      }
    }
    
    cosmic_sample_info <- read.table("/lustre/scratch126/casm/team268im/al28/TwinsUK_Manuscript/COSMIC_v99_WholeGenomeExome_SamplePhenotypeInfo.tsv", sep = "\t", stringsAsFactors = F, header = T, comment = "", quote = "")
    cosmic_squamous_samples <- cosmic_sample_info$COSMIC_SAMPLE_ID[which(cosmic_sample_info$HISTOLOGY_SUBTYPE_1 == "squamous_cell_carcinoma")]
    
    cosmic_squamous_gene_muts <- cosmic_gene_muts[which(cosmic_gene_muts$sampleID %in% cosmic_squamous_samples),]
  }
  
  fixed_binned_muts <- as.data.frame(array(data = NA, dim = c(ceiling(nrow(gene_coverage) / fixed_bin_size) * length(mut_types),3)))
  colnames(fixed_binned_muts) <- c("idx","mut_type","count")
  
  fixed_binned_muts$idx <- rep(seq(from = (fixed_bin_size)/2, to = (ceiling(nrow(gene_coverage) / fixed_bin_size) * fixed_bin_size) - (fixed_bin_size)/2, by = fixed_bin_size), each = length(mut_types))
  fixed_binned_muts$mut_type <- rep(mut_types, times = ceiling(nrow(gene_coverage) / fixed_bin_size))
  
  if(nrow(cosmic_gene_muts) > 0){
    cosmic_fixed_binned_muts <- fixed_binned_muts
    
    for(i in 1:nrow(cosmic_fixed_binned_muts)){
      if(buccal_genes$strand[which(buccal_genes$gene == gene)] == 1){
        cosmic_fixed_binned_muts$count[i] <- length(which(cosmic_gene_muts$impact == cosmic_fixed_binned_muts$mut_type[i] & cosmic_gene_muts$pos >= gene_coverage$pos[cosmic_fixed_binned_muts$idx[i] - ((fixed_bin_size)/2 - 1)] & cosmic_gene_muts$pos <= gene_coverage$pos[min(cosmic_fixed_binned_muts$idx[i] + (fixed_bin_size)/2,nrow(gene_coverage))]))
      }else if(buccal_genes$strand[which(buccal_genes$gene == gene)] == -1){
        cosmic_fixed_binned_muts$count[i] <- length(which(cosmic_gene_muts$impact == cosmic_fixed_binned_muts$mut_type[i] & cosmic_gene_muts$pos <= gene_coverage$pos[cosmic_fixed_binned_muts$idx[i] - ((fixed_bin_size)/2 - 1)] & cosmic_gene_muts$pos >= gene_coverage$pos[min(cosmic_fixed_binned_muts$idx[i] + (fixed_bin_size)/2,nrow(gene_coverage))]))
      }
    }
    
    if(nrow(cosmic_squamous_gene_muts) > 0){
      cosmic_squamous_fixed_binned_muts <- fixed_binned_muts

      for(i in 1:nrow(cosmic_squamous_fixed_binned_muts)){
        if(buccal_genes$strand[which(buccal_genes$gene == gene)] == 1){
          cosmic_squamous_fixed_binned_muts$count[i] <- length(which(cosmic_squamous_gene_muts$impact == cosmic_squamous_fixed_binned_muts$mut_type[i] & cosmic_squamous_gene_muts$pos >= gene_coverage$pos[cosmic_squamous_fixed_binned_muts$idx[i] - ((fixed_bin_size)/2 - 1)] & cosmic_squamous_gene_muts$pos <= gene_coverage$pos[min(cosmic_squamous_fixed_binned_muts$idx[i] + (fixed_bin_size)/2,nrow(gene_coverage))]))
        }else if(buccal_genes$strand[which(buccal_genes$gene == gene)] == -1){
          cosmic_squamous_fixed_binned_muts$count[i] <- length(which(cosmic_squamous_gene_muts$impact == cosmic_squamous_fixed_binned_muts$mut_type[i] & cosmic_squamous_gene_muts$pos <= gene_coverage$pos[cosmic_squamous_fixed_binned_muts$idx[i] - ((fixed_bin_size)/2 - 1)] & cosmic_squamous_gene_muts$pos >= gene_coverage$pos[min(cosmic_squamous_fixed_binned_muts$idx[i] + (fixed_bin_size)/2,nrow(gene_coverage))]))
        }
      }
    }
  }
  
  for(i in 1:nrow(fixed_binned_muts)){
    if(buccal_genes$strand[which(buccal_genes$gene == gene)] == 1){
      fixed_binned_muts$count[i] <- length(which(gene_muts$impact == fixed_binned_muts$mut_type[i] & gene_muts$pos >= gene_coverage$pos[fixed_binned_muts$idx[i] - ((fixed_bin_size)/2 - 1)] & gene_muts$pos <= gene_coverage$pos[min(fixed_binned_muts$idx[i] + (fixed_bin_size)/2,nrow(gene_coverage))]))
    }else if(buccal_genes$strand[which(buccal_genes$gene == gene)] == -1){
      fixed_binned_muts$count[i] <- length(which(gene_muts$impact == fixed_binned_muts$mut_type[i] & gene_muts$pos <= gene_coverage$pos[fixed_binned_muts$idx[i] - ((fixed_bin_size)/2 - 1)] & gene_muts$pos >= gene_coverage$pos[min(fixed_binned_muts$idx[i] + (fixed_bin_size)/2,nrow(gene_coverage))]))
    }
  }
  
  variable_binned_muts <- as.data.frame(array(data = NA, dim = c(ceiling(nrow(gene_coverage) / variable_bin_size) * length(mut_types),3)))
  colnames(variable_binned_muts) <- c("idx","mut_type","count")
  
  variable_binned_muts$idx <- rep(seq(from = (variable_bin_size)/2, to = (ceiling(nrow(gene_coverage) / variable_bin_size) * variable_bin_size) - (variable_bin_size)/2, by = variable_bin_size), each = length(mut_types))
  variable_binned_muts$mut_type <- rep(mut_types, times = ceiling(nrow(gene_coverage) / variable_bin_size))
  
  if(nrow(cosmic_gene_muts) > 0){
    cosmic_variable_binned_muts <- variable_binned_muts
    
    for(i in 1:nrow(cosmic_variable_binned_muts)){
      if(buccal_genes$strand[which(buccal_genes$gene == gene)] == 1){
        cosmic_variable_binned_muts$count[i] <- length(which(cosmic_gene_muts$impact == cosmic_variable_binned_muts$mut_type[i] & cosmic_gene_muts$pos >= gene_coverage$pos[cosmic_variable_binned_muts$idx[i] - ((variable_bin_size)/2 - 1)] & cosmic_gene_muts$pos <= gene_coverage$pos[min(cosmic_variable_binned_muts$idx[i] + (variable_bin_size)/2,nrow(gene_coverage))]))
      }else if(buccal_genes$strand[which(buccal_genes$gene == gene)] == -1){
        cosmic_variable_binned_muts$count[i] <- length(which(cosmic_gene_muts$impact == cosmic_variable_binned_muts$mut_type[i] & cosmic_gene_muts$pos <= gene_coverage$pos[cosmic_variable_binned_muts$idx[i] - ((variable_bin_size)/2 - 1)] & cosmic_gene_muts$pos >= gene_coverage$pos[min(cosmic_variable_binned_muts$idx[i] + (variable_bin_size)/2,nrow(gene_coverage))]))
      }
    }
    
    if(nrow(cosmic_squamous_gene_muts) > 0){
      cosmic_squamous_variable_binned_muts <- variable_binned_muts

      for(i in 1:nrow(cosmic_squamous_variable_binned_muts)){
        if(buccal_genes$strand[which(buccal_genes$gene == gene)] == 1){
          cosmic_squamous_variable_binned_muts$count[i] <- length(which(cosmic_squamous_gene_muts$impact == cosmic_squamous_variable_binned_muts$mut_type[i] & cosmic_squamous_gene_muts$pos >= gene_coverage$pos[cosmic_squamous_variable_binned_muts$idx[i] - ((variable_bin_size)/2 - 1)] & cosmic_squamous_gene_muts$pos <= gene_coverage$pos[min(cosmic_squamous_variable_binned_muts$idx[i] + (variable_bin_size)/2,nrow(gene_coverage))]))
        }else if(buccal_genes$strand[which(buccal_genes$gene == gene)] == -1){
          cosmic_squamous_variable_binned_muts$count[i] <- length(which(cosmic_squamous_gene_muts$impact == cosmic_squamous_variable_binned_muts$mut_type[i] & cosmic_squamous_gene_muts$pos <= gene_coverage$pos[cosmic_squamous_variable_binned_muts$idx[i] - ((variable_bin_size)/2 - 1)] & cosmic_squamous_gene_muts$pos >= gene_coverage$pos[min(cosmic_squamous_variable_binned_muts$idx[i] + (variable_bin_size)/2,nrow(gene_coverage))]))
        }
      }
    }
  }
  
  for(i in 1:nrow(variable_binned_muts)){
    if(buccal_genes$strand[which(buccal_genes$gene == gene)] == 1){
      variable_binned_muts$count[i] <- length(which(gene_muts$impact == variable_binned_muts$mut_type[i] & gene_muts$pos >= gene_coverage$pos[variable_binned_muts$idx[i] - ((variable_bin_size)/2 - 1)] & gene_muts$pos <= gene_coverage$pos[min(variable_binned_muts$idx[i] + (variable_bin_size)/2,nrow(gene_coverage))]))
    }else if(buccal_genes$strand[which(buccal_genes$gene == gene)] == -1){
      variable_binned_muts$count[i] <- length(which(gene_muts$impact == variable_binned_muts$mut_type[i] & gene_muts$pos <= gene_coverage$pos[variable_binned_muts$idx[i] - ((variable_bin_size)/2 - 1)] & gene_muts$pos >= gene_coverage$pos[min(variable_binned_muts$idx[i] + (variable_bin_size)/2,nrow(gene_coverage))]))
    }
  }
  
  interpro_files <- list.files("/lustre/scratch126/casm/team268im/al28/TwinsUK_Manuscript/InterProDomains/")
  domain_file <- interpro_files[which(grepl(paste0(gene,"_"),interpro_files) & !(grepl("DifferentIsoform",interpro_files)))]
  
  if(length(domain_file) == 1){
    domains <- read.table(file = paste0("/lustre/scratch126/casm/team268im/al28/TwinsUK_Manuscript/InterProDomains/",domain_file), sep = "\t", stringsAsFactors = F, header = T, quote = "")
    domains <- domains[which(domains$Source.Database == "interpro" & domains$Type %in% c("domain","repeat")),]
    
    if(nrow(domains) > 0){
      rows_to_split <- which(grepl(",",domains$Matches))
      
      if(length(rows_to_split) >= 1){
        for(i in 1:length(rows_to_split)){
          temp_rows <- as.data.frame(array(data = NA, dim = c(lengths(regmatches(domains$Matches[rows_to_split[i]], gregexpr(",", domains$Matches[rows_to_split[i]]))) + 1,ncol(domains))))
          colnames(temp_rows) <- colnames(domains)
          for(j in 1:nrow(temp_rows)){
            temp_rows[j,which(colnames(temp_rows) != "Matches")] <- domains[rows_to_split[i],which(colnames(domains) != "Matches"),]
          }
          temp_rows$Matches <- unlist(strsplit(domains$Matches[rows_to_split[i]], split = ","))
          
          domains <- rbind(domains,temp_rows)
        }
        
        domains <- domains[-(rows_to_split),]
      }
      
      domains$plot <- NA
      domains$aa_start <- NA
      domains$aa_end <- NA
      
      for(i in 1:nrow(domains)){
        domains$aa_start[i] <- as.numeric(unlist(strsplit(domains$Matches[i], split = "\\.."))[1])
        domains$aa_end[i] <- as.numeric(unlist(strsplit(domains$Matches[i], split = "\\.."))[2])
      }
      
      domains$idx_start <- (3 * (domains$aa_start - 1)) + 1
      domains$idx_end <- (3 * (domains$aa_end))
      
      if(unique(domains$Protein.Length) != ((nrow(gene_coverage) / 3) - 1)){
        domains$plot <- 0
      }else{
        domains$plot <- 1
      }
      
      domains$overlaps <- NA
      domains$encapsulates <- NA
      domains$num_encapsulated <- NA
      domains$high_deg_overlap <- NA
      
      for(i in 1:nrow(domains)){
        overlaps <- NULL
        encapsulates <- NULL
        high_deg_overlap <- NULL
        
        for(j in 1:nrow(domains)){
          if(i != j){
            per_aa_i <- seq(from = domains$aa_start[i], to = domains$aa_end[i], by = 1)
            per_aa_j <- seq(from = domains$aa_start[j], to = domains$aa_end[j], by = 1)
            
            if(length(which(per_aa_j %in% per_aa_i)) > 0){
              overlaps <- c(overlaps,j)
            }
            
            if(length(which(per_aa_j %in% per_aa_i)) == length(per_aa_j)){
              encapsulates <- c(encapsulates,j)
            }
            
            if(length(which(per_aa_j %in% per_aa_i)) > 0.8 * length(per_aa_j) & length(which(per_aa_i %in% per_aa_j)) > 0.8 * length(per_aa_i)){
              high_deg_overlap <- c(high_deg_overlap,j)
            }
          }
          
          domains$overlaps[i] <- paste(overlaps, collapse = ",")
          domains$encapsulates[i] <- paste(encapsulates, collapse = ",")
          domains$num_encapsulated[i] <- length(encapsulates)
          domains$high_deg_overlap[i] <- paste(high_deg_overlap, collapse =  ",")
          
          if(domains$plot[i] != 0){
            domains$plot[high_deg_overlap] <- 0
          }
        }
      }
      
      domains <- domains[order(domains$num_encapsulated, decreasing = T),]
      
      if(length(which(domains$plot ==1)) > 0){
        domains_plot <- ggplot(data = domains[which(domains$plot == 1),]) +
          geom_rect(xmin = min(gene_pos_transition_table$start_idx), xmax = max(gene_pos_transition_table$end_idx), ymin = -0.15, ymax = 0.15) +
          geom_rect(aes(xmin = idx_start, xmax = idx_end, fill = Name), color = "black", ymin = -0.5, ymax = 0.5) +
          scale_fill_brewer(palette = "Set3") +
          scale_x_continuous(expand = c(0,6), limits = c(min(gene_pos_transition_table$start_idx),max(gene_pos_transition_table$end_idx))) +
          scale_y_continuous(expand = c(0,0), limits = c(-0.5,0.5)) +
          theme_bw() +
          labs(fill = "Domain type", y = "Domains") +
          theme(panel.grid = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), panel.border = element_blank(), axis.title.y = element_text(angle = 0, vjust = 0.5), panel.background = element_rect(fill='transparent'))
      }
    }
  }
  
  fixed_binned_muts_limit <- ceiling(max(aggregate(fixed_binned_muts$count ~ fixed_binned_muts$idx, data=fixed_binned_muts, sum)[,2]) * 1.1)
  
  mut_fixed_histo_no_cosmic_plot <- ggplot() +
    geom_area(data = gene_coverage_smooth, aes(x = idx, y = (bin_average * fixed_binned_muts_limit / (ceiling((max(total))/10000) * 10000))), color = "grey70", fill = "grey", alpha = 0.25) +
    geom_bar(data = fixed_binned_muts, aes(x = idx, y = count, fill = mut_type), stat = "identity", width = fixed_bin_size) +
    scale_fill_manual(values = c("Synonymous" = "grey70", "Missense" = "cadetblue", "Nonsense" = "darkorchid4", "Essential_Splice" = "darkorchid2", "no-SNV" = "chocolate3","non_coding" = "burlywood1","Start_loss" = "firebrick","Stop_loss" = "red")) +
    theme_bw() +
    scale_x_continuous(expand = c(0,6)) +
    labs(x = "", fill = "Mutation type", y = "TwinsUK\nmutation\ncount", title = paste0(gene," - ",(fixed_bin_size / 3)," aa bins\n",sum(fixed_binned_muts$count)," TwinsUK mutations")) +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(axis.text.x = element_blank()) +
    theme(axis.title.x = element_blank()) +
    theme(axis.ticks.x = element_blank()) +
    theme(axis.title.y = element_text(angle = 0, vjust = 0.5)) +
    theme(legend.key=element_rect(colour="black")) +
    theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank(), axis.title.y = element_text(angle = 0)) +
    scale_y_continuous(expand = c(0,0), 
                       limits = c(0,fixed_binned_muts_limit),
                       breaks = pretty_breaks(),
                       sec.axis = sec_axis(~ . * ((ceiling(max(gene_coverage_smooth$total)/10000) * 10000)/(ceiling(max(aggregate(fixed_binned_muts$count ~ fixed_binned_muts$idx, data=fixed_binned_muts, sum)[,2]) * 1.1))), name = "Cumulative\nduplex coverage"))
  
  variable_binned_muts_limit <- ceiling(max(aggregate(variable_binned_muts$count ~ variable_binned_muts$idx, data=variable_binned_muts, sum)[,2]) * 1.1)
  
  mut_variable_histo_no_cosmic_plot <- ggplot() +
    geom_area(data = gene_coverage_smooth, aes(x = idx, y = (bin_average * variable_binned_muts_limit / (ceiling((max(total))/10000) * 10000))), color = "grey70", fill = "grey", alpha = 0.25) +
    geom_bar(data = variable_binned_muts, aes(x = idx, y = count, fill = mut_type), stat = "identity", width = variable_bin_size) +
    scale_fill_manual(values = c("Synonymous" = "grey70", "Missense" = "cadetblue", "Nonsense" = "darkorchid4", "Essential_Splice" = "darkorchid2", "no-SNV" = "chocolate3","non_coding" = "burlywood1","Start_loss" = "firebrick","Stop_loss" = "red")) +
    theme_bw() +
    scale_x_continuous(expand = c(0,6)) +
    labs(x = "", fill = "Mutation type", y = "TwinsUK\nmutation\ncount", title = paste0(gene," - ",(variable_bin_size / 3)," aa bins\n",sum(variable_binned_muts$count)," TwinsUK mutations")) +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(axis.text.x = element_blank()) +
    theme(axis.title.x = element_blank()) +
    theme(axis.ticks.x = element_blank()) +
    theme(axis.title.y = element_text(angle = 0, vjust = 0.5)) +
    theme(legend.key=element_rect(colour="black")) +
    theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank(), axis.title.y = element_text(angle = 0)) +
    scale_y_continuous(expand = c(0,0), 
                       limits = c(0,variable_binned_muts_limit),
                       breaks = pretty_breaks(),
                       sec.axis = sec_axis(~ . * ((ceiling(max(gene_coverage_smooth$total)/10000) * 10000)/(ceiling(max(aggregate(variable_binned_muts$count ~ variable_binned_muts$idx, data=variable_binned_muts, sum)[,2]) * 1.1))), name = "Cumulative\nduplex coverage"))
  
  if(nrow(gene_pos_transition_table) >= 2){
    exome_break_plot <- ggplot() +
      geom_rect(data = gene_pos_transition_table[seq(from = 1, to = nrow(gene_pos_transition_table), by = 2),], aes(xmin = start_idx, xmax = end_idx), fill = "firebrick", color = "black", ymin = 0, ymax = 1, alpha = 0.35) +
      geom_rect(data = gene_pos_transition_table[seq(from = 2, to = nrow(gene_pos_transition_table), by = 2),], aes(xmin = start_idx, xmax = end_idx), fill = "firebrick", color = "black", ymin = 0, ymax = 1, alpha = 0.65) +
      scale_x_continuous(expand = c(0,6)) +
      theme_bw() +
      scale_y_continuous(expand = c(0,0)) +
      theme(axis.text = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank()) +
      labs(y = "Exons") +
      theme(axis.title.y = element_text(angle = 0, vjust = 0.5)) +
      theme(panel.border = element_rect(color = "black", fill = NA))
  }else{
    exome_break_plot <- ggplot() +
      geom_rect(data = gene_pos_transition_table[1,], aes(xmin = start_idx, xmax = end_idx), fill = "firebrick", color = "black", ymin = 0, ymax = 1, alpha = 0.35) +
      scale_x_continuous(expand = c(0,6)) +
      theme_bw() +
      scale_y_continuous(expand = c(0,0)) +
      theme(axis.text = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank()) +
      labs(y = "Exons") +
      theme(axis.title.y = element_text(angle = 0, vjust = 0.5)) +
      theme(panel.border = element_rect(color = "black", fill = NA))
  }
  
  if(nrow(cosmic_gene_muts) > 0){
    
    cosmic_fixed_binned_muts_limit <- ceiling(max(aggregate(cosmic_fixed_binned_muts$count ~ cosmic_fixed_binned_muts$idx, data=cosmic_fixed_binned_muts, sum)[,2]) * 1.1)
    cosmic_fixed_breaks <- pretty_breaks()(-cosmic_fixed_binned_muts_limit:0)
    
    mut_fixed_histo_cosmic_plot <- ggplot() +
      geom_bar(data = cosmic_fixed_binned_muts, aes(x = idx, y = -(count), fill = mut_type), stat = "identity", width = fixed_bin_size) +
      scale_fill_manual(values = c("Synonymous" = "grey70", "Missense" = "cadetblue", "Nonsense" = "darkorchid4", "Essential_Splice" = "darkorchid2", "no-SNV" = "chocolate3","non_coding" = "burlywood1","Start_loss" = "firebrick","Stop_loss" = "red")) +
      theme_bw() +
      scale_x_continuous(expand = c(0,6)) +
      labs(x = "", fill = "Mutation type", y = "COSMIC\nmutation\ncount", caption = paste0(sum(cosmic_fixed_binned_muts$count)," COSMIC (all cancer types) mutations")) +
      theme(plot.caption = element_text(hjust=0.5, size=rel(1.2))) +
      theme(axis.text.x = element_blank()) +
      theme(axis.title.x = element_blank()) +
      theme(axis.ticks.x = element_blank()) +
      theme(axis.title.y = element_text(angle = 0, vjust = 0.5)) +
      theme(legend.position = "none") +
      theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank()) +
      scale_y_continuous(expand = c(0,0), 
                         limits = c(-cosmic_fixed_binned_muts_limit,0),
                         breaks = cosmic_fixed_breaks,
                         labels = -(cosmic_fixed_breaks))
    
    cosmic_variable_binned_muts_limit <- ceiling(max(aggregate(cosmic_variable_binned_muts$count ~ cosmic_variable_binned_muts$idx, data=cosmic_variable_binned_muts, sum)[,2]) * 1.1)
    cosmic_variable_breaks <- pretty_breaks()(-cosmic_variable_binned_muts_limit:0)
    
    mut_variable_histo_cosmic_plot <- ggplot() +
      geom_bar(data = cosmic_variable_binned_muts, aes(x = idx, y = -(count), fill = mut_type), stat = "identity", width = variable_bin_size) +
      scale_fill_manual(values = c("Synonymous" = "grey70", "Missense" = "cadetblue", "Nonsense" = "darkorchid4", "Essential_Splice" = "darkorchid2", "no-SNV" = "chocolate3","non_coding" = "burlywood1","Start_loss" = "firebrick","Stop_loss" = "red")) +
      theme_bw() +
      scale_x_continuous(expand = c(0,6)) +
      labs(x = "", fill = "Mutation type", y = "COSMIC\nmutation\ncount", caption = paste0(sum(cosmic_variable_binned_muts$count)," COSMIC (all cancer types) mutations")) +
      theme(plot.caption = element_text(hjust=0.5, size=rel(1.2))) +
      theme(axis.text.x = element_blank()) +
      theme(axis.title.x = element_blank()) +
      theme(axis.ticks.x = element_blank()) +
      theme(axis.title.y = element_text(angle = 0, vjust = 0.5)) +
      theme(legend.position = "none") +
      theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank()) +
      scale_y_continuous(expand = c(0,0), 
                         limits = c(-cosmic_variable_binned_muts_limit,0),
                         breaks = cosmic_variable_breaks,
                         labels = -(cosmic_variable_breaks))
    
    if(nrow(cosmic_squamous_gene_muts) > 0){
      cosmic_squamous_fixed_binned_muts_limit <- ceiling(max(aggregate(cosmic_squamous_fixed_binned_muts$count ~ cosmic_squamous_fixed_binned_muts$idx, data=cosmic_squamous_fixed_binned_muts, sum)[,2]) * 1.1)
      cosmic_squamous_fixed_breaks <- pretty_breaks()(-cosmic_squamous_fixed_binned_muts_limit:0)

      mut_fixed_histo_cosmic_squamous_plot <- ggplot() +
        geom_bar(data = cosmic_squamous_fixed_binned_muts, aes(x = idx, y = -(count), fill = mut_type), stat = "identity", width = fixed_bin_size) +
        scale_fill_manual(values = c("Synonymous" = "grey70", "Missense" = "cadetblue", "Nonsense" = "darkorchid4", "Essential_Splice" = "darkorchid2", "no-SNV" = "chocolate3","non_coding" = "burlywood1","Start_loss" = "firebrick","Stop_loss" = "red")) +
        theme_bw() +
        scale_x_continuous(expand = c(0,6)) +
        labs(x = "", fill = "Mutation type", y = "COSMIC\nmutation\ncount", caption = paste0(sum(cosmic_squamous_fixed_binned_muts$count)," COSMIC (squamous cell carcinoma) mutations")) +
        theme(plot.caption = element_text(hjust=0.5, size=rel(1.2))) +
        theme(axis.text.x = element_blank()) +
        theme(axis.title.x = element_blank()) +
        theme(axis.ticks.x = element_blank()) +
        theme(axis.title.y = element_text(angle = 0, vjust = 0.5)) +
        theme(legend.position = "none") +
        theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank()) +
        scale_y_continuous(expand = c(0,0),
                           limits = c(-cosmic_squamous_fixed_binned_muts_limit,0),
                           breaks = cosmic_squamous_fixed_breaks,
                           labels = -(cosmic_squamous_fixed_breaks))

      cosmic_squamous_variable_binned_muts_limit <- ceiling(max(aggregate(cosmic_squamous_variable_binned_muts$count ~ cosmic_squamous_variable_binned_muts$idx, data=cosmic_squamous_variable_binned_muts, sum)[,2]) * 1.1)
      cosmic_squamous_variable_breaks <- pretty_breaks()(-cosmic_squamous_variable_binned_muts_limit:0)

      mut_variable_histo_cosmic_squamous_plot <- ggplot() +
        geom_bar(data = cosmic_squamous_variable_binned_muts, aes(x = idx, y = -(count), fill = mut_type), stat = "identity", width = variable_bin_size) +
        scale_fill_manual(values = c("Synonymous" = "grey70", "Missense" = "cadetblue", "Nonsense" = "darkorchid4", "Essential_Splice" = "darkorchid2", "no-SNV" = "chocolate3","non_coding" = "burlywood1","Start_loss" = "firebrick","Stop_loss" = "red")) +
        theme_bw() +
        scale_x_continuous(expand = c(0,6)) +
        labs(x = "", fill = "Mutation type", y = "COSMIC\nmutation\ncount", caption = paste0(sum(cosmic_squamous_variable_binned_muts$count)," COSMIC (squamous cell carcinoma) mutations")) +
        theme(plot.caption = element_text(hjust=0.5, size=rel(1.2))) +
        theme(axis.text.x = element_blank()) +
        theme(axis.title.x = element_blank()) +
        theme(axis.ticks.x = element_blank()) +
        theme(axis.title.y = element_text(angle = 0, vjust = 0.5)) +
        theme(legend.position = "none") +
        theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank()) +
        scale_y_continuous(expand = c(0,0),
                           limits = c(-cosmic_squamous_variable_binned_muts_limit,0),
                           breaks = cosmic_squamous_variable_breaks,
                           labels = -(cosmic_squamous_variable_breaks))
    }
  }
  
  if(!(exists("domains_plot"))){
    fixed_no_cosmic_no_domain_plot <- mut_fixed_histo_no_cosmic_plot / plot_spacer() / exome_break_plot + plot_layout(heights = c(20,-1.05,1), guides = "collect")
    ggsave(fixed_no_cosmic_no_domain_plot, file = paste0(gene,"_coding_fixed_no_cosmic_no_domain.pdf"),width = 15, height = 5)
    
    if(exists("mut_fixed_histo_cosmic_plot")){
      fixed_cosmic_all_cancer_no_domain_plot <- mut_fixed_histo_no_cosmic_plot / plot_spacer() / exome_break_plot / plot_spacer() / mut_fixed_histo_cosmic_plot + plot_layout(heights = c(20,-1.05,1,-1.05,20), guides = "collect")
      ggsave(fixed_cosmic_all_cancer_no_domain_plot, file = paste0(gene,"_coding_fixed_cosmic_all_no_domain.pdf"),width = 15, height = 10)
    }
    if(exists("mut_fixed_histo_cosmic_squamous_plot")){
      fixed_cosmic_squamous_cancer_no_domain_plot <- mut_fixed_histo_no_cosmic_plot / plot_spacer() / exome_break_plot / plot_spacer() / mut_fixed_histo_cosmic_squamous_plot + plot_layout(heights = c(20,-1.05,1,-1.05,20), guides = "collect")
      ggsave(fixed_cosmic_squamous_cancer_no_domain_plot, file = paste0(gene,"_coding_fixed_cosmic_squamous_no_domain.pdf"),width = 15, height = 10)
    }
  }else {
    fixed_no_cosmic_domain_plot <- mut_fixed_histo_no_cosmic_plot / plot_spacer() / exome_break_plot / domains_plot + plot_layout(heights = c(20,-1.05,1,1), guides = "collect")
    ggsave(fixed_no_cosmic_domain_plot, file = paste0(gene,"_coding_fixed_no_cosmic_domain.pdf"),width = 15, height = 5)
    
    if(exists("mut_fixed_histo_cosmic_plot")){
      fixed_cosmic_all_cancer_domain_plot <- mut_fixed_histo_no_cosmic_plot / plot_spacer() / exome_break_plot / domains_plot / plot_spacer() / mut_fixed_histo_cosmic_plot + plot_layout(heights = c(20,-1.05,1,1,-1.05,20), guides = "collect")
      ggsave(fixed_cosmic_all_cancer_domain_plot, file = paste0(gene,"_coding_fixed_cosmic_all_domain.pdf"),width = 15, height = 10)
    }
    if(exists("mut_fixed_histo_cosmic_squamous_plot")){
      fixed_cosmic_squamous_cancer_domain_plot <- mut_fixed_histo_no_cosmic_plot / plot_spacer() / exome_break_plot / domains_plot / plot_spacer() / mut_fixed_histo_cosmic_squamous_plot + plot_layout(heights = c(20,-1.05,1,1,-1.05,20), guides = "collect")
      ggsave(fixed_cosmic_squamous_cancer_domain_plot, file = paste0(gene,"_coding_fixed_cosmic_squamous_domain.pdf"),width = 15, height = 10)
    }
  }
  
  if(variable_bin_size != fixed_bin_size){
    if(!(exists("domains_plot"))){
      variable_no_cosmic_no_domain_plot <- mut_variable_histo_no_cosmic_plot / plot_spacer() / exome_break_plot + plot_layout(heights = c(20,-1.05,1), guides = "collect")
      ggsave(variable_no_cosmic_no_domain_plot, file = paste0(gene,"_coding_variable_no_cosmic_no_domain.pdf"),width = 15, height = 5)
      
      if(exists("mut_fixed_histo_cosmic_plot")){
        variable_cosmic_all_cancer_no_domain_plot <- mut_variable_histo_no_cosmic_plot / plot_spacer() / exome_break_plot / plot_spacer() / mut_variable_histo_cosmic_plot + plot_layout(heights = c(20,-1.05,1,-1.05,20), guides = "collect")
        ggsave(variable_cosmic_all_cancer_no_domain_plot, file = paste0(gene,"_coding_variable_cosmic_all_no_domain.pdf"),width = 15, height = 10)
      }
      if(exists("mut_fixed_histo_cosmic_squamous_plot")){
        variable_cosmic_squamous_cancer_no_domain_plot <- mut_variable_histo_no_cosmic_plot / plot_spacer() / exome_break_plot / plot_spacer() / mut_variable_histo_cosmic_squamous_plot + plot_layout(heights = c(20,-1.05,1,-1.05,20), guides = "collect")
        ggsave(variable_cosmic_squamous_cancer_no_domain_plot, file = paste0(gene,"_coding_variable_cosmic_squamous_no_domain.pdf"),width = 15, height = 10)
      }
    }else {
      variable_no_cosmic_domain_plot <- mut_variable_histo_no_cosmic_plot / plot_spacer() / exome_break_plot / domains_plot + plot_layout(heights = c(20,-1.05,1,1), guides = "collect")
      ggsave(variable_no_cosmic_domain_plot, file = paste0(gene,"_coding_variable_no_cosmic_domain.pdf"),width = 15, height = 5)
      
      if(exists("mut_fixed_histo_cosmic_plot")){
        variable_cosmic_all_cancer_domain_plot <- mut_variable_histo_no_cosmic_plot / plot_spacer() / exome_break_plot / domains_plot / plot_spacer() / mut_variable_histo_cosmic_plot + plot_layout(heights = c(20,-1.05,1,1,-1.05,20), guides = "collect")
        ggsave(variable_cosmic_all_cancer_domain_plot, file = paste0(gene,"_coding_variable_cosmic_all_domain.pdf"),width = 15, height = 10)
      }
      if(exists("mut_fixed_histo_cosmic_squamous_plot")){
        variable_cosmic_squamous_cancer_domain_plot <- mut_variable_histo_no_cosmic_plot / plot_spacer() / exome_break_plot / domains_plot / plot_spacer() / mut_variable_histo_cosmic_squamous_plot + plot_layout(heights = c(20,-1.05,1,1,-1.05,20), guides = "collect")
        ggsave(variable_cosmic_squamous_cancer_domain_plot, file = paste0(gene,"_coding_variable_cosmic_squamous_domain.pdf"),width = 15, height = 10)
      }
    }
  }
}

# cosmic_all_cancer_full_plot <- mut_histo_no_cosmic_plot / plot_spacer() / exome_break_plot / domains_plot / plot_spacer() / mut_histo_cosmic_plot + plot_layout(heights = c(20,-1.05,1,1,-1.05,20), guides = "collect")
# ggsave(cosmic_all_cancer_full_plot, file = paste0(gene,"_coding_cosmic_all.pdf"),width = 15, height = 10)
# 
# cosmic_squamous_cancer_full_plot <- mut_histo_no_cosmic_plot / plot_spacer() / exome_break_plot / domains_plot / plot_spacer() / mut_histo_cosmic_squamous_plot + plot_layout(heights = c(20,-1.05,1,1,-1.05,20), guides = "collect")
# ggsave(cosmic_squamous_cancer_full_plot, file = paste0(gene,"_coding_cosmic_squamous_only.pdf"),width = 15, height = 10)
# 
# if(nrow(cosmic_gene_muts) > 0){
#   binned_muts_limit <- ceiling(max(aggregate(binned_muts$count ~ binned_muts$idx, data=binned_muts, sum)[,2]) * 1.1)
#   cosmic_binned_muts_limit <- ceiling(max(aggregate(cosmic_binned_muts$count ~ cosmic_binned_muts$idx, data=cosmic_binned_muts, sum)[,2]) * 1.1)
#   plot_limit <- min(binned_muts_limit,cosmic_binned_muts_limit)
#   binned_muts_ratio <- 1 / (binned_muts_limit / plot_limit)
#   cosmic_binned_muts_ratio <- 1 / (cosmic_binned_muts_limit / plot_limit)
#   cosmic_breaks <- pretty_breaks()((-1 * cosmic_binned_muts_limit):0)
#   cosmic_breaks <- cosmic_breaks[which(cosmic_breaks > (-1 * cosmic_binned_muts_limit))]
#   binned_muts_breaks <- pretty_breaks()(0:binned_muts_limit)
#   binned_muts_breaks <- binned_muts_breaks[which(binned_muts_breaks < binned_muts_limit)]
#   break_list <- c((cosmic_breaks * cosmic_binned_muts_ratio),(binned_muts_breaks * binned_muts_ratio))
#   
#   mut_histo_full_cosmic_plot <- ggplot() +
#     geom_bar(data = binned_muts, aes(x = idx, y = (count * binned_muts_ratio), fill = mut_type), stat = "identity", width = bin_size) +
#     geom_bar(data = cosmic_binned_muts, aes(x = idx, y = (count * -1 * cosmic_binned_muts_ratio), fill = mut_type), stat = "identity", width = bin_size) +
#     scale_fill_manual(values = c("Synonymous" = "grey70", "Missense" = "cadetblue", "Nonsense" = "darkorchid4", "Essential_Splice" = "darkorchid2", "no-SNV" = "chocolate3","non_coding" = "burlywood1","Stop_loss" = "red")) +
#     theme_bw() +
#     scale_x_continuous(expand = c(0,0), breaks = gene_pos_transition_table$start_idx[c(2:nrow(gene_pos_transition_table))]) +
#     labs(x = "", fill = "Type", y = "Mutation count", title = paste0(gene," - ",sum(binned_muts$count)," Targeted NanoSeq mutations; ",sum(cosmic_binned_muts$count)," COSMIC mutations (all cancers)")) +
#     theme(plot.title = element_text(hjust = 0.5)) +
#     theme(axis.text.x = element_blank()) +
#     theme(axis.ticks.length=unit(.25, "cm")) +
#     theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank()) +
#     scale_y_continuous(expand = c(0,0), 
#                        limits = c((-1 * plot_limit),(1 * plot_limit)),
#                        breaks = break_list,
#                        labels = c((-1 * cosmic_breaks), binned_muts_breaks))
#   ggsave(mut_histo_full_cosmic_plot, filename = paste0(gene,"_coding_cosmic_all.pdf"), width = 20, height = 5)
# }
# 
# binned_muts_limit <- ceiling(max(aggregate(binned_muts$count ~ binned_muts$idx, data=binned_muts, sum)[,2]) * 1.1)
# cosmic_squamous_binned_muts_limit <- ceiling(max(aggregate(cosmic_squamous_binned_muts$count ~ cosmic_squamous_binned_muts$idx, data=cosmic_squamous_binned_muts, sum)[,2]) * 1.1)
# plot_limit <- min(binned_muts_limit,cosmic_squamous_binned_muts_limit)
# binned_muts_ratio <- 1 / (binned_muts_limit / plot_limit)
# cosmic_squamous_binned_muts_ratio <- 1 / (cosmic_squamous_binned_muts_limit / plot_limit)
# cosmic_squamous_breaks <- pretty_breaks()((-1 * cosmic_squamous_binned_muts_limit):0)
# cosmic_squamous_breaks <- cosmic_squamous_breaks[which(cosmic_squamous_breaks > (-1 * cosmic_squamous_binned_muts_limit))]
# binned_muts_breaks <- pretty_breaks()(0:binned_muts_limit)
# binned_muts_breaks <- binned_muts_breaks[which(binned_muts_breaks < binned_muts_limit)]
# break_list <- c((cosmic_squamous_breaks * cosmic_squamous_binned_muts_ratio),(binned_muts_breaks * binned_muts_ratio))
# 
# mut_histo_squamous_cosmic_plot <- ggplot() +
#   geom_bar(data = binned_muts, aes(x = idx, y = (count * binned_muts_ratio), fill = mut_type), stat = "identity", width = bin_size) +
#   geom_bar(data = cosmic_squamous_binned_muts, aes(x = idx, y = (count * -1 * cosmic_squamous_binned_muts_ratio), fill = mut_type), stat = "identity", width = bin_size) +
#   scale_fill_manual(values = c("Synonymous" = "grey70", "Missense" = "cadetblue", "Nonsense" = "darkorchid4", "Essential_Splice" = "darkorchid2", "no-SNV" = "chocolate3","non_coding" = "burlywood1","Stop_loss" = "red")) +
#   theme_bw() +
#   scale_x_continuous(expand = c(0,0), breaks = gene_pos_transition_table$start_idx[c(2:nrow(gene_pos_transition_table))]) +
#   labs(x = "", fill = "Type", y = "Mutation count", title = paste0(gene," - ",sum(binned_muts$count)," Targeted NanoSeq mutations; ",sum(cosmic_squamous_binned_muts$count)," COSMIC mutations (squamous carcinoma only)")) +
#   theme(plot.title = element_text(hjust = 0.5)) +
#   theme(axis.text.x = element_blank()) +
#   theme(axis.ticks.length=unit(.25, "cm")) +
#   theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank()) +
#   scale_y_continuous(expand = c(0,0), 
#                      limits = c((-1 * plot_limit),(1 * plot_limit)),
#                      breaks = break_list,
#                      labels = c((-1 * cosmic_squamous_breaks), binned_muts_breaks))
# ggsave(mut_histo_squamous_cosmic_plot, filename = paste0(gene,"_coding_cosmic_squamous_only.pdf"), width = 20, height = 5)
#
# mut_histo_faceted_plot <- ggplot(data = binned_muts[which(binned_muts$mut_type != "Stop_loss"),]) +
#   geom_bar(aes(x = idx, y = count, fill = mut_type), stat = "identity", width = bin_size) +
#   scale_fill_manual(values = c("Synonymous" = "grey70", "Missense" = "cadetblue", "Nonsense" = "darkorchid4", "Essential_Splice" = "darkorchid2", "no-SNV" = "chocolate3")) +
#   theme_bw() +
#   scale_x_continuous(expand = c(0,0), breaks = gene_pos_transition_table$start_idx[c(2:nrow(gene_pos_transition_table))]) +
#   labs(x = "", fill = "Type", y = "Mutation count", title = paste0(gene," - ",sum(binned_muts$count[which(binned_muts$mut_type == "Stop_loss")])," stop loss mutations")) +
#   facet_wrap(~mut_type,ncol=1, scales = "free") +
#   theme(plot.title = element_text(hjust = 0.5)) +
#   theme(axis.text.x = element_blank()) +
#   theme(axis.ticks.length=unit(.25, "cm")) +
#   theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank()) +
#   theme(strip.text = element_blank())
# ggsave(mut_histo_faceted_plot, filename = paste0(gene,"_coding_faceted.pdf"), width = 30, height = 18)

# coverage_plot <- ggplot(data = gene_coverage) +
#   geom_bar(aes(x = idx, y = total), stat = "identity", fill = "gray50") +
#   scale_fill_manual(values = c("1" = "firebrick","0" = "steelblue")) +
#   theme_bw() +
#   scale_x_continuous(expand = c(0,0), breaks = gene_pos_transition_table$start_idx[c(2:nrow(gene_pos_transition_table))]) +
#   labs(fill = "Coding", x = "", y = "Total duplex coverage", title = gene) +
#   theme(plot.title = element_text(hjust = 0.5)) +
#   theme(axis.text.x = element_blank()) +
#   theme(axis.ticks.length=unit(.25, "cm")) +
#   theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank()) +
#   scale_y_continuous(expand = c(0,0), limits = c(0,ceiling((max(gene_coverage$total))/10000) * 10000))
# ggsave(coverage_plot, filename = paste0(gene,"_coding_coverage.pdf"), width = 20, height = 5)


