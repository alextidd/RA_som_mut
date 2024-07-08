#!/usr/bin/env nextflow

// using DSL-2
nextflow.enable.dsl=2

// all of the default parameters are being set in `nextflow.config`

// import functions / modules / subworkflows / workflows
include { validateParameters; paramsHelp; paramsSummaryLog; samplesheetToList } from 'plugin/nf-schema'

// print help message, supply typical command line usage for the pipeline
if (params.help) {
  log.info paramsHelp("nextflow run nf-infercnv --mappings mappings.csv --annotations annotations.tsv --outdir out/")
  exit 0
}

// perform infercnv on each sample
process infercnv {
  container "trinityctat/infercnv"
  tag "${meta.id}"
  label "week16core60gb"
  errorStrategy = "retry"
  publishDir(
    path: "${params.outdir}/${meta.id}", 
    mode: "copy",
    saveAs: { fn -> fn.substring(fn.lastIndexOf("/")+1) }
  ) 
  
  input:
    tuple val(meta), path(raw_counts_matrix)
    path(annotations)
    path(gene_order_file)
    
  output:
    tuple val(meta), path("out/*")

  script:
    """
    #!/usr/bin/env Rscript
    
    # libraries
    library(magrittr)
    library(Matrix)
    
    # generate annotations file
    annotations <- read.delim("${annotations}")
    annotations <- annotations[annotations\$id == "${meta.id}",
                               c("cell", "${params.annotation_col}")]
    write.table(annotations, "formatted_annotations.tsv", col.names = F, 
                quote = F, row.names = F, sep = "\t")
    
    # subset matrix to cells in the sample
    raw_counts_matrix <- readRDS("${raw_counts_matrix}")
    raw_counts_matrix <-
      raw_counts_matrix[, colnames(raw_counts_matrix) %in% annotations\$cell]

    # create infercnv object
    infercnv_obj <-
        infercnv::CreateInfercnvObject(
            raw_counts_matrix = raw_counts_matrix,
            annotations_file = "formatted_annotations.tsv",
            gene_order_file = "${gene_order_file}",
            ref_group_names = NULL)
    
    # run infercnv 
    # (parameters recommended by mp34@sanger.ac.uk)
    # (/lustre/scratch126/casm/team268im/mp34/analysis/synovium_scRNA/infercnv_analysis.R)
    options(scipen = 100)
    infercnv_obj <-
        infercnv::run(
            infercnv_obj,
            out_dir = "out/",
            num_threads = 8, 
            cluster_by_groups = as.logical("${params.cluster_by_groups}"),
            analysis_mode = c("${params.analysis_mode}"),
            denoise = T,
            noise_logistic = F,
            resume_mode = T,
            
            # window of 51 is too small
            window_length = 151,              

            # cutoff of 0.1 recommended for 10x-genomics
            cutoff = 0.1,                     
            
            # turn on to auto-run the HMM prediction of CNV levels
            HMM = TRUE,                       
            HMM_transition_prob = 1e-6,
            HMM_report_by = c("subcluster"),
            
            # # https://github.com/harbourlab/UPhyloplot2
            # tumor_subcluster_partition_method = "random_trees",
            
            # sets midpoint for logistic
            sd_amplifier = 0.65
            )
            
    # apply median filtering
    infercnv_obj <- 
      infercnv_obj %>%
      infercnv::apply_median_filtering()
      
    # plot smoothed output
    infercnv_obj %>%
      infercnv::plot_cnv(
        output_filename = "out/infercnv.median_filtered",
        x.center = 1,
        color_safe_pal = F)
      
    # write metadata table
    infercnv::add_to_seurat(
      infercnv_output_path = "out/")
    """
}

// main workflow
workflow {

  // check inputs
  if (params.mappings == null) {
    error "Please provide a mappings CSV file via --mappings"
  }
  if (params.gene_order_file == null) {
    error "Please provide a gene order file via --gene_order_file"
  }
  if (params.annotations == null) {
    error "Please provide annotations TSV file via --annotations"
  }

  // get annotations file
  Channel.fromPath(params.annotations, checkIfExists: true) 
  | set { ch_annotations }

  // get gene order file
  Channel.fromPath(params.gene_order_file, checkIfExists: true)
  | set { ch_gene_order_file }
  
  // get metadata + bam paths
  Channel.fromPath(params.mappings, checkIfExists: true) 
  | splitCsv(header: true)
  | map { row ->
    meta = row.subMap("id")
    [meta, file(row.raw_counts_matrix, checkIfExists: true)]
  }
  | set { ch_mappings }
  
  // run infercnv
  infercnv(ch_mappings, ch_annotations, ch_gene_order_file)
  
}
