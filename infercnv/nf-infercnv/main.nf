#!/usr/bin/env nextflow

// using DSL-2
nextflow.enable.dsl=2

// all of the default parameters are being set in `nextflow.config`

// import functions / modules / subworkflows / workflows
include { validateParameters; paramsHelp; paramsSummaryLog; samplesheetToList } from 'plugin/nf-schema'
import groovy.json.JsonOutput

// print help message, supply typical command line usage for the pipeline
if (params.help) {
  log.info paramsHelp("nextflow run nf-infercnv --mappings mappings.csv --annotations annotations.tsv --outdir out/")
  exit 0
}

// save params
process save_params {
  publishDir path: "${params.out_dir}/", mode: "copy"

  output:
    path "params.json"

  script:
    "echo '${JsonOutput.toJson(params)}' > params.json"
}

// perform infercnv on each sample
process infercnv {
  container "trinityctat/infercnv"
  tag "${meta.id}"
  label "yday16core60gb"
  errorStrategy = "retry"
  publishDir(
    path: "${params.out_dir}/${meta.id}", 
    mode: "copy",
    pattern: "out/*",
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
                quote = F, row.names = F, sep = "\\t")
    
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
    options(scipen = 100)
    infercnv_obj <-
        infercnv::run(
            infercnv_obj,
            out_dir = "out/",
            num_threads = ${task.cpus},
            resume_mode = TRUE,
            cutoff = ${params.cutoff},
            min_cells_per_gene = ${params.min_cells_per_gene},
            output_format = "${params.output_format}",

            # Smoothing options
            window_length = ${params.window_length},
            smooth_method = "${params.smooth_method}",
            num_ref_groups = ${params.num_ref_groups ?: 'NULL'},
            ref_subtract_use_mean_bounds = as.logical("${params.ref_subtract_use_mean_bounds}"),
            cluster_by_groups = as.logical("${params.cluster_by_groups}"),
            cluster_references = as.logical("${params.cluster_references}"),
            k_obs_groups = ${params.k_obs_groups},
            hclust_method = "${params.hclust_method}",
            max_centered_threshold = ${params.max_centered_threshold},
            scale_data = as.logical("${params.scale_data}"),
            
            # Downstream analyses (HMM or non-DE-masking) based on tumour subclusters
            HMM = as.logical("${params.HMM}"),
            HMM_transition_prob = ${params.HMM_transition_prob},
            HMM_report_by = c("${params.HMM_report_by}"),
            HMM_type = "${params.HMM_type}",
            HMM_i3_pval = ${params.HMM_i3_pval},
            HMM_i3_use_KS = as.logical("${params.HMM_i3_use_KS}"),
            
            # Filtering low confidence HMM predictions via BayesNet P(Normal)
            BayesMaxPNormal = ${params.BayesMaxPNormal},
            reassignCNVs = as.logical("${params.reassignCNVs}"),
            
            # Tumour subclustering
            analysis_mode = c("${params.analysis_mode}"),
            tumor_subcluster_partition_method = "${params.tumor_subcluster_partition_method}",
            tumor_subcluster_pval = ${params.tumor_subcluster_pval},
            k_nn = ${params.k_nn},
            leiden_method = "${params.leiden_method}",
            leiden_function = "${params.leiden_function}",
            leiden_resolution = ${params.leiden_resolution},
            leiden_method_per_chr = "${params.leiden_method_per_chr}",
            leiden_function_per_chr = "${params.leiden_function_per_chr}",
            leiden_resolution_per_chr = ${params.leiden_resolution_per_chr},
            per_chr_hmm_subclusters = as.logical("${params.per_chr_hmm_subclusters}"),
            per_chr_hmm_subclusters_references = as.logical("${params.per_chr_hmm_subclusters_references}"),
            z_score_filter = ${params.z_score_filter},
            
            # Denoising
            denoise = as.logical("${params.denoise}"),
            noise_filter = ${params.noise_filter ?: 'NA'},
            sd_amplifier = ${params.sd_amplifier},
            noise_logistic = as.logical("${params.noise_logistic}"),
            
            # Outlier pruning
            outlier_method_bound = "${params.outlier_method_bound}",
            outlier_lower_bound = ${params.outlier_lower_bound ?: 'NA'},
            outlier_upper_bound = ${params.outlier_upper_bound ?: 'NA'},

            # Miscellandeous
            final_scale_limits = "${params.final_scale_limits ?: 'NULL'}",
            final_center_val = "${params.final_center_val ?: 'NULL'}",
            debug = as.logical("${params.debug}"),
            plot_steps = as.logical("${params.plot_steps}"),
            inspect_subclusters = as.logical("${params.inspect_subclusters}"),
            png_res = "${params.png_res}",
            no_plot = as.logical("${params.no_plot}"),
            no_prelim_plot = as.logical("${params.no_prelim_plot}"),
            write_expr_matrix = as.logical("${params.write_expr_matrix}"),
            write_phylo = as.logical("${params.write_phylo}"),
            plot_chr_scale = as.logical("${params.plot_chr_scale}"),
            useRaster = as.logical("${params.useRaster}"),
            plot_probabilities = as.logical("${params.plot_probabilities}"),
            save_rds = as.logical("${params.save_rds}"),
            save_final_rds = as.logical("${params.save_final_rds}")
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
  
  // save params
  save_params()

  // run infercnv
  infercnv(ch_mappings, ch_annotations, ch_gene_order_file)
  
}
