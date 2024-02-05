nextflow.enable.dsl=2

// expected command line argyments
params.mappings = null
params.gene_order_file = null
params.out_dir = './'

// perform infercnv on each sample
process infercnv {
  tag "${sample_id}"
  label "week16core10gb"

  publishDir (
    path: "${params.out_dir}/${id}/${sample_id}",
    mode: 'copy',
    pattern: "{infercnv*,preliminary.infercnv_obj,run.final.infercnv_obj}"
    )
    
  input:
    tuple val(sample_id), path(raw_counts_matrix), path(annotations), val(id)
    
  output:
    tuple val(sample_id), path('*'), val(id)

  script:
    """
    #!/usr/bin/env /nfs/users/nfs_a/at31/miniforge3/envs/jupy/bin/Rscript
    
    # libraries
    library(magrittr)
    library(Matrix)

    # subset matrix to cells in the sample
    annotations <- 
        readr::read_tsv("${annotations}", col_names = c('cell', 'celltype'))
    raw_counts_matrix <- 
        readRDS("${raw_counts_matrix}")
    raw_counts_matrix <-
      raw_counts_matrix[, colnames(raw_counts_matrix) %in% annotations\$cell]

    # create infercnv object
    infercnv_obj <-
        infercnv::CreateInfercnvObject(
            raw_counts_matrix = raw_counts_matrix,
            annotations_file = "${annotations}",
            delim = '\t',
            gene_order_file = "${params.gene_order_file}",
            ref_group_names = NULL)
    
    # run infercnv 
    # (parameters taken from mp34)
    # (/lustre/scratch126/casm/team268im/mp34/analysis/synovium_scRNA/infercnv_analysis.R)
    options(scipen = 100)
    infercnv_obj <-
        infercnv::run(
            infercnv_obj,
            out_dir = './',
            num_threads = 8, 
            cluster_by_groups = TRUE,
            analysis_mode = c('subclusters'),
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
            
            # sets midpoint for logistic
            sd_amplifier = 0.65,
            
            # https://github.com/harbourlab/UPhyloplot2
            tumor_subcluster_partition_method = 'random_trees'
            )
            
    # apply median filtering
    infercnv_obj <- 
      infercnv_obj %>%
      infercnv::apply_median_filtering()
      
    # plot smoothed output
    infercnv_obj %>%
      infercnv::plot_cnv(
        output_filename = 'infercnv.median_filtered',
        x.center = 1,
        colour_safe_pal = F)
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

    // load the mappings
    mappings = Channel.fromPath(params.mappings, checkIfExists: true) \
        | splitCsv(header:true, strip:true) \
        | map { row -> 
            tuple(row.sample_id, 
                  file(row.raw_counts_matrix, checkIfExists: true), 
                  file(row.annotations, checkIfExists: true), 
                  row.id) } \
        | infercnv

}
