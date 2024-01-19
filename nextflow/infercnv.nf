nextflow.enable.dsl=2

// Expected command line argyments
params.mappings = null
params.gene_order_file = null
params.output_dir = './'

// perform infercnv on each sample
process infercnv {
  tag "${sample_id}"
  label "week16core10gb"
  errorStrategy = {task.attempt <= 1 ? 'retry' : 'ignore'}

  publishDir "${params.output_dir}/${id}/${sample_id}",
    mode: 'copy',
    pattern: "infercnv*"

  input:
    tuple val(sample_id), path(raw_counts_matrix), path(annotations), val(id)
    
  output:
    tuple val(sample_id), path('infercnv*'), val(id)

  script:
    """
    #!/usr/bin/env /nfs/users/nfs_a/at31/miniforge3/envs/jupy/bin/Rscript
    library(magrittr)

    # subset matrix to cells in the sample
    annotations <- 
        readr::read_tsv("${annotations}", col_names = c('cell', 'celltype'))
    raw_counts_matrix <- 
        readRDS("${raw_counts_matrix}")[, annotations\$cell]    

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
    infercnv_obj <-
        infercnv::run(
            infercnv_obj,
            cutoff = 0.1, # 0.1 for 10x-genomics
            out_dir = './',
            num_threads = 8, 
            window_length = 151, # 51 is too small
            cluster_by_groups = TRUE,
            HMM = TRUE, # turn on to auto-run the HMM prediction of CNV levels
            HMM_transition_prob = 1e-6,
            HMM_report_by = c("subcluster"),
            analysis_mode = c('subclusters'),
            denoise = T,
            sd_amplifier = 0.65, # sets midpoint for logistic
            noise_logistic = F,
            resume_mode = T)

    # save output infercnv object
    saveRDS(infercnv_obj, 'infercnv_obj.rds')
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
            tuple(row.sample_id, row.raw_counts_matrix, row.annotations, row.id) } \
        | infercnv

}
