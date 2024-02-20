nextflow.enable.dsl=2

// command line arguments
params.help             = false
params.mappings         = null
params.annotations      = null
params.gene_order_file  = null
params.out_dir          = './'
params.annotation_col   = 'celltype'

// help
if (params.help) {
  help = \
  """
  |infercnv.nf: run inferCNV
  |
  |Required arguments:
  |   --mappings        Path to the mappings CSV file, with columns `id` and 
  |                     `raw_counts_matrix`.
  |   --annotations     Path to the annotations TSV file, with columns `id`,
  |                     `cell`, and the annotation column (specify with the 
  |                     --annotation_col argument, default: 'celltype')
  |   --gene_order_file Path to the headerless gene order TSV file, with columns 
  |                     `gene`, `chr`, `start`, `stop`. See the inferCNV wiki 
  |                     for details.
  |
  |Optional arguments:
  |   --out_dir         Path to output directory. default is `out/`.
  |                     [default: ${params.out_dir}]
  |   --annotation_col  Column of the annotations file to use for defining 
  |                     groups of cells.
  |                     [default: ${params.annotation_col}]  
  """.stripMargin()
  
  // print help and exit
  println(help)
  exit(0)
}

// perform infercnv on each sample
process infercnv {
  tag "${meta.id}"
  label 'week16core60gb'
  errorStrategy = 'retry'
  publishDir(
    path: "${params.out_dir}/${meta.id}", 
    mode: 'copy',
    saveAs: { fn -> fn.substring(fn.lastIndexOf('/')+1) }
  ) 
  
  input:
    tuple val(meta), path(raw_counts_matrix)
    
  output:
    tuple val(meta), path("out/*")

  script:
    """
    #!/usr/bin/env /nfs/users/nfs_a/at31/miniforge3/envs/jupy/bin/Rscript
    
    # libraries
    library(magrittr)
    library(Matrix)
    
    # generate annotations file
    annotations <- 
        readr::read_tsv('${params.annotations}') %>%
        dplyr::filter(id == '${meta.id}') %>%
        dplyr::transmute(cell, annotation = ${params.annotation_col})
    annotations %>% 
      readr::write_tsv('annotations.tsv', col_names = F)
    
    # subset matrix to cells in the sample
    raw_counts_matrix <- 
        readRDS('${raw_counts_matrix}')
    raw_counts_matrix <-
      raw_counts_matrix[, colnames(raw_counts_matrix) %in% annotations\$cell]

    # create infercnv object
    infercnv_obj <-
        infercnv::CreateInfercnvObject(
            raw_counts_matrix = raw_counts_matrix,
            annotations_file = 'annotations.tsv',
            gene_order_file = '${params.gene_order_file}',
            ref_group_names = NULL)
    
    # run infercnv 
    # (parameters taken from mp34)
    # (/lustre/scratch126/casm/team268im/mp34/analysis/synovium_scRNA/infercnv_analysis.R)
    options(scipen = 100)
    infercnv_obj <-
        infercnv::run(
            infercnv_obj,
            out_dir = 'out/',
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
            HMM_report_by = c('subcluster'),
            
            # # https://github.com/harbourlab/UPhyloplot2
            # tumor_subcluster_partition_method = 'random_trees',
            
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
        output_filename = 'out/infercnv.median_filtered',
        x.center = 1,
        color_safe_pal = F)
      
    # write metadata table
    infercnv::add_to_seurat(
      infercnv_output_path = 'out/')
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
  
  // get metadata + bam paths
  Channel.fromPath(params.mappings, checkIfExists: true) 
  | splitCsv(header: true)
  | map { row ->
    meta = row.subMap('id')
    [meta, file(row.raw_counts_matrix, checkIfExists: true)]
  }
  | set { mappings }
  
  // run infercnv
  mappings |
  infercnv
  
}
