nextflow.enable.dsl-2

// Expected command line argyments
params.mappings = null
params.celltypes = null 

// perform infercnv on each sample
process infercnv {
  tag "${sample}"
  memory { 32.GB * task.attempt }
  cpus { 12 * task.attempt }
  time { 24.hour * task.attempt }
  queue { 'week' }

  publishDir "${params.output.dir}/${id}/",
    mode: 'copy',
    pattern: "infercnv.*"

  input:
    tuple val(id), val(rds_file)
    
  output:
    tuple val(id), path('infercnv.*')

  script:
    """
    #!/usr/bin/env Rscript
    seu <- readRDS("${rds_file}")
    
    """
}