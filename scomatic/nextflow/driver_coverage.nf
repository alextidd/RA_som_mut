nextflow.enable.dsl=2

// command line arguments
params.help       = false
params.mappings   = null
params.drivers    = null
params.out_dir    = 'out/'
params.location   = 'irods'
params.min_MQ     = 30

// help
if (params.help) {
  help = \
  """
  |driver_converage.nf: get the coverage of a list of genes in a list of BAMs
  |
  |Required arguments:
  |   --mappings    Path to the mappings CSV file with columns `id`, `celltype`, 
  |                 and `bam`.
  |   --drivers     Path to the drivers TSV file with columns `gene` and `coords`.
  |                 Coords must be in the format 'chr:start-stop' 
  |                 (e.g. 1:12345-12456).
  |
  |Optional arguments:
  |   --out_dir     Path to the output directory. default is `out/`.
  |                 [default: ${params.out_dir}]
  |   --location    Either "irods" or "local".
  |                 [default: ${params.location}]  
  """.stripMargin()
  
  // print help and exit
  println(help)
  exit(0)
}

// Download a given sample's BAM from iRODS
// Then either retrieve the BAI or make one via indexing
// The maxForks of 10 was set after asking jc18 about best iRODS practices
process irods {
  tag "${meta.id}_${meta.celltype}"
  maxForks 10
  label 'normal4core'
  input:
    tuple val(meta), val(bam)
  output:
    tuple val(meta), path("${meta.celltype}.bam"), path("${meta.celltype}.bam.bai"), emit: bams
  script:
    """
    iget -K ${bam} ${meta.celltype}.bam
    if [[ `ils ${bam}.bai | wc -l` == 1 ]]
    then
        iget -K ${bam}.bai ${meta.celltype}.bam.bai
    else
        samtools index -@ ${task.cpus} ${meta.celltype}.bam
    fi
    """
}

// The equivalent of an irods download, but for a local copy of mappings
// Symlink the BAM/BAI appropriately so they're named the right thing for downstream
process local {
  tag "${meta.id}_${meta.celltype}"
  maxForks 10
  label 'normal4core'
  errorStrategy = {task.attempt <= 1 ? 'retry' : 'ignore'}
  input:
    tuple val(meta), val(bam)
  output:
    tuple val(meta), path("${meta.celltype}.bam"), path("${meta.celltype}.bam.bai")
  script:
    """
    # create local symbolic link 
    ln -s ${bam} ${meta.celltype}.bam
    if [ -f "${bam}.bai" ] ; then
        ln -s ${bam}.bai ${meta.celltype}.bam.bai
    else
        samtools index -@ ${task.cpus} ${meta.celltype}.bam
    fi
    """
}

// calculate nucleotidic coverage of each driver position in the BAMs
// (using samtools depth)
process coverage {
  tag "${meta.id}_${meta.celltype}"
  memory = { 5.GB * task.attempt }
  label 'normal'
  input:
    tuple val(meta), path(bam), path(bai)
  output:
    tuple val(meta), path("${meta.celltype}_coverage.tsv"), emit: coverages
  script:
    """
    (
      # header
      echo -e 'chr\\tpos\\tgene\\t${meta.celltype}' ;
      
      # samtools depth to pile up driver regions
      while read -r gene coords ; do
        samtools depth -a -Q ${params.min_MQ} -r \$coords ${bam} |
        awk -v FS='\\t' -v OFS='\\t' -v gene=\$gene '{print \$1,\$2,gene,\$3}'
      done < <(sed 1d ${params.drivers}) ;

    ) | cat > ${meta.celltype}_coverage.tsv
    """
}

// concat the coverage files
process cbind_coverages {
  tag "${meta.id}"
  label 'normal4core'
  publishDir "${params.out_dir}/${meta.id}", mode: 'copy'
  input:
    tuple val(meta), path(coverages)
  output:
    path("${meta.id}_coverage.tsv")
  script:
    """
    #!/usr/bin/env /nfs/users/nfs_a/at31/miniforge3/envs/jupy/bin/Rscript
    
    library(magrittr)
    
    # load in files, join
    cov <-
      list.files(pattern = '_coverage.tsv') %>%
      purrr::set_names(., gsub('_coverage.tsv', '', .)) %>%
      purrr::map(function(file) {
        readr::read_tsv(file, show_col_types = F)
      }) %>%
      purrr::reduce(dplyr::full_join, by = c('chr', 'pos', 'gene'))
      
    # write joined file
    cov %>%
      readr::write_tsv('${meta.id}_coverage.tsv')
    """
}

// knit coverage report
process report {
  tag "${meta.id}"
  label 'long16core10gb'
  publishDir "${params.out_dir}/summary/", mode: 'copy'
  input:
    path(coverages)
    path(rmd)
  output:
    path('driver_coverage.html')
  script:
    """
    #!/usr/bin/env /nfs/users/nfs_a/at31/miniforge3/envs/jupy/bin/Rscript
    rmarkdown::render(
      "${rmd}",
      params = list(
        ref_cds = "${params.ref_cds}",
        cache_dir = "${params.out_dir}/summary/driver_coverage_cache/",
        drivers = "${params.drivers}",
        rerun = T),
      output_file = "driver_coverage.html",
      output_dir = getwd()
    )
    """
}

// main workflow
workflow {

  // check input files have been passed
  if (params.mappings == null) {
    error "Please provide a mappings CSV file via --mappings"
  }
  if (params.drivers == null) {
    error "Please provide a drivers TSV file via --drivers"
  }
  
  // get metadata + bam paths  
  Channel.fromPath(params.mappings, checkIfExists: true)
  | splitCsv(header: true)
  | map { row ->
      meta = row.subMap('id', 'celltype') 
      [meta, file(row.bam, checkIfExists: true)]
  }
  | set { mappings }
  
  // download or locally link bams
  if (params.location == "irods") {
    // download mappings from irods
    sample_bams = irods(mappings)
  }
  else if (params.location == "local") {
    // mappings are locally available
    sample_bams = local(mappings)
  }
  
  // get coverage
  sample_bams 
  | coverage
  
  // group by id, concat coverage files, collect all 
  coverage.out.coverages 
  | map { meta, coverages -> [meta.subMap('id'), coverages] }
  | groupTuple
  | cbind_coverages
  | collect 
  | view
  | set { collected_coverages }
  
  // knit report
  report(collected_coverages, "${baseDir}/../reports/driver_coverage.Rmd")
  
}
