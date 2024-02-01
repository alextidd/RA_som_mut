nextflow.enable.dsl=2

// expected command line arguments
params.mappings = null
params.drivers = null
// we expect the mappings to be on irods by default
params.location = "irods"

// Download a given sample's BAM from iRODS
// Then either retrieve the BAI or make one via indexing
// The maxForks of 10 was set after asking jc18 about best iRODS practices
process irods {
  tag "${meta.id}"
  maxForks 10
  label 'normal4core'
  input:
    tuple val(meta), val(bam)
  output:
    tuple val(meta), path("${meta.id}.bam"), path("${meta.id}.bam.bai")
  script:
    """
    iget -K ${bam} ${meta.id}.bam
    if [[ `ils ${bam}.bai | wc -l` == 1 ]]
    then
        iget -K ${bam}.bai ${meta.id}.bam.bai
    else
        samtools index -@ ${task.cpus} ${meta.id}.bam
    fi
    """
}

// The equivalent of an irods download, but for a local copy of mappings
// Symlink the BAM/BAI appropriately so they're named the right thing for downstream
process local {
  tag "${meta.id}"
  label 'normal4core'
  errorStrategy = {task.attempt <= 1 ? 'retry' : 'ignore'}
  input:
    tuple val(meta), val(bam)
  output:
    tuple val(meta), path("${meta.id}.bam"), path("${meta.id}.bam.bai")
  script:
    """
    ln -s ${bam} ${meta.id}.bam
    if [ -f "${bam}.bai" ]
    then
        ln -s ${bam}.bai ${meta.id}.bam.bai
    else
        samtools index -@ ${task.cpus} ${meta.id}.bam
    fi
    """
}

// calculate nucleotidic coverage of each driver position in the BAMs
process coverage {
  tag "${meta.id}"
  label 'long16core10gb'
  publishDir "${meta.out_dir}", mode: 'copy'
  input:
    tuple val(meta), path(bam), path(bai)
  output:
    tuple val(meta), path('depths.tsv')
  script:
    """
    (
    
      # header
      echo -e 'chr\tpos\tdepth\tgene' ;
      
      # samtools depth to pile up driver regions
      while read -r gene coords ; do
        samtools depth -a -r \$coords -Q 30 ${bam} |
        awk -v gene=\$gene '{print \$0"\t"gene}'
      done < <(sed 1d ${params.drivers}) ;
      
    ) | cat > depths.tsv
    """
}

// knit coverage report
process report {
  tag "${meta.id}"
  label 'long16core10gb'
  publishDir "${meta.out_dir}", mode: 'copy'
  input:
    tuple val(meta), path(depths)
    path(rmd)
  output:
    path('driver_coverage.html')
  script:
    """
    #!/usr/bin/env Rscript
    rmarkdown::render(
      "${rmd}",
      params = list(
        id = "${meta.id}",
        depths = "${depths}",
        cache_dir = "${meta.out_dir}/driver_coverage_cache/",
        rerun = F),
      output_file = "driver_coverage.html",
      output_dir = getwd()
    )
    """
}

// main scomatic workflow
workflow {
  // make sure input files passed
  if (params.mappings == null) {
    error "Please provide a mappings CSV file via --mappings"
  }
  if (params.drivers == null) {
    error "Please provide a drivers TSV file via --drivers"
  }
  
  // get metadata + bam paths  
  mappings = Channel.fromPath(params.mappings, checkIfExists: true)
    | splitCsv(header: true)
    | map { row ->
      meta = [id: row.id, out_dir: row.out_dir]
      [meta,  file(row.bam, checkIfExists: true)]
    }
  
  // download or locally link bams
  if (params.location == "irods") {
    // download mappings from irods
    sample_bams = irods(mappings)
  }
  else if (params.location == "local") {
    // mappings are locally available
    sample_bams = local(mappings)
  }
  
  // get coverage and knit report
  sample_bams \
  | coverage
  | report
}
