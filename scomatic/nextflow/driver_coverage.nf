nextflow.enable.dsl=2

// expected command line arguments
params.mappings = null
params.drivers = null
params.out_dir = null
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
    tuple val(meta), path("${meta.id}.bam"), path("${meta.id}.bam.bai"), emit: bams
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
  maxForks 10
  label 'normal4core'
  input:
    tuple val(meta), val(bam)
  output:
    tuple val(meta), path("${meta.id}.bam"), path("${meta.id}.bam.bai"), emit: bams
  script:
    """
    # create local symbolic link 
    ln -s ${bam} ${meta.id}.bam
    if [ -f "${bam}.bai" ] ; then
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
  input:
    tuple val(meta), path(bam), path(bai)
  output:
    tuple val(meta), path('depths_*.tsv'), emit: depths
  script:
    """
    (
      # header
      echo -e 'chr\\tpos\\tdepth\\tgene\\tcelltype' ;
      
      # samtools depth to pile up driver regions
      while read -r gene coords ; do
        samtools depth -a -r \$coords -Q 30 ${bam} |
        awk -v gene=\$gene -v celltype="${meta.celltype}" '{print \$0"\\t"gene"\\t"celltype}'
      done < <(sed 1d ${params.drivers}) ;

    ) | cat > depths_${meta.celltype}.tsv
    """
}

// concat the depths files
process concat_depths {
  tag "${meta.id}"
  label 'normal4core'
  publishDir "${params.out_dir}/${meta.id}", mode: 'copy'
  input:
    tuple val(meta), path(depths)
  output:
    tuple val(meta), path('depths.tsv'), emit: depths
  script:
    """
    head -n1 depths_NA.tsv > depths.tsv
    tail -n +2 -q depths_*.tsv >> depths.tsv
    """
}

// knit coverage report
process report {
  tag "${meta.id}"
  label 'long16core10gb'
  publishDir "${params.out_dir}/${meta.id}/", mode: 'copy'
  input:
    tuple val(meta), path(depths)
    path(rmd)
  output:
    path('driver_coverage.html')
  script:
    """
    #!/usr/bin/env /nfs/users/nfs_a/at31/miniforge3/envs/jupy/bin/Rscript
    rmarkdown::render(
      "${rmd}",
      params = list(
        id = "${meta.id}",
        depths = "${depths}",
        ref_cds = "${params.ref_cds}",
        cache_dir = "${params.out_dir}/${meta.id}/driver_coverage_cache/",
        drivers = "${params.drivers}",
        rerun = T),
      output_file = "driver_coverage.html",
      output_dir = getwd()
    )
    """
}

// main scomatic workflow
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
  | view 
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
  
  // group by id, concat depths files
  coverage.out.depths 
  | map { meta, depths -> [meta.subMap('id'), depths] }
  | groupTuple
  | view 
  | concat_depths
  
  // knit report
  report(concat_depths.out.depths, "${baseDir}/../reports/driver_coverage.Rmd")
  
}
