nextflow.enable.dsl=2

// command line arguments
params.help         = false
params.mappings     = null
params.out_dir      = 'out/'
params.location     = 'local'
// targeted nanoseq, target genes, fasta to get ref alleles
params.recursites    = '/lustre/scratch125/casm/team268im/at31/RA_som_mut/scomatic/out/Zhang2023/targeted_mutation_calling/summary/recursites.tsv'
params.fasta        = '/lustre/scratch125/casm/team268im/al28/bed_ref/GRCh38_full_analysis_set_plus_decoy_hla_genome.fa'
// bam2r filtering arguments (recommended by Andrew)
params.min_phred    = 30
params.sam_flag     = 3844
params.mapq         = 25
params.window       = 10

// help
if (params.help) {
  help = \
  """
  |driver_converage.nf: get the coverage of a list of genes in a list of BAMs
  |
  |Required arguments:
  |   --mappings    Path to the mappings CSV file with columns `id`, `celltype`, 
  |                 and `bam`.
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

// get windows of interest from dndscv
process get_windows {
  label 'normal'
  output:
    path('refs.tsv')
  script:
    """
    #!/usr/bin/env /nfs/users/nfs_a/at31/miniforge3/envs/jupy/bin/Rscript
    
    #libraries
    library(magrittr)
    library(GenomicRanges)
    
    # get recurrent mutated sites
    recursites <-
      readr::read_tsv('${params.recursites}')
    
    # get each position in a x bp window around each recursite
    recursites_windows <-
      recursites %>%
      dplyr::group_by(chr, pos) %>%
      dplyr::reframe(pos = (pos - ${params.window}):(pos + ${params.window})) %>%
      dplyr::transmute(chr = paste0('chr', chr), start = pos - 1, end = pos)
    
    # save bed 
    recursites_windows %>%
      readr::write_tsv('recursites_windows.bed', col_names = F)
    
    # get reference alleles at each position in the window
    system(
      'bedtools getfasta -fi ${params.fasta} -bed recursites_windows.bed -tab', 
      intern = T) %>%
      paste(collapse = '\\n') %>%
      readr::read_tsv(col_names = c('pos', 'ref')) %>%
      tidyr::separate_wider_delim(
        cols = pos, delim = stringr::regex(':|-'),
        names = c('chr', 'start', 'end')) %>%
      type.convert(as.is = T) %>%
      dplyr::transmute(chr = gsub('chr', '', chr), pos = end, ref) %>%
      dplyr::distinct() %>%
      readr::write_tsv('refs.tsv')
    """
}

// get reads at each window of interest from dndscv
process get_reads {
  tag "${meta.id}_${meta.celltype}"
  memory = { 5.GB * task.attempt }
  label 'normal'
  input:
    tuple val(meta), path(bam), path(bai)
    path(refs)
  output:
    tuple val(meta), path("${meta.id}_${meta.celltype}_reads.tsv")
  script:
    """
    #!/usr/bin/env /nfs/users/nfs_a/at31/miniforge3/envs/jupy/bin/Rscript
    
    # libraries
    library(magrittr)
    
    # all possible nucleotides
    nucleotides <- c("A", "T", "C", "G", "-", "INS", "a", "t", "c", "g", "_", "ins")
    
    # read recursites
    recursites <- readr::read_tsv('${params.recursites}')
    
    # read ref alleles
    refs <- readr::read_tsv('${refs}')
    
    sites <-
      recursites %>%
      purrr::pmap(function(chr, pos, gene, ref, mut, ...) {

        # get all positions
        positions <- (pos - ${params.window}):(pos + ${params.window})
        
        # get reads at each position of the window and total n reads per position
        reads_mat <-
          deepSNV::bam2R(
            '${bam}',
            chr = chr,
            # extend pileup window x bp either side of mutated position
            start = pos - ${params.window},
            stop = pos + ${params.window},
            # quality controls
            q = ${params.min_phred},
            mask = ${params.sam_flag},
            mq = ${params.mapq}
          )[, nucleotides]
        
        # count total reads, return
        reads_mat %>%
          tibble::as_tibble() %>%
          dplyr::mutate(
            celltype = '${meta.celltype}',
            chr = chr, mut_pos = pos, pos = positions,  gene = gene,
            n_reads = rowSums(reads_mat)) %>%
          dplyr::select(
            celltype, 
            chr, pos, mut_pos, gene, n_reads, 
            dplyr::all_of(nucleotides))
        
      }) %>%
      dplyr::bind_rows()
    
    # add refs, get vafs
    vafs <-
      sites %>%
      dplyr::left_join(refs) %>%
      tidyr::pivot_longer(cols = dplyr::all_of(nucleotides),
                          values_to = 'nt_n_reads', names_to = 'nt') %>%
      dplyr::mutate(n_ref = ifelse(ref == toupper(nt), nt_n_reads, 0),
                    n_non_ref = ifelse(ref != toupper(nt), nt_n_reads, 0)) %>%
      dplyr::group_by(celltype, chr, pos, mut_pos, gene, n_reads, ref) %>%
      # get number of ref and non-ref alleles at each position
      dplyr::summarise(n_ref = sum(n_ref), n_non_ref = sum(n_non_ref)) %>%
      dplyr::mutate(vaf = n_non_ref / n_reads)
      
    muts <-
      sites %>%
      dplyr::left_join(vafs)
      
    # save
    muts %>%
      readr::write_tsv("${meta.id}_${meta.celltype}_reads.tsv")
    """
}

// concat reads
process concat_reads {
  tag "${meta.id}"
  label 'normal'
  publishDir "${params.out_dir}/${meta.id}/", mode: 'copy'
  input:
    tuple val(meta), path(reads)
  output:
    path("${meta.id}_reads.tsv")
  script:
    """
    head -1 \$(ls *_reads.tsv | head -1) > ${meta.id}_reads.tsv.tmp
    for file in *_reads.tsv ; do
      echo \$file
      tail -n +2 -q \$file >> ${meta.id}_reads.tsv.tmp
    done
    mv ${meta.id}_reads.tsv.tmp ${meta.id}_reads.tsv
    """
}

// knit coverage report
process report {
  label 'normal10gb'
  publishDir "${params.out_dir}/summary/", mode: 'copy'
  input:
    path(reads)
    path(refs)    
    path(rmd)
  output:
    path('targeted_mutation_calling.html')
    path('targeted_mutation_calling_files/*')
    path('reads.rds')
  script:
    """
    #!/usr/bin/env /nfs/users/nfs_a/at31/miniforge3/envs/jupy/bin/Rscript
    
    # knit
    rmarkdown::render(
      "${rmd}",
      params = list(
        recursites = "${params.recursites}",
        cache_dir = "${params.out_dir}/summary/targeted_mutation_calling_cache/",
        rerun = T),
      output_file = "targeted_mutation_calling.html",
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
  
  // get windows to interrogate
  get_windows()
  
  // get reads
  get_reads(sample_bams, get_windows.out)
  
  // group by id, concat coverage files, collect all 
  get_reads.out 
  | map { meta, reads -> [meta.subMap('id'), reads] }
  | groupTuple
  | concat_reads
  | collect 
  | set { collected_reads }
  
  // knit report
  report(
    collected_reads, 
    get_windows.out,
    "${baseDir}/../reports/targeted_mutation_calling.Rmd")
  
}
