// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

process ITSXPRESS {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:options, publish_dir:getSoftwareName(task.process), publish_id:meta.id) }

    container "nakor/itsxpress:1.8.0"
    conda (params.conda ? "bioconda::itsxpress=1.8.0" : null)

    input:
    tuple val(meta), path(reads)
    val options

    output:
    tuple val(meta), path("*.fastq.gz"), emit: fastq
    path "*.log", emit: log
    path "*.version.txt", emit: version

    script:
    def software = getSoftwareName(task.process)
    def ioptions = initOptions(options)
    def rev_arg = meta.paired ? "--fastq2 ${reads[1]}" : "--single_end"
    """
    itsxpress --region ${params.locus} --threads ${task.cpus} \\
        --outfile ${meta.id}-${params.locus}.fastq.gz --log ITSxpress_${meta.id}.log \\
        --fastq ${reads[0]} $rev_arg

    pip show itsxpress | grep ^Version | sed 's/.*: //' > "${software}.version.txt"
    """
}
