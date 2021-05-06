// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

options = initOptions(params.options)

process ITSXPRESS {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options,
                                        publish_dir:getSoftwareName(task.process),
                                        meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? "bioconda::bbmap=38.69 bioconda::itsxpress=1.8.0" : null)
    container "quay.io/biocontainers/itsxpress:1.8.0--py_1"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.fastq.gz"), emit: fastq
    path "*.log", emit: log
    path "*.version.txt", emit: version

    script:
    def software = getSoftwareName(task.process)
    def prefix   = options.suffix ? "${meta.id}.${options.suffix}" : "${meta.id}"
    def rev_arg = params.paired_end ? "--fastq2 ${reads[1]}" : "--single_end"
    """
    itsxpress --taxa All --region $params.locus --threads $task.cpus \\
        --outfile ${meta.id}-${params.locus}.fastq.gz --log ITSxpress_${meta.id}.log \\
        --fastq ${reads[0]} $rev_arg

    pip show itsxpress | grep ^Version | sed 's/.*: //' > "${software}.version.txt"
    """
}
