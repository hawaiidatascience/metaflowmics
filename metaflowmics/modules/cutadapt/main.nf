// Pulled from https://github.com/nf-core/modules/blob/master/modules/cutadapt/main.nf

include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process CUTADAPT {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? 'bioconda::cutadapt=3.2' : null)
    container 'quay.io/biocontainers/cutadapt:3.2--py38h0213d0e_0'

    input:
    tuple val(meta), path(reads)
    path barcodes

    output:
    tuple val(meta), path('*.fastq.gz'), emit: reads
    path '*.log', emit: log
    path '*.version.txt'                    , emit: version

    script:
    def software = getSoftwareName(task.process)
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    def trimmed  = params.single_end ? "-o ${prefix}_{name}.fastq.gz" : "-o ${prefix}_{name}_R1.fastq.gz -p ${prefix}_{name}_R2.fastq.gz"
    """
    cutadapt --discard-untrimmed --rc \\
        --cores $task.cpus \\
        -e $params.max_error_rate \\
        -a file:$barcodes \\
        $trimmed \\
        $reads \\
        > ${prefix}_cutadapt.log

    echo \$(cutadapt --version) > ${software}.version.txt
    """
}
