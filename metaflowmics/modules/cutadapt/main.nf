// Pulled from https://github.com/nf-core/modules/blob/master/modules/cutadapt/main.nf

include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process CUTADAPT {
    tag "$meta.id"
    label "process_medium"
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process)) }

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
    def trimmed  = params.single_end ?
        "-o {name}.fastq.gz" :
        "-o {name}_R1.fastq.gz -p {name}_R2.fastq.gz"

    barcodes = barcodes instanceof List ? barcodes : [barcodes, barcodes]

    demux_opt = ""
    if (params.paired_end && params.linked_bc) {
        demux_opt += " -a file:${barcodes[0]}"
    } else {
        if (params.forward_bc.matches('5|3')) {
            demux_opt += params.forward_bc == '5' ?
                      " -g file:${barcodes[0]}" :
                      " -a file:${barcodes[0]}"
        }
        if (params.reverse_bc.matches('5|3')) {
            demux_opt += params.reverse_bc == '5' ?
                " -G file:${barcodes[1]}" :
                " -A file:${barcodes[1]}"
        }
    }

    demux_opt += params.single_end ? " --rc" : " --pair-adapters"
    """
    cutadapt --discard-untrimmed $demux_opt \\
        --cores $task.cpus \\
        -e $params.max_error_rate \\
        $trimmed \\
        $reads \\
        > cutadapt.log

    echo \$(cutadapt --version) > ${software}.version.txt
    """
}
