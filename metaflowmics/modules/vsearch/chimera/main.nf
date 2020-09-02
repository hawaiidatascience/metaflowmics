// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

process VSEARCH_CHIMERA {
    tag "$meta.id"
    label 'process_high'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:options, publish_dir:getSoftwareName(task.process), publish_id:meta.id) }

    container "nakor/vsearch:2.15.0"
    conda (params.conda ? "bioconda::vsearch=2.15.0" : null)

    input:
    tuple val(meta), path(fasta)
    val options

    output:
    tuple val(meta), path("*-nochimera.fasta"), emit: fasta
    path "*.version.txt", emit: version

    script:
    def software = getSoftwareName(task.process)
    def ioptions = initOptions(options)
    """
    #!/usr/bin/env bash
    
    vsearch $ioptions.args \\
        --threads $task.cpus \\
        --uchime3_denovo $fasta \\
        --sizein \\
        --sizeout \\
        --fasta_width 0 \\
        --nonchimeras ${meta.id}-nochimera.fasta

    echo \$(vsearch --version 2>&1) | grep "RAM" | sed 's/vsearch v//' | sed 's/, .*//' > ${software}.version.txt
    """
}
