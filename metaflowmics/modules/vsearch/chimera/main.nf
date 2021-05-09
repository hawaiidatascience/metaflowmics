// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from "./functions"

options = initOptions(params.options)

process VSEARCH_CHIMERA {
    tag "$meta.id"
    label "process_high"
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options,
                                        publish_dir:getSoftwareName(task.process),
                                        meta:meta, publish_by_meta:["id"]) }

    container "quay.io/biocontainers/vsearch:2.17.0--h95f258a_1"
    conda (params.enable_conda ? "bioconda::vsearch=2.17.0" : null)

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*-nochimera.fasta"), emit: fasta
    path "*.version.txt", emit: version

    script:
    def software = getSoftwareName(task.process)
    """
    #!/usr/bin/env bash

    vsearch $options.args \\
        --threads $task.cpus \\
        --uchime3_denovo $fasta \\
        --sizein \\
        --sizeout \\
        --fasta_width 0 \\
        --nonchimeras ${meta.id}-nochimera.fasta

    echo \$(vsearch --version 2>&1) | grep "RAM" | sed "s/vsearch v//" | sed "s/, .*//" > ${software}.version.txt
    """
}
