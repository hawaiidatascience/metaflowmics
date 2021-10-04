// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from "./functions"

options = initOptions(params.options)

process CDHIT {
    label "process_high"
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options,
                                        publish_dir:getSoftwareName(task.process)) }

    container "quay.io/biocontainers/cd-hit:4.8.1--h2e03b76_5"
    conda (params.enable_conda ? "bioconda::cd-hit=4.8.1" : null)

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*.clstr"), emit: cluster
    tuple val(meta), path("*.repr.fa"), emit: repr
    path "*.version.txt", emit: version

    script:
    def software = getSoftwareName(task.process)
    """
    #!/usr/bin/env bash
    
    cd-hit \\
      -T $task.cpus \\
      -M ${task.memory.getMega()} \\
      -i $fasta \\
      -o cdhit.repr.fa \\
      -c $params.identity

    cd-hit | grep -i "CD-HIT version" | cut -d' ' -f4 > ${software}.version.txt
    """
}
