// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from "./functions"

options = initOptions(params.options)

process VSEARCH_SINTAX{
    tag "$otu_id"
    label "process_high"
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options,
                                        publish_dir:getSoftwareName(task.process),
                                        meta:meta) }

    container "quay.io/biocontainers/vsearch:2.17.0--h95f258a_1"
    conda (params.enable_conda ? "bioconda::vsearch=2.17.0" : null)

    input:
    tuple val(otu_id), path(fasta)
    path database

    output:
    tuple val(otu_id), path("annotations_*.tsv"), emit: taxonomy
    path "*.version.txt", emit: version

    script:
    def software = getSoftwareName(task.process)
    """
    #!/usr/bin/env bash

    vsearch \\
        --threads $task.cpus \\
        --db $database \\
        --sintax $fasta \\
        --sintax_cutoff $params.tax_confidence \\
        --tabbedout annotations_sintax-${otu_id}.tsv

    echo \$(vsearch --version 2>&1) | grep "RAM" | sed "s/vsearch v//" | sed "s/, .*//" > ${software}.version.txt
    """
}
