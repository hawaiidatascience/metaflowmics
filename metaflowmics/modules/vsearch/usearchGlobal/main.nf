// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from "./functions"

options = initOptions(params.options)

process VSEARCH_USEARCH_GLOBAL {
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

    output:
    tuple val(otu_id), path("matchlist-*.tsv"), emit: tsv
    path "*.version.txt", emit: version

    script:
    def software = getSoftwareName(task.process)
    """
    #!/usr/bin/env bash

    vsearch $options.args \\
        --usearch_global $fasta \\
        --threads $task.cpus \\
        --id ${params.lulu_min_match/100} \\
        --db $fasta --self \\
        --iddef 1 \\
        --userout matchlist-${otu_id}.tsv \\
        -userfields query+target+id \\
        --maxaccepts 0 \\
        --query_cov .9 \\
        --maxhits 10

    echo \$(vsearch --version 2>&1) | grep "RAM" | sed "s/vsearch v//" | sed "s/, .*//" > ${software}.version.txt
    """
}
