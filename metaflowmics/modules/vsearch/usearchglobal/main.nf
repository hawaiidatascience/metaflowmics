// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

process VSEARCH_USEARCHGLOBAL {
    tag "$otu_id"
    label 'process_high'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:options, publish_dir:getSoftwareName(task.process), publish_id:"$otu_id") }

    container "nakor/vsearch:2.15.0"
    conda (params.conda ? "bioconda::vsearch=2.15.0" : null)

    input:
    tuple val(otu_id), path(fasta)
    val options

    output:
    tuple val(otu_id), path("matchlist-*.tsv"), emit: tsv
    path "*.version.txt", emit: version

    script:
    def software = getSoftwareName(task.process)
    def ioptions = initOptions(options)
    """
    #!/usr/bin/env bash
    vsearch $ioptions.args \\
        --usearch_global $fasta \\
        --threads $task.cpus \\
        --id $params.lulu_min_match \\
        --db $fasta --self \\
        --iddef 1 \\
        --userout matchlist-${otu_id}.tsv \\
        -userfields query+target+id \\
        --maxaccepts 0 \\
        --query_cov .9 \\
        --maxhits 10

    echo \$(vsearch --version 2>&1) | grep "RAM" | sed 's/vsearch v//' | sed 's/, .*//' > ${software}.version.txt
    """
}
