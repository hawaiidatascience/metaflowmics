// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

process VSEARCH_CLUSTERING {
    tag "$otu_id"
    label 'process_high'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:options, publish_dir:getSoftwareName(task.process), publish_id:"$otu_id") }

    container "nakor/vsearch:2.15.0"
    conda (params.conda ? "bioconda::vsearch=2.15.0" : null)

    input:
    path(fasta)
    each otu_id
    val options

    output:
    tuple val(otu_id), path("vsearch_OTUs-*.tsv"), emit: abundance
    tuple val(otu_id), path("vsearch_OTUs-*.fasta"), emit: fasta
    path "*.version.txt", emit: version

    script:
    def software = getSoftwareName(task.process)
    def ioptions = initOptions(options)
    def otu_id_pct = otu_id / 100
    """
    #!/usr/bin/env bash
    vsearch $ioptions.args \\
        --threads $task.cpus \\
        --sizein \\
        --sizeout \\
        --fasta_width 0 \\
        --id ${otu_id_pct} \\
        --strand plus \\
        --cluster_size ${fasta} \\
        --uc clusters${otu_id}.uc \\
        --uc clusters${otu_id}.uc \\
        --relabel OTU${otu_id}_ \\
        --centroids vsearch_OTUs-${otu_id}.fasta \\
        --otutabout vsearch_OTUs-${otu_id}.tsv

    echo \$(vsearch --version 2>&1) | grep "RAM" | sed 's/vsearch v//' | sed 's/, .*//' > ${software}.version.txt
    """
}
