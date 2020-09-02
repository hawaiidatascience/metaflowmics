// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

process VSEARCH_TAXONOMY{
    tag "$meta.id"
    label 'process_high'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:options, publish_dir:getSoftwareName(task.process), publish_id:meta.id) }

    container "nakor/vsearch:2.15.0"
    conda (params.conda ? "bioconda::vsearch=2.15.0" : null)

    input:
    tuple val(meta), path(fasta), path(abundance)
    path database
    val options

    output:
    tuple val(meta), path("annotations_*.tsv"), emit: taxonomy
    path "*.version.txt", emit: version

    script:
    def software = getSoftwareName(task.process)
    def ioptions = initOptions(options)
    """
    #!/usr/bin/env bash
    
    vsearch \\
        --threads $task.cpus \\
        --db $database \\
        --sintax $fasta \\
        --sintax_cutoff $params.confidence_threshold \\
        --tabbedout annotations_sintax-${meta.otu_id}.tsv

    echo \$(vsearch --version 2>&1) | grep "RAM" | sed 's/vsearch v//' | sed 's/, .*//' > ${software}.version.txt
    """
}
