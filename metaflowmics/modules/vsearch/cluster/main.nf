// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

options = initOptions(params.options)

process VSEARCH_CLUSTER {
    tag "$otu_id"
    label 'process_high'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options,
                                        publish_dir:getSoftwareName(task.process),
                                        meta:meta) }

    container "quay.io/biocontainers/vsearch:2.17.0--h95f258a_1"
    conda (params.enable_conda ? "bioconda::vsearch=2.17.0" : null)

    input:
    path fasta
    each otu_id

    output:
    tuple val(otu_id), path("vsearch_OTUs-*.tsv"), emit: tsv
    tuple val(otu_id), path("vsearch_OTUs-*.fasta"), emit: fasta
    path "summary.csv", emit: summary
    path "*.version.txt", emit: version

    script:
    def software = getSoftwareName(task.process)
    def otu_id_pct = otu_id / 100
    """
    #!/usr/bin/env bash
    vsearch $options.args \\
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

    cut -f2,3 vsearch_OTUs-${otu_id}.tsv | \\
        awk 'NR==1 { for (i=1; i<=NF; i++) { hd[i]=\$i } } \\
        { for (i=1;i<=NF;i++) if(\$i>0) {c1[hd[i]]+=\$i ; c2[hd[i]]++} }; END \\
        {for (i in c1) print "clustering-${otu_id},"i","c1[i]","c2[i];}' \\
        > summary.csv

    echo \$(vsearch --version 2>&1) | grep "RAM" | sed 's/vsearch v//' | sed 's/, .*//' > ${software}.version.txt
    """
}
