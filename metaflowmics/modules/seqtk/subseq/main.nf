// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from "./functions"

options = initOptions(params.options)

process SEQTK_SUBSEQ {
    tag "$meta"
    label "process_low"
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options,
                                        publish_dir:getSoftwareName(task.process)) }

    container "quay.io/biocontainers/seqtk:1.3--h5bf99c6_3"
    conda (params.enable_conda ? "bioconda::seqtk=1.3" : null)

    input:
    tuple val(meta), path(fasta), path(ids)

    output:
    tuple val(meta), path("*.fasta"), emit: fasta
    path "*.version.txt", emit: version

    script:
    def software = getSoftwareName(task.process)
    outprefix = fasta.getBaseName()
    """
    #!/usr/bin/env bash

    for f in $ids; do
        prefix=${outprefix}_\$(basename \$f .txt)
        seqtk subseq $fasta \$f > \${prefix}.fasta
    done

    seqtk -h 2>&1 | grep -i version | cut -d':' -f2 > ${software}.version.txt
    """
}
