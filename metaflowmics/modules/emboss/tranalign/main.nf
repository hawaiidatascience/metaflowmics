// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from "./functions"

options = initOptions(params.options)

process EMBOSS_TRANALIGN {
    label "process_medium"
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options,
                                        publish_dir:getSoftwareName(task.process)) }

    container "nakor/metaflowmics-biotools:0.0.1"
    conda (params.enable_conda ? "bioconda::emboss=6.6.0" : null)

    input:
    tuple val(meta), path(afa), path(fna)

    output:
    tuple val(meta), path("*.codons.afa"), emit: fna
    path "*.version.txt", emit: version

    script:
    def software = getSoftwareName(task.process)
    def prefix = fna.getSimpleName()
    """
    #!/usr/bin/env bash

    tranalign \\
      -asequence $fna \\
      -bsequence $afa \\
      -outseq stdout \\
      -table 5 |
        seqtk seq -l0 > rev.afa

    [ -s rev.afa ] && echo "tranalign was successfull" || (echo "Something went wrong. Are sequences in the same order?" && exit 1)

    aln_len=\$(grep -v '>' rev.afa | wc -L)

    awk -v L=\$aln_len '{if (NR % 2 == 0 && length(\$1)<1707) \$1=\$1"-"; print}' rev.afa \\
        > ${prefix}.codons.afa

    tranalign -h 2>&1 | grep -i version | cut -d':' -f3 > ${software}.version.txt
    """
}
