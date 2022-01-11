// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from "./functions"

options = initOptions(params.options)

process PEAR {
    label "process_medium"
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options,
                                        publish_dir:getSoftwareName(task.process)) }

    container "quay.io/biocontainers/pear:0.9.6--h36cd882_7"
    conda (params.enable_conda ? "bioconda::pear=0.9.6" : null)

    input:
    tuple val(meta), path(fastq)

    output:
    tuple val(meta_upd), path("*.assembled.fastq"), emit: assembled
    tuple val(meta_upd), path("*.unassembled.*.fastq"), optional: true, emit: unassembled
    tuple val(meta_upd), path("*.discarded.fastq"), optional: true, emit: discarded
    path "*.version.txt", emit: version

    script:
	meta_upd = meta + [paired_end: false]

    def software = getSoftwareName(task.process)
    """
    #!/usr/bin/env bash

    pear -f ${fastq[0]} -r ${fastq[1]} -o ${meta.id}.fastq \\
      -j $task.cpus \\
      -n $params.min_contig_length -m $params.max_contig_length \\
      -v $params.min_overlap \\
      -q $params.quality_threshold

    pear | grep -i "PEAR v[0-9]" | cut -d' ' -f2 > ${software}.version.txt
    """
}
