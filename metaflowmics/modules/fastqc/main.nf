// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from "./functions"

options = initOptions(params.options)

process FASTQC {
    label "process_low"
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options,
                                        publish_dir:getSoftwareName(task.process)) }

    container "quay.io/biocontainers/fastqc:0.11.9--hdfd78af_1"
    conda (params.enable_conda ? "bioconda::fastqc=0.11.9" : null)

    input:
    tuple val(meta), path(fastqs)

    output:
    tuple val(meta), path("*.html"), emit: html
    tuple val(meta), path("*.zip"), emit: zip
    path "*.version.txt", emit: version

    script:
    def software = getSoftwareName(task.process)
    """
    #!/usr/bin/env bash

    fastqc -o . --threads $task.cpus $fastqs
	fastqc --version | sed -e "s/FastQC v//g" > ${software}.version.txt
    """
}
