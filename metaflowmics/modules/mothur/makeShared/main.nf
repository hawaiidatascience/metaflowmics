// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from "./functions"

options = initOptions(params.options)

process MOTHUR_MAKE_SHARED {
    label "process_medium"
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options,
                                        publish_dir:getSoftwareName(task.process)) }

    container "quay.io/biocontainers/mothur:1.44.1--hf0cea05_2"
    conda (params.enable_conda ? "bioconda::mothur:1.44.1" : null)

    input:
    tuple val(otu_id), file(list), file(count)

    output:
    tuple val(otu_id), path("${outprefix}.shared"), emit: shared
    path "*.version.txt", emit: version

    script:
    def software = getSoftwareName(task.process)
    def procname = "${task.process.tokenize(':')[-1].toLowerCase()}"
    outprefix = options.suffix ? "$options.suffix" : "${procname}"
    """
    mothur "#make.shared(count=$count, list=$list)"

    # rename output
    mv *.shared ${outprefix}.shared

    # print version
    mothur -v | tail -n+2 | head -1 | cut -d"=" -f2 > ${software}.version.txt
    """
}
