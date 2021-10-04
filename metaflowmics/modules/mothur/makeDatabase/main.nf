// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from "./functions"

options = initOptions(params.options)

process MOTHUR_MAKE_DATABASE {
    tag "$meta.id"
    label "process_high"
    publishDir "${params.outdir}", mode: params.publish_dir_mode

    container "quay.io/biocontainers/mothur:1.44.1--hf0cea05_2"
    conda (params.enable_conda ? "bioconda::mothur:1.44.1" : null)

    input:
    tuple val(meta), file(shared), file(constax), file(repfasta), file(repcount)

    output:
    tuple val(meta), path("*.database"), emit: database
    path "*.version.txt", emit: version

    script:
    def software = getSoftwareName(task.process)
    def procname = "${task.process.tokenize(':')[-1].toLowerCase()}"
    def outprefix = "${procname}.${meta.id}"
    """
    mothur "#create.database(shared=$shared,repfasta=$repfasta,constaxonomy=$constax,count=$repcount)"
    mv *.database ${outprefix}.database

    # print version
    mothur -v | tail -n+2 | head -1 | cut -d"=" -f2 > ${software}.version.txt
    """
}
