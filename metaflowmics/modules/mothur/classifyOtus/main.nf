// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from "./functions"

options = initOptions(params.options)

process MOTHUR_CLASSIFY_OTUS {
    tag "$meta.id"
    label "process_medium"
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options,
                                        publish_dir:getSoftwareName(task.process)) }

    container "quay.io/biocontainers/mothur:1.47.0--hb64bf22_2"
    conda (params.enable_conda ? "bioconda::mothur=1.47.0" : null)

    input:
    tuple val(meta), file(list), file(count), file(tax)

    output:
    tuple val(meta), path("*.cons.taxonomy"), emit: taxonomy
    path "*.version.txt", emit: version

    script:
    def software = getSoftwareName(task.process)
    def procname = "${task.process.tokenize(':')[-1].toLowerCase()}"
    def outprefix = "${procname}.${meta.id}"
    """
    mothur "#
    list.seqs(list=$list); get.seqs(taxonomy=$tax, accnos=current);
    classify.otu(taxonomy=current, list=current, count=$count, probs=f)"

    # rename output
    mv *.cons.taxonomy ${outprefix}.cons.taxonomy

    # print version
    mothur -v | tail -n+2 | head -1 | cut -d"=" -f2 > ${software}.version.txt
    """
}
