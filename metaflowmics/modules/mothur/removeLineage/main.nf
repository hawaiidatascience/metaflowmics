// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from "./functions"

options = initOptions(params.options)

process MOTHUR_REMOVE_LINEAGE {
    tag "$otu_id"
    label "process_medium"
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options,
                                        publish_dir:getSoftwareName(task.process)) }

    container "quay.io/biocontainers/mothur:1.44.1--hf0cea05_2"
    conda (params.enable_conda ? "bioconda::mothur:1.44.1" : null)

    input:
    tuple val(otu_id), file(count), file(taxonomy), file(list)

    output:
    tuple val(otu_id), path("${outprefix}.count_table"), emit: count_table
    tuple val(otu_id), path("${outprefix}.taxonomy"), emit: taxonomy
    tuple val(otu_id), path("${outprefix}.list"), emit: list
    tuple val(otu_id), path("${outprefix}.shared"), emit: shared
    path "*.version.txt", emit: version

    script:
    def software = getSoftwareName(task.process)
    def procname = "${task.process.tokenize(':')[-1].toLowerCase()}"
    outprefix = options.suffix ? options.suffix : procname
    """
    mothur "#remove.lineage(taxonomy=$taxonomy, count=$count, list=$list, taxon="$params.taxa_to_filter"); make.shared(list=current, count=current)"

    # rename outputs
    mv *.pick.taxonomy ${outprefix}.taxonomy
    mv *.shared ${outprefix}.shared

    # if it exists
    [ -f *.pick.count_table ] && mv *.pick.count_table ${outprefix}.count_table || cp $count ${outprefix}.count_table
    [ -f *.pick.list ] && mv *.pick.list ${outprefix}.list || cp $list ${outprefix}.list

    # print version
    mothur -v | tail -n+2 | head -1 | cut -d"=" -f2 > ${software}.version.txt
    """
}
