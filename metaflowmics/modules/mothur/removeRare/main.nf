// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from "./functions"

options = initOptions(params.options)

process MOTHUR_REMOVE_RARE {
    tag "$meta.id"
    label "process_medium"
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options,
                                        publish_dir:getSoftwareName(task.process)) }

    container "quay.io/biocontainers/mothur:1.47.0--hb64bf22_2"
    conda (params.enable_conda ? "bioconda::mothur=1.47.0" : null)

    input:
    tuple val(meta), file(list), file(count)

    output:
    tuple val(meta), path("${outprefix}.list"), emit: list
    tuple val(meta), path("${outprefix}.shared"), emit: shared
    tuple val(meta), path("${outprefix}.count_table"), emit: count_table
    path "*.version.txt", emit: version

    script:
    def software = getSoftwareName(task.process)
    def procname = "${task.process.tokenize(':')[-1].toLowerCase()}"
    outprefix = "${procname}.${meta.id}"
    """
    mothur "#
    remove.rare(count=$count, list=$list, nseqs=${params.min_abundance});
    make.shared(count=current, list=current);
    count.seqs(count=current, compress=f)"

    # rename outputs
    mv *.pick.shared ${outprefix}.shared
    mv *.pick.list ${outprefix}.list
    mv *.pick.full.count_table ${outprefix}.count_table

    # print version
    mothur -v | tail -n+2 | head -1 | cut -d"=" -f2 > ${software}.version.txt
    """
}
