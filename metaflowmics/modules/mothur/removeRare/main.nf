// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from "./functions"

options = initOptions(params.options)

process MOTHUR_REMOVE_RARE {
    tag "$otu_id"
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
    tuple val(otu_id), path("${outprefix}.list"), emit: list
    tuple val(otu_id), path("${outprefix}.shared"), emit: shared
    tuple val(otu_id), path("${outprefix}.count_table"), emit: count_table
    path "*.version.txt", emit: version

    script:
    def software = getSoftwareName(task.process)
    def procname = "${task.process.tokenize(':')[-1].toLowerCase()}"
    outprefix = options.suffix ? "$options.suffix" : "${procname}.${otu_id}"
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
