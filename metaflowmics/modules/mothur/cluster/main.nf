// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from "./functions"

options = initOptions(params.options)

process MOTHUR_CLUSTER {
    tag "${otu_id}"
    label "process_high"
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options,
                                        publish_dir:getSoftwareName(task.process)) }

    container "quay.io/biocontainers/mothur:1.47.0--hb64bf22_2"
    conda (params.enable_conda ? "bioconda::mothur=1.47.0" : null)

    input:
    tuple val(meta), path(fasta), path(count), path(tax)
    each otu_id

    output:
    tuple val(meta_upd), path("*.list"), emit: list
    tuple val(meta_upd), path("*.shared"), emit: shared
    tuple val(meta_upd), path(fasta), emit:fasta
    tuple val(meta_upd), path(count), emit:count_table
    tuple val(meta_upd), path(tax), emit:taxonomy
    path "*.version.txt", emit: version

    script:
    meta_upd = meta + [id: "${otu_id}", otu_id: otu_id]
    def ext = ["rep.fasta", "cons.taxonomy", "shared", "list", "database"]
    def software = getSoftwareName(task.process)
    def procname = "${task.process.tokenize(':')[-1].toLowerCase()}"
    def outprefix = "${procname}.${meta_upd.id}"
    """
    mothur "#
    cluster(count=$count, fasta=$fasta, method=${otu_id == 100 ? "unique" : "dgc"}, cutoff=${(100-meta_upd.otu_id) / 100});
    make.shared(count=current, list=current)"

    mv *.list ${outprefix}.list
    mv *.shared ${outprefix}.shared

    # print version
    mothur -v | tail -n+2 | head -1 | cut -d"=" -f2 > ${software}.version.txt
    """
}
