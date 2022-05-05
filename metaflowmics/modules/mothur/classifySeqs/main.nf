// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from "./functions"

options = initOptions(params.options)

process MOTHUR_CLASSIFY_SEQS {
    tag "$meta.id"
    label "process_high"
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options,
                                        publish_dir:getSoftwareName(task.process)) }        

    container "quay.io/biocontainers/mothur:1.47.0--hb64bf22_2"
    conda (params.enable_conda ? "bioconda::mothur=1.47.0" : null)

    input:
    tuple val(meta), path(fasta), path(count)
    path(db_aln)
    path(db_tax)

    output:
    tuple val(meta), path("*.taxonomy"), emit: taxonomy
    path "*.version.txt", emit: version

    script:
    def software = getSoftwareName(task.process)
    def procname = "${task.process.tokenize(':')[-1].toLowerCase()}"
    def outprefix = "${procname}.${meta.id}"
    """
    mothur "#classify.seqs(fasta=$fasta, count=$count, template=$db_aln, taxonomy=$db_tax, cutoff=$params.classification_consensus)"
    # rename output
    mv *.wang.taxonomy ${outprefix}.taxonomy

    # print version
    mothur -v | tail -n+2 | head -1 | cut -d"=" -f2 > ${software}.version.txt
    """
}
