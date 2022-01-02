// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from "./functions"

options = initOptions(params.options)

process MOTHUR_CLASSIFY_SEQS {
	tag "$meta.id"
    label "process_high"
    publishDir "${params.outdir}", mode: params.publish_dir_mode

    container "quay.io/biocontainers/mothur:1.44.1--hf0cea05_2"
    conda (params.enable_conda ? "bioconda::mothur:1.44.1" : null)

    input:
    tuple val(meta), path(fasta), path(count), file(db)
    // tuple val(db_meta), path(db_aln), path(db_tax)

    output:
    tuple val(meta), path("*.taxonomy"), emit: taxonomy
    path "*.version.txt", emit: version

    script:
	def db = db[0].getExtension() == "tax" ? db.reverse() : db
    def software = getSoftwareName(task.process)
    def procname = "${task.process.tokenize(':')[-1].toLowerCase()}"
    def outprefix = "${procname}.${meta.id}"
    """
    mothur "#classify.seqs(fasta=$fasta, count=$count, template=${db[0]}, taxonomy=${db[1]}, cutoff=$params.classification_consensus)"
    # rename output
    mv *.wang.taxonomy ${outprefix}.taxonomy

    # print version
    mothur -v | tail -n+2 | head -1 | cut -d"=" -f2 > ${software}.version.txt
    """
}
