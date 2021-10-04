// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from "./functions"

options = initOptions(params.options)

process MOTHUR_ALIGN_SEQS {
	tag "$meta.id"
    label "process_high"
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options,
                                        publish_dir:getSoftwareName(task.process)) }

    container "quay.io/biocontainers/mothur:1.44.1--hf0cea05_2"
    conda (params.enable_conda ? "bioconda::mothur:1.44.1" : null)

    input:
    tuple val(meta), path(fasta), val(db_meta), file(db_aln)

    output:
    tuple val(meta_upd), path("*.fasta"), emit: fasta
    path "*.version.txt", emit: version

    script:
	meta_upd = meta + db_meta
    def software = getSoftwareName(task.process)
    def procname = "${task.process.tokenize(':')[-1].toLowerCase()}"
    def outprefix = "${procname}.${meta.id}"
    """
    mothur "#
	align.seqs(fasta=$fasta, reference=$db_aln);
	filter.seqs(fasta=current, vertical=T);"

	mv *.filter.fasta ${outprefix}.fasta

    # print version
    mothur -v | tail -n+2 | head -1 | cut -d"=" -f2 > ${software}.version.txt
    """
}
