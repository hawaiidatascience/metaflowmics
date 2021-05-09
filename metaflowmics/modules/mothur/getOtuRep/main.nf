// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from "./functions"

options = initOptions(params.options)

process MOTHUR_GET_OTU_REP {
    tag "$otu_id"
    label "process_high"
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options,
                                        publish_dir:getSoftwareName(task.process)) }

    container "quay.io/biocontainers/mothur:1.44.1--hf0cea05_2"
    conda (params.enable_conda ? "bioconda::mothur:1.44.1" : null)

    input:
    tuple val(otu_id), file(list), file(fasta), file(count)

    output:
    tuple val(otu_id), path("${outprefix}.rep.fasta"), emit: fasta
    tuple val(otu_id), path("${outprefix}.rep.count_table"), emit: count_table
    path "*.version.txt", emit: version

    script:
    def software = getSoftwareName(task.process)
    def procname = "${task.process.tokenize(':')[-1].toLowerCase()}"
    outprefix = options.suffix ? "$options.suffix" : "${procname}.${otu_id}"
    """
    mothur "#get.oturep(fasta=$fasta, list=$list, count=$count, method=abundance, rename=t)"

    mv *.rep.fasta ${outprefix}.rep.fasta
    mv *.rep.count_table ${outprefix}.rep.count_table

    # print version
    mothur -v | tail -n+2 | head -1 | cut -d"=" -f2 > ${software}.version.txt
    """
}
