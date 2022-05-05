// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from "./functions"

options = initOptions(params.options)

process MOTHUR_SCREEN_SEQS {
    tag "$meta.id"
    label "process_low"
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options,
                                        publish_dir:getSoftwareName(task.process)) }

    container "quay.io/biocontainers/mothur:1.47.0--hb64bf22_2"
    conda (params.enable_conda ? "bioconda::mothur=1.47.0" : null)

    input:
    tuple val(meta), file(fasta), file(count)

    output:
    tuple val(meta), path("${outprefix}.fasta"), emit: fasta
    tuple val(meta), path("${outprefix}.count_table"), emit: count_table
    path "*.version.txt", emit: version

    script:
    def software = getSoftwareName(task.process)
    def procname = "${task.process.tokenize(':')[-1].toLowerCase()}"
    outprefix = "${procname}.${meta.id}"
    """
    mothur "#
    screen.seqs(fasta=$fasta, count=$count, minlength=${params.min_aln_len});
    screen.seqs(fasta=current, count=current, optimize=start-end, criteria=${params.criteria});
    summary.seqs(fasta=current)"

    mv *.good.good.${fasta.getExtension()} ${outprefix}.fasta

    # if it exists
    if [ -f *.good.good.count_table ]; then
        mv *.good.good.count_table ${outprefix}.count_table
    elif [ -f *.good.count_table ]; then
        mv *.good.count_table ${outprefix}.count_table
    else
        cp $count ${outprefix}.count_table
    fi

    # print version
    mothur -v | tail -n+2 | head -1 | cut -d"=" -f2 > ${software}.version.txt
    """
}
