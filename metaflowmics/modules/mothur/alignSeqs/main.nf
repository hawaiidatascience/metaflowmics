// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from "./functions"

options = initOptions(params.options)

process MOTHUR_ALIGN_SEQS {
    label "process_high"
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options,
                                        publish_dir:getSoftwareName(task.process)) }

    container "quay.io/biocontainers/mothur:1.44.1--hf0cea05_2"
    conda (params.enable_conda ? "bioconda::mothur:1.44.1" : null)

    input:
    tuple file(fasta), file(count)
    path db_aln

    output:
    path "${outprefix}.fasta", emit: fasta
    path "${outprefix}.count_table", emit: count_table
    path "*.version.txt", emit: version

    script:
    def software = getSoftwareName(task.process)
    def procname = "${task.process.tokenize(':')[-1].toLowerCase()}"
    outprefix = options.suffix ? "$options.suffix" : "${procname}"
    """
    mothur "#
    align.seqs(fasta=$fasta, reference=$db_aln);
    filter.seqs(fasta=current, vertical=T);
	screen.seqs(fasta=current, count=$count, minlength=${params.min_aln_len});
	screen.seqs(fasta=current, count=current, optimize=start-end, criteria=${params.criteria});
	summary.seqs(fasta=current)"

    mv *.filter.good.good.fasta ${outprefix}.fasta

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
