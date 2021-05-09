// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from "./functions"

options = initOptions(params.options)

process MOTHUR_CHIMERA {
    label "process_high"
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options,
                                        publish_dir:getSoftwareName(task.process)) }

    container "quay.io/biocontainers/mothur:1.44.1--hf0cea05_2"
    conda (params.enable_conda ? "bioconda::mothur:1.44.1" : null)

    input:
    tuple file(fasta), file(count)

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
    chimera.${params.chimera_tool}(fasta=$fasta, count=$count, dereplicate=t);
    remove.seqs(fasta=current, accnos=current, dups=f)"

    # rename outputs
    mv *.pick.fasta ${outprefix}.fasta
    mv *.denovo.${params.chimera_tool}.pick.count_table ${outprefix}.count_table

    # print version
    mothur -v | tail -n+2 | head -1 | cut -d"=" -f2 > ${software}.version.txt
    """
}
