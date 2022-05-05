// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from "./functions"

options = initOptions(params.options)

process MOTHUR_CHIMERA {
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

    output:
    tuple val(meta), path("${outprefix}.fasta"), emit: fasta
    tuple val(meta), path("${outprefix}.count_table"), emit: count_table
    path "*.version.txt", emit: version

    script:
    def software = getSoftwareName(task.process)
    def procname = "${task.process.tokenize(':')[-1].toLowerCase()}"
    outprefix = "${procname}.${meta.id}"
    """
    mothur "#chimera.${params.chimera_tool}(fasta=$fasta, count=$count, dereplicate=t)"

    # rename outputs
	if compgen -G "*.${params.chimera_tool}*.fasta" > /dev/null; then
        mv *.${params.chimera_tool}*.fasta ${outprefix}.fasta
	else
		cp $fasta ${outprefix}.fasta
	fi

    mv *.${params.chimera_tool}*.count_table ${outprefix}.count_table

    # print version
    mothur -v | tail -n+2 | head -1 | cut -d"=" -f2 > ${software}.version.txt
    """
}
