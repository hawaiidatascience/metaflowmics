// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

options = initOptions(params.options)

process MOTHUR_REMOVE_LINEAGE {
    label "process_medium"
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options,
                                        publish_dir:getSoftwareName(task.process)) }

    container "quay.io/biocontainers/mothur:1.44.1--hf0cea05_2"
    conda (params.enable_conda ? "bioconda::mothur:1.44.1" : null)

    input:
    tuple file(fasta), file(count), file(taxonomy)

    output:
    path "${outprefix}.fasta", emit: fasta
    path "${outprefix}.count_table", emit: count_table
    path "${outprefix}.taxonomy", emit: taxonomy    
    path "*.version.txt", emit: version

    script:
    def software = getSoftwareName(task.process)
    def procname = "${task.process.tokenize(':')[-1].toLowerCase()}"
    outprefix = options.suffix ? options.suffix : procname
    """
    mothur "#remove.lineage(taxonomy=$taxonomy, fasta=$fasta, count=$count, taxon='$params.taxa_to_filter')"
    # rename outputs
    mv *.pick.fasta ${outprefix}.fasta
    mv *.pick.count_table ${outprefix}.count_table
    mv *.pick.taxonomy ${outprefix}.taxonomy

    # print version
    mothur -v | tail -n+2 | head -1 | cut -d'=' -f2 > ${software}.version.txt
    """
}
