// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

options = initOptions(params.options)

process MOTHUR_REMOVE_LINEAGE {
    tag "$otu_id"
    label "process_high"
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options,
                                        publish_dir:getSoftwareName(task.process)) }

    container "quay.io/biocontainers/mothur:1.44.1--hf0cea05_2"
    conda (params.enable_conda ? "bioconda::mothur:1.44.1" : null)

    input:
    tuple val(otu_id), file(list), file(constax), file(shared), file(count_table)

    output:
    tuple val(otu_id), path("${outprefix}.cons.taxonomy"), emit: constaxonomy
    tuple val(otu_id), path("${outprefix}.list"), emit: list
    tuple val(otu_id), path("${outprefix}.shared"), emit: shared
    tuple val(otu_id), path("${outprefix}.count_table"), emit: count_table    
    path "*.version.txt", emit: version

    script:
    def software = getSoftwareName(task.process)
    def procname = "${task.process.tokenize(':')[-1].toLowerCase()}"
    outprefix = options.suffix ? "${options.suffix}.${otu_id}" : "${procname}.${otu_id}"
    """
    mothur "#remove.lineage(constaxonomy=$constax, list=$list, shared=$shared, taxon='$params.taxa_to_filter');list.seqs(list=current);get.seqs(accnos=current,count=$count)"
    # rename outputs
    mv *.pick.list ${outprefix}.list
    mv *.pick.cons.taxonomy ${outprefix}.cons.taxonomy
    mv *.pick.shared ${outprefix}.shared
    mv *.pick.count_table ${outprefix}.count_table

    # print version
    mothur -v | tail -n+2 | head -1 | cut -d'=' -f2 > ${software}.version.txt
    """
}
