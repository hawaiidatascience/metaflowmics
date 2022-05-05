// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from "./functions"

options = initOptions(params.options)

process MOTHUR_CONSENSUS {
    tag "$meta.id"
    label "process_high"
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options,
                                        publish_dir:getSoftwareName(task.process)) }

    container "quay.io/biocontainers/mothur:1.47.0--hb64bf22_2"
    conda (params.enable_conda ? "bioconda::mothur=1.47.0" : null)

    input:
    tuple val(meta), file(shared), file(list), file(fasta), file(count), file(tax)

    output:
    tuple val(meta), path("*.rep.fasta"), emit: repfasta
    tuple val(meta), path("*.rep.count_table"), emit: repcount_table
    tuple val(meta), path("*.cons.taxonomy"), emit: constaxonomy
    path "*.version.txt", emit: version

    script:
    def software = getSoftwareName(task.process)
    outprefix = "OTUs.${meta.id}"
    """
    mothur "#
    classify.otu(taxonomy=$tax, list=$list, count=$count, probs=f);
    get.oturep(fasta=$fasta, list=current, count=current, method=abundance, rename=t)"

    mv *.rep.fasta ${outprefix}.rep.fasta
    mv *.rep.count_table ${outprefix}.rep.count_table
    mv *.cons.taxonomy ${outprefix}.cons.taxonomy

    # print version
    mothur -v | tail -n+2 | head -1 | cut -d"=" -f2 > ${software}.version.txt
    """
}
