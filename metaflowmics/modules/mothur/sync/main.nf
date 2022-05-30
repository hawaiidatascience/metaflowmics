// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from "./functions"

options = initOptions(params.options)

process MOTHUR_SYNC {
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
    tuple val(meta), path("*.fasta"), emit: fasta
    tuple val(meta), path("*.count_table"), emit: count_table
    tuple val(meta), path("*.taxonomy"), emit: taxonomy    
    tuple val(meta), path("*.list"), emit: list
    path "*.version.txt", emit: version

    script:
    def software = getSoftwareName(task.process)
    def outprefix = "OTUs.${meta.id}"
    """
    mothur "#
    list.otulabels(shared=$shared);
    get.otus(accnos=current, list=$list);
    list.seqs(list=current);
    get.seqs(accnos=current, fasta=$fasta);    
    get.seqs(accnos=current, count=$count);    
    get.seqs(accnos=current, taxonomy=$tax)"
 
    mv *.pick.list ${outprefix}.list
    mv *.pick.fasta ${outprefix}.fasta
    mv *.pick.count_table ${outprefix}.count_table
    mv *.pick.taxonomy ${outprefix}.taxonomy

    # print version
    mothur -v | tail -n+2 | head -1 | cut -d"=" -f2 > ${software}.version.txt
    """
}
