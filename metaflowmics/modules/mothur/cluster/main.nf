// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from "./functions"

options = initOptions(params.options)

process MOTHUR_CLUSTER {
    tag "$otu_id"
    label "process_high"
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options,
                                        publish_dir:getSoftwareName(task.process)) }

    container "quay.io/biocontainers/mothur:1.44.1--hf0cea05_2"
    conda (params.enable_conda ? "bioconda::mothur:1.44.1" : null)

    input:
    tuple file(fasta), file(count)
    each otu_id

    output:
    tuple val(otu_id), path("${outprefix}.list"), emit: list
    tuple val(otu_id), path("${outprefix}.shared"), emit: shared
    path "*.version.txt", emit: version

    script:
    def software = getSoftwareName(task.process)
    def method = otu_id == 100 ? "unique" : "dgc"
    def otu_diff = (100-otu_id) / 100
    def procname = "${task.process.tokenize(':')[-1].toLowerCase()}"
    def ext = ["rep.fasta", "cons.taxonomy", "shared", "list", "database"]
    outprefix = options.suffix ? "$options.suffix" : "${procname}.${otu_id}"
    """
    mothur "#
    cluster(count=$count, fasta=$fasta, method=$method, cutoff=$otu_diff);
    make.shared(count=current, list=current)"

    mv *.list ${outprefix}.list
    mv *.shared ${outprefix}.shared

    # print version
    mothur -v | tail -n+2 | head -1 | cut -d"=" -f2 > ${software}.version.txt
    """
}
