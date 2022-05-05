// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from "./functions"

options = initOptions(params.options)

process MOTHUR_SUMMARY_SHARED {
    tag "$meta.id"
    label "process_low"
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options,
                                        publish_dir:getSoftwareName(task.process)) }

    container "quay.io/biocontainers/mothur:1.47.0--hb64bf22_2"
    conda (params.enable_conda ? "bioconda::mothur=1.47.0" : null)

    input:
    tuple val(meta), file(shared)

    output:
    path "*.summary", emit: summary
    path "*.version.txt", emit: version

    script:
    def software = getSoftwareName(task.process)
    def ext = shared.getBaseName()
    """
    mothur "#summary.shared(shared=$shared, calc=${params.calc}, distance=lt)"
    mv *.summary beta-diversity.${meta.otu_id}.summary

    # print version
    mothur -v | tail -n+2 | head -1 | cut -d"=" -f2 > ${software}.version.txt
    """
}
