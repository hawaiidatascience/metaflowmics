// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from "./functions"

options = initOptions(params.options)

process MOTHUR_SUMMARY_SHARED {
    tag "$otu_id"
    label "process_low"

    container "quay.io/biocontainers/mothur:1.44.1--hf0cea05_2"
    conda (params.enable_conda ? "bioconda::mothur:1.44.1" : null)

    input:
    tuple val(otu_id), file(shared)

    output:
    path "*.csv", emit: csv
    path "*.version.txt", emit: version

    script:
    def software = getSoftwareName(task.process)
    def ext = shared.getBaseName()
    """
    mothur "#summary.shared(shared=$shared, calc=${params.calc}, distance=lt)"

    mv *.summary beta-diversity_${otu_id}.csv

    # print version
    mothur -v | tail -n+2 | head -1 | cut -d"=" -f2 > ${software}.version.txt
    """
}
