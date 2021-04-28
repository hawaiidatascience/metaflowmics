// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

options = initOptions(params.options)

process MOTHUR_SUMMARY_SINGLE {
    tag "$step"
    label "process_high"

    container "quay.io/biocontainers/mothur:1.44.1--hf0cea05_2"
    conda (params.enable_conda ? "bioconda::mothur:1.44.1" : null)

    input:
    tuple val(step), val(otu_id), file(shared)

    output:
    path "*.summary", emit: summary
    path "*.version.txt", emit: version

    script:
    def software = getSoftwareName(task.process)
    def ext = shared.getBaseName()
    """
    mothur '#summary.single(shared=$shared, calc=nseqs-sobs)'

    # summary
    cut -f2- *.groups.summary | tail -n+2 \\
    | awk '{OFS=","}{print "$step","$otu_id",\$1,int(\$2),int(\$3)}' \\
    | tr '\\t' ',' \\
    > ${ext}.summary

    rm -f *.groups.summary

    # print version
    mothur -v | tail -n+2 | head -1 | cut -d'=' -f2 > ${software}.version.txt
    """
}
