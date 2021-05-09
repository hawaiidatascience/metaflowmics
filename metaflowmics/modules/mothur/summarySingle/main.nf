// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from "./functions"

options = initOptions(params.options)

process MOTHUR_SUMMARY_SINGLE {
    label "process_low"

    container "quay.io/biocontainers/mothur:1.44.1--hf0cea05_2"
    conda (params.enable_conda ? "bioconda::mothur:1.44.1" : null)

    input:
    tuple val(meta), file(shared)

    output:
    path "*.csv", emit: summary
    path "*.version.txt", emit: version

    script:
    def software = getSoftwareName(task.process)
    def ext = shared.getBaseName()
    """
    mothur "#summary.single(shared=$shared, calc=$params.calc)"

    # summary
    if [ "$params.calc" = "nseqs-sobs" ]; then
        cut -f2- *.groups.summary | tail -n+2 \\
        | awk '{OFS=","}{print "$meta.step","$meta.otu_id",\$1,int(\$2),int(\$3)}' \\
        | tr "\\t" "," \\
        > ${ext}.csv
        rm -f *.groups.summary
    else
        mv *.summary alpha-diversity_${meta.otu_id}.csv
    fi

    # print version
    mothur -v | tail -n+2 | head -1 | cut -d"=" -f2 > ${software}.version.txt
    """
}
