// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from "./functions"

options = initOptions(params.options)

process MOTHUR_SUMMARY_SINGLE {
    tag "$meta.id"
    label "process_low"
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        pattern: "*.summary",
        saveAs: { filename -> saveFiles(filename:filename, options:params.options,
                                        publish_dir:getSoftwareName(task.process)) }

    container "quay.io/biocontainers/mothur:1.47.0--hb64bf22_2"
    conda (params.enable_conda ? "bioconda::mothur=1.47.0" : null)

    input:
    tuple val(meta), val(step), file(shared)

    output:
    tuple val(meta), path("*.{csv,summary}"), emit: summary
    path "*.version.txt", emit: version

    script:
    def software = getSoftwareName(task.process)
    def ext = shared.getBaseName()
    """
    mothur "#make.shared(shared=$shared); summary.single(shared=current, calc=$params.calc)"

    # summary
    if [ "$params.calc" = "nseqs-sobs" ]; then
        cut -f2- *.groups.summary | tail -n+2 \\
        | awk '{OFS=","}{print "$step","$meta.otu_id",\$1,int(\$2),int(\$3)}' \\
        | tr "\\t" "," \\
        > ${ext}.csv
        rm -f *.groups.summary
    else
        mv *.summary alpha-diversity.${meta.otu_id}.summary
    fi

    # print version
    mothur -v | tail -n+2 | head -1 | cut -d"=" -f2 > ${software}.version.txt
    """
}
