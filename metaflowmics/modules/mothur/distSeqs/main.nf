// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from "./functions"

options = initOptions(params.options)

process MOTHUR_DIST_SEQS {
    tag "$meta.id"
    label "process_high"

    container "quay.io/biocontainers/mothur:1.47.0--hb64bf22_2"
    conda (params.enable_conda ? "bioconda::mothur=1.47.0" : null)

    input:
    tuple val(meta), file(fasta)

    output:
    tuple val(meta), path("$outprefix*.dist"), emit: dist
    path "*.version.txt", emit: version

    script:
    def software = getSoftwareName(task.process)
    def cutoff = params.cutoff ?: 1
    def procname = "${task.process.tokenize(':')[-1].toLowerCase()}"
    outprefix = "${procname}.${meta.id}"    
    """
    # Manually rename sequences
    # sed "/^>/s/.*\\(Otu[0-9]*\\)\\(.*\\)/>\\1\\t\\1\\2/" $fasta > renamed.fasta

    mothur "#filter.seqs(fasta=$fasta, trump=.); dist.seqs(fasta=current, cutoff=$cutoff)"

    if [ "${params.format.toLowerCase()}" == "vsearch" ]; then
        awk '{OFS="\\t"}{print \$1,\$2,100*(1-\$3)}' *.dist > ${outprefix}.dist
    else
        mv *.dist ${outprefix}.dist
    fi

    # print version
    mothur -v | tail -n+2 | head -1 | cut -d"=" -f2 > ${software}.version.txt
    """
}
