// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process HMMER_HMMBUILD {
	tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::hmmer=3.3.2" : null)
    container "quay.io/biocontainers/hmmer:3.3.2--h1b792b2_1"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("hmmdb"), emit: hmm
    path "*.version.txt", emit: version

    script:
    def software = getSoftwareName(task.process)
    """
    hmmbuild $options.args --cpu $task.cpus hmmdb $fasta

    echo \$(hmmbuild -h | grep -o '^# HMMER [0-9.]*') | sed 's/^# HMMER *//' > ${software}.version.txt
    """
}
