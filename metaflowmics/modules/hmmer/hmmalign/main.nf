// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process HMMER_HMMALIGN {
    label "process_medium"
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options,
                                        publish_dir:getSoftwareName(task.process)) }

    conda (params.enable_conda ? "bioconda::hmmer=3.3.2" : null)
    container "quay.io/biocontainers/hmmer:3.3.2--h1b792b2_1"

    input:
    tuple val(meta), path(fasta)
    tuple val(db_meta), path(hmm)

    output:
    tuple val(meta_upd), path("*.afa"), emit: afa
    path "*.version.txt", emit: version

    script:
	meta_upd = meta + db_meta
    def software = getSoftwareName(task.process)
    def prefix   = fasta.getBaseName()
    def fastacmd = fasta.getExtension() == 'gz' ? "gunzip -c $fasta" : "cat $fasta"
    """
    $fastacmd |
        sed 's/[[:space:]]/|/g' |
        hmmalign --outformat Pfam $options.args $hmm - |
        egrep -v "^#|^\$|^//\$" |
        awk '{print ">"\$1"\\n"\$2}' |
        sed 's/|/ /g' > ${prefix}.afa

    echo \$(hmmalign -h | grep -o '^# HMMER [0-9.]*') | sed 's/^# HMMER *//' > ${software}.version.txt
    """
}
