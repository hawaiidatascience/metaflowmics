// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process RDP_CLASSIFY {
    label "process_high"
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options,
                                        publish_dir:getSoftwareName(task.process)) }

    conda (params.enable_conda ? "bioconda::rdptools=2.0.3" : null)
    container "quay.io/biocontainers/rdptools:2.0.3--hdfd78af_1"

    input:
    tuple val(meta), path(fasta)
    path rdp_db_files

    output:
    tuple val(meta), path("*.tsv"), emit: tsv
    tuple val(meta), path("*.taxonomy"), emit: taxonomy

    script:
    def software = getSoftwareName(task.process)
    def prefix   = fasta.getBaseName()
    """
    classifier -Xmx${task.memory.getGiga()}g classify -t rRNAClassifier.properties -o ${prefix}.tsv $fasta -c $params.rdp_confidence
	awk -F'\\t' '{printf \$1""FS; for(i=6;i<=NF;i=i+3) printf \$i";"; print ""}' ${prefix}.tsv > ${prefix}.taxonomy
    """
}
