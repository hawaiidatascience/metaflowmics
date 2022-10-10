// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process RDP_TRAIN {
    label "process_high"
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options,
                                        publish_dir:getSoftwareName(task.process)) }

    conda (params.enable_conda ? "bioconda::rdptools=2.0.3" : null)
    container "quay.io/biocontainers/rdptools:2.0.3--hdfd78af_1"

    input:
	path fasta
	path tax

    output:
    path "*"

    script:
    def software = getSoftwareName(task.process)
    def prefix   = fasta.getBaseName()
    """
    classifier -Xmx${task.memory.getGiga()}g train -o db -s $fasta -t $tax 
    """
}


