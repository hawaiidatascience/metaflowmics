// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options = initOptions(params.options)

process DIAMOND_BLASTP {
    tag "$meta"
    label "process_medium"
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options,
                                        publish_dir:getSoftwareName(task.process),
                                        meta:meta) }    

    conda (params.enable_conda ? "bioconda::diamond=2.0.11" : null)
    container "quay.io/biocontainers/diamond:2.0.11--hdcc8f71_0"

    input:
    tuple val(meta), path(query), path(ref)
    
    output:
    tuple val(meta), path("*.tsv"), emit: tsv

    script:
    """
    diamond blastp --header $options.args \\
    --threads $task.cpus \\
    --query $query \\
    --db $ref \\
    --out blastp.tsv
    """
}


