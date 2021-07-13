process DIAMOND_MAKEDB {
    label "process_low"

    conda (params.enable_conda ? "bioconda::diamond=2.0.11" : null)
    container "quay.io/biocontainers/diamond:2.0.11--hdcc8f71_0"

    input:
    path fasta
    
    output:
    path "*.dmnd", emit: db

    script:
    """
    diamond makedb \\
    --threads $task.cpus \\
    --db ${fasta.getBaseName()} \\
    --in $fasta
    """
}


