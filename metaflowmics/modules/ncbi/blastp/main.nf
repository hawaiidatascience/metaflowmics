process BLASTP {
    label "process_low"

    conda (params.enable_conda ? "bioconda::blast=2.11.0" : null)
    container "quay.io/biocontainers/blast:2.11.0--pl5262h3289130_1"

    input:
    tuple val(meta), path(query), path(ref)
    
    output:
    path "blastp.tsv", emit: tsv

    script:
    """
    makeblastdb -in $ref -dbtype prot

    blastp $params.args -db $ref -query $query \\
        -outfmt "6 qseqid sseqid qlen slen length nident mismatch gapopen pident qstart qend sstart send evalue bitscore" \\
        -out blastp.tsv
    """    
}


