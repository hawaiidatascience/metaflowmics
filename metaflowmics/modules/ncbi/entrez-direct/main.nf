process GET_TRANSLATION_TABLES {
    label "process_low"

    conda (params.enable_conda ? "bioconda::entrez-direct=13.9" : null)
    container "quay.io/biocontainers/entrez-direct:13.9--pl5262he881be0_2"

    input:
    path labels
    
    output:
    path "tables.tsv", emit: txt

    script:
    """
    ids=\$(sed ':a;N;\$!ba;s/\\n/ [$params.level] OR /g' $labels)

    esearch -db taxonomy -query "\$ids" | 
        efetch -format xml |
        xtract -pattern Taxon -element ScientificName,TaxId,MGCId >> tables.tsv
    """    
}


