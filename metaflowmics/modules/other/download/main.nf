process DOWNLOAD_UNITE {
    tag "download_unite_db"
    label "process_low"

    conda (params.enable_conda ? "bash:5.0.018" : null)
    container "nakor/bash:5.1.4"

    output:
    path "*.fasta", emit: fasta

    shell:
    url_base = "https://files.plutof.ut.ee/public/orig/BF/49"
    file = "BF49FBF4B47314A1CC5238B280FC58BFB8CEBD44A8D45F4A2BF5B8A466715693"
    """
    wget $url_base/${file}.gz && gunzip ${file}.gz
    mv $file unite_db.fasta
    """
}

process DOWNLOAD_SILVA_FOR_MOTHUR {
    tag "download_silva_db"
    label "process_low"

    conda (params.enable_conda ? "bash:5.0.018" : null)
    container "nakor/bash:5.1.4"

    output:
    path "*.tax", emit: tax
    path "*.align", emit: align

    shell:
    url_base = "https://mothur.s3.us-east-2.amazonaws.com/wiki"
    file = "silva.nr_v138_1.tgz"
    """
    wget -qO- $url_base/${file}.tgz | tar xz
    """
}

