process DOWNLOAD_UNITE {
    tag "download_unite_db"
    label "process_low"

    conda (params.enable_conda ? "bash:5.0.018" : null)
    container "nakor/bash:5.1.4"

    output:
    path "*.fasta", emit: fasta

    script:
    root_url = "https://files.plutof.ut.ee/public/orig"
    if (params.release == 'fungi') {
        url_base = "${root_url}/1E/66"
        file = "1E662B6EB320312A61E7E3218327F34C7DB09CFF8E4686A89EF47886822DA6AB.gz"
        """
        wget -qO- $url_base/$file | tar xz
        iconv -f utf-8 -t ascii sh_general_release*/*.fasta \
            > unite_fungi.fasta
        rm -rf sh_general_release*
        """
    } else {
        url_base = "${root_url}/BF/49"
        file = "BF49FBF4B47314A1CC5238B280FC58BFB8CEBD44A8D45F4A2BF5B8A466715693.gz"
        """
        wget -qO- $url_base/$file \\
            | gunzip \\
            | iconv -f utf-8 -t ascii \\
            > unite_all_eukaryotes.fasta
        rm -f $file
        """
    }
}

process DOWNLOAD_SILVA_FOR_MOTHUR {
    tag "download_silva_db"
    label "process_low"

    conda (params.enable_conda ? "bash:5.0.018" : null)
    container "nakor/bash:5.1.4"

    output:
    path "*.tax", emit: tax
    path "*.align", emit: align

    script:
    url_base = "https://mothur.s3.us-east-2.amazonaws.com/wiki"
    file = "silva.nr_v138_1.tgz"
    """
    wget -qO- $url_base/${file}.tgz | tar xz
    """
}

