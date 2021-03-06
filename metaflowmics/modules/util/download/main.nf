process DOWNLOAD_IBOL {
    tag "v$version"
    label "process_low"

    conda (params.enable_conda ? "bash:5.0.018" : null)
    container "nakor/bash:5.1.4"

    input:
    each version
    
    output:
    tuple val(version), path("iBOL_COI*.fna"), emit: fna
    path "iBOL_COI*.tsv", emit: tsv

    script:
    url = "https://v3.boldsystems.org/data/datarelease/NewPackages"
    name = "iBOL_phase_${version}_COI.tsv"
    prefix = "iBOL_COI"
    ranks = ['phylum', 'class', 'order', 'family', 'subfamily', 'genus', 'species']
    tax_fmt = ranks.withIndex().collect{ r,i ->
        symbol = r[0]
        if (r=='subfamily') {symbol='sf'}
        "${symbol}__\"\$${i+9}\""
    }.join(',')
    regex = "[^A-Za-z0-9_-]"
    """
    wget $url/${name}.zip && unzip ${name}.zip && rm -f ${name}.zip
    
    awk -F'\\t' '\$31 != ""' $name |
      awk '{FS=OFS="\\t"} { 
        if (\$14=="") \$14=gensub(/ .*/,"","g",\$15);
        \$15=gensub(/sp\\./,"sp","g",\$15);
        \$15=gensub(/$regex/,"_","g",\$15);
        print
      }' > ${prefix}.tsv && rm -f $name
    
    tail -n+2 ${prefix}.tsv | 
      awk -F'\\t' '{
        print ">"\$28";"\$36";k__Animalia,$tax_fmt\\n"\$31
      }' > ${prefix}.fna
    """
}

process DOWNLOAD_UNITE {
    tag "$params.db_release"
    label "process_low"

    conda (params.enable_conda ? "bash:5.0.018" : null)
    container "nakor/bash:5.1.4"

    output:
    path "*.fasta", emit: fasta

    script:
    root_url = "https://files.plutof.ut.ee/public/orig"
    if (params.db_release == 'fungi') {
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
    tag "$params.db_release"
    label "process_low"

    conda (params.enable_conda ? "bash:5.0.018" : null)
    container "nakor/bash:5.1.4"

    output:
    path "*.tax", emit: tax
    path "*.align", emit: align

    script:
    url_base = "https://mothur.s3.us-east-2.amazonaws.com/wiki"
    file = "silva.${params.db_release}_v138.tgz"
    """
    wget -qO- $url_base/$file | tar xz
    """
}

process DOWNLOAD_RDP_FOR_DADA2 {
    tag "$params.db_release"
    label "process_low"

    conda (params.enable_conda ? "bash:5.0.018" : null)
    container "nakor/bash:5.1.4"

    output:
    path "*.fa", emit: fasta

    script:
    url_base = "https://zenodo.org/record/4587955/files"
    file = "silva_nr99_v138.1_wSpecies_train_set.fa.gz?download=1"
    """
    wget -qO- "$url_base/$file" | gunzip \\
        > silva_nr99_v138.1_wSpecies_train_set.fa
    """
}

