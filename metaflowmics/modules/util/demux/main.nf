// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

options = initOptions(params.options)

process GUESS_MATCH_ORDER {
    label "process_low"

    container "nakor/metaflowmics-python"
    conda (params.enable_conda ? "conda-forge::" : null)

    input:
    tuple val(meta), path(fastqs)
    path metadata // rename it to make sure we don't have already the same name

    output:
    path "barcodes_meta.csv"

    script:
    z = fastqs[0].getExtension() == 'gz' ? 'z' : ''
    fwd = fastqs[0]
    rev = fastqs.size() == 2 ? fastqs[1] : ''
    """
    #!/usr/bin/env bash

    # Remove carriage returns at the end of file
    sed -e 's/\r//g' $metadata > metadata.csv

    # Reverse complement if set
    if [ "$params.rc" == true ]; then
        cat metadata.csv | rev | cut -d, -f1 | tr "ATGC" "TACG" > rev_idx_rc.csv
        cat metadata.csv | rev | cut -d, -f2- | rev > samp_idx1.csv
        paste -d, samp_idx1.csv rev_idx_rc.csv | sed -e 's/\r//g' > metadata.csv
    fi

    if [ ${fastqs.size()} -eq 2 ]; then

        reversed=[ "$params.matching" == reversed ] && echo true || echo false)

        # if "auto", we check which matching order looks best
        if [ "$params.matching" == auto ]; then
            symbol=\$(${z}cat $fwd | head -c 2)
            
            # extract the most frequent fwd/rev index reads pairs and compare it to the metadata 
            paste --delimiters=' ' \\
                <( ${z}grep --no-group-separator "^\$symbol" $fwd -A1 | grep -v \$symbol ) \\
                <( ${z}grep --no-group-separator "^\$symbol" $rev -A1 | grep -v \$symbol ) \\
                | grep -v N \\
                | sort | uniq -c | sort -rnk1 | head -n 20 \\
                | sed 's/^[[:space:]]*//g' | sed 's/ /,/g' \\
                | cut -d, -f2,3 > freqs.txt
            )

            # take the order that matches best
            n1=\$(comm -12 <(sort <(cut -d, -f2,3 metadata.csv)) <(sort freqs.txt) | wc -l)
            n2=\$(comm -12 <(sort <(cut -d, -f2,3 metadata.csv)) <(sort freqs_rc.txt) | wc -l)

            reversed=\$([ \$n2 -ge \$n1 ] && echo true || echo false)

        [ \$reversed==true ] \\
            && awk -F"," '{OFS=","}{print \$1,\$3,\$2}' metadata.csv > barcodes_meta.csv \\
            || mv metadata.csv barcodes_meta.csv
    else
        awk -F"," '{OFS=","}{print \$1,\$2,NaN}' metadata.csv > barcodes_meta.csv
    fi
    """
}

process TO_H5 {
    tag "$split"
    label "process_medium"

    container "nakor/metaflowmics-python"
    conda (params.enable_conda ? "conda-forge::" : null)

    input:
    tuple val(split), file(fastqs)

    output:
    

    script:
    """
    #!/usr/bin/env python

    
    """
}
        
process ERROR_MODEL {
    tag "$meta.id"
    label "process_low"

    container "nakor/metaflowmics-python"
    conda (params.enable_conda ? "conda-forge::" : null)

    input:
    tuple val(meta), path()

    output:


    script:
    """
    #!/usr/bin/env python

    """
    
}

process IndexMapping {
    tag "$meta.id"
    label "process_low"

    container "nakor/metaflowmics-python"
    conda (params.enable_conda ? "conda-forge::" : null)

    input:
    tuple val(meta), path()

    output:


    script:
    """
    #!/usr/bin/env python

    """
    
}

process SampleSizeDistribution {
    tag "$meta.id"
    label "process_low"

    container "nakor/metaflowmics-python"
    conda (params.enable_conda ? "conda-forge::" : null)

    input:
    tuple val(meta), path()

    output:


    script:
    """
    #!/usr/bin/env python

    """
    
}

process WRITE_SAMPLE_FASTQ {
    tag "$meta.id"
    label "process_low"

    container "nakor/metaflowmics-python"
    conda (params.enable_conda ? "conda-forge::" : null)

    input:
    tuple val(meta), path()

    output:


    script:
    """
    #!/usr/bin/env python

    """
    
}

process GZIP {
    tag "$meta.id"
    label "process_low"

    container "nakor/metaflowmics-python"
    conda (params.enable_conda ? "conda-forge::" : null)

    input:
    tuple val(meta), path()

    output:


    script:
    """
    #!/usr/bin/env python

    """
    
}
