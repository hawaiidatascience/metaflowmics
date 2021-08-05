// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from "./functions"

options = initOptions(params.options)

process MUSCLE {
    tag "${fasta.getBaseName()}"
    label "process_high"
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options,
                                        publish_dir:getSoftwareName(task.process)) }

    container "quay.io/biocontainers/muscle:3.8.1551--h7d875b9_6"
    conda (params.enable_conda ? "bioconda::muscle=3.8" : null)

    input:
    path fasta

    output:
    path "*.afa", emit: afa
    path "*.version.txt", emit: version

    script:
    def software = getSoftwareName(task.process)
    def arg = fasta.getExtension() == "afa" ? "-profile" : ""
    """
    #!/usr/bin/env bash

    n=\$(grep -c "^>" $fasta)

    if [ "$arg" == "" ] && [ \$n -ge 500 ] ; then
        opt="-maxiters 1 -diags -sv -distance1 kbit20_3"
    else
        opt=""
    fi

    # Run muscle and convert to 2-line fasta
    cat $fasta | muscle $arg \$opt \\
        | awk '/^>/ {printf("\\n%s\\n",\$0);next; } { printf("%s",\$0);}  END {printf("\\n");}' \\
        | tail -n+2 > ${fasta.getBaseName()}.afa

    echo \$(muscle --version 2>&1) | cut -d" " -f2  > ${software}.version.txt
    """
}
