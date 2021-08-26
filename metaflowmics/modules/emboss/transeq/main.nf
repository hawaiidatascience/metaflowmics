// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from "./functions"

options = initOptions(params.options)

process EMBOSS_TRANSEQ {
    label "process_high"
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options,
                                        publish_dir:getSoftwareName(task.process)) }

    container "nakor/metaflowmics-biotools:0.0.1"
    conda (params.enable_conda ? "bioconda::emboss=6.6.0" : null)

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*.faa"), emit: faa // raw output from transeq
    tuple val(meta), path("*_single.faa"), emit: single
    tuple val(meta), path("*_multiple.faa"), emit: multiple
    path "*.version.txt", emit: version

    script:
    def software = getSoftwareName(task.process)
    outprefix = fasta.getBaseName()
    """
    #!/usr/bin/env bash

    ## 1) Translate the sequence
    ## 2) Replace title suffix pattern from "_{frame}" to " {frame}"
    ## 3) Reformat to two-line fasta and remove short sequences
    ## 4) Put the commas back in the lineage header
    ## 5) Remove trailing stop codons

    transeq -auto -trim \\
      -outseq stdout \\
      -sequence $fasta \\
      -table $params.table \\
      -frame $params.frames \\
    | sed -E 's/_([1-6])\$/ \\1/g' \\
    | seqtk seq -l0 \\
    | sed -E '/^>/ s/_([a-z][a-z]?)__/,\\1__/g' \\
    > ${outprefix}.faa

    ## 6) Discard sequence with stop codons in the middle of the sequence
    grep -B1 -v '[>*]' --no-group-separator ${outprefix}.faa > ${outprefix}_no_stop.faa

    ## 7) Find the sequence with single or multiple possible ORFs
    grep '^>' ${outprefix}_no_stop.faa | sed 's/ [1-6]\$//g' | uniq -c > freqs.txt
    awk '\$1==1' freqs.txt | cut -d'>' -f2 > single_frame.txt
    awk '\$1>1' freqs.txt | cut -d'>' -f2 > multiple_frames.txt

    ## Subset the translated sequences with those ids
    seqtk subseq ${outprefix}_no_stop.faa single_frame.txt > ${outprefix}_single.faa
    seqtk subseq ${outprefix}_no_stop.faa multiple_frames.txt > ${outprefix}_multiple.faa

    transeq -h 2>&1 | grep -i version | cut -d':' -f3 > ${software}.version.txt
    """
}
