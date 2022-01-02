// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from "./functions"

options = initOptions(params.options)

process EMBOSS_TRANSEQ {
    label "process_medium"
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
		pattern: "*.single.faa",
        saveAs: { filename -> saveFiles(filename:filename, options:params.options,
                                        publish_dir:getSoftwareName(task.process)) }

    container "nakor/metaflowmics-biotools:0.0.1"
    conda (params.enable_conda ? "bioconda::emboss=6.6.0" : null)

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*.faa"), emit: faa // raw output from transeq
    tuple val(meta), path("*.single.faa"), emit: single
    tuple val(meta), path("*.multiple.faa"), emit: multiple
    path "*.version.txt", emit: version

    script:
    def software = getSoftwareName(task.process)
    outprefix = fasta.getSimpleName()
    """
    #!/usr/bin/env bash

    ## 1) Translate the sequence
    ## 2) Replace title suffix pattern from "_{frame}" to " {frame}"
    ## 3) Reformat to two-line fasta and remove short sequences

    transeq -auto -trim \\
      -outseq stdout \\
      -sequence $fasta \\
      -table $params.table \\
      -frame $params.frames \\
    | sed -E '/^>/s/([^ ]+)_([1-6])( .*)?/\\1\\3 frame=\\2/g' \\
    | seqtk seq -l0 \\
    > ${outprefix}.faa

    ## 4) Discard sequences with stop codons in the middle of the sequence
    grep -B1 -v '[>*]' --no-group-separator ${outprefix}.faa > ${outprefix}.nostop.faa

    ## 5) Find the sequence with single or multiple possible ORFs
    grep '^>' ${outprefix}.nostop.faa | cut -d' ' -f1 | uniq -c > freqs.txt
    awk '\$1==1' freqs.txt | cut -d'>' -f2 > single_frame.txt
    awk '\$1>1' freqs.txt | cut -d'>' -f2 > multiple_frames.txt

    ## 6) Subset the translated sequences with those ids
    seqtk subseq ${outprefix}.nostop.faa single_frame.txt > ${outprefix}.single.faa
    seqtk subseq ${outprefix}.nostop.faa multiple_frames.txt > ${outprefix}.multiple.faa

    transeq -h 2>&1 | grep -i version | cut -d':' -f3 > ${software}.version.txt
    """
}
