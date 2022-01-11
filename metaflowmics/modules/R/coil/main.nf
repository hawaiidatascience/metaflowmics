// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

options = initOptions(params.options)

process COIL {
    tag "$meta.id"
    label "process_high"
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options,
                                        publish_dir:getSoftwareName(task.process),
                                        meta:meta) }

    container "nakor/coil:0.0.1"
    // no conda repository for coil

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*.fasta"), emit: afasta
    path "*.version.txt", emit: version

    script:
    def software = getSoftwareName(task.process)
    def ext = abundance.getExtension()
    """
    #!/usr/bin/env Rscript

    library(seqinr)
    library(coil)

    seq = read.fasta("$fasta", as.string=TRUE)

    #parse the data in the header line by splitting the name on the | character

    #step 1: build the coi5p object
    dat = coi5p(seq, name="example_sequence_1")

    #step 2: frame the sequence
    dat = frame(dat)

    #step 3: by default censored translation is performed - see vignette for details
    dat = translate(dat)

    ##step 3a: if taxonomy is known, but the translation table is not, a helper function
    #can be used to look up the proper translation table.
    which_trans_table("Scyliorhinidae")

    #step 3a: the proper translation table can be passed to the translation function
    dat = translate(dat, trans_table = 2)

    writeLines(paste0(packageVersion('coil')), "${software}.version.txt")
    """
}
