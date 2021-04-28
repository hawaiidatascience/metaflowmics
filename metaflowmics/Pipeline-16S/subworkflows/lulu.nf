nextflow.enable.dsl=2

moduledir = "../../modules"

include { VSEARCH_USEARCH_GLOBAL } from "$moduledir/vsearch/usearchGlobal/main.nf"
include { LULU } from "$moduledir/lulu/main.nf" \
    addParams( options: [publish_dir: "interm/13-Lulu"] )

workflow lulu {
    take:
    fasta
    abundance

    main:
    matchlist = VSEARCH_USEARCH_GLOBAL(
        fasta
    )
    lulu = LULU(
        matchlist.tsv.join(abundance).join(fasta)
    )

    emit:
    fasta = lulu.fasta
    abundance = lulu.abundance
    tracking = lulu.summary
}
