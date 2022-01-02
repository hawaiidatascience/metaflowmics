module_dir = "../modules"

// `compile` sub-workflow
include { MOTHUR_CLASSIFY_OTUS } from "$module_dir/mothur/classifyOtus/main.nf" \
    addParams( options: [publish_dir: "consensus_taxonomy"] )
include { MOTHUR_GET_OTU_REP } from "$module_dir/mothur/getOtuRep/main.nf" \
    addParams( options: [publish_dir: "representative_sequences"] )
include { MOTHUR_MAKE_DATABASE } from "$module_dir/mothur/makeDatabase/main.nf" \
    addParams( options: [publish_dir: "mothur_db_files"] )

// `sync` sub-workflow
include { MOTHUR_GET_SEQS } from "$module_dir/mothur/getSeqs/main.nf"
include { MOTHUR_GET_OTUS } from "$module_dir/mothur/getOtus/main.nf"

// read tracking
include { MOTHUR_SUMMARY_SINGLE } from "$module_dir/mothur/summarySingle/main.nf" \
    addParams( calc: "nseqs-sobs" )
include{ SUMMARIZE_TABLE } from "$module_dir/util/misc/main.nf"


workflow SYNC {
    // Subset all files to make them coherent
    take:
    files
    list
    shared

    main:
    // The shared file is the reference
    // Use it to fix the list file
    list = MOTHUR_GET_OTUS(
        shared.join(list)
    )

    // Fix the other files with the list
    files = MOTHUR_GET_SEQS(
        list.combine(files, by: 0)
    ).branch{
		fasta: it[-1].getExtension() == "fasta"
		count_table: it[-1].getExtension() == "count_table"
		taxonomy: it[-1].getExtension() == "taxonomy"
		shared: it[-1].getExtension() == "shared"
		summary: it.getExtension() == "summary"		
	}

    emit:
    fasta=files.fasta
    count_table=files.count_table
    taxonomy=files.taxonomy
    list=list
}


workflow CONSENSUS {
    take:
    fasta
    count_table
    taxonomy
    list
    shared

    main:
    // 1) Consensus taxonomy
    cons = MOTHUR_CLASSIFY_OTUS(
        list.join(count_table).join(taxonomy)
    )

    // 2) Representative sequences
    rep = MOTHUR_GET_OTU_REP(
        list.join(fasta).join(count_table)
    )

    // 3) Database
    if (params.compute_mothur_db) {
        db = MOTHUR_MAKE_DATABASE(
            shared.join(cons.taxonomy).join(rep.fasta).join(rep.count_table)
        ).database
    }

    shared.map{it[1].copyTo("$params.outdir/abundance_tables/${it[0].id}.shared")}

    emit:
    constaxonomy=cons.taxonomy
    repfasta=rep.fasta
    repcount=rep.count_table
}
