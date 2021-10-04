#!/usr/bin/env nextflow

nextflow.enable.dsl=2
params.options = [:]
params.is_coding = false

module_dir = "../modules"

include{ HMMER_HMMBUILD } from "$module_dir/hmmer/hmmbuild/main.nf" \
    addParams( options: [args: "--amino"] )
include{ HMMER_HMMALIGN } from "$module_dir/hmmer/hmmalign/main.nf" \
    addParams( options: [args: "--amino", publish_dir: "hmmalign"] )
include{ EMBOSS_TRANALIGN } from "$module_dir/emboss/tranalign/main.nf"
include{ SYNC_SEQIDS } from "$module_dir/python/biopython/sync_seqids/main.nf"
include { MOTHUR_ALIGN_SEQS } from "$module_dir/mothur/alignSeqs/main.nf"

include { TRANSLATE } from "./translation.nf"



workflow MSA {
    take:
    contigs // unaligned DNA sequence
	count_table // corresponding abundance
    dbs // reference alignments

	main:

	/*
	 ========================================================================================
     Translation (if coding DNA --> skipped for 16S and ITS)
	 ========================================================================================
	 */
	if (params.is_coding) {
		proteins = TRANSLATE(contigs).faa
		synced = SYNC_SEQIDS(
			proteins,
			contigs.mix(count_table).collect{it[1]}
		)
		contigs = synced.fna
		count_table = synced.count_table
	}

	// Split dbs in nucl or prot
	dbs = dbs.branch{
		nucl: it[0].db_type == "nucl"
		prot: it[0].db_type == "prot"
		tax: true
	}

	/*
	 ========================================================================================
	 Protein alignment with hmmer (if coding DNA and protein ref is provided)
	 ========================================================================================
	 */
	msa_nucl_hmmer = Channel.empty()
	if (params.is_coding) {
		hmmdb = HMMER_HMMBUILD( dbs.prot ).hmm
		msa_prot_hmmer = HMMER_HMMALIGN(proteins, hmmdb).afa
	
		// codon reverse translation for compatibility with the remainder of the pipelines
		msa_nucl_hmmer = EMBOSS_TRANALIGN(
			msa_prot_hmmer.combine(contigs.map{it[1]})
		).fna
	}

	/*
	 ========================================================================================
	 DNA alignment (if dna ref is provided)
	 ========================================================================================
	 */

	msa_nucl_mother = MOTHUR_ALIGN_SEQS(
		contigs.combine(dbs.nucl)
	).fasta

	emit:
	aln = msa_nucl_hmmer.mix(msa_nucl_mother)
	count_table = count_table.map{it[1]} // needed in case some contigs were filtered during translation
}
