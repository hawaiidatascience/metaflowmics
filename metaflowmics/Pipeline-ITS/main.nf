#!/usr/bin/env nextflow

nextflow.enable.dsl=2
params.options = [:]

def helpMessage() {
    log.info"""
    ===================================
     Metaflow|mics ITS pipeline
    ===================================
    Usage:
    nextflow run ITS-pipeline --reads '*_R1.fastq.gz' --locus ITS1 -profile local

    ---------------------------------- Mandatory arguments ----------------------------------------

      --reads         Path to input data

    ---------------------------------- Optional arguments ----------------------------------------

      --paired_end     If your data is paired-end
      --locus         Sequenced ITS region (ITS1,ITS2 or ALL). Default: ITS1
      --outdir        The output directory where the results will be saved. 
                      Default: ./ITS-pipeline_outputs

    [Quality filtering]
      --max_expected_errors  Discard reads with more than <max_expected_errors> errors. Default: 3
      --min_read_count       Discard samples with less than <minReads> reads. Default: 50
      --min_derep_count      Discard samples with less than <minDerep> unique contigs. Default: 5

    [Taxonomy assignment]
      --tax_confidence     The minimum consensus for calling a taxonomy. Default: 50%

    [Co-occurence pattern correction]
	--skip_lulu                 Skip the Lulu step (can save a lot of computation time)
    --lulu_min_ratio_type       Function to compare the abundance of a parent OTU and its daughter.
                                  Default: min
    --lulu_min_ratio            Minimum abundance ratio between parent and daughter
                                  (across all samples). Default: 1
    --lulu_min_rel_cooccurence  Proportion of the parent samples where daughter occurs. Default: 1
    """.stripIndent()
}

// Show help message
params.help = false
if (params.help){
    helpMessage()
    exit 0
}

def summary = [:]
summary['paired end'] = params.paired_end
summary['locus'] = params.locus
summary['max expected errors'] = params.max_expected_error
summary['min reads per sample'] = params.min_read_count
summary['min unique sequence per sample'] = params.min_derep_count
summary['clustering similarity thresholds'] = params.clustering_thresholds
summary['skip lulu'] = params.skip_lulu
summary['lulu ratio type'] = params.lulu_min_ratio_type
summary['lulu parent/daughter min similarity'] = params.lulu_min_match
summary['lulu parent/daughter min abundance ratio'] = params.lulu_min_ratio
summary['lulu parent/daughter min co-occurrence'] = params.lulu_min_rel_cooccurence
summary['consensus taxonomy confidence threshold'] = params.tax_confidence

file(params.outdir).mkdir()
File f = new File("${params.outdir}/parameters_summary.log")
f.write("====== Parameter summary =====\n")
summary.each { key, val -> f.append("\n${key.padRight(50)}: ${val}") }

include{ DOWNLOAD_UNITE } from '../modules/util/download/main.nf' \
    addParams( db_release: 'fungi' )
include{ ITSXPRESS } from '../modules/itsxpress/main.nf' \
    addParams( options: [publish_dir: '1-ITSxpress'] )
include{ DADA2_FILTERANDTRIM } from '../modules/dada2/filterAndTrim/main.nf' \
    addParams( options: [publish_dir: '2-Quality_filtering'], trunc_len: 0, trunc_quality: 2, min_read_len: 20 )
include{ DADA2_DEREPFASTQ } from '../modules/dada2/derepFastq/main.nf' \
    addParams( options: [publish_dir: '3-Chimera'] )
include{ VSEARCH_CHIMERA } from '../modules/vsearch/chimera/main.nf' \
    addParams( options: [publish_dir: '3-Chimera'] )
include{ DADA2_LEARNERRORS } from '../modules/dada2/learnErrors/main.nf' \
    addParams( options: [publish_dir: '4-Denoising'] )
include{ DADA2_DADA } from '../modules/dada2/dada/main.nf' \
    addParams( options: [publish_dir: '4-Denoising'] )
include{ SUBSET_READS_RDS; BUILD_ASV_TABLE } from '../modules/util/dada2/main.nf' \
    addParams( options: [publish_dir: '5-OTU_clustering'], format: 'VSEARCH' )
include{ VSEARCH_CLUSTER } from '../modules/vsearch/cluster/main.nf' \
    addParams( options: [publish_dir: '5-OTU_clustering'] )
include{ VSEARCH_USEARCHGLOBAL } from '../modules/vsearch/usearchglobal/main.nf' \
    addParams( options: [publish_dir: '6-LULU'] )
include{ LULU } from '../modules/lulu/main.nf' \
    addParams( options: [publish_dir: '6-LULU'] )
include{ DADA2_ASSIGN_TAXONOMY } from '../modules/dada2/assignTaxonomy/main.nf' \
    addParams( options: [publish_dir: '7-Taxonomy'] )
include{ SUMMARIZE_TABLE; READ_TRACKING } from '../modules/util/misc/main.nf' \
    addParams( options: [publish_dir: '8-Read_tracking'] )

workflow pipeline_ITS {
    take:
    reads

    main:
    tracked_files = Channel.empty()
    raw_counts = reads.map{['raw', it[1][0].countFastq(), it]}

    // Extract ITS marker
    its = ITSXPRESS(
        raw_counts.filter{it[1] > params.min_read_count}.map{it[2]}
    )
    its_counts = its.fastq.map{['ITS', it[1].countFastq(), it]}

    // Remove low quality reads
    qc = DADA2_FILTERANDTRIM(
        its_counts.filter{it[1] > params.min_read_count/2}.map{it[2]}
    )
    qc_counts = qc.fastq.map{['QC', it[1].countFastq(), it]}

    // Dereplicate before chimera removal and denoising
    derep = DADA2_DEREPFASTQ(
        qc_counts.filter{it[1] > params.min_read_count/4}.map{it[2]}
    )

    // Flag chimera
    chimera = VSEARCH_CHIMERA(
        derep.fasta
    )

    // Remove chimera from DADA2 derep object
    nochim = SUBSET_READS_RDS(
        derep.rds.join(chimera.fasta)
    )

    // Build Illumina reads error model
    error_model = DADA2_LEARNERRORS(
        nochim.rds
    )

    // Denoise reads
    dada = DADA2_DADA(
        nochim.rds.join(error_model.rds)
    )

    // Make raw ASV table
    asvs = BUILD_ASV_TABLE(
        dada.denoised.collect{it[1]}
    )

    // Cluster ASVs into OTUs
    otus = VSEARCH_CLUSTER(
        asvs.fasta_dup,
        params.clustering_thresholds.split(',').findAll({it!='100'}).collect{it as int}
    )
    otus_summary = SUMMARIZE_TABLE(otus.tsv.map{["clustering-${it[0]}", it[1]]})

    otu_fastas = asvs.fasta.map{[100, it]}.mix(otus.fasta)
    otu_tables = asvs.count_table.map{[100, it]}.mix(otus.tsv)

    if (!params.skip_lulu) {
        // Co-occurrence pattern detection
        matchlist = VSEARCH_USEARCHGLOBAL(
            otu_fastas
        )
        lulu = LULU(
            matchlist.tsv.join(otu_tables).join(otu_fastas)
        )
        tracked_files = tracked_files.mix(lulu.summary)
        otu_fastas = lulu.fasta
    }

    // Taxonomy assignment with Sintax
    unite_db = DOWNLOAD_UNITE()
    taxonomy = DADA2_ASSIGN_TAXONOMY(
        otu_fastas,
        unite_db
    )

    // Track read in the pipeline
    fq_counts = raw_counts.mix(its_counts).mix(qc_counts)
        .collectFile(name: 'summary.csv', newLine: true){"${it[0]},${it[2][0].id},${it[1]},"}
    
    tracked_files = tracked_files
        .concat(fq_counts, derep.summary, nochim.summary, dada.summary, otus_summary)

    summary = READ_TRACKING(
        tracked_files.collectFile()
    )
}

workflow {
    reads = Channel.fromFilePairs(params.reads, size: params.paired_end ? 2 : 1)
        .map{[ [id: it[0], paired: params.paired_end], it[1] ]}
    pipeline_ITS(reads)
}
