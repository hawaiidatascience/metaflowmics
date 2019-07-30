#!/usr/bin/env nextflow

def helpMessage() {
    log.info"""
    ===================================
    16S-rDNA-pipeline
    ===================================
    Usage:
    nextflow run 16S-pipeline --reads '*_R{1,2}.fastq.gz' -profile manoa_hpc
    """.stripIndent()
}

/*
 * SET UP CONFIGURATION VARIABLES
 */

// Show help message
params.help = false
if (params.help){
    helpMessage()
    exit 0
}

//script_dir = workflow.projectDir+"/scripts"
script_dir = '/workspace'

def clusteringThresholds = params.clusteringThresholds.split(',').collect{it as int}

/*
 *
 Beginning of the pipeline
 *
 */

if ( params.singleEnd ) {
    read_path = params.reads.replaceAll("\\{1,2\\}","1")
} else {
    read_path = params.reads
}

Channel
    .fromFilePairs( read_path, size: params.singleEnd ? 1 : 2 )
    .ifEmpty { error "Cannot find any reads matching: ${params.reads}" }
    .map { it -> tuple(it[0].replaceAll("-","_"), it[1]) }
    .into { INPUT_FASTQ ; INPUT_FASTQ_FOR_COUNT }

RAW_COUNTS = INPUT_FASTQ_FOR_COUNT.collect{ ["'${it[0]}': ${it[1][0].countFastq()}"] }

/*
 *
 * Dada2's filterAndTrim function (in scripts/util.R)
 * Parameters: {params.maxEE} (expected errors), {params.truncLen}, {params.truncQ}
 *
 */

process FilterAndTrim {
    // Quality filter and trimming using dada2. 
    tag { "FilterAndTrim.${pairId}" }
    publishDir params.outdir+"1-FilterAndTrim", mode: "copy", pattern: "*.{fastq,png}"
    label "low_computation"
    label "r_script"
    
    input:
        set val(pairId), file(fastq) from INPUT_FASTQ
    output:
        set val(pairId), file("${pairId}*_trimmed.fastq") into FASTQ_TRIMMED_RAW
        file("*.png") optional true

    script:
    """
    #!/usr/bin/env Rscript

    source("${script_dir}/util.R")  

    fastqs <- c("${fastq.join('","')}")
    rev <- NULL

    if ( "${params.singleEnd}" == "false" ) {
        rev <- fastqs[2]
    }

    tryCatch(
    {
        filterReads("${pairId}", fastqs[1], rev=rev, minLen=${params.minLength},maxEE=${params.maxEE},truncLen=${params.truncLen},rm.phix=${params.rmphix}, truncQ=${params.truncQ})
    },
    error=function(cond) {
        file.create("${pairId}_R1_trimmed.fastq")
    })
    """
}

FASTQ_TRIMMED_RAW
    .filter{ (it[1].size() > 0)  }
    .into{ FASTQ_TRIMMED ; FASTQ_TRIMMED_FOR_MODEL ; FASTQ_TRIMMED_FOR_COUNT }

FILTERED_COUNTS = FASTQ_TRIMMED_FOR_COUNT.collect{ ["'${it[0]}': ${it[1][0].countFastq()}"] }

/*
 *
 * Dada2 Error model (in scripts/util.R)
 *
 */

process LearnErrors {
    tag { "LearnErrors.${pairId}" }
    publishDir params.outdir+"2-ErrorModel", mode: "copy", pattern: "*.{RDS,png}"
    label "medium_computation"
    label "r_script"
    
    input:
	set val(pairId), file(fastq) from FASTQ_TRIMMED_FOR_MODEL
    output:
	set val(pairId), file("${pairId}*.RDS") into ERROR_MODEL
        file "*.png"
    
    script:
    """
    #!/usr/bin/env Rscript
    source("${script_dir}/util.R") 

    fastqs <- c("${fastq.join('","')}")
    learnErrorRates(fastqs,"${pairId}")
    """
}

/*
 *
 * Dada2 main denoising algorithm (in scripts/util.R)
 *
 */

process Denoise {
    tag { "Denoising.${pairId}" }
    publishDir params.outdir+"3-Denoising", mode: "copy", pattern: "*.RDS"
    label "medium_computation"
    label "r_script"
    
    input:
        set val(pairId), file(err), file(fastq) from ERROR_MODEL.join(FASTQ_TRIMMED)
    output:
        file("*.RDS") into DADA_RDS

    script:
    """
    #!/usr/bin/env Rscript
    source("${script_dir}/util.R")

    errors <- c("${err.join('","')}")
    fastqs <- c("${fastq.join('","')}")

    dadaDenoise(errors[1], fastqs[1], "${pairId}_R1")

    if ("${params.singleEnd}" == "false") {
        dadaDenoise(errors[2], fastqs[2], "${pairId}_R2")
    }
    """
}

/*
 *
 * Dada2 main reads merging algorithm (in scripts/util.R)
 *
 */

process Esv {
    tag { "Esv" }

    publishDir params.outdir+"4-ESV", mode: "copy"
    
    label "high_computation"
    label "r_script"
    
    input:
	file dadas from DADA_RDS.collect()
    output:
        file("all.esv.{count_table,fasta}")  into DEREP_CONTIGS
        file("*.RDS")
        file("count_summary.tsv") into COUNT_SUMMARIES
    
    script:
    """
    #!/usr/bin/env Rscript
    source("${script_dir}/util.R")

    esvTable(${params.minOverlap},${params.maxMismatch},"${params.singleEnd}")       
    """
}

/*
 *
 * Mothur MSA and filtering (in scripts/mothur.sh)
 * Parameters: {params.criteria}
 *
 */

process MultipleSequenceAlignment {
    tag { "MSA" }
    publishDir params.outdir+"5-MultipleSequenceAlignment", mode: "copy"
    label "high_computation"
    label "mothur_script"
    
    input:
        set file(count), file(fasta) from DEREP_CONTIGS
        file refAln from Channel.fromPath(params.referenceAln)
    output:
        file("all_MSA.{count_table,fasta}") into DEREP_CONTIGS_ALN
        file("all_MSA.count_table") into MSA_TO_COUNT
        file("*.summary") optional true
    
    script:
    """
    #!/usr/bin/env bash

    ${script_dir}/mothur.sh \
       --step=MSA \
       --refAln=${refAln} \
       --criteria=${params.criteria} \
       --optimize=start-end
    """
}

/*
 *
 * Chimera removal using mothur's unoise algorithm (in scripts/mothur.sh)
 *
 */

process ChimeraRemoval {
    tag { "chimeraRemoval" }
    publishDir params.outdir+"6-ChimeraRemoval", mode: "copy", pattern: "*.{fasta,count_table}"
    label "high_computation"
    label "mothur_script"
    
    input:
	set file(fasta), file(count) from DEREP_CONTIGS_ALN
    output:
        file("all_chimera.{fasta,count_table}") into (NO_CHIMERA_FASTA, FASTA_FOR_REPR)
        file("all_chimera.count_table") into CHIMERA_TO_COUNT
    script:
    """
    ${script_dir}/mothur.sh --step=chimera
    """
}

process PreClassification {
    tag { "preClassification" }
    publishDir params.outdir+"7-PreClassification", mode: "copy", pattern: "*.taxonomy"
    label "high_computation"
    label "mothur_script"
    
    input:
	set file(count), file(fasta) from NO_CHIMERA_FASTA
        file refAln from Channel.fromPath(params.referenceAln)
        file refTax from Channel.fromPath(params.referenceTax)
    output:
        set file(count), file(fasta), file("all_preClassification.taxonomy") into PRE_CLASSIFIED_CONTIGS

    script:
    """
    ${script_dir}/mothur.sh \
        --step=preClassification \
 	--refAln=${refAln} \
 	--refTax=${refTax}
    """
}

/*
 *
 * Clustering with VSEARCH, dgc method (in scripts/mothur.sh)
 *
 */

process Clustering {
    tag { "clustering.${idThreshold}" }
    publishDir params.outdir+"8-Clustering", mode: "copy", pattern: "all_clustering*.{fasta,shared,list}"
    label "high_computation"
    label "mothur_script"
    
    input:
	set file(count), file(fasta), file(tax) from PRE_CLASSIFIED_CONTIGS
        each idThreshold from clusteringThresholds
    output:
        set val(idThreshold), file(count), file(tax), file("all_clustering_*.list"), file("all_clustering_*.shared") into CONTIGS_FOR_CLASSIFICATION
        file("all_clustering_*.shared") into CLUSTERING_TO_COUNT
    script:
    """
    ${script_dir}/mothur.sh --step=clustering --idThreshold=${idThreshold}
    """
}

/*
 *
 * Consensus classification (in scripts/mothur.sh)
 *
 */

process ConsensusClassification {
    tag { "consensusClassification.${idThreshold}" }
    publishDir params.outdir+"9-ConsensusClassification", mode: "copy", pattern: "all_consensusClassification*.{summary,taxonomy}"
    label "medium_computation"    
    label "mothur_script"
    
    input:
	set val(idThreshold), file(count), file(tax), file(list), file(shared) from CONTIGS_FOR_CLASSIFICATION
    output:
	file("all_consensusClassification_*.summary") into CLASSIFICATION_SUMMARY
        set val(idThreshold), file("all_consensusClassification_*.taxonomy"), file(list), file(shared) into CONSTAXONOMY_CONTIGS
    script:
    """
    ${script_dir}/mothur.sh --step=consensusClassification --idThreshold=${idThreshold}
    """
}

/*
** Taxa Filter
*/

process TaxaFilter {
    tag { "convertToMothur.${idThreshold}" }
    publishDir params.outdir+"10-TaxaFilter", mode: "copy"
    errorStrategy "${params.errorsHandling}"
    label "medium_computation"    
    label "mothur_script"
    
    input:
	set val(idThreshold), file(tax), file(list), file(shared) from CONSTAXONOMY_CONTIGS
    output:
        set val(idThreshold), file("all_taxaFilter*.taxonomy"), file("all_taxaFilter*.list"), file("all_taxaFilter*.shared") \
    into TAXA_FILTERED
        file("all_taxaFilter*.shared") into TAXA_FILTER_TO_COUNT
    script:
    """
    ${script_dir}/mothur.sh \
	--step=taxaFilter \
        --idThreshold=${idThreshold} \
        --taxaToFilter='${params.taxaToFilter}' \
	--refAln=${params.referenceAln} \
	--refTax=${params.referenceTax}
    """
}

process OtuRepresentative {
    tag { "OtuRepresentative.${idThreshold}" }
    errorStrategy "${params.errorsHandling}"
    label "medium_computation"    
    label "mothur_script"
    
    input:
	set val(idThreshold), file(tax), file(list), file(shared), file(fasta), file(count) \
        from TAXA_FILTERED.combine(FASTA_FOR_REPR)
    output:
        set val(idThreshold), file("*.fasta"), file(tax), file(shared) into FOR_SUBSAMPLING
        set val(idThreshold), file(shared) into SUBSAMPLING_EST    
    script:
    """
    ${script_dir}/mothur.sh --step=otuRepr --idThreshold=${idThreshold}
    """    
}

process GetSubsamlingValue {
    errorStrategy "${params.errorsHandling}"
    label "low_computation"    
    label "python_script"
    
    input:
	set val(idThreshold), file(shared) from SUBSAMPLING_EST
    output:
	set val(idThreshold), stdout into SUBSAMPLING_THRESHOLDS
    
    script:
    """
    #!/usr/bin/env python3
     
    from util import getSubsamplingThreshold

    getSubsamplingThreshold("${shared}",${params.subsamplingQuantile},${params.minSubsampling})   
    """    
}


(SUBSAMPLING_IN, ALT_CHANNEL) = ( params.skipSubsampling
				 ? [Channel.empty(), FOR_SUBSAMPLING]
				 : [FOR_SUBSAMPLING, Channel.empty()] )

/*
 *
 * Subsampling the samples to the 10th percentile or 1k/sample if the 10th pct is below 1k (in scripts/mothur.sh)
 *
 */

process Subsampling {
    tag { "subsampling" }
    publishDir params.outdir+"11-Subsampling", mode: "copy"
    label "medium_computation"
    label "mothur_script"
    
    input:
	set val(idThreshold), file(fasta), file(tax), file(shared), val(subSampThresh) \
        from SUBSAMPLING_IN.join(SUBSAMPLING_THRESHOLDS)

    output:
	set val(idThreshold), file("all_subsampling*.fasta"), file("all_subsampling*.taxonomy"), file("all_subsampling*.shared") \
    into SUBSAMPLED_OUT
        file("all_subsampling*.shared") into SUBSAMPLING_TO_COUNT
    
    script:
    """
    #!/usr/bin/env bash

    ${script_dir}/mothur.sh --step=subsampling --idThreshold=${idThreshold} --subsamplingNb=${subSampThresh} 
    """
}

(CONTIGS_FOR_PRELULU,FASTA_TO_FILTER,SUBSAMPLED_TAX,ABUNDANCE_TABLES_FOR_LULU) = SUBSAMPLED_OUT
    .mix(ALT_CHANNEL)
    .separate(4) { x -> [ tuple(x[0],x[1]), tuple(x[0],x[1]), tuple(x[0],x[2]), tuple(x[0],x[3]) ] }

/*
 *
 * Pre-lulu step
 *    - Blast each contig against each other
 *    - Keep top 10 hits with at least 84% similarity and 90% query coverage
 *
 */

process PreLulu {
    tag { "preLulus.${idThreshold}" }
    publishDir params.outdir+"12-Lulu", mode: "copy"
    label "medium_computation"
    label "vsearch"
    
    input:
	set val(idThreshold),file(fasta) from CONTIGS_FOR_PRELULU
    output:
	set val(idThreshold),file("match_list_${idThreshold}.txt") into MATCH_LISTS
    script:
	
    """
    fasta_noGap="contigs_${idThreshold}_nogap.fasta"

    sed '/^>/! s/[\\.-]//g' ${fasta} > \$fasta_noGap

    vsearch --usearch_global \$fasta_noGap \
            --db \$fasta_noGap --self \
            --id 0.${params.min_match} \
            --iddef 1 \
            --userout match_list_${idThreshold}.txt \
            -userfields query+target+id \
            --maxaccepts 0 \
            --query_cov .9 \
            --maxhits 10
    """
}

/*
 *
 * Lulu
 *
 */

process Lulu {
    tag { "Lulu.${idThreshold}" }
    publishDir params.outdir+"12-Lulu", mode: "copy"
    errorStrategy "${params.errorsHandling}"
    label "high_computation"    
    label "r_script"
    
    input:
	set val(idThreshold),file(matchlist),file(table) from MATCH_LISTS.join(ABUNDANCE_TABLES_FOR_LULU)
    output:
        set val(idThreshold), file("lulu_table_${idThreshold}.csv") into TABLE_TO_FILTER
        file("lulu_ids_${idThreshold}.csv") into IDS_LULU
        file("lulu*.log_*") optional true
        file("lulu_table_${idThreshold}.csv") into LULU_TO_COUNT
    script:
	
    """
    #!/usr/bin/env Rscript
    source("${script_dir}/util.R")

    luluCurate("${table}","${matchlist}","${idThreshold}","${params.min_ratio_type}","${params.min_ratio}","${params.min_match}","${params.min_rel_cooccurence}")
    """
}

/*
 *
 * Filter out fasta sequences that LULU merged with the most abundant sequence
 * Remove OTUs with an abundance lower than {minAbundance}
 * Convert the abundance table to .shared format
 *
 */

process SingletonFilter {
    tag { "SingletonFilter.${idThreshold}" }
    publishDir params.outdir+"13-SingletonFilter", mode:"copy", pattern:"*.shared"    
    publishDir params.outdir+"14-Results", mode:"copy", pattern:"*.{fasta,shared,taxonomy}"
    errorStrategy "${params.errorsHandling}"
    label "low_computation"    
    label "python_script"
    
    input:
	set idThreshold, file(fasta), file(abundance), file(tax) from FASTA_TO_FILTER.join(TABLE_TO_FILTER).join(SUBSAMPLED_TAX)
    output:
        set idThreshold, file("*.{fasta,shared,taxonomy}") into MOTHUR_TO_PROCESS
        file("*.shared") into SINGLETON_FILTER_TO_COUNT
    script:
    """
    #!/usr/bin/env python3

    from util import filterIds, filterAbundance, csvToShared

    filterAbundance("${abundance}",minAbundance=${params.minAbundance})
    filterIds("curated_${abundance}","${fasta}","${tax}","${idThreshold}")
    csvToShared("curated_${abundance}","${idThreshold}")
    """
}

/*
 * File summary (in scripts/generate_step_summary.py)
 */
	 
process SummaryFile {
    tag { "SummaryFile" }
    publishDir params.outdir+"14-Results", mode: "copy"
    errorStrategy "${params.errorsHandling}"
    label "low_computation"
    label "python_script"
    
    input:
        file f from COUNT_SUMMARIES
        val(raw_counts) from RAW_COUNTS
        val(filtered_counts) from FILTERED_COUNTS
        file f1 from MSA_TO_COUNT.collect()
        file f2 from CHIMERA_TO_COUNT.collect()
	file f3 from CLUSTERING_TO_COUNT.collect()
	file f4 from SUBSAMPLING_TO_COUNT.collect()
	file f5 from TAXA_FILTER_TO_COUNT.collect()
	file f6 from LULU_TO_COUNT.collect()
	file f7 from SINGLETON_FILTER_TO_COUNT.collect()
    output:
        file("sequences_per_sample_per_step_*.tsv") into STEPS_SUMMARY
    script:
    """
    #!/usr/bin/env python3
    from generate_step_summary import write_summary

    clustering_threshold = ${clusteringThresholds}
    counts = { '0-RawData': {${raw_counts.join(', ')}},
               '1-FilterAndTrim': {${filtered_counts.join(', ')}} }

    steps = [("0-RawData",-1),
             ("1-FilterAndTrim",-1),
             ("2-ErrorModel",None),
             ("3.1-Dereplication",-1),
             ("3.2-Denoising",-1),
             ("4-ESV",-1),
             ("5-MultipleSequenceAlignment","all_MSA.count_table"),
             ("6-ChimeraRemoval","all_chimera.count_table"),
             ("7-PreClassification",None),
             ("8-ConsensusClassification",None),             
             ("9-Clustering","all_clustering_*.shared"),
             ("10-TaxaFilter","all_taxaFilter_*.shared"),
             ("11-Subsampling","all_subsampling_*.shared"),
             ("12-Lulu","lulu_table_*.csv"),
             ("13-SingletonFilter","abundance_table_*.shared")
    ]

    write_summary(steps,counts,clustering_threshold)
    """
}

/*
 *
 * Generates some results with mothur 
 *
 */

process Postprocessing {
    tag { "mothurResults" }
    publishDir params.outdir+"15-Postprocessing", mode: "copy"
    errorStrategy "${params.errorsHandling}"
    label "medium_computation"
    label "mothur_script"
    
    input:
	set val(idThreshold), file(f) from MOTHUR_TO_PROCESS
	// set val(idThreshold), file(fasta), file(shared), file(tax) from MOTHUR_TO_PROCESS
    output:
	set file("*.relabund"), file("*summary"), file("*.tre") into RESULTS
    script:
    """
    ${script_dir}/mothur.sh --step=postprocessing --idThreshold=${idThreshold}
    """    
}


