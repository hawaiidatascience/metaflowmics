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

script_dir = workflow.projectDir+"/scripts"

// Add the step number
output_dirs = params.results_folders.indexed().collect { idx,name -> "${params.outdir}/${1+idx}-${name}" }

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
    .set { INPUT_FASTQ }

/*
 *
 * Dada2's filterAndTrim function (in scripts/util.R)
 * Parameters: {params.maxEE} (expected errors), {params.truncLen}, {params.truncQ}
 *
 */

process FilterAndTrim {
    // Quality filter and trimming using dada2. 
    tag { "FilterAndTrim.${pairId}" }
    publishDir output_dirs[0], mode: "copy", pattern: "*.{fastq,png}"
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

    source("${workflow.projectDir}/scripts/util.R")  

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
    .into{ FASTQ_TRIMMED ; FASTQ_TRIMMED_FOR_MODEL }

/*
 *
 * Dada2 Error model (in scripts/util.R)
 *
 */

process LearnErrors {
    tag { "LearnErrors.${pairId}" }
    publishDir output_dirs[1], mode: "copy", pattern: "*.{RDS,png}"
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
    source("${workflow.projectDir}/scripts/util.R") 

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
    publishDir output_dirs[2], mode: "copy", pattern: "*.RDS"
    label "medium_computation"
    label "r_script"
    
    input:
        set val(pairId), file(err), file(fastq) from ERROR_MODEL.join(FASTQ_TRIMMED)
    output:
        file("*.RDS") into DADA_RDS

    script:
    """
    #!/usr/bin/env Rscript
    source("${workflow.projectDir}/scripts/util.R")

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

    publishDir output_dirs[3], mode: "copy"
    
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
    source("${workflow.projectDir}/scripts/util.R")

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
    publishDir output_dirs[4], mode: "copy"
    label "high_computation"
    label "mothur_script"
    
    input:
        set file(count), file(fasta) from DEREP_CONTIGS
    output:
        file("all_MSA.{count_table,fasta}") into DEREP_CONTIGS_ALN
        file("*.summary") optional true
    
    script:
    """
    #!/usr/bin/env bash

    ${script_dir}/mothur.sh \
       --step=MSA \
       --refAln=${params.referenceAln} \
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
    publishDir output_dirs[5], mode: "copy", pattern: "*.{fasta,count_table}"
    label "high_computation"
    label "mothur_script"
    
    input:
	set file(fasta), file(count) from DEREP_CONTIGS_ALN

    output:
        file("all_chimera.{fasta,count_table}") into (NO_CHIMERA_FASTA, FASTA_FOR_REPR)
    
    script:
    """
    ${script_dir}/mothur.sh --step=chimera
    """
}

process PreClassification {
    tag { "preClassification" }
    publishDir output_dirs[6], mode: "copy", pattern: "*.taxonomy"
    label "high_computation"
    label "mothur_script"
    
    input:
	set file(count), file(fasta) from NO_CHIMERA_FASTA

    output:
        set file(count), file(fasta), file("all_preClassification.taxonomy") into PRE_CLASSIFIED_CONTIGS
    
    script:
    """
    ${script_dir}/mothur.sh \
        --step=preClassification \
 	--refAln=${params.referenceAln} \
 	--refTax=${params.referenceTax}
    """
}

/*
 *
 * Clustering with VSEARCH, dgc method (in scripts/mothur.sh)
 *
 */

process Clustering {
    tag { "clustering.${idThreshold}" }
    publishDir output_dirs[7], mode: "copy", pattern: "all_clustering*.{fasta,shared,list}"
    label "high_computation"
    label "mothur_script"
    
    input:
	set file(count), file(fasta), file(tax) from PRE_CLASSIFIED_CONTIGS
        each idThreshold from clusteringThresholds
    output:
        set val(idThreshold), file(count), file(tax), file("all_clustering_*.list"), file("all_clustering_*.shared") \
    into CONTIGS_FOR_CLASSIFICATION
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
    publishDir output_dirs[8], mode: "copy", pattern: "all_consensusClassification*.{summary,taxonomy}"
    label "medium_computation"    
    label "mothur_script"
    
    input:
	set val(idThreshold), file(count), file(tax), file(list), file(shared) from CONTIGS_FOR_CLASSIFICATION
    output:
	file("all_consensusClassification_*.summary") into CLASSIFICATION_SUMMARY
        set val(idThreshold), file("all_consensusClassification_*.taxonomy"), file(list), file(shared) into CONSTAXONOMY_CONTIGS
    script:
    """
    ${script_dir}/mothur.sh --step=consensusClassification \
        --idThreshold=${idThreshold} \
 	--refAln=${params.referenceAln} \
 	--refTax=${params.referenceTax}
    """
}

/*
Positioning of taxa filtering (here or at the end)
--> Redirect channel depending on the choice
*/

(CLASSIFIED_CONTIGS_1, CLASSIFIED_CONTIGS_2) = ( params.taxaFilterAtEnd
						? [Channel.empty(), CONSTAXONOMY_CONTIGS]
						: [CONSTAXONOMY_CONTIGS, Channel.empty()] )

process TaxaFilterInterm {
    tag { "convertToMothur.${idThreshold}" }
    publishDir output_dirs[9], mode: "copy"
    errorStrategy "${params.errorsHandling}"
    label "medium_computation"    
    label "mothur_script"
    
    input:
	set val(idThreshold), file(tax), file(list), file(shared) from CLASSIFIED_CONTIGS_1
    output:
        set val(idThreshold), file("all_taxaFilter*.taxonomy"), file("all_taxaFilter*.list"), file("all_taxaFilter*.shared") \
        into INTERM_TAXA_FILTER
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

OTU_REPR_IN = INTERM_TAXA_FILTER.mix(CLASSIFIED_CONTIGS_2)

process OtuRepresentative {
    tag { "OtuRepresentative.${idThreshold}" }
    errorStrategy "${params.errorsHandling}"
    label "medium_computation"    
    label "mothur_script"
    
    input:
	set val(idThreshold), file(tax), file(list), file(shared), file(fasta), file(count) \
        from OTU_REPR_IN.combine(FASTA_FOR_REPR)
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
    publishDir output_dirs[10], mode: "copy"
    label "medium_computation"
    label "mothur_script"
    
    input:
	set val(idThreshold), file(fasta), file(tax), file(shared), val(subSampThresh) \
        from SUBSAMPLING_IN.join(SUBSAMPLING_THRESHOLDS)

    output:
	set val(idThreshold), file("all_subsampling*.fasta"), file("all_subsampling*.taxonomy"), file("all_subsampling*.shared") \
        into SUBSAMPLED_OUT
    
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
    publishDir output_dirs[11], mode: "copy"
    label "medium_computation"
    
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
    publishDir output_dirs[11], mode: "copy"
    errorStrategy "${params.errorsHandling}"
    label "high_computation"    
    label "r_script"
    
    input:
	set val(idThreshold),file(matchlist),file(table) from MATCH_LISTS.join(ABUNDANCE_TABLES_FOR_LULU)
    output:
        set val(idThreshold), file("lulu_table_${idThreshold}.csv") into TABLE_TO_FILTER
        file("lulu_ids_${idThreshold}.csv") into IDS_LULU
        file("lulu*.log_*") optional true
    script:
	
    """
    #!/usr/bin/env Rscript
    source("${workflow.projectDir}/scripts/util.R")

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
    publishDir output_dirs[12], mode: "copy"
    errorStrategy "${params.errorsHandling}"
    label "low_computation"    
    label "python_script"
    
    input:
	set idThreshold, file(fasta), file(abundance), file(tax) from FASTA_TO_FILTER.join(TABLE_TO_FILTER).join(SUBSAMPLED_TAX)
    output:
        set idThreshold, file("singleton_filtered_*.{fasta,shared,taxonomy}") into FOR_TAXA_FILTER
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
 * Taxa filter in the end
 * Removing contigs that match the ${params.taxaToFilter} list (in scripts/mothur.sh)
 * Default: unknown
 *
 */

(TAXA_FILT_IN, TAX_FILT_ALT) = ( params.taxaFilterAtEnd
				? [FOR_TAXA_FILTER, Channel.empty()]
				: [Channel.empty(), FOR_TAXA_FILTER] )

process TaxaFilterEnd {
    tag { "taxaFilter.${idThreshold}" }
    publishDir output_dirs[13], mode: "copy"
    errorStrategy "${params.errorsHandling}"
    label "medium_computation"
    label "mothur_script"
    
    input:
	set val(idThreshold), file(f) from TAXA_FILT_IN
    output:
        set val(idThreshold), file("all_taxaFilter_*.{fasta,shared,taxonomy}") into TAXA_FILT_OUT
    
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

TAXA_FILT_OUT
    .mix(TAX_FILT_ALT)
    .into {FINAL_SUMMARY ; MOTHUR_TO_PROCESS}

/*
 * File summary (in scripts/generate_step_summary.py)
 */

process SummaryFile {
    tag { "SumaryFile" }
    publishDir output_dirs[14], mode: "copy"
    errorStrategy "${params.errorsHandling}"
    label "low_computation"
    label "python_script"
    
    input:
	file f1 from COUNT_SUMMARIES
        set val(idThreshold), file(f2) from FINAL_SUMMARY
    output:
        file("sequences_per_sample_per_step_*.tsv") into STEPS_SUMMARY
        file("*.{fasta,shared,taxonomy}") into FINAL_RESULTS
    script:
    """
    #!/usr/bin/env python3

    from shutil import copyfile
    from generate_step_summary import write_summary

    write_summary("${params.outdir}","${params.reads}","${idThreshold}")

    ## Rename files

    for f in ["${f2.join('","')}"]:
        if f.endswith("fasta"):
            copyfile(f,"sequences_${idThreshold}.fasta")
        elif f.endswith("shared"):
            copyfile(f,"abundance_table_${idThreshold}.shared")
        else:
            copyfile(f,"annotations_${idThreshold}.taxonomy")
    """

}

/*
 *
 * Generates some results with mothur 
 *
 */

process Postprocessing {
    tag { "mothurResults" }
    publishDir output_dirs[15], mode: "copy"
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


