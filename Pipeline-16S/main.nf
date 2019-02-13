#!/usr/bin/env nextflow

def helpMessage() {
    log.info"""
    ===================================
     16S-rDNA-pipeline
    ===================================
    Usage:
    nextflow run 16S-pipeline --reads '*_R{1,2}.fastq.gz' --reference 'uniteDB_01-12-2017.fasta' -profile manoa      
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

params.script_dir = workflow.projectDir+"/scripts"

// Header log info
def summary = [:]
summary['Run Name'] = workflow.runName
summary['Reads'] = params.reads
summary['Output dir'] = params.outdir
summary['Script dir'] = params.script_dir
summary['Working dir'] = workflow.workDir
summary['Current home'] = "$HOME"
summary['Current user'] = "$USER"
summary['Current path'] = "$PWD"
summary['Config Profile'] = workflow.profile
log.info summary.collect { k,v -> "${k.padRight(15)}: $v" }.join("\n")
log.info "========================================="

Channel
    .fromFilePairs( params.reads, size: params.revRead ? 2 : 1 )
    .ifEmpty { error "Cannot find any reads matching: ${params.reads}" }
    .into { INPUT_FASTQ ; INPUT_FASTQ_TO_QC }

/*
 *
 * Dada2's filterAndTrim function (in scripts/util.R)
 * Parameters: {params.maxEE} (expected errors), {params.truncLen}, {params.truncQ}
 *
 */

process FilterAndTrim {
    // Quality filter and trimming using dada2. 
    tag { "FilterAndTrim.${pairId}" }
    publishDir "${params.outdir}/1-filterAndTrim", mode: "copy"
    
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

    if ( ${params.revRead} == 1 ) {
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
    publishDir "${params.outdir}/2-errorModel", mode: "copy"
    label "medium_computation"

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
    publishDir "${params.outdir}/3-denoising", mode: "copy"
    label "medium_computation"
    
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

    if (${params.revRead} == 1) {
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

    publishDir "${params.outdir}/4-esv", mode: "copy"
    
    label "high_computation"
    
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

    esvTable(${params.minOverlap},${params.maxMismatch},${params.revRead})       
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
    publishDir "${params.outdir}/5-multipleSequenceAlignment", mode: "copy"
    label "high_computation"
    
    input:
        set file(count), file(fasta) from DEREP_CONTIGS
    output:
        file("all_MSA.{fasta,count_table}") into DEREP_CONTIGS_ALN
        file("*.summary")
    
    script:
    """
    ${params.script_dir}/mothur.sh \
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
    publishDir "${params.outdir}/6-chimeraRemoval", mode: "copy"
    label "high_computation"
    
    input:
	set file(count), file(fasta) from DEREP_CONTIGS_ALN

    output:
        file("all_chimera.{fasta,count_table}") into NO_CHIMERA_FASTA
    
    script:
    """
    ${params.script_dir}/mothur.sh --step=chimera
    """
}

/*
 *
 * Subsampling the samples to the 10th percentile or 1k/sample if the 10th pct is below 1k (in scripts/mothur.sh)
 *
 */

process Subsampling {
    tag { "subsampling" }
    publishDir "${params.outdir}/7-subsampling", mode: "copy"
    
    input:
	set file(fasta), file(count) from NO_CHIMERA_FASTA
    output:
	file("all_subsampling.{fasta,count_table}") into SUBSAMPLED_CONTIGS

    script:
    """

    awk '{for (i=3;i<=NF;i++) sum[i]+=\$i;}; END{for (i in sum) print sum[i]}' ${count} | tail -n +2 |sort -n > sample_size.txt

    percentile_value=`awk '{all[NR] = \$1} END{print all[int(NR*${params.subsamplingQuantile})]}' sample_size.txt`

    if [ \$percentile_value -lt ${params.minSubsampling} ] || [ -z \$percentile_value ]
        then percentile_value=${params.minSubsampling}
    fi

    ${params.script_dir}/mothur.sh --step=subsampling --subsamplingNb=\$percentile_value
    """
}

/*
 *
 * Clustering with VSEARCH, dgc method (in scripts/mothur.sh)
 *
 */

process Clustering {
    tag { "clustering.${idThreshold}" }
    publishDir "${params.outdir}/8-clustering", mode: "copy", pattern: "*.{fasta,shared,list}"
    label "high_computation"
    
    input:
	set file(fasta), file(count) from SUBSAMPLED_CONTIGS
        each idThreshold from (0,0.03)
    output:
        set val(idThreshold), file("all_clustering_*.fasta") into PRELULU_FASTA, FASTA_TO_FILTER
        set val(idThreshold), file("all_clustering_*.shared") into ABUNDANCE_TABLES
        set val(idThreshold), file(fasta), file("all_clustering_*.list"), file(count) into CONTIGS_FOR_CLASSIFICATION
    script:
    """
    ${params.script_dir}/mothur.sh --step=clustering --idThreshold=${idThreshold}
    """
}

/*
 *
 * Consensus classification (in scripts/mothur.sh)
 *
 */

process ConsensusClassification {
    tag { "consensusClassification.${idThreshold}" }
    publishDir "${params.outdir}/9-consensusClassification", mode: "copy"
    
    input:
	set val(idThreshold), file(fasta), file(list), file(count) from CONTIGS_FOR_CLASSIFICATION
    output:
	file("all_consensusClassification_*.summary") into CLASSIFICATION_SUMMARY
	set val(idThreshold), file("all_consensusClassification_*.taxonomy") into CONSENSUS_TAXONOMY
    script:
    """
    ${params.script_dir}/mothur.sh --step=consensusClassification \
        --idThreshold=${idThreshold} \
 	--refAln=${params.referenceAln} \
 	--refTax=${params.referenceTax}
    """
}

/*
 *
 * Pre-lulu step
 *    - Blast each contig against each other
 *    - Keep top 10 hits with at least 84% similarity and 90% query coverage
 *
 */

process PreLulu {
    tag { "preLulus.${idThreshold}" }
    publishDir "${params.outdir}/10-lulu", mode: "copy"

    input:
	set val(idThreshold),file(fasta) from PRELULU_FASTA
    output:
	set val(idThreshold),file("match_list_${idThreshold}.txt") into MATCH_LISTS
    script:
	
    """
    fasta_clean="contigs_${idThreshold}_no_gap.fasta"
    sed 's/[\\.]//g' ${fasta} > \$fasta_clean

    vsearch --usearch_global \$fasta_clean \
            --db \$fasta_clean --self \
            --id .84 \
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
 * Lulu (in scripts/util.R)
 *
 */

process Lulu {
    tag { "Lulu.${idThreshold}" }
    publishDir "${params.outdir}/10-lulu", mode: "copy"
    errorStrategy "${params.errorsHandling}"

    input:
	set val(idThreshold),file(matchlist),file(table) from MATCH_LISTS.join(ABUNDANCE_TABLES)
    output:
        set val(idThreshold), file("lulu_table_${idThreshold}.csv") into LULU_TO_FILTER
    script:
	
    """
    #!/usr/bin/env Rscript
    source("${workflow.projectDir}/scripts/util.R")

    luluCurate("${table}","${matchlist}","${idThreshold}")
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
    publishDir "${params.outdir}/11-singletonFilter", mode: "copy"
    errorStrategy "${params.errorsHandling}"
    
    input:
	set idThreshold,file(abundance),file(tax),file(fasta) from LULU_TO_FILTER.join(CONSENSUS_TAXONOMY).join(FASTA_TO_FILTER)
    output:
	// set idThreshold, file("OTU_${idThreshold}.fasta"), file("abundance_${idThreshold}.shared") into MOTHUR_INPUTS
        set idThreshold, file("OTU_${idThreshold}.fasta"), file("abundance_${idThreshold}.shared"), file("consensus_${idThreshold}.taxonomy") into FOR_TAXA_FILTER
    script:
    """
    #!/usr/bin/env python

    import sys
    sys.path.append("${workflow.projectDir}/scripts")

    from util import filterIds, filterAbundance, csvToShared

    filterAbundance("${abundance}",minAbundance=${params.minAbundance})
    filterIds("curated_${abundance}","${fasta}","${tax}","${idThreshold}")
    csvToShared("curated_${abundance}","${idThreshold}")
    """
}

/*
 *
 * Removing contigs that match the ${params.taxaToFilter} list (in scripts/mothur.sh)
 * Default: unknown
 *
 */

process TaxaFilter {
    tag { "convertToMothur.${idThreshold}" }
    publishDir "${params.outdir}/12-taxaFilter", mode: "copy", pattern: "*.{fasta,shared,taxonomy}"
    errorStrategy "${params.errorsHandling}"
    
    input:
	set val(idThreshold), file(fasta), file(shared), file(tax) from FOR_TAXA_FILTER
    output:
	file("all_taxaFilter_${idThreshold}.{fasta,shared,taxonomy}") into MOTHUR_TO_PROCESS
    script:
    """
    ${params.script_dir}/mothur.sh \
	--step=taxaFilter \
        --idThreshold=${idThreshold} \
        --taxaToFilter='${params.taxaToFilter}' \
	--refAln=${params.referenceAln} \
	--refTax=${params.referenceTax}

    python ${params.script_dir}/filter_shared.py ${shared} ${idThreshold}
    """
}

/*
 *
 * Generates some results with mothur 
 *
 */

process Results {
    tag { "mothurResults.${idThreshold}" }
    publishDir "${params.outdir}/13-Postprocessing", mode: "copy"
    errorStrategy "${params.errorsHandling}"
    
    input:
	set file(fasta), file(shared), file(tax) from MOTHUR_TO_PROCESS
    output:
	set file("*.relabund"), file("*.wsummary"), file("*.tre") into RESULTS
    script:
    """
    mothur "#get.relabund(shared=${shared})"
    mothur "#clearcut(fasta=${fasta}, DNA=T) ; count.seqs(shared=${shared}) ; unifrac.weighted(tree=current,count=current)"
    """    
}

/*
 *
 * File summary (in scripts/generate_step_summary.py)
 *
 */

process SummaryFile {
    tag { "mothurResults" }
    publishDir "${params.outdir}/13-Postprocessing", mode: "copy"
    errorStrategy "${params.errorsHandling}"
    
    input:
	file f1 from COUNT_SUMMARIES
	file f2 from RESULTS.collect()
    output:
        file("sequences_per_sample_per_step.tsv") into STEPS_SUMMARY
    script:
    """
    #!/usr/bin/env python
    import sys
    sys.path.append("${workflow.projectDir}/scripts")
    from generate_step_summary import write_summary

    write_summary("${params.outdir}","${params.reads}")
    """

}
