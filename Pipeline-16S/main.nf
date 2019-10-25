#!/usr/bin/env nextflow

def helpMessage() {
    log.info"""
    ===================================
    16S-rDNA-pipeline
    ===================================
    Usage:
    nextflow run 16S-pipeline --reads "*_R{1,2}.fastq.gz" -profile local
    
    --------------- Mandatory arguments ---------------
    --reads         Path to input data (glob pattern)
    --referenceAln  Path to the SILVA reference database (fasta format)
    --referenceTax  Path to the SILVA reference database (taxonomy file)

    ---------------- Optional arguments ---------------
    -profile        Select a configuration from the conf/ folder. Default is "local"
    --outdir        Path to output directory. Default: "./16S-pipeline_outputs"
    
    --singleEnd   If your data is single end

    [Quality filtering]
    --minReads    Sample with less than [minRead] are discarded. Default: 50
    --truncLen    Where to truncate the forward and reverse reads. Default: "220,190"
    --minLength   Reads short than [minLength] are discarded. Default: 20
    --maxEE       Reads with an expected number of errrors above [maxEE] are discarded. Default: 3
    --truncQ      Read truncation at the first occurence of a base of quality below [truncQ]. Default: 2
    --keepPhix    Keep reads matching phiX genome.

    [Read merging]
    --minOverlap    Minimum overlap between forward and reverse read. Default: 20
    --maxMismatch   Maximum number of mismatches in the overlap between forward and reverse read. Default: 1

    [Contig filtering]
    --criteria      Optimization criteria when aligning sequences against reference database. 
                    Discard any sequence starting after where [criteria]% of the sequences start, or end before [criteria]% of the sequences end. 
                    Default: 95
    --minAlnLen     Minimum alignment length in MSA. Default: 50
    --taxaToFilter  Set of taxa to exclude from the analysis. 
                    Default: "Bacteria;Proteobacteria;Alphaproteobacteria;Rickettsiales;Mitochondria;-Bacteria;Cyanobacteria;Oxyphotobacteria;Chloroplast;-unknown;"
    
    [Subsampling]
    --customSubsamplingLevel  User defined subsampling level. Ignored if <= 0
    --subsamplingQuantile     Automatic subsampling level is at quantile [subsamplingQuantile] of the sample sizes. 
                              Ignored if customSubsamplingLevel or skipSubsampling are set. Default: 0.1
    --minSubsamplingLevel     Minimum subsampling level used if the automatic level falls below [minSubsamplingLevel]. Default: 5000
    --skipSubsampling         Skip the subssampling step

    [OTU clustering]
    --clusteringThresholds   Percent similarity threshold for OTU clustering (100 is no clustering). Default: 100,97

    [Co-occurence pattern correction]
    --min_ratio_type         Function to compare the abundance of a parent OTU and its daughter (default: min)
    --min_ratio              Minimum abundance ratio between parent and daughter (across all samples). Default: 1
    --min_rel_cooccurence    Proportion of the parent samples where daughter occurs (max=1)
   
    [Other]
    --minAbundance    Remove OTUs with a total abundance equal or below [minAbundance]. Default: 2
    """.stripIndent()
}

// Show help message
params.help = false
if (params.help){
    helpMessage()
    exit 0
}

def clusteringThresholds = params.clusteringThresholds.toString().split(',').collect{it as int}

if ( params.singleEnd ) {
    read_path = params.reads.replaceAll("\\{1,2\\}","1")
} else {
    read_path = params.reads
}

/*
 *
 Beginning of the pipeline
 *
 */

Channel
    .fromFilePairs( read_path, size: params.singleEnd ? 1 : 2 )
    .ifEmpty { error "Cannot find any reads matching: ${params.reads}" }
    .map { it -> tuple(it[0].replaceAll("-","_"), it[1]) }
    .into { INPUT_FASTQ ; INPUT_FASTQ_FOR_COUNT }

RAW_COUNTS = INPUT_FASTQ_FOR_COUNT.collect{ ["'${it[0]}': ${it[1][0].countFastq()}"] }

INPUT_FASTQ_PREFILT = INPUT_FASTQ.filter{it[1][0].countFastq()>params.minReads}

/*
 *
 * Dada2's filterAndTrim function (in scripts/util.R)
 * Parameters: {params.maxEE} (expected errors), {params.truncLen}, {params.truncQ}
 *
 */

process FilterAndTrim {
    // Quality filter and trimming using dada2
    tag { "FilterAndTrim.${pairId}" }
    label "low_computation"
    label "r_script"
    publishDir params.outdir+"Misc/1-FilterAndTrim", mode: "copy", pattern: "*.{fastq.gz,png}"

    input:
    set val(pairId), file(fastq) from INPUT_FASTQ_PREFILT

    output:
    set val(pairId), file("${pairId}*_trimmed.fastq.gz") optional true into FASTQ_TRIMMED_RAW
    file("*.png") optional true

    script:
    """
    #!/usr/bin/env Rscript

    source("${params.script_dir}/util.R")

    fastqs <- c("${fastq.join('","')}")
    fwd <- fastqs[1]
    rev <- NULL

    if ( "${params.singleEnd}" == "false" ) {
        rev <- fastqs[2]
    }

    rmphix <- !("${params.keepPhix}" == "true")

    filterReads("${pairId}", fwd, rev=rev, minLen=${params.minLength}, maxEE=${params.maxEE}, truncLen=c(${params.truncLen}),rm.phix=rmphix, truncQ=${params.truncQ})
    """
}

FASTQ_TRIMMED_RAW
	.groupTuple() // For single-end, to convert singleton to list for file
    .filter{ (it[1][0].countFastq() > params.minReads)  }
    .ifEmpty { error "All reads were filtered out. Consider relaxing your filtering parameters" }
    .into{ FASTQ_TRIMMED ; FASTQ_TRIMMED_FOR_MODEL ; FASTQ_TRIMMED_FOR_COUNT }

FILTERED_COUNTS = FASTQ_TRIMMED_FOR_COUNT.collect{ ["'${it[0]}': ${it[1][0].countFastq()}"] }

/*
 *
 * Dada2 Error model (in scripts/util.R)
 *
 */

process LearnErrors {
    tag { "LearnErrors.${pairId}" }
    label "medium_computation"
    label "r_script"
    publishDir params.outdir+"Misc/2-ErrorModel", mode: "copy", pattern: "*.{RDS,png}"

    input:
    set val(pairId), file(fastq) from FASTQ_TRIMMED_FOR_MODEL

    output:
    set val(pairId), file("*.RDS") into ERROR_MODEL
    file "*.png"

    script:
    """
    #!/usr/bin/env Rscript
    source("${params.script_dir}/util.R")

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
    label "medium_computation"
    label "r_script"
    publishDir params.outdir+"Misc/3-Denoising", mode: "copy", pattern: "*.RDS"

    input:
    set val(pairId), file(err), file(fastq) from ERROR_MODEL.join(FASTQ_TRIMMED)

    output:
    file("*.RDS") into DADA_RDS

    script:
    """
    #!/usr/bin/env Rscript
    source("${params.script_dir}/util.R")

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
    label "high_computation"
    label "r_script"
    publishDir params.outdir+"Misc/4-ESV", mode: "copy"

    input:
    file dadas from DADA_RDS.collect()

    output:
    file("all_esv.{count_table,fasta}")  into DEREP_CONTIGS
    file("*.RDS")
    file("count_summary.tsv") into COUNT_SUMMARIES

    script:
    """
    #!/usr/bin/env Rscript
    source("${params.script_dir}/util.R")

    esvTable(${params.minOverlap},${params.maxMismatch},"${params.singleEnd}"=="false")
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
    label "high_computation"
    label "mothur_script"
	label "python_script"
    publishDir params.outdir+"Misc/5-MultipleSequenceAlignment", mode: "copy"

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

    ${params.script_dir}/mothur.sh \
    --step=MSA \
    --refAln=${refAln} \
    --criteria=${params.criteria} \
    --minAlnLen=${params.minAlnLen} \
    --optimize=start-end

    python ${params.script_dir}/patch_mothur.py all_MSA.count_table
    """
}

/*
 *
 * Chimera removal using mothur's unoise algorithm (in scripts/mothur.sh)
 *
 */

process ChimeraRemoval {
    tag { "chimeraRemoval" }
    label "high_computation"
    label "mothur_script"
    publishDir params.outdir+"Misc/6-ChimeraRemoval", mode: "copy", pattern: "*.{fasta,count_table}"
    publishDir params.outdir+"Results/raw", mode: "copy", pattern: "*.fasta"

    input:
    set file(fasta), file(count) from DEREP_CONTIGS_ALN

    output:
    file("all_chimera.{fasta,count_table}") into (NO_CHIMERA_FASTA, FASTA_FOR_REPR)
    file("all_chimera.count_table") into CHIMERA_TO_COUNT

    script:
    """
    ${params.script_dir}/mothur.sh --step=chimera
    """
}

process PreClassification {
    tag { "preClassification" }
    label "high_computation"
    label "mothur_script"
    publishDir params.outdir+"Misc/7-PreClassification", mode: "copy", pattern: "*.taxonomy"

    input:
    set file(count), file(fasta) from NO_CHIMERA_FASTA
    file refAln from Channel.fromPath(params.referenceAln)
    file refTax from Channel.fromPath(params.referenceTax)

    output:
    set file(count), file(fasta), file("all_preClassification.taxonomy") into PRE_CLASSIFIED_CONTIGS

    script:
    """
    ${params.script_dir}/mothur.sh \
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
    label "high_computation"
    label "mothur_script"
    publishDir params.outdir+"Misc/8-Clustering", mode: "copy", pattern: "all_clustering*.shared"
    publishDir params.outdir+"Results/raw", mode: "copy", pattern: "all_clustering*.shared"

    input:
    set file(count), file(fasta), file(tax) from PRE_CLASSIFIED_CONTIGS
    each idThreshold from clusteringThresholds

    output:
    set val(idThreshold), file(count), file(tax), file("all_clustering_*.list"), file("all_clustering_*.shared") into CONTIGS_FOR_CLASSIFICATION
    file("all_clustering_*.shared") into CLUSTERING_TO_COUNT

    script:
    """
    ${params.script_dir}/mothur.sh \
    --step=clustering \
    --idThreshold=${idThreshold}
    """
}

/*
 *
 * Consensus classification (in scripts/mothur.sh)
 *
 */

process ConsensusClassification {
    tag { "consensusClassification.${idThreshold}" }
    label "medium_computation"
    label "mothur_script"
    publishDir params.outdir+"Misc/9-ConsensusClassification", mode: "copy", pattern: "*.{summary,taxonomy}"
    publishDir params.outdir+"Results/raw", mode: "copy", pattern: "*.taxonomy"

    input:
    set val(idThreshold), file(count), file(tax), file(list), file(shared) from CONTIGS_FOR_CLASSIFICATION

    output:
    file("*.summary") into CLASSIFICATION_SUMMARY
    set val(idThreshold), file("*.taxonomy"), file(list), file(shared) into CONSTAXONOMY_CONTIGS

    script:
    """
    ${params.script_dir}/mothur.sh \
    --step=consensusClassification \
    --idThreshold=${idThreshold}
    """
}

/*
** Taxa Filter
*/

process TaxaFilter {
    tag { "taxaFilter.${idThreshold}" }
    label "medium_computation"
    label "mothur_script"
    publishDir params.outdir+"Misc/10-TaxaFilter", mode: "copy"

    input:
    set val(idThreshold), file(tax), file(list), file(shared) from CONSTAXONOMY_CONTIGS

    output:
    file("all_taxaFilter*.shared") into TAXA_FILTER_TO_COUNT
    set val(idThreshold), file("all_taxaFilter*.shared") into FOR_MULTIPLETONS_FILTER
    set val(idThreshold), file("all_taxaFilter*.taxonomy"), file("all_taxaFilter*.list") into TAXA_FILTERED

    script:
    """
    ${params.script_dir}/mothur.sh \
    --step=taxaFilter \
    --idThreshold=${idThreshold} \
    --taxaToFilter='${params.taxaToFilter}' \
    --refAln=${params.referenceAln} \
    --refTax=${params.referenceTax}
    """
}

process MultipletonsFilter {
    tag { "MultipletonsFilter.${idThreshold}" }
    label "medium_computation"
    label "mothur_script"
    publishDir params.outdir+"Misc/11-MultipletonsFilter", mode: "copy"

    input:
    set val(idThreshold), file(shared) from FOR_MULTIPLETONS_FILTER

    output:
    set val(idThreshold), file("all_multipletonsFilter*.shared") into (SUBSAMPLING_EST, SUBSAMPLING_OFF_TO_COUNT, MULTIPLETONS_FILTERED)
    file("all_multipletonsFilter*.shared") into MULTIPLETONS_FILTER_TO_COUNT

    script:
    """
    ${params.script_dir}/mothur.sh \
    --step=multipletonsFilter \
    --idThreshold=${idThreshold} \
    --minAbundance=${params.minAbundance}
    """
}

TAXA_AND_MULTIPLETONS_FILTERED = TAXA_FILTERED.join(MULTIPLETONS_FILTERED).combine(FASTA_FOR_REPR)

process OtuRepresentative {
    tag { "OtuRepresentative.${idThreshold}" }
    label "medium_computation"
    label "mothur_script"

    input:
    set val(idThreshold), file(tax), file(list), file(shared), file(fasta), file(count) from TAXA_AND_MULTIPLETONS_FILTERED

    output:
    set val(idThreshold), file("*.fasta"), file(tax), file(shared) into FOR_SUBSAMPLING

    script:
    """
    ${params.script_dir}/mothur.sh \
    --step=otuRepr \
    --idThreshold=${idThreshold}
    """
}

process GetSubsamlingValue {
    tag "getSubsamplingValue_${idThreshold}"
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

    getSubsamplingThreshold("${shared}", ${params.subsamplingQuantile}, ${params.minSubsamplingLevel}, ${params.customSubsamplingLevel})
    """
}


(SUBSAMPLING_IN, ALT_CHANNEL) = ( params.skipSubsampling
                  ? [Channel.empty(), FOR_SUBSAMPLING ]
                  : [FOR_SUBSAMPLING, Channel.empty()] )
/*
 *
 * Subsampling the samples to the 10th percentile or 1k/sample if the 10th pct is below 1k (in scripts/mothur.sh)
 *
 */

process Subsampling {
    tag { "subsampling" }
    label "medium_computation"
    label "mothur_script"
    publishDir params.outdir+"Misc/12-Subsampling", mode: "copy"

    input:
    set val(idThreshold), file(fasta), file(tax), file(shared), val(subSampThresh) \
    from SUBSAMPLING_IN.join(SUBSAMPLING_THRESHOLDS)

    output:
    set val(idThreshold), file("all_subsampling*.fasta"), file("all_subsampling*.taxonomy"), file("all_subsampling*.shared") into SUBSAMPLED_OUT
    file("all_subsampling*.shared") into SUBSAMPLING_TO_COUNT

    script:
    """
    #!/usr/bin/env bash

    ${params.script_dir}/mothur.sh \
    --step=subsampling \
    --idThreshold=${idThreshold} \
    --subsamplingNb=${subSampThresh}
    """
}

(CONTIGS_FOR_PRELULU,FASTA_TO_FILTER,SUBSAMPLED_TAX,ABUNDANCE_TABLES_FOR_LULU) = SUBSAMPLED_OUT
    .mix(ALT_CHANNEL)
    .separate(4) { x -> [tuple(x[0],x[1]), tuple(x[0],x[1]), tuple(x[0],x[2]), tuple(x[0],x[3])] }

/*
 *
 * Pre-lulu step
 *    - Blast each contig against each other
 *    - Keep top 10 hits with at least 84% similarity and 90% query coverage
 *
 */

process PreLulu {
    tag { "preLulus.${idThreshold}" }
    label "medium_computation"
    label "require_vsearch"
    publishDir params.outdir+"Misc/13-Lulu", mode: "copy"

    input:
    set val(idThreshold),file(fasta) from CONTIGS_FOR_PRELULU

    output:
    set val(idThreshold),file("match_list_${idThreshold}.txt") into MATCH_LISTS

    script:

    """
    fasta_noGap="contigs_${idThreshold}_nogap.fasta"

    sed '/^>/! s/[\\.-]//g' ${fasta} > \$fasta_noGap

    vsearch --usearch_global \$fasta_noGap \
            --threads ${task.cpus} \
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
    label "high_computation"
    label "r_script"
    publishDir params.outdir+"Misc/13-Lulu", mode: "copy"

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
    source("${params.script_dir}/util.R")

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

process RareSeqsFilter {
    tag { "RareSeqsFilter.${idThreshold}" }
    label "low_computation"
    label "python_script"
    publishDir params.outdir+"Misc/14-RareSeqsFilter", mode:"copy", pattern:"*.shared"
    publishDir params.outdir+"Results/main", mode:"copy", pattern:"*.{fasta,shared,taxonomy}"

    input:
    set idThreshold, file(fasta), file(abundance), file(tax) from FASTA_TO_FILTER.join(TABLE_TO_FILTER).join(SUBSAMPLED_TAX)

    output:
    set idThreshold, file("*.fasta"), file("*.shared"), file("*.taxonomy") into (FOR_DATABASE, FOR_PLOT)
    file("*.shared") into RARE_SEQS_FILTER_TO_COUNT

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
 *
 * Filter out fasta sequences that LULU merged with the most abundant sequence
 * Remove OTUs with an abundance lower than {minAbundance}
 * Convert the abundance table to .shared format
 *
 */

process SummaryPlot {
    tag { "SummaryPlot.${idThreshold}" }
    label "medium_computation"
    label "python_script"
    publishDir params.outdir+"Results/figures", mode:"copy", pattern:"*.png"

    input:
    set idThreshold, file(fasta), file(shared), file(taxonomy) from FOR_PLOT

    output:
    file("*.png")

    script:
    """
    python ${params.script_dir}/visualization.py --thresh ${idThreshold}
    """
}


/*
 * File summary (in scripts/generate_step_summary.py)
 */

process SummaryFile {
    tag { "SummaryFile" }
    label "medium_computation"
    label "python_script"
    publishDir params.outdir+"Results", mode: "copy"

    input:
    file f from COUNT_SUMMARIES
    val(raw_counts) from RAW_COUNTS
    val(filtered_counts) from FILTERED_COUNTS
    file f from MSA_TO_COUNT
        .mix(CHIMERA_TO_COUNT)
        .mix(CLUSTERING_TO_COUNT)
        .mix(MULTIPLETONS_FILTER_TO_COUNT)
        .mix(SUBSAMPLING_TO_COUNT)
        .mix(TAXA_FILTER_TO_COUNT)
        .mix(LULU_TO_COUNT)
        .mix(RARE_SEQS_FILTER_TO_COUNT)
        .collect()

    output:
    file("sequences_per_sample_per_step_*.tsv") into STEPS_SUMMARY

    script:
    """
    #!/usr/bin/env python3
    from generate_step_summary import write_summary

    clustering_thresholds = ${clusteringThresholds}
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
             ("11-MultipletonsFilter","all_multipletonsFilter_*.shared"),
             ("12-Subsampling","all_subsampling_*.shared"),
             ("13-Lulu","lulu_table_*.csv"),
             ("14-RareSeqsFilter","abundance_table_*.shared")
    ]

    write_summary(steps,counts,clustering_thresholds)
    """
}

/*
 *
 * Collect all results into mothur database before moving on
 *
 */

process Database {
    tag { "database" }
    label "medium_computation"
    label "mothur_script"
    // publishDir params.outdir+"Results/main", mode: "copy", pattern: "*.db"
    publishDir params.outdir+"Results/postprocessing", mode: "copy", pattern: "*.relabund"

    input:
    set val(idThreshold), file(fasta), file(shared), file(taxa) from FOR_DATABASE
    file f from STEPS_SUMMARY

    output:
    set val(idThreshold), file(fasta), file("*.count_table") into FOR_UNIFRAC
    set val(idThreshold), file(shared) into FOR_ALPHADIV, FOR_BETADIV
    file("*.relabund")
    // file("*.db")

    script:
    """
    ${params.script_dir}/mothur.sh \
    --step=postprocessing \
    --idThreshold=${idThreshold}
    """
}

/*
 *
 * Generates some results with mothur
 *
 */

process UnifracDist {
    tag { "Unifrac_${idThreshold}" }
    label "high_computation"
    label "mothur_script"
    publishDir params.outdir+"Results/postprocessing/unifrac", mode: "copy", pattern: "*.{summary,dist}"

    input:
    set val(idThreshold), file(fasta), file(count) from FOR_UNIFRAC

    output:
    set file("*.summary"), file("*.dist")

    script:
    """
    ${params.script_dir}/mothur.sh \
    --step=unifrac \
    --idThreshold=${idThreshold}
    """
}

process AlphaDiv {
    tag { "alphaDiv_${idThreshold}" }
    label "high_computation"
    label "mothur_script"
    publishDir params.outdir+"Results/postprocessing/alpha_diversity", mode: "copy", pattern: "*.summary"

    input:
    set val(idThreshold), file(shared) from FOR_ALPHADIV

    output:
    file("*.summary")

    script:
    """
    ${params.script_dir}/mothur.sh \
    --step=alphaDiversity \
    --idThreshold=${idThreshold}
    """
}

process BetaDiv {
    tag { "betaDiv_${idThreshold}" }
    label "high_computation"
    label "mothur_script"
    publishDir params.outdir+"Results/postprocessing/beta_diversity", mode: "copy", pattern: "*.summary"

    input:
    set val(idThreshold), file(shared) from FOR_BETADIV

    output:
    file("*.summary")

    script:
    """
    ${params.script_dir}/mothur.sh \
    --step=betaDiversity \
    --idThreshold=${idThreshold}
    """
}
