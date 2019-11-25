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
    label "medium_computation"
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

    filter_reads("${pairId}", fwd, rev=rev, minLen=${params.minLength}, maxEE=${params.maxEE}, truncLen=c(${params.truncLen}),rm.phix=rmphix, truncQ=${params.truncQ})
    """
}

FASTQ_TRIMMED_RAW
	.map{ params.singleEnd ? [it[0], [it[1]]] : it }
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
    learn_error_rates(fastqs,"${pairId}")
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

    dada_denoise(errors[1], fastqs[1], "${pairId}_R1")

    if ("${params.singleEnd}" == "false") {
        dada_denoise(errors[2], fastqs[2], "${pairId}_R2")
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

    esv_table(${params.minOverlap},${params.maxMismatch},"${params.singleEnd}"=="false")
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

    mothur "#make.table(count=all_MSA.count_table, compress=f)"
    mv all_MSA.full.count_table all_MSA.count_table
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
    publishDir params.outdir+"Results/raw/details", mode: "copy", pattern: "raw*.{shared,list,fasta,count_table}"
    publishDir params.outdir+"Results/raw", mode: "copy", pattern: "*.database"	

    input:
    set file(count), file(fasta), file(tax) from PRE_CLASSIFIED_CONTIGS
    each idThreshold from clusteringThresholds

    output:
	set val(idThreshold), file(fasta), file(count), file(tax), file("raw_clustering_*.list") into FOR_TAXA_FILTER
    file("raw_*.summary") into CLUSTERING_TO_COUNT
	file("raw*.{fasta,shared,taxonomy,taxonomy,database}")

    script:
    """
    ${params.script_dir}/mothur.sh \
    --step=clustering \
    --idThreshold=${idThreshold} \
    --prefix=raw

    ${params.script_dir}/mothur.sh \
    --step=consensusClassification \
    --idThreshold=${idThreshold} \
    --prefix=raw

    ${params.script_dir}/mothur.sh \
    --step=abundanceTable \
    --idThreshold=${idThreshold} \
    --prefix=raw

    ${params.script_dir}/mothur.sh \
    --step=otuRepr \
    --idThreshold=${idThreshold} \
    --rename=f \
    --prefix=raw

    mv raw_otuRepr_${idThreshold}.fasta raw_sequences_${idThreshold}.fasta

    ${params.script_dir}/mothur.sh \
    --step=postprocessing \
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
    publishDir params.outdir+"Misc/9-TaxaFilter", mode: "copy"

    input:
    set val(idThreshold), file(fasta), file(count), file(tax), file(list) from FOR_TAXA_FILTER

    output:
    file("all_taxaFilter*.count_table") into TAXA_FILTER_TO_COUNT
    set val(idThreshold), file("all_taxaFilter*.{list,taxonomy,fasta,count_table}") into FOR_MULTIPLETONS_FILTER

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
    publishDir params.outdir+"Misc/10-MultipletonsFilter", mode: "copy"

    input:
    set val(idThreshold), file(f) from FOR_MULTIPLETONS_FILTER

    output:
    set val(idThreshold), file("all_multipletonsFilter*.count_table") into SUBSAMPLING_EST
	set val(idThreshold), file("all_multipletonsFilter*.{count_table,fasta,list,taxonomy}") into FOR_SUBSAMPLING
    file("all_multipletonsFilter*.count_table") into MULTIPLETONS_FILTER_TO_COUNT

    script:
    """
    ${params.script_dir}/mothur.sh \
    --step=multipletonsFilter \
    --idThreshold=${idThreshold} \
    --minAbundance=${params.minAbundance}
    """
}

process GetSubsamlingValue {
    tag "getSubsamplingValue_${idThreshold}"
    label "low_computation"
    label "python_script"

    input:
    set val(idThreshold), file(count) from SUBSAMPLING_EST

    output:
    set val(idThreshold), stdout into SUBSAMPLING_THRESHOLDS

    script:
    """
    #!/usr/bin/env python3

    from util import get_subsampling_threshold

    get_subsampling_threshold("${count}", ${params.subsamplingQuantile}, ${params.minSubsamplingLevel}, ${params.customSubsamplingLevel})
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
    publishDir params.outdir+"Misc/11-Subsampling", mode: "copy"

    input:
    set val(idThreshold), file(f), val(subSampThresh) from SUBSAMPLING_IN.join(SUBSAMPLING_THRESHOLDS)

    output:
    set val(idThreshold), file("all_subsampling*.{count_table,fasta,list,taxonomy}") into SUBSAMPLED_OUT
    file("all_subsampling*.count_table") into SUBSAMPLING_TO_COUNT

    script:
    """
    #!/usr/bin/env bash

    ${params.script_dir}/mothur.sh \
    --step=subsampling \
    --idThreshold=${idThreshold} \
    --subsamplingNb=${subSampThresh}
    """
}

SUBSAMPLED_OUT.mix(ALT_CHANNEL).into{FOR_PRELULU ; SUBSAMPLED_ALL}
SUBSAMPLED_ALL.map{it -> [it[0], it[1][0], it[1][1], it[1][3]]}.set{SUBSAMPLED_NO_LIST}

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
    publishDir params.outdir+"Misc/12-Lulu", mode: "copy"

    input:
    set val(idThreshold), file(f) from FOR_PRELULU

    output:
    set val(idThreshold), file("match_list_*.txt"), file("all_abundanceTable_*.shared"), file("${f[2]}") into FOR_LULU

    script:
    """
    ${params.script_dir}/mothur.sh --step=abundanceTable --idThreshold=${idThreshold}
    ${params.script_dir}/mothur.sh --step=otuRepr --idThreshold=${idThreshold} --rename=t

    sed -i '/^>/! s/[\\.-]//g' all_otuRepr_${idThreshold}.fasta

    vsearch --usearch_global all_otuRepr_${idThreshold}.fasta \
            --threads ${task.cpus} \
            --self --db all_otuRepr_${idThreshold}.fasta \
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
    publishDir params.outdir+"Misc/12-Lulu", mode: "copy"

    input:
    set val(idThreshold), file(matchlist), file(table), file(list) from FOR_LULU

    output:
    set val(idThreshold), file("all_lulu*.list") into MERGED_OTUS
	file("all_lulu_*.csv") into LULU_TO_COUNT

    script:
    """
    #!/usr/bin/env Rscript
    source("${params.script_dir}/util.R")

    lulu_curate("${table}","${matchlist}","${idThreshold}",
                "${params.min_ratio_type}","${params.min_ratio}","${params.min_match}","${params.min_rel_cooccurence}",
                cores=${task.cpus})
    merge_otu_list("${list}", "mapping_discarded_${idThreshold}.txt", ${idThreshold}, cores=${task.cpus})
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
    label "mothur_script"
    publishDir params.outdir+"Misc/13-RareSeqsFilter", mode:"copy", pattern:"*.shared"

    input:
	set idThreshold, file(count), file(fasta), file(taxonomy), file(list) from SUBSAMPLED_NO_LIST.join(MERGED_OTUS) 

    output:
	set idThreshold, file("all_rareSeqFilter*.{count_table,fasta,list,taxonomy}") into FOR_DB
	file("*.count_table") into RARE_SEQS_FILTER_TO_COUNT

    script:
    """
    ${params.script_dir}/mothur.sh \
    --step=multipletonsFilter \
    --idThreshold=${idThreshold} \
    --minAbundance=${params.minAbundance}

    for f in \$(ls all_multipletonsFilter*); do
        mv \$f \${f/multipletonsFilter/rareSeqFilter}
    done
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
    publishDir params.outdir+"Results/main/details", mode: "copy", pattern: "*.{taxonomy,shared,fasta}"
    publishDir params.outdir+"Results/main", mode: "copy", pattern: "*.database"	
    publishDir params.outdir+"Results/postprocessing", mode: "copy", pattern: "*.relabund"

    input:
    set val(idThreshold), file(f) from FOR_DB

    output:
    set val(idThreshold), file("sequences_*.fasta"), file("abundance_table_*.shared") into FOR_TREE
    set val(idThreshold), file("abundance_table_*.shared") into FOR_ALPHADIV, FOR_BETADIV
    set val(idThreshold), file("*.database") into FOR_PLOT
    file("*.relabund")

    script:
    """
    ${params.script_dir}/mothur.sh \
    --step=consensusClassification \
    --idThreshold=${idThreshold}

    mv all_consensusClassification_${idThreshold}.taxonomy annotations_${idThreshold}.taxonomy

    ${params.script_dir}/mothur.sh \
    --step=abundanceTable \
    --idThreshold=${idThreshold}

    mv all_abundanceTable_${idThreshold}.shared abundance_table_${idThreshold}.shared

    ${params.script_dir}/mothur.sh \
    --step=otuRepr \
    --idThreshold=${idThreshold} \
    --rename=f

    mv all_otuRepr_${idThreshold}.fasta sequences_${idThreshold}.fasta

    ${params.script_dir}/mothur.sh \
    --step=postprocessing \
    --idThreshold=${idThreshold}
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
    publishDir params.outdir+"Results/figures", mode:"copy", pattern:"*.pdf"

    input:
    set idThreshold, file(db) from FOR_PLOT

    output:
    file("*.pdf")

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
             ("8-Clustering","raw_abundanceTable_*.groups.summary"),
             ("9-TaxaFilter","all_taxaFilter_*.count_table"),
             ("10-MultipletonsFilter","all_multipletonsFilter_*.count_table"),
             ("11-Subsampling","all_subsampling_*.count_table"),
             ("12-Lulu","all_lulu_*.csv"),
             ("13-RareSeqsFilter","all_rareSeqFilter_*.count_table")
    ]

    write_summary(steps,counts,clustering_thresholds)
    """
}
/*
 *
 * Generates some results with mothur
 *
 */

process FastTree {
    tag { "FastTree_${idThreshold}" }
    label "high_computation"
    // label "mothur_script"
    publishDir params.outdir+"Results/postprocessing/unifrac", mode: "copy", pattern: "*.tre"

    input:
    set val(idThreshold), file(fasta), file(shared) from FOR_TREE

    output:
    set val(idThreshold), file(shared), file("*.tre") into CLEARCUT_TREE

    script:
    """
    sed -r 's/.*(Otu[0-9]+)\\|.*/\\>\\1/' ${fasta} > relabeled_${fasta}
    FastTree relabeled_${fasta} > FastTree_${idThreshold}.tre
    """
}

process UnifracDist {
    tag { "Unifrac_${idThreshold}" }
    label "high_computation"
    label "r_script"
    publishDir params.outdir+"Results/postprocessing/unifrac", mode: "copy", pattern: "*.csv"

    input:
    set val(idThreshold), file(shared), file(tree) from CLEARCUT_TREE
	each mode from ('weighted', 'unweighted')

    output:
    file("unifrac*.csv")

    script:
    """
    #!/usr/bin/env Rscript
    source("${params.script_dir}/util.R")
 
    calculate_unifrac("${tree}", "${shared}", method="${mode}", otu_thresh=${idThreshold}, cores=${task.cpus})
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
