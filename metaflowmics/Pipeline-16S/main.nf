#!/usr/bin/env nextflow

def helpMessage() {
    log.info"""
    ===================================
    16S-rDNA-pipeline
    ===================================

    Usage:
    nextflow run 16S-pipeline -profile local --reads "*_R{1,2}.fastq.gz"

    For more detailed information, see https://metagenomics-pipelines.readthedocs.io/en/latest/
    
    ---------------------------------- Mandatory arguments ----------------------------------------

    --reads         Path to input data (glob pattern)
    --referenceAln  Path to the SILVA reference database (fasta format)
    --referenceTax  Path to the SILVA reference database (taxonomy file)

	---------------------------------- Optional arguments ----------------------------------------

    --outdir        Path to output directory. Default: "./16S-pipeline_outputs"
    --singleEnd     If your data is single end

    [Quality filtering]
    --minReads    Sample with less than <minRead> are discarded. Default: 50
    --truncLen    Where to truncate the forward and reverse reads. Default: "220,190"
    --minLength   Reads short than <minLength> are discarded. Default: 20
    --maxEE       Reads with an expected number of errrors above <maxEE> are discarded. Default: 3
    --truncQ      Read truncation at the 1st occurence of a base of quality <= <truncQ>. Default: 2
    --keepPhix    Keep reads matching phiX genome.

    [Read merging]
    --minOverlap    Minimum overlap between forward and reverse read. Default: 20
    --maxMismatch   Maximum number of mismatches in the overlap between forward and reverse read. 
                    Default: 1

    [Contig filtering]
    --criteria      Optimization criteria when aligning sequences against reference database. 
                    Discard any sequence starting after where <criteria>% of the sequences 
                    start or end before <criteria>% of the sequences end. 
                    Default: 95
    --minAlnLen     Minimum alignment length in MSA. Default: 50
    --taxaToFilter  Set of taxa to exclude from the analysis. 
                    Default: "Bacteria;Proteobacteria;Alphaproteobacteria;Rickettsiales;Mitochondria;-Bacteria;Cyanobacteria;Oxyphotobacteria;Chloroplast;-unknown;"
    
    [Subsampling]
    --customSubsamplingLevel  User defined subsampling level. Ignored if <= 0
    --subsamplingQuantile     Automatic subsampling level is at quantile <subsamplingQuantile>
                              of the sample sizes. 
                              Ignored if customSubsamplingLevel or skipSubsampling are set.
                              Default: 0.1
    --minSubsamplingLevel     Minimum subsampling level used if the automatic level falls below 
                              this value. Default: 5000
    --skipSubsampling         Skip the subsampling step

    [OTU clustering]
    --clusteringThresholds   Percent similarity threshold for OTU clustering (100 is no clustering)
                             Default: 100,97

    [Co-occurence pattern correction]
	--skipLulu               Skip the Lulu step
    --min_ratio_type         Function to compare the abundance of a parent OTU and its daughter.
                             Default: min
    --min_ratio              Minimum abundance ratio between parent and daughter
                             (across all samples). Default: 1
    --min_rel_cooccurence    Proportion of the parent samples where daughter occurs. Default: 1

    [Other]
    --minAbundance    Remove OTUs with a total abundance equal or below <minAbundance>. Default: 2
	--mothurDb        Compute mothur database summary file. Can be memory intensive. Default: false
	--biom            Compute biom summary file. Default: false
	--skipBetaDiv     Skip beta diversity calculation. Recommended if your memory is limited and
	                  your have a significant amount of OTUs. Default: false
    """.stripIndent()
}

// Show help message
params.help = false
if (params.help){
    helpMessage()
    exit 0
}

def summary = [:]
summary['database aln'] = params.referenceAln
summary['database tax'] = params.referenceTax
summary['single end'] = params.singleEnd
summary['min reads per sample'] = params.minReads
summary['min read length'] = params.minLength
summary['read truncation'] = params.truncLen
summary['max number of expected errors'] = params.maxEE
summary['quality truncation'] = params.truncQ
summary['keep phix genome'] = params.keepPhix
summary['min overlap for merging'] = params.minOverlap
summary['max mismatches for merging'] = params.maxMismatch
summary['Percentile for start/end contig filtering'] = params.criteria
summary['Minimum contig alignment length against db'] = params.minAlnLen
summary['Filtered taxa'] = params.taxaToFilter
summary['Skip subsampling'] = params.skipSubsampling

if (params.customSubsamplingLevel > 0){
	summary['Subsampling'] = params.customSubsamplingLevel
} else if (!params.skipSubsampling){
	summary['Percentile for automatic subsampling'] = params.subsamplingQuantile * 100
	summary['Hard threshold for automatic subsampling'] = params.minSubsamplingLevel
}
summary['clustering similarity thresholds'] = params.clusteringThresholds
summary['Min OTU abundace filter'] = params.minAbundance

summary['Skip LULU'] = params.skipLulu
if (!params.skipLulu){
	summary['Lulu ratio type'] = params.min_ratio_type
	summary['Lulu parent/daughter min similarity'] = params.min_match
	summary['Lulu parent/daughter min abundance ratio'] = params.min_ratio
	summary['Lulu parent/daughter min co-occurrence'] = params.min_rel_cooccurence
}

file(params.outdir).mkdir()
summary_handle = file("${params.outdir}/parameters_summary.log")
summary_handle << summary.collect { k,v -> "${k.padRight(50)}: $v" }.join("\n")

/*
 *
 Beginning of the pipeline
 *
 */

def clusteringThresholds = params.clusteringThresholds.toString().split(',').collect{it as int}

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
    publishDir params.outdir+"/Misc/1-FilterAndTrim", mode: "copy", pattern: "*.{fastq.gz,png}"

    input:
    tuple val(pairId), file(fastq) from INPUT_FASTQ_PREFILT

    output:
    tuple val(pairId), file("${pairId}*_trimmed.fastq.gz") optional true into FASTQ_TRIMMED_RAW
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

    filter_reads("${pairId}", fwd, rev=rev, minLen=c(${params.minLength}), maxEE=c(${params.maxEE}), truncLen=c(${params.truncLen}),rm.phix=rmphix, truncQ=c(${params.truncQ}))
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
    publishDir params.outdir+"/Misc/2-ErrorModel", mode: "copy", pattern: "*.{RDS,png}"

    input:
    tuple val(pairId), file(fastq) from FASTQ_TRIMMED_FOR_MODEL

    output:
    tuple val(pairId), file("*.RDS") into ERROR_MODEL
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
    publishDir params.outdir+"/Misc/3-Denoising", mode: "copy", pattern: "*.RDS"

    input:
    tuple val(pairId), file(err), file(fastq) from ERROR_MODEL.join(FASTQ_TRIMMED)

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
    publishDir params.outdir+"/Misc/4-ESV", mode: "copy"

    input:
    file dadas from DADA_RDS.collect()

    output:
    file("all_esv.{count_table,fasta}")  into DEREP_CONTIGS
    file("count_summary.tsv") into COUNT_SUMMARIES
    file("*.RDS")
    file("{nmatch,nmismatch}_summary.tsv") optional true

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
    publishDir params.outdir+"/Misc/5-MultipleSequenceAlignment", mode: "copy"

    input:
    tuple file(count), file(fasta) from DEREP_CONTIGS
    file(refAln) from Channel.fromPath(params.referenceAln)

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
    publishDir params.outdir+"/Misc/6-ChimeraRemoval", mode: "copy", pattern: "*.{fasta,count_table}"
    publishDir params.outdir+"/Results/raw/details", mode: "copy", pattern: "raw*.fasta"

    input:
    tuple file(fasta), file(count) from DEREP_CONTIGS_ALN

    output:
    file("all_chimera.{fasta,count_table}") into (NO_CHIMERA_FASTA, FASTA_FOR_REPR)
    file("all_chimera.count_table") into CHIMERA_TO_COUNT
    file("raw*.fasta")

    script:
    """
    ${params.script_dir}/mothur.sh --step=chimera
    cp all_chimera*.fasta raw_sequences.fasta
    """
}

process PreClassification {
    tag { "preClassification" }
    label "high_computation"
    label "mothur_script"
    publishDir params.outdir+"/Misc/7-PreClassification", mode: "copy", pattern: "*.taxonomy"

    input:
    tuple file(count), file(fasta) from NO_CHIMERA_FASTA
    file(refAln) from Channel.fromPath(params.referenceAln)
    file(refTax) from Channel.fromPath(params.referenceTax)

    output:
    tuple file(count), file(fasta), file("all_preClassification.taxonomy") into PRE_CLASSIFIED_CONTIGS

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
    tag { "clustering.${otuId}" }
    label "high_computation"
    label "mothur_script"
    publishDir params.outdir+"/Misc/8-Clustering", mode: "copy", pattern: "raw.*"
    publishDir params.outdir+"/Results/raw/details", mode: "copy", pattern: "raw*.{shared,fasta,taxonomy}"
    publishDir params.outdir+"/Results/raw", mode: "copy", pattern: "*.biom"	

    input:
    tuple file(count), file(fasta), file(tax) from PRE_CLASSIFIED_CONTIGS
    each otuId from clusteringThresholds

    output:
	tuple val(otuId), file(fasta), file(count), file(tax), file("raw_clustering_*.list") into FOR_TAXA_FILTER
    file("raw_*.summary") into CLUSTERING_TO_COUNT
	file("raw*.{fasta,shared,taxonomy}")
	file("*.{biom,database}") optional true

    script:
    """
    ${params.script_dir}/mothur.sh \
    --step=clustering \
    --otuId=${otuId} \
    --prefix=raw

    ${params.script_dir}/mothur.sh \
    --step=consensusClassification \
    --otuId=${otuId} \
    --prefix=raw

    ${params.script_dir}/mothur.sh \
    --step=abundanceTable \
    --otuId=${otuId} \
    --prefix=raw

    ${params.script_dir}/mothur.sh \
    --step=otuRepr \
    --otuId=${otuId} \
    --rename=f \
    --prefix=raw

    ${params.script_dir}/mothur.sh \
    --step=relabund \
    --otuId=${otuId}

	[ "${params.biom}" == "true" ] &&
    ${params.script_dir}/mothur.sh \
    --step=biom \
    --otuId=${otuId} || echo "Skipping biom"

	[ "${params.mothurDb}" == "true" ] &&
    ${params.script_dir}/mothur.sh \
    --step=database \
    --otuId=${otuId} || echo "Skipping mothur db"

    """
}

/*
** Taxa Filter
*/

process TaxaFilter {
    tag { "taxaFilter.${otuId}" }
    label "medium_computation"
    label "mothur_script"
    publishDir params.outdir+"/Misc/9-TaxaFilter", mode: "copy"

    input:
    tuple val(otuId), file(fasta), file(count), file(tax), file(list) from FOR_TAXA_FILTER

    output:
    file("all_taxaFilter*.summary") into TAXA_FILTER_TO_COUNT
    tuple val(otuId), file("all_taxaFilter*.{list,taxonomy,fasta,count_table}") into FOR_MULTIPLETONS_FILTER

    script:
    """
    ${params.script_dir}/mothur.sh \
    --step=taxaFilter \
    --otuId=${otuId} \
    --taxaToFilter='${params.taxaToFilter}' \
    --refAln=${params.referenceAln} \
    --refTax=${params.referenceTax}
    """
}

process MultipletonsFilter {
    tag { "MultipletonsFilter.${otuId}" }
    label "medium_computation"
    label "mothur_script"
    publishDir params.outdir+"/Misc/10-MultipletonsFilter", mode: "copy"

    input:
    tuple val(otuId), file(f) from FOR_MULTIPLETONS_FILTER

    output:
    tuple val(otuId), file("all_multipletons*.count_table") into FILTERED_TABLE
	tuple val(otuId), file("all_multipletons*.count_table"), file("all_multipletons*.fasta"), file("all_multipletons*.taxonomy"), file("all_multipletons*.list") into FILTERED_OTUS
    file("all_multipletonsFilter*.summary") into MULTIPLETONS_FILTER_TO_COUNT

    script:
    """
    ${params.script_dir}/mothur.sh \
    --step=multipletonsFilter \
    --otuId=${otuId} \
    --minAbundance=${params.minAbundance}
    """
}

// FILTERED_TABLE.branch{off: (params.customSubsamplingLevel>0) || params.skipSubsampling
// 					  on: true}
// 	.set{TABLE_FOR_SUBSAMPLING}

if (params.customSubsamplingLevel<=0 && !params.skipSubsampling) {
	process GetSubsamlingValue {
		tag "getSubsamplingValue_${otuId}"
		label "low_computation"
		label "python_script"

		input:
		tuple val(otuId), file(count) from FILTERED_TABLE

		output:
		tuple val(otuId), stdout into SUBSAMPLING_THRESHOLDS

		script:
		"""
		#!/usr/bin/env python3

		from util import get_subsampling_threshold
		
		get_subsampling_threshold("${count}", ${params.subsamplingQuantile}, ${params.minSubsamplingLevel}, ${params.customSubsamplingLevel})
		"""
	}
}

// FILTERED_OTU_DATA
// 	.branch {off: params.skipSubsampling
// 			 on: true}
// 	.set { FOR_SUBSAMPLING }

/*
 *
 * Subsampling the samples to the 10th percentile or 1k/sample if the 10th pct is below 1k (in scripts/mothur.sh)
 *
 */

if (!params.skipSubsampling) {
	process Subsampling {
		tag { "subsampling" }
		label "medium_computation"
		label "mothur_script"
		publishDir params.outdir+"/Misc/11-Subsampling", mode: "copy"
		
		input:
		tuple val(otuId), file(count), file(fasta), file(tax), file(list), val(subSampThresh) from FILTERED_OTUS.join(SUBSAMPLING_THRESHOLDS, remainder: true)

		output:
		tuple val(otuId), file("all_subsampling*.count_table"), file("all_subsampling*.fasta"), file("all_subsampling*.taxonomy"), file("all_subsampling*.list") into SUBSAMPLED_DATA
		file("all_subsampling*.summary") into SUBSAMPLING_TO_COUNT

		script:
		"""
		#!/usr/bin/env bash

		[ "${subSampThresh}" == "null" ] && subsampling_level=${params.customSubsamplingLevel} || subsampling_level=${subSampThresh}

		${params.script_dir}/mothur.sh \
			--step=subsampling \
			--otuId=${otuId} \
			--subsamplingNb=\$subsampling_level
		"""
	}
} else {
	FILTERED_OTUS.set{SUBSAMPLED_DATA}
	Channel.empty().set{SUBSAMPLING_TO_COUNT}
}

/*
 *
 * Pre-lulu step
 *    - Blast each contig against each other
 *    - Keep top 10 hits with at least 84% similarity and 90% query coverage
 *
 */

if (!params.skipLulu) {
	process PreLulu {
		tag { "preLulus.${otuId}" }
		label "medium_computation"
		label "mothur_script"
		publishDir params.outdir+"/Misc/12-Lulu", mode: "copy"

		input:
    	tuple val(otuId), file(count), file(fasta), file(tax), file(list) from SUBSAMPLED_DATA

		output:
		tuple val(otuId), file("match_list_*.txt"), file("all_abundanceTable_*.shared"), file(list), file(count), file(fasta), file(tax) into FOR_LULU

		script:
		"""
		${params.script_dir}/mothur.sh --step=abundanceTable --otuId=${otuId}
		${params.script_dir}/mothur.sh --step=otuRepr --otuId=${otuId} --rename=t

		sed -i '/^>/! s/[\\.-]//g' all_otuRepr_${otuId}.fasta

		vsearch --usearch_global all_otuRepr_${otuId}.fasta \
            --threads ${task.cpus} \
            --self --db all_otuRepr_${otuId}.fasta \
            --id 0.${params.min_match} \
            --iddef 1 \
            --userout match_list_${otuId}.txt \
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
		tag { "Lulu.${otuId}" }
		label "high_computation"
		label "r_script"
		publishDir params.outdir+"/Misc/12-Lulu", mode: "copy"
		
		input:
		tuple val(otuId), file(matchlist), file(table), file(list), file(count), file(fasta), file(tax) from FOR_LULU

		output:
 		tuple val(otuId), file(count), file(fasta), file(tax), file("all_lulu*.list") into LULU_FILTERED
    	file("all_lulu_*.csv") into LULU_TO_COUNT

		script:
		"""
		#!/usr/bin/env Rscript
		source("${params.script_dir}/util.R")

		lulu_curate("${table}","${matchlist}","${otuId}",
            "${params.min_ratio_type}","${params.min_ratio}","${params.min_match}","${params.min_rel_cooccurence}",
            cores=${task.cpus})
		merge_otu_list("${list}", "mapping_discarded_${otuId}.txt", ${otuId}, cores=${task.cpus})
		"""
	}
} else {
	SUBSAMPLED_DATA.set{LULU_FILTERED}
	Channel.empty().set{LULU_TO_COUNT}
}

/*
 *
 * Collect all results into mothur database before moving on
 *
 */

process Database {
    tag { "database_${otuId}" }
    label "medium_computation"
    label "mothur_script"
    publishDir params.outdir+"/Results/main/details", mode: "copy", pattern: "*.{taxonomy,shared,fasta,relabund}"
    publishDir params.outdir+"/Results/main", mode: "copy", pattern: "*.biom"	

    input:
    tuple val(otuId), file(count), file(fasta), file(tax), file(list) from LULU_FILTERED

    output:
	tuple val(otuId), file("otu_repr_*.fasta"), file("abundance_*.shared"), file("*.taxonomy") into DB_OUT
    file("all_esv.fasta")
	file("*.relabund")
	file("*.{biom,database}") optional true

    script:
    """
    ${params.script_dir}/mothur.sh \
    --step=consensusClassification \
    --otuId=${otuId}

    mv all_consensusClassification_${otuId}.taxonomy annotations_${otuId}.taxonomy

    ${params.script_dir}/mothur.sh \
    --step=abundanceTable \
    --otuId=${otuId}

    mv all_abundanceTable_${otuId}.shared abundance_table_${otuId}.shared

    ${params.script_dir}/mothur.sh \
    --step=otuRepr \
    --otuId=${otuId} \
    --rename=t

    mv all_otuRepr_${otuId}.fasta otu_repr_${otuId}.fasta

    ${params.script_dir}/mothur.sh \
    --step=relabund \
    --otuId=${otuId}

	[ "${params.biom}" == "true" ] &&
    ${params.script_dir}/mothur.sh \
    --step=biom \
    --otuId=${otuId} || echo "Skipping biom"

	[ "${params.mothurDb}" == "true" ] &&
    ${params.script_dir}/mothur.sh \
    --step=database \
    --otuId=${otuId} || echo "Skipping mothur db"

    cp ${fasta} all_esv.fasta
    """
}

// order: [thresh, fasta, shared, tax]
DB_OUT
	.multiMap{ it ->
	plot: [it[0], it[2], it[3]]
	alpha_div: [it[0], it[2]]
	beta_div: [it[0], it[2]]
	tree: [it[0], it[1], it[2]]
	species_assign: [it[0], it[1], it[3]] }
	.set{MAIN_OUTPUTS}

/*
 *
 * Filter out fasta sequences that LULU merged with the most abundant sequence
 * Remove OTUs with an abundance lower than {minAbundance}
 * Convert the abundance table to .shared format
 *
 */

process SummaryPlot {
    tag { "SummaryPlot.${otuId}" }
    label "medium_computation"
    label "python_script"
    publishDir params.outdir+"/Results/figures", mode:"copy", pattern:"*.html"
	errorStrategy 'ignore'

    input:
    tuple otuId, file(shared), file(tax) from MAIN_OUTPUTS.plot

    output:
    file("*.html")

    script:
    """
    python3 ${params.script_dir}/bokeh_viz.py --shared ${shared} --tax ${tax} --thresh ${otuId}
    """
}


/*
 * File summary (in scripts/generate_step_summary.py)
 */

process SummaryFile {
    tag { "SummaryFile" }
    label "medium_computation"
    label "python_script"
    publishDir params.outdir+"/Results", mode: "copy"

    input:
    file(f) from COUNT_SUMMARIES
    val(raw_counts) from RAW_COUNTS
    val(filtered_counts) from FILTERED_COUNTS
    file(f) from MSA_TO_COUNT
        .mix(CHIMERA_TO_COUNT)
        .mix(CLUSTERING_TO_COUNT)
        .mix(MULTIPLETONS_FILTER_TO_COUNT)
        .mix(SUBSAMPLING_TO_COUNT)
        .mix(TAXA_FILTER_TO_COUNT)
        .mix(LULU_TO_COUNT)
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
             ("9-TaxaFilter","all_taxaFilter_*.groups.summary"),
             ("10-MultipletonsFilter","all_multipletonsFilter_*.groups.summary"),
             ("11-Subsampling","all_subsampling_*.groups.summary"),
             ("12-Lulu","all_lulu_*.csv")
    ]

    write_summary(steps,counts,clustering_thresholds)
    """
}

/*
 *
 * Postprocessing
 *
 */


process AlphaDiversity {
    tag { "alphaDiv_${otuId}" }
    label "high_computation"
    label "mothur_script"
    publishDir params.outdir+"/Results/postprocessing/alpha_diversity", mode: "copy", pattern: "*.summary"

    input:
    tuple val(otuId), file(shared) from MAIN_OUTPUTS.alpha_div

    output:
    file("*.summary")

    script:
    """
    ${params.script_dir}/mothur.sh \
    --step=alphaDiversity \
    --otuId=${otuId}
    """
}

process FastTree {
	tag { "FastTree_${otuId}" }
	label "high_computation"
	// label "mothur_script"
	publishDir params.outdir+"/Results/postprocessing/unifrac", mode: "copy", pattern: "*.tre"

	input:
	tuple val(otuId), file(fasta), file(shared) from MAIN_OUTPUTS.tree

	output:
	tuple val(otuId), file(shared), file("*.tre") into FOR_UNIFRAC

	script:
	"""
	sed -r 's/.*(Otu[0-9]+)\\|.*/\\>\\1/' ${fasta} > relabeled_${fasta}
	FastTree -nt relabeled_${fasta} > FastTree_${otuId}.tre
	"""
}

if (!params.skipBetaDiversity) {
	process BetaDiv {
		tag { "betaDiv_${idThreshold}" }
		label "high_computation"
		label "mothur_script"
		publishDir params.outdir+"/Results/postprocessing/beta_diversity", mode: "copy", pattern: "*.summary"

		input:
		tuple val(idThreshold), file(shared) from MAIN_OUTPUTS.beta_div

		output:
		file("*.summary")

		script:
		"""
		${params.script_dir}/mothur.sh \
			--step=betaDiversity \
			--idThreshold=${idThreshold}
		"""
	}

	process UnifracDistPhylo {
		tag { "Unifrac_${otuId}_${mode}" }
		label "high_computation"
		label "r_script"
		publishDir params.outdir+"/Results/postprocessing/unifrac", mode: "copy", pattern: "*.csv"

		input:
		tuple val(otuId), file(shared), file(tree) from FOR_UNIFRAC
		each mode from ('weighted', 'unweighted')

		output:
		file("unifrac*.csv")

		script:
		"""
		#!/usr/bin/env Rscript
		source("${params.script_dir}/util.R")
		
		calculate_unifrac("${tree}", "${shared}", method="${mode}", otu_thresh=${otuId}, cores=${task.cpus})
		"""
	}
}

process SpeciesAssignment {
    tag { "species_${otuId}" }
    label "high_computation"
    label "r_script"
    publishDir params.outdir+"/Results/postprocessing/", mode: "copy", pattern: "*.csv"
	errorStrategy 'ignore'

    input:
    tuple val(otuId), file(fasta), file(tax) from MAIN_OUTPUTS.species_assign

    output:
    file("*.csv")

    script:
    """
    #!/usr/bin/env Rscript
    source("${params.script_dir}/util.R")

    res <- get_species("${fasta}", "${tax}", "${params.speciesDB}")

    write.csv(res, "species_${otuId}.csv")
    """
}
