#!/usr/bin/env nextflow

def helpMessage() {
    log.info"""
    ===================================
     ITS-rDNA-pipeline
    ===================================
    Usage:
    nextflow run ITS-pipeline --reads '*_R1.fastq.gz' --locus ITS1 --pairedEnd 0 -profile local

    ---------------------------------- Mandatory arguments ----------------------------------------

      --reads         Path to input data
      --uniteDB       Path to taxonomic database to be used for annotation (e.g. /path/unite.fa.gz)

    ---------------------------------- Optional arguments ----------------------------------------

      --pairedEnd     If your data is paired-end
      --locus         Sequenced ITS region (ITS1,ITS2 or ALL). Default: ITS1
      --outdir        The output directory where the results will be saved. 
                      Default: ./ITS-pipeline_outputs

    [Trimming]
      --minLen        Remove reads with length less than minLen. 
                      minLen is checked after trimming and truncation. Default: 20
      --minQuality    Quality filtering. Keep contigs with more than <minQuality> bases in 
                      <minPercentHighQ>% bases. Default: 25

    [Quality filtering]
      --minPercentHighQ  Quality filtering. Keep contigs with more than <minQuality> bases in 
                         <minPercentHighQ>% bases. Default: 90
      --minReads         Discard samples with less than <minReads> reads. Default: 50
      --minDerep         Discard samples with less than <minDerep> unique contigs. Default: 5

    [Taxonomy assignment]
      --taxaMinId        The minimum identity for calling a taxonomy. Default 50%

    [Co-occurence pattern correction]
	--skipLulu               Skip the Lulu step
    --min_ratio_type         Function to compare the abundance of a parent OTU and its daughter.
                             Default: min
    --min_ratio              Minimum abundance ratio between parent and daughter
                             (across all samples). Default: 1
    --min_rel_cooccurence    Proportion of the parent samples where daughter occurs. Default: 1

    """.stripIndent()
}

// Show help message
params.help = false
if (params.help){
    helpMessage()
    exit 0
}

def summary = [:]
summary['database'] = params.uniteDB
summary['paired end'] = params.pairedEnd
summary['locus'] = params.locus
summary['min read quality'] = params.minQuality
summary['% of bases with {read quality}'] = params.minPercentHighQ
summary['min read length'] = params.minLen
summary['min reads per sample'] = params.minReads
summary['min unique sequence per sample'] = params.minDerep
summary['clustering similarity thresholds'] = params.clusteringThresholds
summary['Skip LULU'] = params.skipLulu
summary['Lulu ratio type'] = params.min_ratio_type
summary['Lulu parent/daughter min similarity'] = params.min_match
summary['Lulu parent/daughter min abundance ratio'] = params.min_ratio
summary['Lulu parent/daughter min co-occurrence'] = params.min_rel_cooccurence

file(params.outdir).mkdir()
summary_handle = file("${params.outdir}/parameters_summary.log")
summary_handle << summary.collect { k,v -> "${k.padRight(50)}: $v" }.join("\n")


/*
 *
 Beginning of the pipeline
 *
 */

def clusteringThresholds = params.clusteringThresholds.toString().split(',').collect{it as int}

if ( !params.pairedEnd ) {
    read_path = params.reads.replaceAll("\\{1,2\\}","1")
} else {
    read_path = params.reads
}

Channel
    .fromFilePairs( read_path, size: params.pairedEnd ? 2 : 1, flat: true )
    .ifEmpty { error "Cannot find any reads matching: ${params.reads}" }
    .into { INPUT_FASTQ ; INPUT_FASTQ_FOR_COUNT }

RAW_COUNTS = INPUT_FASTQ_FOR_COUNT.collect{ ["'${it[0]}': ${it[1].countFastq()}"] }

/*
 *
 * Remove small files (< 50 {params.minReads}) to prevent odd behaviors 
 *
 */


INPUT_FASTQ_OK = INPUT_FASTQ
    .filter{ it[1].countFastq() > params.minReads }
    .map{ it -> tuple(it[0],it[1..-1]) }

/*
 *
 * Extract ITS(1,2 or both) region with pair-end (or not)
 *
 */

process ExtractITS {
    tag { "ITSextraction.${pairId}" }
    label "medium_computation"
	label "python_script"
    
    publishDir params.outdir+"/Misc/1-ITSxpress", mode: "copy"
    
    input:
        tuple val(pairId), file(reads) from INPUT_FASTQ_OK    
    output:
        tuple val(pairId), file("${pairId}_${params.locus}.fastq") \
    into (ITS_FASTQ, ITS_FASTQ_FOR_COUNT)
    	
    script:
    """
    #!/usr/bin/env bash

    read_pair=(${reads})
    
    rev_arg=\$([ ${params.pairedEnd} == true ] && echo "--fastq2 \${read_pair[1]}" || echo "--single_end")

    itsxpress --region ${params.locus} --threads ${task.cpus} \
    	      --outfile ${pairId}_${params.locus}.fastq --log ITSxpress_${pairId}.log \
              --fastq \${read_pair[0]} \${rev_arg}
    """
}

ITS_FASTQ_COUNTS = ITS_FASTQ_FOR_COUNT.collect{ ["'${it[0]}': ${it[1].countFastq()}"] }

/*
 *
 * Remove any sequence with a "N" nucleotide (for dada2) or with less than {params.minLen} nucleotides
 *
 */

process RemoveSmallAndNseqs {
    tag { "removeN.${pairId}" }
    label "low_computation"        
    label "python_script"

    publishDir params.outdir+"/Misc/2-NandSmallSeqsFiltering", mode: "copy"

    input:
    tuple val(pairId), file(itsFile) from ITS_FASTQ
    
    output:
    tuple val(pairId), file("${pairId}_noN.fastq") into (NO_N_FASTQ,NO_N_FASTQ_FOR_COUNT)

    script:
    """
    #!/usr/bin/env python3
    from Bio import SeqIO

    seqs = [ seq for seq in SeqIO.parse("${itsFile}","fastq")
             if "N" not in seq.seq and len(seq) > ${params.minLen} ]

    SeqIO.write(seqs, "${pairId}_noN.fastq", "fastq")

    """        
}

NO_N_COUNTS = NO_N_FASTQ_FOR_COUNT.collect{ ["'${it[0]}': ${it[1].countFastq()}"] }

/*
 *
 * Keep contigs with more than {params.minPercentLowQ}% bases with quality {params.minQuality}
 *
 */


process QcFilter {
    tag { "filteredITS.${pairId}" }
    label "low_computation"

    publishDir params.outdir+"/Misc/3-QualityFiltering", mode: "copy"

    input:
    tuple val(pairId), file(itsFile) from NO_N_FASTQ.filter{ it[1].size() > 0 }
    
    output:
	tuple val(pairId), file("${pairId}_filtered.fastq.gz") into FILTERED_OUT

    script:
    """
    fastq_quality_filter -z -q ${params.minQuality} \
                            -p ${params.minPercentHighQ} \
                            -i ${itsFile} \
                            -o ${pairId}_filtered.fastq.gz \
    """    
}

// Remove files with less than params.minReads
FILTERED_OUT
	.map { [it[1].countFastq(), it[0], it[1]] }
	.filter{ it[0] > params.minReads }
	.multiMap { it ->
	    all: [it[1], it[2]]
        counts: [it[1],it[0]] }
	.set{FILTERED_FASTQ}

FILTERED_FASTQ.counts
	.collect{ ["'${it[0]}': ${it[1]}"]}
	.set{FILTERED_COUNTS}

/*
 *
 * Dereplicates fastq and store in an abundance vector ("RDS" R object)
 *
 */

process Dereplication {
    tag { "dereplication.${pairId}" }
    label "low_computation"        
    label "r_script"

    publishDir params.outdir+"/Misc/4-Dereplication", mode: "copy"

    input:
    tuple val(pairId), file(filtRead) from FILTERED_FASTQ.all
    
    output:
	tuple val(pairId), file("*_derep.RDS") into DEREP_RDS
    tuple val(pairId), file("*_derep.fasta") into DEREP_FASTA
    file("*derep.RDS") into DEREP_FOR_COUNT_SUMMARY
	
    script:
	"""
	#!/usr/bin/env Rscript 
	library(dada2)

	derep <- derepFastq("${filtRead}")
	saveRDS(derep,"${pairId}_R12_derep.RDS")
	uniquesToFasta(derep,"${pairId}_R12_derep.fasta")
	"""
}

/*
 *
 * Remove chimeric sequences using VSEARCH
 *
 */

process ChimeraRemoval {
    tag { "ChimeraRemoval.${pairId}" }
    label "low_computation"        
    label "require_vsearch"

    publishDir params.outdir+"/Misc/5-ChimeraRemoval", mode: "copy"

    input:
    tuple val(pairId),file(derepFa) from DEREP_FASTA.filter{ it[1].size() > 0 }
	
    output:
    tuple val(pairId),file("*_noChimeras.fasta") into DEREP_FASTA_NO_CHIMERA
	
    script:
    """
	vsearch --threads ${task.cpus} \
		--uchime3_denovo ${derepFa} \
		--sizein \
		--sizeout \
		--fasta_width 0 \
		--nonchimeras ${pairId}_R12_noChimeras.fasta \
    """
}

DEREP_FASTA_NO_CHIMERA
    .filter{ it[1].countFasta() > params.minDerep }
    .set{ DEREP_FASTA_NO_CHIMERA_FILT }

/*
 *
 * Build error model using dada2
 *
 */

process LearnErrors {
    tag { "LearnErrors.${pairId}" }
    label "medium_computation"
    label "r_script"
    
    publishDir params.outdir+"/Misc/6-Denoising", mode: "copy", pattern: "*{.RDS,.png}"
    
    input:
    tuple val(pairId), file(derepRDS), file(sample_nochim) \
    from DEREP_RDS.join(DEREP_FASTA_NO_CHIMERA_FILT)
    
    output:
    tuple val(pairId), file("*_errors.RDS"), file("*_nochimera.RDS") into ERRORS_AND_DEREP

    script:
    """
    #!/usr/bin/env Rscript
    source("${params.script_dir}/util.R")

    # Extract no chimera sequences and store into an R object for dada2
    extract_chimeras("${derepRDS}","${sample_nochim}","${pairId}")

    # Build error model
    learn_error_rates("${pairId}_R12_nochimera.RDS","${pairId}",rds=TRUE)
    """
}

/*
 *
 * Run dada2 denoising algorithm with the error models in the previous step
 *
 */

process Denoise {
    tag { "Denoising.${pairId}" }
    label "medium_computation"
    label "r_script"
    
    publishDir params.outdir+"/Misc/6-Denoising", mode: "copy"

    input:
    tuple pairId, file(err), file(nochimera) from ERRORS_AND_DEREP
	
    output:
    tuple file(nochimera),file("*_dada.RDS") into DADA_RDS
	
    script:
    """
    #!/usr/bin/env Rscript
    source("${params.script_dir}/util.R")

    dada_denoise("${err}","${nochimera}","${pairId}_R12",derep.rds=TRUE,paired=FALSE)
    """
}

/*
 *
 * Construct ESV tables with denoised contigs
 *
 */

process MakeEsvTable {
    tag { "esvTable" }
    label "high_computation"
    label "r_script"
    
    publishDir params.outdir+"/Misc/7-Clustering", mode: "copy"
    
    input:
	file(denoised) from DADA_RDS.collect()
	file(derep) from DEREP_FOR_COUNT_SUMMARY.collect()
	
    output:
    file("esv_merged_to_cluster.fasta") into ESV_ALL_SAMPLES_TO_CLUSTER
    tuple val(100), file("all_esv.fasta"), file("clustering_100.csv") into ESV_OUT
    file("count_summary.tsv") into DENOISING_SUMMARY
	file('clustering_100.csv') into ESV_TO_COUNT
	
    script:
    """
    #!/usr/bin/env Rscript
    source("${params.script_dir}/util.R")

    esv_table(paired=FALSE)
    merge_fasta_for_vsearch("raw_sequence_table.csv")
    """    
}

/*
 *
 * VSEARCH OTU clustering at several thresholds 
 *
 */

process Clustering {
    tag { "clustering.${otuId}" }
    label "high_computation"
    label "require_vsearch"
    
    publishDir params.outdir+"/Misc/7-Clustering", mode: "copy"
    
    input:
    each otuId from clusteringThresholds.findAll{ it != 100}
    file(esvFasta) from ESV_ALL_SAMPLES_TO_CLUSTER
	
    output:
	tuple val(otuId), file("OTUs_*.fasta"), file("clustering*.csv") into CLUSTERING_OUT
	file("clustering*.csv") into CLUSTERING_TO_COUNT
    
    script:
    """
    vsearch --threads ${task.cpus} \
            --cluster_size ${esvFasta} \
            --id 0.${otuId} \
            --strand plus \
            --sizein \
            --sizeout \
            --fasta_width 0 \
            --uc clusters${otuId}.uc \
            --relabel OTU${otuId}_ \
            --centroids OTUs_${otuId}.fasta \
            --otutabout clustering_${otuId}.tsv

    cat clustering_${otuId}.tsv | tr "\\t" "," > clustering_${otuId}.csv
    """
}

/*
 *
 * Preliminary step before LULU: blast each contig against the others to look for sequence similarities
 *
 */

if (!params.skipLulu) {
	process PreLulu {
		tag { "preLulus.${otuId}" }
		label "high_computation"
		label "require_vsearch"

		publishDir params.outdir+"/Misc/8-LULU_correction", mode: "copy"

		input:
		tuple val(otuId), file(fasta), file(abundance) from CLUSTERING_OUT.mix(ESV_OUT)

		output:
		tuple val(otuId), file("match_list_${otuId}.txt"), file(fasta), file(abundance) into FOR_LULU
		
		script:	
		"""
		vsearch --usearch_global ${fasta} \
            --threads ${task.cpus} \
            --db ${fasta} --self \
            --id .84 \
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
	 * LULU
	 *
	 */

	process Lulu {
		tag { "Lulu.${otuId}" }
		label "medium_computation"
		label "r_script"

		publishDir params.outdir+"/Misc/8-LULU_correction", mode: "copy", pattern: "{abundance_table_*.csv}"

		input:
		tuple val(otuId), file(matchlist), file(fasta), file(abundance) from FOR_LULU

		output:
		tuple val(otuId), file("lulu_ids_${otuId}.csv"), file(fasta), file("all_lulu_${otuId}.csv") into FOR_LULU_FILTER
		file("all_lulu_*.csv") into ABUNDANCE_LULU_TO_COUNT
		file("lulu*.log_*") optional true

		script:
		"""
		#!/usr/bin/env Rscript
		source("${params.script_dir}/util.R")

		lulu_curate("${abundance}","${matchlist}","${otuId}","${params.min_ratio_type}","${params.min_ratio}","${params.min_match}","${params.min_rel_cooccurence}")
		"""
	}

	/*
	 *
	 * Filter out FASTA sequences that LULU merged with a more abundant sequence
	 *
	 */

	process ExtractFastaLulu {
		tag { "extractFastaLulu.${otuId}" }
		label "high_computation"
		label "python_script"

		input:
		tuple otuId, file(ids), file(fasta), file(abundance) from FOR_LULU_FILTER

		output:
		tuple otuId, file("lulu_${otuId}.fasta"), file(abundance) into LULU_FILTERED

		script:

		"""
		#!/usr/bin/env python3
		from Bio import SeqIO
		import pandas as pd

		ids = pd.read_csv("${ids}", header=None)[0].values
		sequences = [ seq for seq in SeqIO.parse("${fasta}","fasta") if seq.id.split(";")[0] in ids ]
		SeqIO.write(sequences,"lulu_${otuId}.fasta","fasta")
		"""
	}

} else {
	CLUSTERING_OUT.mix(ESV_OUT).set{LULU_FILTERED}
	Channel.empty().set{ABUNDANCE_LULU_TO_COUNT}
}

Channel.fromPath(params.uniteDB)
	.combine(Channel.from(clusteringThresholds))
	.map{it[0]}
	.set{DB_CHANNEL}

process ClassificationSintax {
    tag { "classification.${otuId}" }
    label "high_computation"
    label "require_vsearch"
	publishDir "${params.outdir}/Results", mode: "copy"

    input:
	tuple otuId, file(fasta), file(abundance) from LULU_FILTERED
	file(db) from DB_CHANNEL

    output:
    tuple file("sequences*.fasta"), file("abundance_table*.tsv"), file("annotations_*.tsv")

    script:

    """
    vsearch --threads ${task.cpus} --db ${db} \
            --sintax ${fasta} --sintax_cutoff ${params.confidenceThresh} \
            --tabbedout annotations_sintax_${otuId}.tsv

	cp ${fasta} sequences_${otuId}.fasta
	cp ${abundance} abundance_table_${otuId}.tsv
    """
}

/*
 *
 * Generates a summary file with the sequence count at each step for each sample (sample x step table) 
 *
 */
// RAW_COUNTS,ITS_FASTQ_COUNTS,NO_N_COUNTS,FILTERED_COUNTS,NO_CHIMERA_COUNTS

process SummaryFile {
    tag { "summary" }
    label "medium_computation"
    label "python_script"

    publishDir "${params.outdir}/Results", mode: "copy", pattern: "*.tsv"
    
    input:
    file(f) from DENOISING_SUMMARY
    val(raw_counts) from RAW_COUNTS
	val(its_counts) from ITS_FASTQ_COUNTS
	val(no_n_counts) from NO_N_COUNTS
	val(filtered_counts) from FILTERED_COUNTS
    file(f) from CLUSTERING_TO_COUNT
		.mix(ESV_TO_COUNT)
		.mix(ABUNDANCE_LULU_TO_COUNT)
		.collect()
	
    output:
    file("sequences_per_sample_per_step_*.tsv") into STEPS_SUMMARY
	
    script:
    """
    #!/usr/bin/env python3
    from generate_step_summary import write_summary

    clustering_threshold = ${clusteringThresholds}
    counts = { '0-RawData': {${raw_counts.join(', ')}},
               '1-ITSxpress': {${its_counts.join(', ')}},
               '2-NandSmallSeqsFiltering': {${no_n_counts.join(', ')}},
               '3-QualityFiltering': {${filtered_counts.join(', ')}}
             }

    steps = [("0-RawData", -1),
             ("1-ITSxpress",-1),
             ("2-NandSmallSeqsFiltering",-1),
             ("3-QualityFiltering",-1),
             ("7-Clustering","clustering_*.csv"),
             ("8-LULU_correction","abundance_table_*.csv")                                         
    ]
    
    write_summary(steps,counts,clustering_threshold)
    """
}

