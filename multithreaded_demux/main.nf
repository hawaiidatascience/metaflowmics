#!/usr/bin/env nextflow

def helpMessage() {
    log.info"""
    ===================================
    nf-multithreaded-demultiplexing
    ===================================

    Usage:
    nextflow run nextflow-demux -profile local --inputdir test/
    
    --------------- Mandatory arguments ---------------
    --inputdir    Path to data folder. It should include:
                      - 2 index fastq files (unzipped) matching the glob pattern "*_I{1,2}*.fastq.gz"
                      - 2 read fastq files (unzipped) matching the glob pattern "*_R{1,2}*.fastq.gz"
                      - 1 barcode file (extension: ".csv"), comma separated, with no header and 3 columns: 
                         (sample name, forward barcode, reverse complement of reverse barcode)

    ---------------- Optional arguments ---------------
    -profile          Select a configuration from the conf/ folder. Default is "local"
    --outdir       Path to output folder. Default: "./demultiplexed"
    --max_mismatches  Maximum number of allowed mismatches between index and barcode. Default: 1
    --nsplits         Number of file chunks to create for multithreading. Default: 10
    --matching        By default, the order in which the barcodes pair match the index pair is inferred from the data.
                      To change this behavior, set this parameter to either "ordered" or "reversed". Default: "auto"
    --multimap        Strategy for handling index pairs matching multiple samples. 
                      Default consists in assigning to the sample with the least mismatches and discarding the pair if 
                      multiple sample achieve the minimum score. To change this behavior, set this parameter to:
					  - "discard" to discard any multimapping. This strategy is strongly affected by 
                        the number of allowed mismatches, and can result in discarding most of the reads.
					  - "min_all" to assign to all the samples with the minimum score.
					  - "all" to assign to all matching samples.
    """.stripIndent()
}

// Show help message
params.help = false
if (params.help){
    helpMessage()
    exit 0
}

Channel
    .fromFilePairs("${params.inputdir}/*_R{1,2}*.fastq*", size: params.singleEnd ? 1 : 2, flat: true)
    .ifEmpty { error "Cannot find any sequencing read in ${params.inputdir}/" }
    .into { INPUT_SEQ ; INPUT_SEQ_FOR_QC }

Channel
    .fromFilePairs("${params.inputdir}/*_I{1,2}*.fastq*", size: params.singleBarcoded ? 1 : 2 , flat: true)
    .ifEmpty { error "Cannot find any indexing read in ${params.inputdir}/" }
    .into { INPUT_INDEX ; INPUT_INDEX_FOR_GUESS }

if (params.singleBarcoded) {
	INPUT_INDEX
		.splitFastq(by: params.n_per_file.toInteger(), file: true)
		.map{[it[0], it[1..-1]]}
		.into{IDX_SPLIT; SPLITS}
} else {
	INPUT_INDEX
		.splitFastq(by: params.n_per_file.toInteger(), file: true, pe: !params.singleBarcoded)
		.map{[it[0], it[1..-1]]}
		.into{IDX_SPLIT; SPLITS}
}

INPUT_SEQ
	.splitFastq(by: params.n_per_file.toInteger(), file: true, pe: !params.singleEnd)
	.map{[it[0], it[1..-1]]}
	.set{SEQ_SPLIT}

SPLITS.count().map{1..it}.flatten().into{COUNTER_IDX ; COUNTER_SEQ}

COUNTER_IDX.merge(IDX_SPLIT).set{INPUT_IDX_SPLIT}
COUNTER_SEQ.merge(SEQ_SPLIT).set{INPUT_SEQ_SPLIT}

process Fastqc {
    tag { "Fastqc" }
    publishDir  params.outdir, mode: "copy", pattern: "*.html"
    label "high_computation"

    input:
    set val(v), file(fastqs) from INPUT_SEQ_FOR_QC.map{[it[0], it[1..-1]]}

    output:
    file("*.html")

    script:
    """
    fastqc -o . --threads ${task.cpus} ${fastqs}
    """
}

process GuessMatchOrder {
    tag { "GuessMatchOrder" }
    label "low_computation"
	publishDir params.outdir, pattern: "barcodes_ok.csv"
    
    input:
    set (val(v), file(fastqs)) from INPUT_INDEX_FOR_GUESS.map{[it[0], it[1..-1]]}
	file(meta) from Channel.fromPath("${params.inputdir}/*.csv")

    output:
	file("barcodes_ok.csv") into (BARCODES_FILE, BARCODES_FILE_FOR_MODEL)

    script:
    """
    #!/usr/bin/env bash

    bash ${workflow.projectDir}/guess_matching_order.sh ${meta} ${params.matching} ${fastqs} > barcodes_ok.csv

    """	
}

process To_h5 {
    tag { "to_h5_${split}" }
    //publishDir "${params.outdir}/h5_data", mode: "copy", pattern: "*.h5"
    label "python_script"
    label "medium_computation"
    
    input:
    set (val(split), val(v), file(fastqs)) from INPUT_IDX_SPLIT

    output:
    set (val(split), file("*.h5")) into H5_FILES
	file "*.h5" into H5_FILES_FOR_ERROR_MODEL

    script:
    """
    python3 ${workflow.projectDir}/load.py --fastqs ${fastqs} --split ${split}
	"""
}

process ErrorModel {
    tag { "ErrorModel" }
    publishDir "${params.outdir}/other", mode: "copy", pattern: "*.{h5,pdf}"
    label "python_script"
    label "high_computation"
    
    input:
	file h5 from H5_FILES_FOR_ERROR_MODEL.collect()
	file meta from BARCODES_FILE_FOR_MODEL

    output:
    file("transition_probs.h5") into ERROR_MODEL
	file("*.pdf")

    script:
    """
    python3 ${workflow.projectDir}/error_model.py --max_dist ${params.max_mismatches} --n_bases ${params.n_bases.toInteger()} --meta ${meta}
	"""
}

process IndexMapping {
    tag { "IndexMapping_${split}" }
    publishDir "${params.outdir}/sample_idx_mapping", mode: "copy", pattern: "*.tsv"
    label "python_script"
    label "medium_computation"
	stageInMode "copy"
    
    input:
    set val(split), file(h5), file(model), file(meta) from H5_FILES.combine(ERROR_MODEL.combine(BARCODES_FILE))

    output:
    set (val(split), file("demux_info*.tsv")) into DEMUX_CHANNEL
	file("sample_counts*.csv") into SAMPLE_COUNTS

    script:
    """
    python3 ${workflow.projectDir}/dada_demux_index.py --data ${h5} --meta ${meta} --error-model ${model} --split ${split} --max-mismatches ${params.max_mismatches}
	"""
}

process SampleSizeDistribution {
    tag { "SampleSizeDistribution" }
    publishDir params.outdir, mode: "copy", pattern: "*.pdf"

    label "python_script"
    label "medium_computation"

    input:
	file f from SAMPLE_COUNTS.collect()
	
    output:
    file("*.pdf")

    script:
    """
    python3 ${workflow.projectDir}/plot_sample_size_distribution.py
    """
}

process WriteSampleFastq {
    tag { "WriteSampleFastq_${split}" }
    stageOutMode "move"
    label "python_script"
    label "low_computation"
    
    input:
    set (val(split), val(v), file(fastqs), file(demux_info)) from INPUT_SEQ_SPLIT.join(DEMUX_CHANNEL)

    output:
    file("*.fastq") optional true into DEMUX_SEQ_SPLIT
    
    script:
    """
    python3 ${workflow.projectDir}/demux_fastq.py --fastqs ${fastqs} --mapping ${demux_info}
    """
}

DEMUX_SEQ_SPLIT.flatten().collectFile(storeDir: "${params.outdir}/reads").set{DEMUX_SEQ}

process Gzip {
	tag { "Gzip" }
    label "high_computation"

    input:
	file f from DEMUX_SEQ.collect()
	
    script:
    """
    ls ${params.outdir}/reads/*.fastq | xargs -P ${task.cpus} gzip
    """	
}
