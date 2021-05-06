#!/usr/bin/env nextflow

nextflow.enable.dsl=2
params.options = [:]

module_dir = "../modules"

// Functions
include { helpMessage; saveParams } from "./util.nf"


/*
 *
 Beginning of the pipeline
 *
 */

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
	label "python_script"

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
	label "python_script"
	publishDir params.outdir, pattern: "barcodes_ok.csv"
    
    input:
    set (val(v), file(fastqs)) from INPUT_INDEX_FOR_GUESS.map{[it[0], it[1..-1]]}
	file(meta) from Channel.fromPath("${params.inputdir}/*.csv")

    output:
	file("barcodes_ok.csv") into (BARCODES_FILE, BARCODES_FILE_FOR_MODEL)

    script:
    """
    #!/usr/bin/env bash

    bash ${params.script_dir}/demux/guess_matching_order.sh ${meta} ${params.matching} ${params.reverseComplement} ${fastqs} > barcodes_ok.csv
    """	
}

process To_h5 {
    tag { "to_h5_${split}" }
    label "python_script"
    label "medium_computation"
    
    input:
    set (val(split), val(v), file(fastqs)) from INPUT_IDX_SPLIT

    output:
    set (val(split), file("*.h5")) into H5_FILES
	file "*.h5" into H5_FILES_FOR_ERROR_MODEL

    script:
    """
    python3 ${params.script_dir}/demux/load.py \
        --fastqs ${fastqs} \
        --split ${split} \
        \$([ "${params.reverseComplement}" == true ] && echo "--rc" || echo "")
	"""
}

process ErrorModel {
    tag { "ErrorModel" }
    publishDir "${params.outdir}/other", mode: "copy", pattern: "*.{h5,html}"
    label "python_script"
    label "high_computation"
    
    input:
	file h5 from H5_FILES_FOR_ERROR_MODEL.collect()
	file meta from BARCODES_FILE_FOR_MODEL

    output:
    file("transition_probs.h5") into ERROR_MODEL
	file("*.html")

    script:
    """
    python3 ${params.script_dir}/demux/error_model.py --max_dist ${params.max_mismatches} --n_bases ${params.n_bases.toInteger()} --meta ${meta}
	"""
}

process IndexMapping {
    tag { "IndexMapping_${split}" }
    publishDir "${params.outdir}/sample_idx_mapping", mode: "copy", pattern: "*.tsv"
    label "python_script"
    label "high_computation"

	stageInMode "copy"
    
    input:
    set val(split), file(h5), file(model), file(meta) from H5_FILES.combine(ERROR_MODEL.combine(BARCODES_FILE))

    output:
    set (val(split), file("demux_info*.tsv")) into DEMUX_CHANNEL
	file("sample_counts*.csv") into SAMPLE_COUNTS

    script:
    """
    python3 ${params.script_dir}/demux/dada_demux_index.py --data ${h5} --meta ${meta} --error-model ${model} --split ${split} --max-mismatches ${params.max_mismatches}
	"""
}

process SampleSizeDistribution {
    tag { "SampleSizeDistribution" }
    publishDir params.outdir, mode: "copy", pattern: "*.{html,csv}"

    label "python_script"
    label "medium_computation"

    input:
	file f from SAMPLE_COUNTS.collect()
	
    output:
    file("*.html")
	file("sample_sizes.csv")

    script:
    """
    python3 ${params.script_dir}/demux/plot_sample_size_distribution.py
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
    python3 ${params.script_dir}/demux/demux_fastq.py --fastqs ${fastqs} --mapping ${demux_info}
    """
}

DEMUX_SEQ_SPLIT.flatten().collectFile(storeDir: "${params.outdir}/reads", sort: true).set{DEMUX_SEQ}

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
