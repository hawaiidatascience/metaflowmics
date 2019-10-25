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
                      - 1 barcode file named "barcodes.csv", comma separated, with no header and 3 columns: 
                         (sample name, forward barcode, reverse complement of reverse barcode)

    ---------------- Optional arguments ---------------
    -profile          Select a configuration from the conf/ folder. Default is "local"
    --outputdir       Path to output folder. Default: "./demultiplexed"
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
    .fromFilePairs("${params.inputdir}/*_R{1,2}*.fastq.gz", flat: true)
    .ifEmpty { error "Cannot find any sequencing read matching in ${params.inputdir}" }
    .into { INPUT_SEQ ; INPUT_SEQ_FOR_QC }

Channel
    .fromFilePairs("${params.inputdir}/*_I{1,2}*.fastq.gz", flat: true)
    .ifEmpty { error "Cannot find any indexing read matching in ${params.inputdir}" }
    .into { INPUT_INDEX ; INPUT_INDEX_FOR_GUESS }

meta = "${params.inputdir}/barcodes.csv"

n_per_file = (file("${params.inputdir}/*_I1*.fastq.gz").get(0).countFastq() / params.nsplits).toInteger()

Channel.from(1..params.nsplits).into{COUNTER_IDX ; COUNTER_SEQ}

INPUT_IDX_SPLIT = COUNTER_IDX.merge(INPUT_INDEX.splitFastq(by: n_per_file, file: true, pe: true))
INPUT_SEQ_SPLIT = COUNTER_SEQ.merge(INPUT_SEQ.splitFastq(by: n_per_file, file: true, pe: true))

process Fastqc {
    tag { "Fastqc" }
    publishDir  params.outputdir, mode: "copy", pattern: "*.html"
    label "high_computation"

    input:
    set val(v), file(fwd), file(rev) from INPUT_SEQ_FOR_QC

    output:
    file("*.html")

    script:
    """
    fastqc -o . --threads ${task.cpus} ${fwd} ${rev}
    """
}

process GuessMatchOrder {
    tag { "GuessMatchOrder" }
    label "low_computation"
    
    input:
    set (val(v), file(fwd), file(rev)) from INPUT_INDEX_FOR_GUESS

    output:
	stdout into GUESS

    script:
    """
    #!/usr/bin/env bash

    if [ ${params.matching} == "auto" ]; then
		symbol=\$(zcat ${fwd} | head -c 2)

		paste --delimiters=' ' <(\
                zgrep --no-group-separator "^\${symbol}" ${fwd} -A1 | grep -v \${symbol}) <( \
                grep --no-group-separator "^\${symbol}" ${rev} -A1 | grep -v \${symbol}) \
			   | grep -v N \
			   | sort | uniq -c | sort -rnk1 | head -n20 \
			   | sed 's/^[[:space:]]*//g' | sed 's/ /,/g' | cut -d, -f2,3 \
			   > pairs_freqs.txt

		awk -F, '{OFS=","} {print \$2,\$1}' pairs_freqs.txt > pairs_freqs_rev.txt

		same1=\$(comm -12 <(sort <(cut -d, -f2,3 ${meta})) <(sort pairs_freqs.txt) | wc -l)
		same2=\$(comm -12 <(sort <(cut -d, -f2,3 ${meta})) <(sort pairs_freqs_rev.txt) | wc -l)

		[ \$same2 -ge \$same1 ] && echo -n "--revMapping" || echo -n " "
    else
        [ ${params.matching} == "reversed" ] && echo -n "--revMapping" || echo -n " "
    fi
    """	
}

process IndexMapping {
    tag { "IndexMapping_${split}" }
    publishDir "${params.outputdir}/other", mode: "copy", pattern: "*.tsv"
    label "python_script"
    label "low_computation"
    
    input:
    set (val(split), val(v), file(fwd), file(rev), val(guess)) from INPUT_IDX_SPLIT.combine(GUESS)

    output:
    set (val(split), file("demux_info*.tsv")) into DEMUX_CHANNEL

    script:
    """
    python3 ${workflow.projectDir}/demux_index.py --fwd ${fwd} --rev ${rev} --meta ${meta} --max_mismatches ${params.max_mismatches} ${guess} --strategy ${params.multimap} > demux_info_${split}.tsv
	"""
}

process WriteSampleFastq {
    tag { "WriteSampleFastq_${split}" }
    stageOutMode "move"
    label "python_script"
    label "low_computation"
    
    input:
    set (val(split), val(v), file(fwd), file(rev), file(demux_info)) from INPUT_SEQ_SPLIT.join(DEMUX_CHANNEL)

    output:
    file("*.fastq") optional true into DEMUX_SEQ_SPLIT
    
    script:
    """
    python3 ${workflow.projectDir}/demux_fastq.py --fwd ${fwd} --rev ${rev} --mapping ${demux_info}
    """
	
}

DEMUX_SEQ_SPLIT.flatten().collectFile(storeDir: "${params.outputdir}/reads").set{DEMUX_SEQ}

process SampleSizeDistribution {
    tag { "SampleSizeDistribution" }
    publishDir params.outputdir, mode: "copy", pattern: "*.png"

    label "python_script"
    label "high_computation"

    input:
    file f from DEMUX_SEQ.collect()
	
    output:
    file("*.png")

    script:
    """
    python3 ${workflow.projectDir}/plot_sample_size_distribution.py ${params.outputdir}/reads
    ls ${params.outputdir}/reads/*.fastq | xargs -P ${task.cpus} gzip
    """
}
