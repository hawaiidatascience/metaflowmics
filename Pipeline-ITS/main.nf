#!/usr/bin/env nextflow

def helpMessage() {
    log.info"""
    ===================================
     ITS-rDNA-pipeline
    ===================================
    Usage:
    nextflow run ITS-pipeline --reads '*_R1.fastq.gz' --locus ITS1 --pairEnd 0 -profile manoa

    Mandatory arguments:
      -profile                      Hardware config to use. local / manoa
      --reads                       Path to input data

    Other arguments:
      --pairEnd                     To set if the data is paired-end. Default: single-end
      --reference                   Path to taxonomic database to be used for annotation (e.g. /path/unite.fa.gz)
      --locus                       Sequenced ITS region (ITS1,ITS2 or ALL)
      --outdir                      The output directory where the results will be saved

    Trimming arguments (optional):
      --minLen                      Remove reads with length less than minLen. minLen is enforced after trimming and truncation; default=20
      --minQuality                  Quality filtering. Keep contigs with more than {minQuality} bases in {minPercentHighQ} % bases. Default: 25
      --minPercentHighQ              Quality filtering. Keep contigs with more than {minQuality} bases in {minPercentHighQ} % bases. Default: 90
      --minReads                    Discard samples with less than {minReads} reads. Default: 50
      --minDerep                    Discard samples with less than {minDerep} unique contigs. Default: 5

    Taxonomy assignment (optional):
      --taxaMinId                   The minimum identity for calling a taxonomy. Default 97%
      
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

// Defining input channels
if( params.pairEnd ){
    Channel
        .fromFilePairs( "${params.reads}" )
        .ifEmpty { error "Cannot find any reads matching: ${params.reads}" }
        .into { INPUT_FASTQ ; INPUT_FASTQ_TO_QC }
} else {
    Channel
	.fromPath( "${params.reads}" )
	.ifEmpty { error "Cannot find any reads matching: ${params.reads}" }
    	.map { filename -> tuple(filename.baseName,
				 tuple(file(filename),filename.baseName))}
	.into { INPUT_FASTQ ; INPUT_FASTQ_TO_QC }
}

// Header log info
def summary = [:]
summary['Reads'] = params.reads
summary['locus'] = params.locus
summary['pairEnd'] = params.pairEnd
summary['Output dir'] = params.outdir
summary['Working dir'] = workflow.workDir
summary['Current home'] = "$HOME"
summary['Current user'] = "$USER"
summary['Current path'] = "$PWD"
summary['Config Profile'] = workflow.profile
log.info summary.collect { k,v -> "${k.padRight(15)}: $v" }.join("\n")
log.info "========================================="


/*
 *
 * Remove small files (< 50 {params.minReads}) to prevent odd behaviors 
 *
 */

process PreFiltering {
    tag { "pre_filter.${pairId}" }
    label "low_computation"
    publishDir "${params.outdir}/1-preFiltering", mode: "copy", pattern: "{*.fastq*, *.txt}"
    errorStrategy "${params.errorsHandling}"

    input:
        set val(pairId), file(fastq) from INPUT_FASTQ    
    output:
	set stdout, val(pairId), file(fastq) into INPUT_FASTQ_LEN
	file "${pairId}.removed.txt" optional true into BAD_SAMPLES

    script:
    """
    #!/usr/bin/env bash

    # Count the number of reads in the forward fastq
    filename="${fastq.get(0)}"

    if [ "\${filename##*.}" == "gz"]; then
        nreads=`zcat \$filename | wc -l`
    else
        nreads=`cat \$filename | wc -l`
    fi

    nreads=\$((\$nreads/4))

    # Write sample name if it does not pass the threshold
    if [ nreads -lt ${params.minReads} ]; then
       echo \$filename > ${pairId}.removed.txt;
    fi

    echo \$nreads
    """    
}

// Filter out files with less than params.minReads
INPUT_FASTQ_LEN
    .filter{ (it[0] as int) > params.minReads }
    .map{ nreads, pairId, filename -> tuple(pairId,filename) }
    .set{ INPUT_FASTQ_OK }

/*
 *
 * FastQC 
 *
 */

process runFastQC {
    tag { "fastQC.${pairId}" }
    label "low_computation"
    publishDir "${params.outdir}/0-fastQC", mode: "copy"
    errorStrategy 'ignore'

    input:
        set val(pairId), file(in_fastq) from INPUT_FASTQ_TO_QC

    output:
        file("${pairId}_fastqc/*.zip") into FASTQC_FILES

    """
    mkdir ${pairId}_fastqc

    if [ ${params.pairEnd} == true ]; then
        fastqc --outdir ${pairId}_fastqc ${in_fastq.get(0)} ${in_fastq.get(1)}
    else
        fastqc --outdir ${pairId}_fastqc ${in_fastq.get(0)}
    fi
    """
}


/*
 *
 * Extract ITS(1,2 or both) region with pair-end (or not)
 *
 */

process ExtractITS {
    tag { "ITSextraction.${pairId}" }
    label "low_computation"    
    publishDir "${params.outdir}/2-ITSxpress", mode: "copy"
    errorStrategy "${params.errorsHandling}"
    
    input:
        set val(pairId), file(reads) from INPUT_FASTQ_OK
    
    output:
        set val(pairId), file("${pairId}_${params.locus}.fastq") \
          into ITS_FASTQ
    	
    script:
    """
    if [ ${params.pairEnd} == true ]; then
        itsxpress --fastq ${reads[0]} \
                  --fastq2 ${reads[1]} \
                  --region ${params.locus} \
                  --log ITSxpress.log \
                  --outfile ${pairId}_${params.locus}.fastq \
                  --threads ${task.cpus}    
    else 
        itsxpress --fastq ${reads[0]} \
                  --single_end \
                  --region ${params.locus} \
                  --log ITSxpress.log \
                  --outfile ${pairId}_${params.locus}.fastq \
                  --threads ${task.cpus}    
    fi

    """
}

/*
 *
 * Remove any sequence with a "N" nucleotide (for dada2) or with less than {params.minLen} nucleotides
 *
 */


process RemoveSmallAndNseqs {
    tag { "removeN.${pairId}" }
    label "low_computation"        
    publishDir "${params.outdir}/3-NandSmallSeqsFiltering", mode: "copy"
    errorStrategy "${params.errorsHandling}"

    input:
        set val(pairId), file(itsFile) from ITS_FASTQ
    
    output:
	set val(pairId), file("${pairId}_noN.fastq") into NO_N_FASTQ

    script:
    """
    #!/usr/bin/env python3
    from Bio import SeqIO

    seqs = [ seq for seq in SeqIO.parse("${itsFile}","fastq")
             if "N" not in seq.seq and len(seq) > ${params.minLen} ]

    SeqIO.write(seqs, "${pairId}_noN.fastq", "fastq")

    """        
}

/*
 *
 * Keep contigs with more than {params.minPercentLowQ}% bases with quality {params.minQuality}
 *
 */


process QcFilter {
    tag { "filteredITS.${pairId}" }
    label "low_computation"
    publishDir "${params.outdir}/4-qualityFiltering", mode: "copy"
    errorStrategy "${params.errorsHandling}"

    input:
        set val(pairId), file(itsFile) from NO_N_FASTQ
    
    output:
	set stdout, val(pairId), file("${pairId}_filtered.fastq.gz") into FILTERED_FASTQ_ALL

    script:
    """
    #!/usr/bin/env bash

    nreads=`cat ${itsFile} | wc -l`

    if [ \$nreads -gt 0 ]; then
        fastq_quality_filter -z -q ${params.minQuality} \
                                -p ${params.minPercentHighQ} \
                                -i ${itsFile} \
                                -o ${pairId}_filtered.fastq.gz \
        > /dev/null
    
        nreads=`zcat ${pairId}_filtered.fastq.gz | wc -l`
        nreads=\$((\$nreads/4))        
    else
        touch "${pairId}_filtered.fastq.gz"
    fi

    echo \$nreads
    """    
}

// Remove files with less than params.minReads
FILTERED_FASTQ_ALL
    .filter{ (it[0] as int) > params.minReads }
    .map{ nreads, pairId, filename -> tuple(pairId,filename) }
    .set{ FILTERED_FASTQ }

/*
 *
 * Dereplicates fastq and store in an abundance vector ("RDS" R object)
 *
 */

process Dereplication {
    tag { "dereplication.${pairId}" }
    label "low_computation"        
    publishDir "${params.outdir}/5-Dereplication", mode: "copy"

    errorStrategy "${params.errorsHandling}"
    
    input:
        set val(pairId), file(filtRead) from FILTERED_FASTQ
    
    output:
	set val(pairId), file("${pairId}_derep.RDS") into DEREP_RDS
        set val(pairId), file("${pairId}_derep.fasta") into DEREP_FASTA
    
    script:
	"""
        #!/usr/bin/env Rscript 
        library(dada2)

        derep <- derepFastq("${filtRead}")
        saveRDS(derep,"${pairId}_derep.RDS")
        uniquesToFasta(derep,"${pairId}_derep.fasta")
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
    publishDir "${params.outdir}/6-ChimeraRemoval", mode: "copy"
    errorStrategy "${params.errorsHandling}"

    input:
        set val(pairId),file(derepFa) from DEREP_FASTA
    output:
        set stdout,val(pairId),file("${pairId}_noChimeras.fasta") into DEREP_FASTA_NO_CHIMERA
    script:
        """
        ncontigs=`grep -c ">" ${derepFa}`

        if [ \$ncontigs -gt 0 ]; then
	    vsearch --threads ${task.cpus} \
		    --uchime3_denovo ${derepFa} \
		    --sizein \
		    --sizeout \
		    --fasta_width 0 \
		    --nonchimeras ${pairId}_noChimeras.fasta \
            > /dev/null
             
        else
            touch "${pairId}_noChimeras.fasta"
        fi

        echo \$ncontigs

        """
}

// Remove files with less than minDerep (unique) contigs
DEREP_FASTA_NO_CHIMERA
    .filter{ (it[0] as int) > params.minDerep }
    .map{ ncontigs, pairId, filename -> tuple(pairId,filename) }
    .set{ DEREP_FASTA_NO_CHIMERA_FILT }

/*
 *
 * Build error model using dada2
 *
 */

process LearnErrors {
    tag { "LearnErrors.${pairId}" }
    label "low_computation"
    publishDir "${params.outdir}/7-Denoising", mode: "copy", pattern: "*{_errors.RDS,.png}"
    errorStrategy "${params.errorsHandling}"
    
    input:
        set val(pairId), file(derepRDS), file(sample_nochim) \
    from DEREP_RDS.join(DEREP_FASTA_NO_CHIMERA_FILT)
    
    output:
        set val(pairId), file("${pairId}_errors.RDS"), file("${pairId}_nochimera.RDS") into ERRORS_AND_DEREP

    script:
    """
    #!/usr/bin/env Rscript
    source("${workflow.projectDir}/scripts/util.R")

    # Extract no chimera sequences and store into an R object for dada2
    extractChimeras("${derepRDS}","${sample_nochim}","${pairId}")

    # Build error model
    learnErrorRates("${pairId}_nochimera.RDS","${pairId}")
    """
}

/*
 *
 * Run dada2 denoising algorithm with the error models in the previous step
 *
 */

process Denoise {
    tag { "Denoising.${pairId}" }
    label "low_computation"        
    publishDir "${params.outdir}/7-Denoising", mode: "copy"
    errorStrategy "${params.errorsHandling}"

    input:
        set pairId, file(err), file(derep) from ERRORS_AND_DEREP

    output:
        set pairId,file("${pairId}.dada.fasta") into DADA_FASTA
        set pairId,file("${pairId}.dada.RDS") into DADA_RDS
    script:
    """
    #!/usr/bin/env Rscript
    source("${workflow.projectDir}/scripts/util.R")

    dadaDenoise("${err}","${derep}","${pairId}")
    """
}

/*
 *
 * Construct ESV tables with denoised contigs
 *
 */

process MakeEsvTable {
    tag { "esvTable.${pairId}" }
    label "medium_computation"        
    publishDir "${params.outdir}/8-Clustering", mode: "copy"
    errorStrategy "${params.errorsHandling}"
    
    input:
	file mr from DADA_RDS.collect()
    output:
	file("merged.RDS")
        set val(100),file("esv_table.csv") into ABUNDANCE_TABLES_ESV
        set val(100),file("esv_seq.fasta") into ESV_ALL_SAMPLES,ESV_ALL_SAMPLES_LULU
        
    script:
    """
    #!/usr/bin/env Rscript
    source("${workflow.projectDir}/scripts/util.R")

    esvTable(pattern="*.dada.RDS")
    """    
}

/*
 *
 * Merge fastas (and keep sample id and abundance in name) for clustering
 *
 */

process MergeFastas {
    tag { "mergeFastas" }
    label "medium_computation"
    publishDir "${params.outdir}/8-Clustering", mode: "copy"
    errorStrategy "${params.errorsHandling}"
    
    input:
	file fastas from DADA_FASTA.collect()
    output:
	set val(100), file("all_samples_merged.fasta") into ESV_ALL_SAMPLES_TO_CLUSTER
    script:
    """
    #!/usr/bin/env python3

    import sys
    sys.path.append("${workflow.projectDir}/scripts")
    from util import mergeSamplesFa

    mergeSamplesFa()
    """    
}

/*
 *
 * VSEARCH OTU clustering at thresholds 95,97,99
 *
 */

process Clustering {
    tag { "clustering.${idThreshold}" }
    label "medium_computation"    
    publishDir "${params.outdir}/8-Clustering", mode: "copy"
    errorStrategy "${params.errorsHandling}"
    
    input:
	each idThreshold from (97)
        set val(esvId),file(esvFasta) from ESV_ALL_SAMPLES_TO_CLUSTER
    output:
	set val(idThreshold),file("otus${idThreshold}_seq.fasta") into OTU_ALL_SAMPLES,OTU_ALL_SAMPLES_LULU
        set val(idThreshold),file("otus${idThreshold}_table.csv") into ABUNDANCE_TABLES_OTU
    
    script:
    """
    vsearch --threads ${task.cpus} \
            --cluster_size ${esvFasta} \
            --id 0.${idThreshold} \
            --strand plus \
            --sizein \
            --sizeout \
            --fasta_width 0 \
            --uc clusters${idThreshold}.uc \
            --relabel OTU${idThreshold}_ \
            --centroids otus${idThreshold}_seq.fasta \
            --otutabout otus${idThreshold}_table.tsv

    cat otus${idThreshold}_table.tsv | tr "\\t" "," > otus${idThreshold}_table.csv
    """
}

// Mix channels with OTU and ESVs
ABUNDANCE_TABLES = ABUNDANCE_TABLES_OTU.mix(ABUNDANCE_TABLES_ESV)
LULU_ALL_SAMPLES = OTU_ALL_SAMPLES_LULU.mix(ESV_ALL_SAMPLES_LULU)
ALL_SAMPLES = OTU_ALL_SAMPLES.mix(ESV_ALL_SAMPLES)

/*
 *
 * Preliminary step before LULU: blast each contig against the others to look for sequence similarities
 *
 */

process PreLulu {
    tag { "preLulus.${idThreshold}" }
    label "medium_computation"    
    publishDir "${params.outdir}/9-LULU_correction", mode: "copy"
    errorStrategy "${params.errorsHandling}"

    input:
	set val(idThreshold),file(fasta) from LULU_ALL_SAMPLES
    output:
	set val(idThreshold),file("match_list_${idThreshold}.txt") into MATCH_LISTS
    script:
	
    """
    vsearch --usearch_global ${fasta} \
            --db ${fasta} --self \
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
 * LULU
 *
 */

process Lulu {
    tag { "Lulu.${idThreshold}" }
    label "medium_computation"    
    publishDir "${params.outdir}/9-LULU_correction", mode: "copy"
    errorStrategy "${params.errorsHandling}"

    input:
	set val(idThreshold),file(matchlist),file(table) from MATCH_LISTS.join(ABUNDANCE_TABLES)
    output:
	set val(idThreshold),file("curated_table_${idThreshold}.csv") into ABUNDANCE_LULU
        set val(idThreshold),file("curated_ids_${idThreshold}.csv") into IDS_LULU
        file("lulu*.log_*") optional true
    script:
	
    """
    #!/usr/bin/env Rscript
    source("${workflow.projectDir}/scripts/util.R")

    luluCurate("${table}","${matchlist}","${idThreshold}")
    """
}

/*
 *
 * Filter out FASTA sequences that LULU merged with a more abundant sequence
 *
 */

process ExtractFastaLulu {
    label "medium_computation"
    tag { "extractFastaLulu.${idThreshold}" }
    publishDir "${params.outdir}/10-taxonomy", mode: "copy"
    errorStrategy "${params.errorsHandling}"
    
    input:
	set idThreshold,file(ids),file(fasta) from IDS_LULU.join(ALL_SAMPLES)
    output:
	set idThreshold,file("lulu_fasta${idThreshold}.fasta") into FASTA_LULU, FASTA_LULU_FOR_SUMMARY
    script:
	
    """
    #!/usr/bin/env python3
    from Bio import SeqIO
    import pandas as pd

    ids = pd.read_csv("${ids}", header=None)[0].values
    sequences = [ seq for seq in SeqIO.parse("${fasta}","fasta") if seq.id.split(";")[0] in ids ]
    SeqIO.write(sequences,"lulu_fasta${idThreshold}.fasta","fasta")
    """
}    

/*
 *
 * Classification using CONSTAX. The algorithm was modified in order to take inputs from the commandline and generate outputs in a specified directory
 *
 */

process ClassificationCONSTAX  {
    tag { "classification.${idThreshold}" }
    label "medium_computation"
    publishDir "${params.outdir}/10-taxonomy", mode: "copy"
    errorStrategy "${params.errorsHandling}"

    input:
	set idThreshold,file(fasta) from FASTA_LULU
    output:
        set idThreshold,file("annotations_${idThreshold}/*") into TAXONOMY
    script:
	
    """

    bash ${workflow.projectDir}/friendly_CONSTAX/constax.sh ${fasta} "${workflow.projectDir}/friendly_CONSTAX/RDPTools" ${params.usearch8} ${params.usearch10}

    mv CONSTAX_outputs/outputs annotations_${idThreshold}

    """
}

/*
 *
 * Generates a summary file with the sequence count at each step for each sample (sample x step table) 
 *
 */

process SummaryFile {
    tag { "summary" }
    publishDir "${params.outdir}/11-Postprocessing", mode: "copy", pattern: "sequences_per_sample_per_step.tsv"
    errorStrategy "${params.errorsHandling}"
    label "medium_computation"
    
    input:
	file f from FASTA_LULU_FOR_SUMMARY.collect()
    output:
        file("sequences_per_sample_per_step.tsv")
    script:
    """
    #!/usr/bin/env python3
    import sys
    sys.path.append("${workflow.projectDir}/scripts")
    from generate_step_summary import write_summary

    write_summary("${params.outdir}","${params.reads}")
    """

}

