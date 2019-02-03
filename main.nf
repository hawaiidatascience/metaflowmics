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

if( params.revRead == 0 ){
    Channel
	.fromPath( params.reads )
	.ifEmpty { error "Cannot find any reads matching: ${params.reads}" }
    	.map { filename -> tuple(filename.baseName.substring(0, filename.baseName.length()-1),
				 file(filename))}
	.into { INPUT_FASTQ ; INPUT_FASTQ_TO_QC }

} else {
    Channel
	.fromFilePairs( params.reads )
	.ifEmpty { error "Cannot find any reads matching: ${params.reads}" }
	.into { INPUT_FASTQ ; INPUT_FASTQ_TO_QC }
}

// Header log info
def summary = [:]
summary['Run Name'] = workflow.runName
summary['Reads'] = params.reads
summary['Output dir'] = params.outdir
summary['Working dir'] = workflow.workDir
summary['Current home'] = "$HOME"
summary['Current user'] = "$USER"
summary['Current path'] = "$PWD"
summary['Config Profile'] = workflow.profile
log.info summary.collect { k,v -> "${k.padRight(15)}: $v" }.join("\n")
log.info "========================================="

/* Pipeline summary

- fastqc

****** dada2 ****** 
- filterAndTrim
- learnErrors
- dada (denoising)
******************

- Deunique (convert unique reads to redundant)

****** Mothur ******
- make.contigs (if paired-end)
- screen.seqs (optimization of minoverlap-mismatches-minlength-maxlength)
- unique.seqs
- align.seqs; filter.seqs; screen.seqs  --> alignment against reference for further filtering
- chimera.vsearch; remove.seqs             --> chimera removal
- classify.seqs; (remove.lineage)       --> classification, optional taxa filtering
- sub.sample
- cluster (at 95,97,99 and 100% identity)
- classify.otu
*******************

- cleanTables (postprocessing before LULU)
- preLulu (create matchlists for LULU)
- LULU
- FilterFasta (remove sequences from FASTA and taxonomy file that Lulu removed)

*/

// process runFastQC {
//     tag { "FastQC.${pairId}" }
//     publishDir "${params.outdir}/0-qualityControl", mode: "copy", overwrite: false
    
//     input:
//         set val(pairId), file(in_fastq) from INPUT_FASTQ_TO_QC
//     output:
//         file("${pairId}_fastqc/*.zip") into FASTQC_FILES

//     script:
//     """
//     mkdir ${pairId}_fastqc

//     fastqc --outdir ${pairId}_fastqc ${in_fastq.join(' ')}

//     """
// }

/*
 *
 * Step 0: Demultiplexing
 *
 */

process FilterAndTrim {
    // Quality filter and trimming using dada2. 
    tag { "FilterAndTrim.${pairId}" }
    publishDir "${params.outdir}/1-filterAndTrim", mode: "copy"
    
    input:
        set val(pairId), file(fastq) from INPUT_FASTQ
    output:
        set val(pairId), file("${pairId}*_trimmed.fastq") into FASTQ_TRIMMED, FASTQ_TRIMMED_FOR_MODEL
        set val(pairId), file("${pairId}*.ids") into FILTERED_READS_IDS
        file "*.pdf" into QUALITY_PROFILE

    script:
    """
    #!/usr/bin/env Rscript
    source("${workflow.projectDir}/scripts/util.R")  

    fastqs <- c("${fastq.join('","')}")

    filterReads("${pairId}", fastqs[1], rev=fastqs[2],
                minLen=${params.minLength}, 
		maxEE=${params.maxEE}, 
		truncLen=${params.truncLen}, 
		rm.phix=${params.rmphix}, 
		truncQ=${params.truncQ})
    """
}

process LearnErrors {
    // Build error model using dada2. 
    tag { "LearnErrors.${pairId}" }
    publishDir "${params.outdir}/2-errorModel", mode: "copy"
    label "medium_computation"

    input:
	set val(pairId), file(fastq) from FASTQ_TRIMMED_FOR_MODEL
    output:
	set val(pairId), file("${pairId}*.RDS") into ERROR_MODEL
        file("*.pdf") into ERROR_PROFILE
    

    script:
    """
    #!/usr/bin/env Rscript
    source("${workflow.projectDir}/scripts/util.R") 

    fastqs <- c("${fastq.join('","')}")
    learnErrorRates(fastqs,"${pairId}")
    """
}

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

process Esv {
    tag { "Esv" }
    publishDir "${params.outdir}/4-esv", mode: "copy"
    label "high_computation"
    
    input:
	file dadas from DADA_RDS.collect()
    output:
        set file("all.esv.count_table"), file("all.esv.fasta")  into DEREP_CONTIGS
        file("dada_merged.RDS")
        file("count_summary.tsv") into COUNT_SUMMARIES
    
    script:
    """
    #!/usr/bin/env Rscript
    source("${workflow.projectDir}/scripts/util.R")

    esvTable(${params.minOverlap},${params.maxMismatch},${params.revRead})       
    """
}

process MultipleSequenceAlignment {
    tag { "MSA" }
    publishDir "${params.outdir}/5-multipleSequenceAlignment", mode: "copy"
    label "high_computation"
    
    input:
        set file(count), file(fasta) from DEREP_CONTIGS
    output:
        set file("all.screening.start_end.count_table"), file("all.screening.start_end.fasta") into DEREP_CONTIGS_ALN
        file("all.screening.start_end.summary") into ALN_SUMMARY
    
    script:
    """
    ${params.scripts}/mothur.sh \
       --step=MSA \
       --rad=all.esv \
       --refAln=${params.referenceAln}

    ${params.scripts}/mothur.sh \
       --step=screening \
       --rad=all.MSA \
       --criteria=${params.criteria} \
       --optimize=start-end

    ${params.scripts}/mothur.sh --step=summary --rad=all.screening.start_end
    """
}

process ChimeraRemoval {
    tag { "chimeraRemoval" }
    publishDir "${params.outdir}/6-chimeraRemoval", mode: "copy"
    label "high_computation"
    
    input:
	set file(count), file(fasta) from DEREP_CONTIGS_ALN

    output:
        set file("all.chimera.count_table"), file("all.chimera.fasta") into NO_CHIMERA_FASTA
    
    script:
    """
    ${params.scripts}/mothur.sh --step=chimera --rad=all.screening.start_end
    """
}

process TaxaFiltering {
    tag { "taxaFilter" }
    publishDir "${params.outdir}/7-taxaFiltering", mode: "copy"
    
    input:
	set file(count), file(fasta) from NO_CHIMERA_FASTA

    output:
	set file("all.taxaFilter.fasta"), file("all.taxaFilter.count_table"), file("all.taxaFilter.taxonomy") into TAXA_FILTERED_CONTIGS

    script:
    """
    ${params.scripts}/mothur.sh \
	--step=taxaFilter \
	--rad=all.chimera \
        --taxaToFilter=${params.taxaToFilter} \
	--refAln=${params.referenceAln} \
	--refTax=${params.referenceTax}

    if [ -z ${params.taxaToFilter} ]; then
        mv ${count} all.taxaFilter.count_table
        mv ${fasta} all.taxaFilter.fasta
    fi

    """
}

process Subsampling {
    tag { "subsampling" }
    publishDir "${params.outdir}/8-subsampling", mode: "copy"
    
    input:
	set file(fasta), file(count), file(tax) from TAXA_FILTERED_CONTIGS
    output:
	set file("all.subsampling.fasta"), file("all.subsampling.count_table"), file("all.subsampling.taxonomy") into SUBSAMPLED_CONTIGS

    script:
    """

    awk '{for (i=3;i<=NF;i++) sum[i]+=\$i;}; END{for (i in sum) print sum[i]}' ${count} | tail -n +2 |sort -n > sample_size.txt

    percentile_value=`awk '{all[NR] = \$1} END{print all[int(NR*${params.subsamplingQuantile})]}' sample_size.txt`

    if [ \$percentile_value < ${params.minSubsampling} ] || [ -z \$percentile_value ]
        then percentile_value=${params.minSubsampling}
    fi

    ${params.scripts}/mothur.sh --step=subsampling --rad=all.taxaFilter --subsamplingNb=\$percentile_value
    """
}

process Clustering {
    tag { "clustering" }
    publishDir "${params.outdir}/9-clustering", mode: "copy", pattern: "*.{fasta,shared,list}"
    label "high_computation"
    
    input:
	set file(fasta), file(count), file(tax) from SUBSAMPLED_CONTIGS
        each idThreshold from (0,0.01,0.03,0.05)
    output:
        set val(idThreshold), file("all.clustering.*.fasta") into PRELULU_FASTA, FASTA_TO_FILTER
        set val(idThreshold), file("all.clustering.*.shared") into ABUNDANCE_TABLES
	set val(idThreshold), file("all.clustering.*.list"), file(count), file(tax) into CONTIGS_FOR_CLASSIFICATION
    script:
    """
    ${params.scripts}/mothur.sh --step=clustering --rad=all.subsampling --idThreshold=${idThreshold}
    """
}

process ConsensusClassification {
    tag { "consensusClassification" }
    publishDir "${params.outdir}/10-consensusClassification", mode: "copy"
    
    input:
	set val(idThreshold), file(list), file(count), file(tax) from CONTIGS_FOR_CLASSIFICATION
    output:
	file("all.classification.*.summary") into CLASSIFICATION_SUMMARY
	set val(idThreshold), file("all.classification.*.taxonomy") into CONSENSUS_TAXONOMY
    script:
    """
    rad=`basename all.clustering.*.list .list`
    ${params.scripts}/mothur.sh --step=consensusClassification --rad=\${rad} --idThreshold=${idThreshold} --refTax=${params.referenceTax} --refAln=${params.referenceAln}
    """
}

/*
 *
 * Pre-lulu step
 *    - Blast each contig against each other
 *
 */

process PreLulu {
    tag { "preLulus" }
    publishDir "${params.outdir}/11-lulu", mode: "copy"

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

process Lulu {
    tag { "Lulu" }
    publishDir "${params.outdir}/11-lulu", mode: "copy"
    errorStrategy "${params.errorsHandling}"

    input:
	set val(idThreshold),file(matchlist),file(table) from MATCH_LISTS.join(ABUNDANCE_TABLES)
    output:
	set val(idThreshold),file("curated_table_${idThreshold}.csv") into ABUNDANCE_LULU,ABUNDANCE_MOTHUR
	set val(idThreshold),file("curated_ids_${idThreshold}.csv") into IDS_LULU
    script:
	
    """
    #!/usr/bin/env Rscript
    source("${workflow.projectDir}/scripts/util.R")

    luluCurate("${table}","${matchlist}","${idThreshold}")
    """
}

process FilterFasta {
    tag { "filterFasta" }
    publishDir "${params.outdir}/11-lulu", mode: "copy"
    errorStrategy "${params.errorsHandling}"
    
    input:
	set idThreshold,file(ids),file(tax),file(fasta) from IDS_LULU.join(CONSENSUS_TAXONOMY).join(FASTA_TO_FILTER)
    output:
	set idThreshold,file("OTU_${idThreshold}.fasta") into FASTA_FOR_MOTHUR
        set file("OTU_${idThreshold}.fasta"), file("taxonomy_${idThreshold}.csv") into OUTPUT_FILES
    script:
    """
    #!/usr/bin/env python

    import pandas as pd
    from Bio import SeqIO

    ids = [ name.strip() for name in open("${ids}","r").readlines() ]
    fasta = [ seq for seq in SeqIO.parse("${fasta}","fasta") if seq.id in ids ]
    SeqIO.write(fasta,"OTU_${idThreshold}.fasta","fasta")

    taxonomyTable = pd.read_table("${tax}",index_col=0)
    taxonomyTable.loc[ids].to_csv("taxonomy_${idThreshold}.csv")
    """
}

process ConvertToMothur {
    tag { "convertToMothur" }
    publishDir "${params.outdir}/12-mothurFmtOutputs", mode: "copy"
    errorStrategy "${params.errorsHandling}"
    
    input:
	set val(idThreshold),file(fasta),file(abundanceTable) from FASTA_FOR_MOTHUR.join(ABUNDANCE_MOTHUR)
    output:
	set val(idThreshold),file(fasta),file("abundance_${idThreshold}.shared") into MOTHUR_INPUTS
    script:
    """
    #!/usr/bin/env python

    import pandas as pd

    shared = pd.read_csv("${abundanceTable}", index_col=0).T
    shared.index.name = "Group"
    shared.insert(loc=0, column="numOtus", value=[shared.shape[1]]*shared.shape[0])
    shared.reset_index(inplace=True)
    shared.insert(loc=0, column="label", value=[${idThreshold}]*shared.shape[0])
    shared.to_csv("abundance_${idThreshold}.shared", sep="\t", index=False)
    """
}

process Results {
    tag { "mothurResults" }
    publishDir "${params.outdir}/13-Postprocessing", mode: "copy"
    errorStrategy "${params.errorsHandling}"
    
    input:
	set val(idThreshold), file(fasta), file(shared) from MOTHUR_INPUTS
    output:
	set file("*.relabund"), file("*.wsummary"), file("*.tre") into RESULTS
    script:
    """
    mothur "#get.relabund(shared=${shared})"
    mothur "#clearcut(fasta=${fasta}, DNA=T) ; count.seqs(shared=${shared}) ; unifrac.weighted(tree=current,count=current)"
    """    
}

process SummaryFile {
    tag { "mothurResults" }
    publishDir "${params.outdir}/13-Postprocessing", mode: "copy"
    errorStrategy "${params.errorsHandling}"
    
    input:
	file f1 from COUNT_SUMMARIES
	file f2 from RESULTS.collect()
        file f3 from OUTPUT_FILES.collect()    
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
