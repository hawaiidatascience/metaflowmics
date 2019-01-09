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
- chimera.seqs; remove.seqs             --> chimera removal
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

process runFastQC {
    tag { "FastQC.${pairId}" }
    publishDir "${params.outdir}/0-qualityControl", mode: "copy", overwrite: false
    
    input:
        set val(pairId), file(in_fastq) from INPUT_FASTQ_TO_QC
    output:
        file("${pairId}_fastqc/*.zip") into FASTQC_FILES

    script:
    """
    mkdir ${pairId}_fastqc

    fastqc --outdir ${pairId}_fastqc ${in_fastq.join(' ')}

    """
}

/*
 *
 * Step 0: Demultiplexing
 *
 */

process FilterAndTrim {
    // Quality filter and trimming using dada2. 
    tag { "FilterAndTrim.${pairId}" }
    publishDir "${params.outdir}/1-filterAndTrim", mode: "copy", overwrite: false
    
    input:
        set val(pairId), file(fastq) from INPUT_FASTQ
    output:
        set val(pairId), file("${pairId}*_trimmed.fastq") into FASTQ_TRIMMED, FASTQ_TRIMMED_FOR_MODEL
        set val(pairId), file("${pairId}*.ids") into FILTERED_READS_IDS

    script:
    """
    #!/usr/bin/env Rscript
    source("${workflow.projectDir}/scripts/util.R")  

    fastqs <- c("${fastq.join('","')}")

    filterReads("${pairId}", "${fastq.join('","')}", 
                minLen=${params.minLength}, maxEE=${params.maxEE})
    """
}

process LearnErrors {
    // Build error model using dada2. 
    tag { "LearnErrors.${pairId}" }
    publishDir "${params.outdir}/2-errorModel", mode: "copy", overwrite: false

    input:
	set val(pairId), file(fastq) from FASTQ_TRIMMED_FOR_MODEL
    output:
	set val(pairId), file("${pairId}*.RDS") into ERROR_MODEL

    script:
    """
    #!/usr/bin/env Rscript
    source("${workflow.projectDir}/scripts/util.R") 

    fastqs <- c("${fastq.join('","')}")
    learnErrorRates(fastqs[1],"${pairId}1")

    if (${params.revRead} == 1) {
        learnErrorRates(fastqs[2],"${pairId}2")
    }
    """
}

process Denoise {
    tag { "Denoising.${pairId}" }
    publishDir "${params.outdir}/3-denoising", mode: "copy", overwrite: false
    
    input:
        set val(pairId), file(err), file(fastq) from ERROR_MODEL.join(FASTQ_TRIMMED)
    output:
        set val(pairId), file("${pairId}*.derep.RDS"), file("${pairId}*.dada.RDS") into DADA_RDS
    script:
    """
    #!/usr/bin/env Rscript
    source("${workflow.projectDir}/scripts/util.R")

    errors <- c("${err.join('","')}")
    fastqs <- c("${fastq.join('","')}")

    dadaDenoise(errors[1], fastqs[1], "${pairId}1")

    if (${params.revRead} == 1) {
        dadaDenoise(errors[2], fastqs[2], "${pairId}2")
    }
    """
}

process Deunique {
    tag { "Deunique.${pairId}" }
    publishDir "${params.outdir}/4-denoisedFasta", mode: "copy", overwrite: false
    
    input:
        set val(pairId), file(derep_rds), file(dada_rds), file(ids) from DADA_RDS.join(FILTERED_READS_IDS)
    
    output:
        set val(pairId),file("${pairId}*.denoised.fasta") into DADA_FASTA
    script:
    """
    #!/usr/bin/env Rscript
    source("${workflow.projectDir}/scripts/util.R")  

    derep_rds <- c("${derep_rds.join('","')}")
    dada_rds <- c("${dada_rds.join('","')}")
    ids <- c("${ids.join('","')}")

    readsFromDenoised(dada_rds[1],derep_rds[1],"${pairId}1",ids[1])

    if (${params.revRead} == 1) {
        readsFromDenoised(dada_rds[2],derep_rds[2],"${pairId}2",ids[2])
    }
    """
}

process ReadsMerging {

    tag { "ContigsMerging.${pairId}" }
    publishDir "${params.outdir}/5-readsMerging", mode: "copy", overwrite: false

    input:
        set val(pairId), file(fasta) from DADA_FASTA
    output:
	set val(pairId), file("${pairId}.merging.fasta") into CONTIGS_FASTA
	
    script:
    """
    fastas=(${fasta.join(' ')})

    if [ ${params.revRead} -eq 1 ]; then
       ${params.scripts}/mothur.sh --step=merging --pairId=${pairId} --fwdFasta=\${fastas[0]} --revFasta=\${fastas[1]} ;
    else 
       # no need for merging
       cp ${fasta} ${pairId}.merging.fasta
    fi
    """
}

/*
 *
 * Outlier removal:
 *   1) Overlap and mismatches optimization at 95%
 *   2) Maxlength optimization at 95%, minLength, maxAmbig, maxHomop
 *
 */

process Screening {

    tag { "screening.${pairId}" }
    publishDir "${params.outdir}/6-screening", mode: "copy", overwrite: false
    
    input:
        set val(pairId), file(fasta) from CONTIGS_FASTA
    output:
	set val(pairId), file("${pairId}.screening.minoverlap_mismatches_minlength_maxlength.fasta") into FILTERED_CONTIGS
        file("${pairId}.screening.minoverlap_mismatches_minlength_maxlength.summary") into SUMMARIES

    script:
    """
    ${params.scripts}/mothur.sh \
       --step=screening \
       --pairId=${pairId} \
       --rad=${pairId}.merging \
       --optimize=minoverlap-mismatches-minlength-maxlength \
       --criteria=${params.criteria}

    ${params.scripts}/mothur.sh --step=summary --rad=${pairId}.screening.minoverlap_mismatches_minlength_maxlength
    """
}


process Dereplication {
    tag { "dereplication" }
    publishDir "${params.outdir}/7-dereplication", mode: "copy", overwrite: false
    
    input:
        file(fasta) from FILTERED_CONTIGS.collect()
    output:
	set file("all.dereplication.names"), file("all.dereplication.fasta"), file("all.dereplication.groups") into DEREP_CONTIGS

    script:
    """
    ${params.scripts}/mothur.sh --step=dereplication

    # different .groups name depending on the number of samples
    mv *.groups all.dereplication.groups

    """
}

process MultipleSequenceAlignment {
    tag { "MSA" }
    publishDir "${params.outdir}/8-multipleSequenceALignment", mode: "copy", overwrite: false
    
    input:
        set file(names), file(fasta), file(groups) from DEREP_CONTIGS
    output:
        set file("all.screening.start_end.names"), file("all.screening.start_end.fasta"), file("all.screening.start_end.groups") into DEREP_CONTIGS_ALN
        file("all.screening.start_end.summary") into ALN_SUMMARY
    
    script:
    """
    ${params.scripts}/mothur.sh \
       --step=MSA \
       --rad=all.dereplication \
       --refAln=${params.referenceAln}

    ${params.scripts}/mothur.sh \
       --step=screening \
       --rad=all.dereplication \
       --criteria=${params.criteria} \
       --optimize=start-end

    ${params.scripts}/mothur.sh --step=summary --rad=all.screening.start_end
    """
}

process ChimeraRemoval {
    tag { "chimeraRemoval" }
    publishDir "${params.outdir}/9-chimeraRemoval", mode: "copy", overwrite: false
    
    input:
	set file(names), file(fasta), file(groups) from DEREP_CONTIGS_ALN

    output:
        set file("all.chimera.names"), file("all.chimera.fasta"), file("all.chimera.groups") into NO_CHIMERA_FASTA
    
    script:
    """
    ${params.scripts}/mothur.sh --step=chimera --rad=all.screening.start_end
    """
}

process TaxaFiltering {
    tag { "taxaFilter" }
    publishDir "${params.outdir}/10-taxaFiltering", mode: "copy", overwrite: false
    
    input:
	set file(names), file(fasta), file(groups) from NO_CHIMERA_FASTA

    output:
	set file("all.taxaFilter.fasta"), file("all.taxaFilter.names"), file("all.taxaFilter.groups"), file("all.taxaFilter.taxonomy") into TAXA_FILTERED_CONTIGS

    script:
    """
    ${params.scripts}/mothur.sh \
	--step=taxaFilter \
	--rad=all.chimera \
	--refAln=${params.referenceAln} \
	--refTax=${params.referenceTax}

    if [ ${params.taxaFilter} -eq 0 ]; then
        mv ${names} all.taxaFilter.names
        mv ${fasta} all.taxaFilter.fasta
        mv ${groups} all.taxaFilter.groups
    fi

    """
}

process Subsampling {
    tag { "subsampling" }
    publishDir "${params.outdir}/11-subsampling", mode: "copy", overwrite: false
    
    input:
	set file(fasta), file(names), file(groups), file(tax) from TAXA_FILTERED_CONTIGS
    output:
	set file("all.subsampling.fasta"), file("all.subsampling.names"), file("all.subsampling.groups"), file("all.subsampling.taxonomy") into SUBSAMPLED_CONTIGS

    script:
    """
    ${params.scripts}/mothur.sh --step=subsampling --rad=all.taxaFilter
    """
}

process Clustering {
    tag { "clustering" }
    publishDir "${params.outdir}/12-clustering", mode: "copy", overwrite: false
    
    input:
	set file(fasta), file(names), file(groups), file(tax) from SUBSAMPLED_CONTIGS
        each idThreshold from (0,0.01,0.03,0.05)
    output:
	set val(idThreshold),file("all.clustering.*.list"), file(tax) into CONTIGS_FOR_CLASSIFICATION
        set val(idThreshold), file("all.clustering.*.shared"), file("all.clustering.*.fasta") into ABUNDANCE_TABLES
        set val(idThreshold), file("all.clustering.*.fasta") into PRELULU_FASTA, FASTA_TO_FILTER

    script:
    """
    ${params.scripts}/mothur.sh --step=clustering --rad=all.subsampling --idThreshold=${idThreshold}
    """
}

process PreClassification {
    tag { "pre-classification" }
    publishDir "${params.outdir}/13-preclassification", mode: "copy", overwrite: false
    
    input:
	set val(idThreshold), file(list), file(tax) from CONTIGS_FOR_CLASSIFICATION
    output:
	set val(idThreshold), file("all.classification.*.taxonomy"), file("all.classification.*.summary") into CLASSIFIED_CONTIGS
    script:
    """
    rad=`basename all.clustering.*.list .list`
    ${params.scripts}/mothur.sh --step=classification --rad=\${rad} --idThreshold=${idThreshold} --refTax=${params.referenceTax} --refAln=${params.referenceAln}
    """
}

process CleanTables {
    tag "cleanTable"
    publishDir "${params.outdir}/14-cleanTables", mode: "copy", overwrite: false

    input:
	set val(idThreshold), file(table), file(fasta), file(taxonomy), file(taxSummary) from ABUNDANCE_TABLES.join(CLASSIFIED_CONTIGS)
    output:
	set val(idThreshold), file("all.abundanceTable.*.csv"), file("all.taxonomyTable.*.csv") into ABUNDANCE_TABLES_CLEAN
        set val(idThreshold), file("all.taxonomyTable.*.csv") into TAXONOMY_TO_FILTER
    script:
    """
    #!/usr/bin/env python
 
    from Bio import SeqIO
    import pandas as pd
    import re

    taxonomyTable = pd.read_table("${taxonomy}", index_col=0)
    abundanceTable = pd.read_table("${table}", index_col=1).drop(["label","numOtus"],axis=1)
    mapping = pd.Series({ 
                      re.split('[\\t|\\|]',seq.description)[1]: seq.id
                      for seq in SeqIO.parse("${fasta}","fasta") 
                  })

    abundanceTable.columns = mapping[abundanceTable.columns]
    taxonomyTable.index = mapping[taxonomyTable.index]

    abundanceTable.T.to_csv("all.abundanceTable.${idThreshold}.csv")
    taxonomyTable.to_csv("all.taxonomyTable.${idThreshold}.csv")
    
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
    publishDir "${params.outdir}/15-lulu", mode: "copy", overwrite: false

    input:
	set val(idThreshold),file(fasta) from PRELULU_FASTA
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

process Lulu {
    tag { "Lulu" }
    publishDir "${params.outdir}/15-lulu", mode: "copy", overwrite: false
    errorStrategy "${params.errorsHandling}"

    input:
	set val(idThreshold),file(matchlist),file(table),file(taxTable) from MATCH_LISTS.join(ABUNDANCE_TABLES_CLEAN)
    output:
	set val(idThreshold),file("curated_table_${idThreshold}.csv") into ABUNDANCE_LULU
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
    publishDir "${params.outdir}/15-lulu", mode: "copy", overwrite: false
    errorStrategy "${params.errorsHandling}"
    
    input:
	set idThreshold,file(ids),file(tax),file(fasta) from IDS_LULU.join(TAXONOMY_TO_FILTER).join(FASTA_TO_FILTER)
    output:
	set file("OTU_${idThreshold}.fasta"), file("taxonomy_${idThreshold}.csv")
    script:
    """
    #!/usr/bin/env python3

    import pandas as pd
    from Bio import SeqIO

    ids = [ name.strip() for name in open("${ids}","r").readlines() ]
    fasta = [ seq for seq in SeqIO.parse("${fasta}","fasta") if seq.id in ids ]
    SeqIO.write(fasta,"OTU_${idThreshold}.fasta","fasta")

    taxonomyTable = pd.read_csv("${tax}",index_col=0)
    taxonomyTable.loc[ids].to_csv("taxonomy_${idThreshold}.csv")    
    """
}

