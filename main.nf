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

- make.contigs
- make.fastq                            --> we'll need them later for dada2
- screen.seqs (x2)                      --> quality control
- unique.seqs; count.seqs               --> dereplication
- align.seqs; filter.seqs; screen.seqs  --> alignment against reference for further filtering
- chimera.seqs; remove.seqs             --> chimera removal
- classify.seqs; remove.lineage         --> remove custom taxa (ex: unknown)
- list.seqs; get.seqs                   --> conversion to dada2 format
- dada2 denoising
- abundance table
- lulu             --> co-occurence pattern correction
- taxonomy

Not implemented yet:
- subsampling      --> depth normalization across samples

*/



/*
 *
 * Step 0: FastQC 
 *
 */

process runFastQC {
    tag { "FastQC.${pairId}" }
    publishDir "${params.outdir}/qualityControl", mode: "copy", overwrite: false
    errorStrategy "${params.errorsHandling}"
    
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

/*
 *
 * Step 1: Filter and Trim
 *
 */

process FilterAndTrim {
    // Quality filter and trimming using dada2. 
    tag { "FilterAndTrim.${pairId}" }
    publishDir "${params.outdir}/FilterAndTrim", mode: "copy", overwrite: false
    errorStrategy "${params.errorsHandling}"
    
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

/*
 *
 * Step 2: Error model
 *
 */

process LearnErrors {
    // Build error model using dada2. 
    tag { "LearnErrors.${pairId}" }
    publishDir "${params.outdir}/ErrorModel", mode: "copy", overwrite: false
    errorStrategy "${params.errorsHandling}"

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

/*
 *
 * Step 3: Denoising
 *
 */

process Denoise {
    tag { "Denoising.${pairId}" }
    publishDir "${params.outdir}/Denoising", mode: "copy", overwrite: false
    errorStrategy "${params.errorsHandling}"
    
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

/*
 *
 * Step 4: DeUnique reads
 *
 */

process Deunique {
    tag { "Deunique.${pairId}" }
    publishDir "${params.outdir}/denoisedFasta", mode: "copy", overwrite: false
    errorStrategy "${params.errorsHandling}"
    
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

/*
 *
 * Step 5: Merge reads
 *
 */

process ReadsMerging {

    tag { "ContigsMerging.${pairId}" }
    publishDir "${params.outdir}/ReadsMerging", mode: "copy", overwrite: false
    errorStrategy "${params.errorsHandling}"

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
       # skip this step
       cp ${fasta} ${pairId}.merging.fasta
    fi
    """
}

/*
 *
 * Step 6: Outlier removal:
 *   1) Overlap and mismatches optimization at 95%
 *   2) Maxlength optimization at 95%, minLength, maxAmbig, maxHomop
 *
 */

process OutliersRemoval {

    tag { "outlierRemoval.${pairId}" }
    publishDir "${params.outdir}/OutlierRemoval", mode: "copy", overwrite: false
    errorStrategy "${params.errorsHandling}"
    
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
       --inputFasta=${fasta} \
       --optimize=minoverlap-mismatches-minlength-maxlength \
       --criteria=${params.criteria}

    ${params.scripts}/mothur.sh --step=summary --pairId=${pairId} --inputFasta=${pairId}.screening.minoverlap_mismatches_minlength_maxlength.fasta
    """
}


process Dereplication {
    tag { "dereplication" }
    publishDir "${params.outdir}/OutlierRemoval", mode: "copy", overwrite: false
    errorStrategy "${params.errorsHandling}"
    
    input:
        file f from FILTERED_CONTIGS.collect()
    output:
	set file("all.dereplication.groups"), file("all.dereplication.count_table"), file("all.dereplication.names"), file("all.dereplication.fasta") \
    into DEREP_CONTIGS_MSA

    script:
    """
    fastas=\$(echo \$(ls *.fasta) | sed 's/\\ /-/g')
    groups=`ls -l *.fasta | awk '{print \$9}' | awk 'BEGIN{FS="."} {\$NF=FILENAME; print \$1}' | paste -sd "-" -`
    cat *.fasta > samples_merged.fasta
    mothur "#make.group(fasta=\$fastas, groups=\$groups)"
    mv groups all.dereplication.groups

    ${params.scripts}/mothur.sh --step=dereplication  --inputFasta=samples_merged.fasta --inputGroups=all.dereplication.groups
    """
}

/*
 *
 * Step 8: Multiple sequence alignment
 *   and filter out of bad alignments 
 *   using a start-end optimization at 95%
 *
 */

process MultipleSequenceAlignment {
    tag { "MSA" }
    publishDir "${params.outdir}/MultipleSequenceALignment", mode: "copy", overwrite: false
    errorStrategy "${params.errorsHandling}"
    
    input:
        set file(groups), file(counts), file(names), file(fasta) from DEREP_CONTIGS_MSA
    output:
        set file("all.screening.start_end.groups"), file("all.screening.start_end.names"), file("all.screening.start_end.fasta") into DEREP_CONTIGS_ALN
        file("all.screening.start_end.summary") into ALN_SUMMARY
    
    script:
    """
    ${params.scripts}/mothur.sh \
       --step=MSA \
       --inputFasta=${fasta} \
       --refAln=${params.referenceAln}

    ${params.scripts}/mothur.sh \
       --step=screening \
       --inputFasta=all.MSA.fasta \
       --inputNames=${names} \
       --inputGroups=${groups} \
       --criteria=${params.criteria} \
       --optimize=start-end

    ${params.scripts}/mothur.sh --step=summary --inputFasta=all.screening.start_end.fasta
    """
}

/*
 *
 * Step 9: Chimera removal
 *
 */


process ChimeraRemoval {
    tag { "chimeraRemoval" }
    publishDir "${params.outdir}/chimeraRemoval", mode: "copy", overwrite: false
    errorStrategy "${params.errorsHandling}"
    
    input:
	set file(groups), file(names), file(fasta) from DEREP_CONTIGS_ALN

    output:
        set file("all.chimera.groups"), file("all.chimera.names"), file("all.chimera.fasta") \
        into NO_CHIMERA_FASTA
    
    script:
    """
    ${params.scripts}/mothur.sh --step=chimera --inputFasta=${fasta} --inputNames=${names} --inputGroups=${groups} 
    """
}

/*
 *
 * Step 10 (optional?): Taxa filtering
 *   user-defined taxa list to exclude
 *   (parameter to add) 
 *
 */


process TaxaFiltering {
    tag { "taxaFilter" }
    publishDir "${params.outdir}/taxaFiltering", mode: "copy", overwrite: false
    errorStrategy "${params.errorsHandling}"
    
    input:
	set file(groups), file(names), file(fasta) \
    from NO_CHIMERA_FASTA

    output:
	set file("all.taxaFilter.names"), file("all.taxaFilter.fasta"), file("all.taxaFilter.taxonomy"), file("all.taxaFilter.groups") \
    into TAXA_FILTERED_CONTIGS

    when:
	params.taxaFilter == 1
    
    script:
    """
    ${params.scripts}/mothur.sh \
        --step=taxaFilter \
        --inputFasta=${fasta} \
        --inputNames=${names} \
        --inputGroups=${groups} \
        --refAln=${params.referenceAln} \
        --refTax=${params.referenceTax}    
    """
}


process Subsampling {
    tag { "subsampling" }
    publishDir "${params.outdir}/subsampling", mode: "copy", overwrite: false
    errorStrategy "${params.errorsHandling}"
    
    input:
	set file(names),file(fasta),file(tax),file(groups) from TAXA_FILTERED_CONTIGS
    output:
	set file("all.subsampling.names"), file("all.subsampling.fasta"), file("all.subsampling.taxonomy"), file("all.subsampling.groups") into SUBSAMPLED_CONTIGS

    script:
    """
    ${params.scripts}/mothur.sh --step=subsampling --inputFasta=${fasta} --inputNames=${names} --inputGroups=${groups} --inputTax=${tax}
    """
}


process Clustering {
    tag { "clustering" }
    publishDir "${params.outdir}/clustering", mode: "copy", overwrite: false
    errorStrategy "${params.errorsHandling}"
    
    input:
	set file(names), file(fasta), file(tax), file(groups) from SUBSAMPLED_CONTIGS
        each idThreshold from (0,0.01,0.03,0.05)
    output:
	set val(idThreshold), file("all.clustering.${idThreshold}.names"), file("all.clustering.${idThreshold}.fasta"), file("all.clustering.${idThreshold}.list"), file(tax), file(groups) \
    into CLUSTERED_CONTIGS

    script:
    """
    if [ ${idThreshold} -eq 0 ]; then
        ${params.scripts}/mothur.sh --step=clustering --otuMethod=unique --inputFasta=${fasta} --inputNames=${names} --clustId=${idThreshold}
    else
        ${params.scripts}/mothur.sh --step=clustering --otuMethod=dgc --inputFasta=${fasta} --inputNames=${names} --clustId=${idThreshold}
    fi
    """
}

process Classification {
    tag { "classification" }
    publishDir "${params.outdir}/classification", mode: "copy", overwrite: false
    errorStrategy "${params.errorsHandling}"
    
    input:
	set val(idThreshold), file(names), file(fasta), file(list), file(tax), file(groups) from CLUSTERED_CONTIGS
    output:
	set val(idThreshold), file(names), file(fasta), file(list), file("all.classification.${idThreshold}.taxonomy"), file(groups) \
    into CLASSIFIED_CONTIGS
        file("all.classification.${idThreshold}.summary") into CLASSIFICATION_SUMMARY
        set val(idThreshold), file("all.classification.${idThreshold}.shared") into ABUNDANCE_TABLE

    script:
    """
    ${params.scripts}/mothur.sh --step=classification --inputNames=${names} --inputTax=${tax} --inputList=${list} --clustId=${idThreshold} --inputGroups=${groups}
    """
}

/*
 *
 * Step 15: Pre-lulu step
 *   - Blast each contig against each other
 *
 */

// process PreLulu {
//     tag { "preLulus" }
//     label "medium_computation"    
//     publishDir "${params.outdir}/lulu", mode: "copy", overwrite: false
//     errorStrategy "${params.errorsHandling}"

//     input:
// 	set val(idThreshold),file(fasta) from FASTA_FOR_LULU
//     output:
// 	set val(idThreshold),file("match_list_${idThreshold}.txt") into MATCH_LISTS
//     script:
	
//     """
//     vsearch --usearch_global ${fasta} \
//             --db ${fasta} --self \
//             --id .84 \
//             --iddef 1 \
//             --userout match_list_${idThreshold}.txt \
//             -userfields query+target+id \
//             --maxaccepts 0 \
//             --query_cov .9 \
//             --maxhits 10
//     """
// }

// /*
//  *
//  * Step 16: Lulu
//  *
//  */

// process Lulu {
//     tag { "Lulu" }
//     label "high_computation"    
//     publishDir "${params.outdir}/lulu", mode: "copy", overwrite: false
//     errorStrategy "${params.errorsHandling}"

//     input:
// 	set val(idThreshold),file(matchlist),file(table) from MATCH_LISTS.join(ABUNDANCE_TABLES)
//     output:
// 	set val(idThreshold),file("curated_table_${idThreshold}.csv") into ABUNDANCE_LULU
// 	set val(idThreshold),file("curated_ids_${idThreshold}.csv") into IDS_LULU
//     script:
	
//     """
//     #!/usr/bin/env Rscript
//     source("${workflow.projectDir}/scripts/util.R")

//     luluCurate("${table}","${matchlist}","${idThreshold}")
//     """
// }

// /*
//  *
//  * Step 17: Retrieve contigs from lulu table for annotation
//  *
//  */

// process ExtractFastaLulu {
//     label "medium_computation"
//     tag { "extractFastaLulu" }
//     publishDir "${params.outdir}/taxonomy", mode: "copy", overwrite: false
//     errorStrategy "${params.errorsHandling}"
    
//     input:
// 	set idThreshold,file(ids),file(fasta) from IDS_LULU.join(ALL_SAMPLES)
//     output:
// 	set idThreshold,file("lulu_fasta${idThreshold}.fasta") into FASTAS_LULU
//     script:
	
//     """
//     #!/usr/bin/env python3

//     import sys
//     sys.path.append("${workflow.projectDir}/scripts")
//     from util import extractFastaLulu

//     extractFastaLulu("${fasta}","${ids}","${idThreshold}")
//     """
// }    

// /*
//  *
//  * Step 19: Annotate abunance tables from LULU
//  *
//  */

// process FillAbundanceTable {
//     tag { "fillAbundanceTable" }
//     label "medium_computation"    
//     publishDir "${params.outdir}/taxonomy", mode: "copy", overwrite: false
//     errorStrategy "${params.errorsHandling}"
    
//     input:
// 	set idThreshold,file(table),file(tax) from ABUNDANCE_LULU.join(TAXONOMY)
//     output:
// 	file("abundance_table_annotated_ID=${idThreshold}.tsv") into TAXONOMY_ABUNDANCE
//     script:
//     """
//     #!/usr/bin/env python3

//     import sys
//     sys.path.append("${workflow.projectDir}/scripts")
//     from util import fillAbundances

//     fillAbundances("${table}","${tax}",${idThreshold})
//     """
// }

