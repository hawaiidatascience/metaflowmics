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

// Has the run name been specified by the user?
custom_runName = params.name
if( !(workflow.runName ==~ /[a-z]+_[a-z]+/) ){
  custom_runName = workflow.runName
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
summary['Run Name'] = custom_runName ?: workflow.runName
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
    tag { "rFQC.${pairId}" }
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
    tag { "FilterAndTrim" }
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
	set val(pairId), file("${pairId}.filter.minoverlap_mismatches_minlength_maxlength.fasta") into FILTERED_CONTIGS
        file("${pairId}.filter.minoverlap_mismatches_minlength_maxlength.summary") into SUMMARIES

    script:
    """
    ${params.scripts}/mothur.sh \
       --step=filter \
       --pairId=${pairId} \
       --inputFasta=${fasta} \
       --optimize=minoverlap-mismatches-minlength-maxlength \
       --criteria=${params.criteria}

    ${params.scripts}/mothur.sh --step=summary --pairId=${pairId} --inputFasta=${pairId}.filter.minoverlap_mismatches_minlength_maxlength.fasta
    """
}

/*
 *
 * Step 7: Dereplication
 *
 */

process Dereplication {
    tag { "dereplication.${pairId}" }
    publishDir "${params.outdir}/OutlierRemoval", mode: "copy", overwrite: false
    errorStrategy "${params.errorsHandling}"
    
    input:
        set val(pairId), file(fasta) from FILTERED_CONTIGS
    output:
	set val(pairId), file("${pairId}.dereplication.names"), file("${pairId}.dereplication.fasta") into DEREP_CONTIGS_MSA

    script:
    """
    ${params.scripts}/mothur.sh --step=dereplication --pairId=${pairId} --inputFasta=${fasta}
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
        set val(pairId), file(names), file(fasta) from DEREP_CONTIGS_MSA
    output:
        set val(pairId), file("${pairId}.filter.start_end.names"), file("${pairId}.filter.start_end.fasta") into DEREP_CONTIGS_ALN
        file("${pairId}.filter.start_end.summary") into ALN_SUMMARY
    
    script:
    """
    ${params.scripts}/mothur.sh \
       --step=MSA \
       --pairId=${pairId} \
       --inputFasta=${fasta} \
       --refAln=${params.referenceAln}

    ${params.scripts}/mothur.sh \
       --step=filter \
       --pairId=${pairId} \
       --inputFasta=${fasta} \
       --inputNames=${names} \
       --criteria=${params.criteria} \
       --optimize=start-end

    ${params.scripts}/mothur.sh --step=summary --pairId=${pairId} --inputFasta=${pairId}.filter.start_end.fasta
    """
}

/*
 *
 * Step 9: Chimera removal
 *
 */


process ChimeraRemoval {
    tag { "chimeraRemoval.${pairId}" }
    publishDir "${params.outdir}/chimeraRemoval", mode: "copy", overwrite: false
    errorStrategy "${params.errorsHandling}"
    
    input:
	set val(pairId), file(names), file(fasta) from DEREP_CONTIGS_ALN

    output:
        set val(pairId), file("${pairId}.chimera.names"), file("${pairId}.chimera.fasta") \
        into NO_CHIMERA_FASTA
    
    script:
    """
    ${params.scripts}/mothur.sh --step=chimera --pairId=${pairId} --inputFasta=${fasta} --inputNames=${names} 
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
    tag { "taxaFilter.${pairId}" }
    publishDir "${params.outdir}/taxaFiltering", mode: "copy", overwrite: false
    errorStrategy "${params.errorsHandling}"
    
    input:
	set val(pairId), file(names), file(fasta) \
    from NO_CHIMERA_FASTA

    output:
	set val(pairId), file("${pairId}.taxaFilter.fasta"), file("${pairId}.taxaFilter.names") \
    into TAXA_FILTERED_CONTIGS

    when:
	params.taxaFilter == 1
    
    script:
    """
    ${params.scripts}/mothur.sh \
        --step=taxaFilter \
        --pairId=${pairId} \
        --inputFasta=${fasta} \
        --inputNames=${names} \
        --refAln=${params.referenceAln} \
        --refTax=${params.referenceTax}
    
    """
}

// /*
//  *
//  * Step 13: Merge individual Fasta files for OTU clustering
//  *
//  */

// process MergeFastas {
//     tag { "mergeFastas" }
//     label "medium_computation"
//     publishDir "${params.outdir}/MergedFastas", mode: "copy", overwrite: false
//     errorStrategy "${params.errorsHandling}"
    
//     input:
// 	file fastas from DADA_FASTA.collect()
//     output:
// 	set val(100), file("all_samples_merged.fasta") into ESV_ALL_SAMPLES_TO_CLUSTER
//     script:
//     """
//     #!/usr/bin/env python3

//     import sys
//     sys.path.append("${workflow.projectDir}/scripts")
//     from util import mergeSamplesFa

//     mergeSamplesFa()
    
//     ## previous behaviour :  cat *dada.fasta > all_sample_merged.fasta
//     """    
// }

// /*
//  *
//  * Step 14: Clustering
//  *
//  */

// process Clustering {
//     tag { "clustering" }
//     label "medium_computation"    
//     publishDir "${params.outdir}/MergedFastas", mode: "copy", overwrite: false
//     errorStrategy "${params.errorsHandling}"
    
//     input:
// 	each idThreshold from (95,97,99)
//         set val(esvId),file(esvFasta) from ESV_ALL_SAMPLES_TO_CLUSTER
//     output:
// 	set val(idThreshold),file("otus${idThreshold}_seq.fasta") into OTU_ALL_SAMPLES,OTU_ALL_SAMPLES_LULU
//         set val(idThreshold),file("otus${idThreshold}_table.csv") into ABUNDANCE_TABLES_OTU
    
//     script:
//     """
//     vsearch --threads ${task.cpus} \
//             --cluster_size ${esvFasta} \
//             --id 0.${idThreshold} \
//             --strand plus \
//             --sizein \
//             --sizeout \
//             --fasta_width 0 \
//             --uc clusters${idThreshold}.uc \
//             --relabel OTU${idThreshold}_ \
//             --centroids otus${idThreshold}_seq.fasta \
//             --otutabout otus${idThreshold}_table.tsv

//     cat otus${idThreshold}_table.tsv | tr "\\t" "," > otus${idThreshold}_table.csv
//     """
// }

// // Mix ESV and OTU channels

// ABUNDANCE_TABLES = ABUNDANCE_TABLES_OTU.mix(ABUNDANCE_TABLES_ESV)
// LULU_ALL_SAMPLES = OTU_ALL_SAMPLES_LULU.mix(ESV_ALL_SAMPLES_LULU)
// ALL_SAMPLES = OTU_ALL_SAMPLES.mix(ESV_ALL_SAMPLES)

// /*
//  *
//  * Step 15: Pre-lulu step
//  *    - Blast each contig against each other
//  *
//  */

// process PreLulu {
//     tag { "preLulus" }
//     label "medium_computation"    
//     publishDir "${params.outdir}/lulu", mode: "copy", overwrite: false
//     errorStrategy "${params.errorsHandling}"

//     input:
// 	set val(idThreshold),file(fasta) from LULU_ALL_SAMPLES
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
//  * Step 18: Annotation
//  *
//  */

// process ClassificationVsearch  {
//     tag { "classification" }
//     label "medium_computation"
//     publishDir "${params.outdir}/taxonomy", mode: "copy", overwrite: false
//     errorStrategy "${params.errorsHandling}"

//     input:
// 	set idThreshold,file(fasta) from FASTAS_LULU
//     output:
//         set idThreshold,file("taxonomy_${idThreshold}.tsv") into TAXONOMY
//     script:
	
//     """
//     vsearch --threads ${task.cpus} \
// 	    --db ${params.referenceVSEARCH} \
// 	    --userfields query+target+id+alnlen+qcov+qstrand \
// 	    --userout taxonomy_${idThreshold}.tsv \
// 	    --alnout aln${idThreshold}.tsv \
// 	    --usearch_global ${fasta} \
// 	    --id ${params.taxaMinId}
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

