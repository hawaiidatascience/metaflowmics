// Help function for 16S pipeline

def helpMessage() {
    log.info"""
    ===================================
    Metaflow|mics 16S-pipeline
    ===================================
    Usage:
    nextflow run 16S-pipeline -profile local --reads "*_R{1,2}.fastq.gz"
    For more detailed information, see https://metagenomics-pipelines.readthedocs.io/en/latest/
    
    ---------------------------------- Mandatory arguments ----------------------------------------
    --reads         Path to input data (glob pattern)

    ---------------------------------- Optional arguments ----------------------------------------
    --outdir        Path to output directory. Default: "./16S-pipeline_outputs"
    --single_end    If your data is single end
                      Default: false
    [Quality filtering]
    --min_read_count      Sample with less than <minRead> are discarded. Default: 50
    --trunc_len           Where to truncate the forward and reverse reads. Default: "220,190"
    --min_read_len        Reads short than <minLength> are discarded. Default: 20
    --max_expected_error  Reads with an expected number of errrors above <maxEE> are discarded. 
                            Default: 3
    --trunc_quality       Read truncation at the 1st occurence of a base of quality <= <truncQ>. 
                            Default: 2
    --keep_phix           Keep reads matching phiX genome.
                            Default: false
    [Denoising]
    --pool                Pooling method for denoising with dada2. Choices are TRUE, FALSE or pseudo
                            Default: FALSE
    [Read merging]
    --min_overlap    Minimum overlap between forward and reverse read. Default: 20
    --max_mismatch   Maximum number of mismatches in the overlap between forward and reverse read. 
                       Default: 1
    [Contig filtering]
    --criteria        Optimization criteria when aligning sequences against reference database. 
                        Discard any sequence starting after where <criteria>% of the sequences 
                        start or end before <criteria>% of the sequences end. 
                        Default: 95
    --min_aln_len     Minimum alignment length in MSA. Default: 50
    --chimera_tool    Software to filter chimera (see mothur documentation for all the available 
                        softwares). Default: vsearch
    --taxa_to_filter  Set of taxa to exclude from the analysis. 
                        Default: "Chloroplast-Mitochondria-unknown"
    
    [Subsampling]
    --custom_subsampling_level  User defined subsampling level. Ignored if <= 0
    --subsampling_quantile      Automatic subsampling level is at quantile <subsamplingQuantile>
                                  of the sample sizes. Ignored if customSubsamplingLevel or 
                                  skipSubsampling are set.
                                  Default: 0.1
    --min_subsampling           Minimum subsampling level used if the automatic level falls below 
                                  this value. Default: 5000
    --skip_subsampling          Skip the subsampling step.
                                  Default: false
    [OTU clustering]
    --clustering_thresholds     Percent similarity threshold for OTU clustering. Use 100 for ASVs.
                                  Default: 100,97
    [Co-occurence pattern correction]
    --skip_lulu                  Skip the Lulu step
    --lulu_min_ratio_type        Function to compare the abundance of a parent OTU and its daughter.
                                   Default: min
    --lulu_min_ratio             Minimum abundance ratio between parent and daughter
                                   (across all samples). Default: 1
    --lulu_min_rel_cooccurence   Proportion of the parent samples where daughter occurs. 
                                   Default: 1
    [Other]
    --silva_release      Name of SILVA reference database (seed or nr)
    --silva_version      Version of silva (e.g. 138_1)
    --min_abundance      Remove OTUs with a total abundance equal or below <minAbundance>. Default: 2
    --compute_mothur_db  Compute mothur database summary file. Can be memory intensive. Default: false
    --skip_unifrac       Skip unifrac calculation. Recommended if your memory is limited and
                           your have a significant amount of OTUs. Default: false
    """.stripIndent()
}

// Show help message
params.help = false
if (params.help){
    helpMessage()
    exit 0
}

def saveParams() {
    def summary = [:]
    summary['silva database'] = "${params.silva_release}_v${params.silva_version}"
    summary['paired end'] = params.paired_end
    summary['min reads per sample'] = params.min_read_count
    summary['min read length'] = params.min_read_len
    summary['read truncation'] = params.trunc_len
    summary['max number of expected errors'] = params.max_expected_error
    summary['quality truncation'] = params.trunc_quality
    summary['keep phix genome'] = params.keep_phix
    summary['Pooling method for denoising'] = params.pool
    summary['min overlap for merging'] = params.min_overlap
    summary['max mismatches for merging'] = params.max_mismatch
    summary['Percentile for start/end contig filtering'] = params.criteria
    summary['Minimum contig alignment length against db'] = params.min_aln_len
    summary['Filtered taxa'] = params.taxa_to_filter
    summary['Skip subsampling'] = params.skip_subsampling

    if (!params.skip_subsampling) {
        if (params.custom_subsampling_level > 0){
	        summary['Subsampling'] = params.custom_subsampling_level
        } else {
	        summary['Percentile for automatic subsampling'] = params.subsampling_quantile * 100
	        summary['Hard threshold for automatic subsampling'] = params.min_subsampling
        }
    }
    summary['clustering similarity thresholds'] = params.clustering_thresholds
    summary['Min OTU abundace filter'] = params.min_subsampling

    summary['Skip LULU'] = params.skip_lulu
    if (!params.skip_lulu){
	    summary['Lulu ratio type'] = params.lulu_min_ratio_type
	    summary['Lulu parent/daughter min similarity'] = params.lulu_min_match
	    summary['Lulu parent/daughter min abundance ratio'] = params.lulu_min_ratio
	    summary['Lulu parent/daughter min co-occurrence'] = params.lulu_min_rel_cooccurence
    }
    summary['Skip mothur database computation'] = params.skip_unifrac
    summary['Skip unifrac computation'] = params.skip_unifrac

    summary['Alpha diversity metrics'] = params.alpha_diversity
    summary['Beta diversity metrics'] = params.beta_diversity

    file(params.outdir).mkdir()
    File f = new File("${params.outdir}/parameters_summary.log")
    f.write("====== Parameter summary =====\n")
    summary.each { key, val -> f.append("\n${key.padRight(50)}: ${val}") }
}
