// Help function for ITS pipeline

def helpMessage() {
    log.info"""
    ===================================
     Metaflow|mics ITS pipeline
    ===================================
    Usage:
    nextflow run ITS-pipeline --reads '*_R1.fastq.gz' --locus ITS1 -profile local

    ---------------------------------- Mandatory arguments ----------------------------------------
      --reads         Path to input data

    ---------------------------------- Optional arguments ----------------------------------------
      --paired_end     If your data is paired-end
      --locus         Sequenced ITS region (ITS1,ITS2 or ALL). Default: ITS1
      --outdir        The output directory where the results will be saved. 
                      Default: ./ITS-pipeline_outputs

    [Quality filtering]
      --max_expected_errors  Discard reads with more than <max_expected_errors> errors. Default: 3
      --min_read_count       Discard samples with less than <minReads> reads. 
                               Default: 50
    [Taxonomy assignment]
      --tax_confidence     The minimum consensus for calling a taxonomy. 
                             Default: 50%
    [Co-occurence pattern correction]
    --skip_lulu                 Skip the Lulu step (can save a lot of computation time)
    --lulu_min_ratio_type       Function to compare the abundance of a parent OTU and its daughter.
                                  Default: min
    --lulu_min_ratio            Minimum abundance ratio between parent and daughter
                                  (across all samples). Default: 1
    --lulu_min_rel_cooccurence  Proportion of the parent samples where daughter occurs. Default: 1
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
    summary['paired end'] = params.paired_end
    summary['locus'] = params.locus
    summary['max expected errors'] = params.max_expected_error
    summary['min reads per sample'] = params.min_read_count
    summary['clustering similarity thresholds'] = params.clustering_thresholds
    summary['skip lulu'] = params.skip_lulu

    if (!params.skip_lulu) {
        summary['lulu ratio type'] = params.lulu_min_ratio_type
        summary['lulu parent/daughter min similarity'] = params.lulu_min_match
        summary['lulu parent/daughter min abundance ratio'] = params.lulu_min_ratio
        summary['lulu parent/daughter min co-occurrence'] = params.lulu_min_rel_cooccurence
    }
    summary['consensus taxonomy confidence threshold'] = params.tax_confidence

    file(params.outdir).mkdir()
    File f = new File("${params.outdir}/parameters_summary.log")
    f.write("====== Parameter summary =====\n")
    summary.each { key, val -> f.append("\n${key.padRight(50)}: ${val}") }
}
