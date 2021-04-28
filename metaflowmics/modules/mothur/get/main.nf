// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

options = initOptions(params.options)

process MOTHUR_GET {
    tag "$otu_id"
    label "process_high"

    container "quay.io/biocontainers/mothur:1.44.1--hf0cea05_2"
    conda (params.enable_conda ? "bioconda::mothur:1.44.1" : null)

    input:
    tuple val(otu_id), file(ref), file(filt)

    output:
    tuple val(otu_id), path("*.pick.*")

    script:
    def ref_ext = ref.getExtension()
    def filt_ext = filt.getExtension().replaceAll('_table', '')
    def fn = filt_ext ==~ /.*\.(shared|list|cons\.taxonomy|rep\.fasta)/ ?
        "otus" : "seqs"
    """
    mothur "#list.${fn}(${ref_ext}=$ref);get.${fn}(accnos=current,$filt_ext=$filt)"
    """
}

process MOTHUR_GET_SEQS {
    tag "$otu_id"
    label "process_high"

    container "quay.io/biocontainers/mothur:1.44.1--hf0cea05_2"
    conda (params.enable_conda ? "bioconda::mothur:1.44.1" : null)

    input:
    tuple val(otu_id), file(ref), file(filt)

    output:
    tuple val(otu_id), path("*.pick.*")

    script:
    def ref_ext = ref.getExtension()
    def filt_ext = filt.getExtension().replaceAll('_table', '')
    """
    mothur "#list.seqs(${ref_ext}=$ref);get.seqs(accnos=current,$filt_ext=$filt)"
    """
}

process MOTHUR_GET_OTUS {
    tag "$otu_id"
    label "process_high"

    container "quay.io/biocontainers/mothur:1.44.1--hf0cea05_2"
    conda (params.enable_conda ? "bioconda::mothur:1.44.1" : null)

    input:
    tuple val(otu_id), file(ref), file(filt)

    output:
    tuple val(otu_id), path("*.pick.*")

    script:
    def ref_ext = ref.getExtension()
    def filt_ext = filt.getExtension().replaceAll('_table', '')
    """
    mothur "#list.otus(${ref_ext}=$ref);get.otus(accnos=current,$filt_ext=$filt)"
    """
}
