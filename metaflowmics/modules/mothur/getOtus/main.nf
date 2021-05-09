// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from "./functions"

options = initOptions(params.options)

process MOTHUR_GET_OTUS {
    tag "$otu_id"
    label "process_low"

    container "quay.io/biocontainers/mothur:1.44.1--hf0cea05_2"
    conda (params.enable_conda ? "bioconda::mothur:1.44.1" : null)

    input:
    tuple val(otu_id), file(ref), file(filt)

    output:
    tuple val(otu_id), path("*.pick.${filt.getExtension()}")

    script:
    def ref_ext = ref.getExtension()
    def arg = filt.getExtension()
    """
    mothur "#list.otus(${ref_ext}=$ref);get.otus(accnos=current,$arg=$filt)"
    """
}
