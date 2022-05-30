// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from "./functions"

options = initOptions(params.options)

process MOTHUR_GET_OTUS {
    tag "$meta.id"
    label "process_low"

    container "quay.io/biocontainers/mothur:1.47.0--hb64bf22_2"
    conda (params.enable_conda ? "bioconda::mothur=1.47.0" : null)

    input:
    tuple val(meta), file(ref), file(filt)

    output:
    tuple val(meta), path("*.pick.${filt.getExtension()}")

    script:
    def ref_ext = ref.getExtension()
    def arg = filt.getExtension()
    """
    mothur "#list.otulabels(${ref_ext}=$ref);get.otus(accnos=current,$arg=$filt)"
    """
}
