#!/usr/bin/env nextflow

nextflow.enable.dsl=2
params.options = [:]

module_dir = "../modules"

include{ FASTTREE } from "$module_dir/fasttree/main.nf"
include{ PHYLOSEQ_UNIFRAC } from "$module_dir/phyloseq/unifrac/main.nf" \
    addParams( method: params.unifrac )
include{ MOTHUR_SUMMARY_SINGLE } from "$module_dir/mothur/summarySingle/main.nf" \
    addParams( calc: params.alpha_diversity)
include{ MOTHUR_SUMMARY_SHARED } from "$module_dir/mothur/summaryShared/main.nf" \
    addParams( calc: params.beta_diversity)

// Main workflow
workflow diversity {
    take:
    repfasta
    shared

    main:
    if (params.alpha_diversity != "") {
        MOTHUR_SUMMARY_SINGLE(
            shared.map{[[otu_id: it[0], step: ''], it[1]]}
        )
    }
    if (params.beta_diversity != "") {
        MOTHUR_SUMMARY_SHARED(
            shared
        )
    }
    if (!params.skip_unifrac) {
        tree = FASTTREE(
            repfasta
        ).nwk

        PHYLOSEQ_UNIFRAC(shared.join(tree))
    }
}
