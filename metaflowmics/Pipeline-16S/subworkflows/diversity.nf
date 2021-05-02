#!/usr/bin/env nextflow

nextflow.enable.dsl=2
params.options = [:]

module_dir = "../../modules"

include{ FASTTREE } from "$module_dir/fasttree/main.nf" \
    addParams( options: [publish_dir: "postprocessing/beta-diversity"] )
include{ PHYLOSEQ_UNIFRAC } from "$module_dir/phyloseq/unifrac/main.nf" \
    addParams( options: [publish_dir: "postprocessing/beta-diversity"], method: params.unifrac )
include{ MOTHUR_SUMMARY_SINGLE } from "$module_dir/mothur/summarySingle/main.nf" \
    addParams( options: [publish_dir: "postprocessing/beta-diversity"], calc: params.alpha_diversity)
include{ MOTHUR_SUMMARY_SHARED } from "$module_dir/mothur/summaryShared/main.nf" \
    addParams( options: [publish_dir: "postprocessing/beta-diversity"], calc: params.beta_diversity)

// Main workflow
workflow diversity {
    take:
    repfasta
    shared

    main:
    tree = FASTTREE(
        repfasta
    ).nwk

    if (params.alpha_diversity != "") {
        MOTHUR_SUMMARY_SINGLE(
            shared.map{[[otu_id: it[0], step: ''], it[1]]}
        )
    }
    if (!params.skip_unifrac) {
        PHYLOSEQ_UNIFRAC(shared.join(tree))
    }

    if (params.beta_diversity != "") {
        MOTHUR_SUMMARY_SHARED(
            shared
        )
    }
    // alpha diversity(shared)
    // FastTree(repfasta) > UniFrac
    // Braycurtis(shared)
}
