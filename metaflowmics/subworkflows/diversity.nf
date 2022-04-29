#!/usr/bin/env nextflow

nextflow.enable.dsl=2
params.options = [:]

module_dir = "../modules"

include{ FASTTREE } from "$module_dir/fasttree/main.nf"
include{ PHYLOSEQ_UNIFRAC as PHYLOSEQ_UNIFRAC_WEIGHTED } from "$module_dir/R/phyloseq/unifrac/main.nf" \
    addParams( unifrac: "weighted" )
include{ PHYLOSEQ_UNIFRAC as PHYLOSEQ_UNIFRAC_UNWEIGHTED } from "$module_dir/R/phyloseq/unifrac/main.nf" \
    addParams( unifrac: "unweighted" )
include{ MOTHUR_SUMMARY_SINGLE } from "$module_dir/mothur/summarySingle/main.nf" \
    addParams( calc: params.alpha_diversity )
include{ MOTHUR_SUMMARY_SHARED } from "$module_dir/mothur/summaryShared/main.nf" \
    addParams( calc: params.beta_diversity )


workflow DIVERSITY {
    take:
    repfasta
    shared

    main:
    if (params.alpha_diversity != "") {
        MOTHUR_SUMMARY_SINGLE(
            shared.map{ [it[0], "", it[1]] }
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

        PHYLOSEQ_UNIFRAC_WEIGHTED(shared.join(tree))
        PHYLOSEQ_UNIFRAC_UNWEIGHTED(shared.join(tree))
    }
}
