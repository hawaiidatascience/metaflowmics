#!/usr/bin/env nextflow

params.options = [:]

moduledir = "../modules/python/holoviews"

include{ HOLOVIEWS_PREPARE } from "$moduledir/util/main.nf"
include{ HOLOVIEWS_BARS } from "$moduledir/bars/main.nf" \
    addParams( min_abund: 0.01 )
include{ HOLOVIEWS_SCATTER } from "$moduledir/scatter/main.nf"
include{ HOLOVIEWS_CLUSTERMAP } from "$moduledir/heatmap/main.nf" \
    addParams( min_abund: 0.05 )


workflow HOLOVIEWS {
    take:
    shared
    constaxonomy

    main:
    data = HOLOVIEWS_PREPARE(
        shared.join(constaxonomy)
    )
    
    HOLOVIEWS_BARS( data )
    HOLOVIEWS_SCATTER( data )
    HOLOVIEWS_CLUSTERMAP( data )
}
