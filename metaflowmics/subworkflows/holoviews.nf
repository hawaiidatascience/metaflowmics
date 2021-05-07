#!/usr/bin/env nextflow

params.options = [:]

moduledir = "../modules"

include{ HOLOVIEWS_PREPARE } from "$moduledir/holoviews/util/main.nf"
include{ HOLOVIEWS_BARS } from "$moduledir/holoviews/bars/main.nf" \
    addParams( min_abund: 0.01 )
include{ HOLOVIEWS_SCATTER } from "$moduledir/holoviews/scatter/main.nf"
include{ HOLOVIEWS_CLUSTERMAP } from "$moduledir/holoviews/heatmap/main.nf" \
    addParams( min_abund: 0.05 )


workflow holoviews {
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
