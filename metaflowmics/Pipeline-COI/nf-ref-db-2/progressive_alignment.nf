nextflow.enable.dsl = 2

params.options = [:]
module_dir = "../../modules"
subworkflow_dir = ".subworkflows"

include {align_tax_groups as align_species} from "$subworkflow_dir/muscle_per_taxa.nf" \
    addParams( field: 6, outdir: "$params.outdir/up/1-species" )
include {align_tax_groups as align_genus} from "$subworkflow_dir/muscle_per_taxa.nf" \
    addParams( field: 5, outdir: "$params.outdir/up/2-genus" )
// include {align_tax_groups as align_subfamily} from "$subworkflow_dir/muscle_per_taxa.nf" \
//     addParams( field: 5, outdir: "$params.outdir/up/3-subfamily" )
include {align_tax_groups as align_family} from "$subworkflow_dir/muscle_per_taxa.nf" \
    addParams( field: 4, outdir: "$params.outdir/up/3-family" )
include {align_tax_groups as align_order} from "$subworkflow_dir/muscle_per_taxa.nf" \
    addParams( field: 3, outdir: "$params.outdir/up/4-order" )
include {align_tax_groups as align_class} from "$subworkflow_dir/muscle_per_taxa.nf" \
    addParams( field: 2, outdir: "$params.outdir/up/5-class" )
include {align_tax_groups as align_phylum} from "$subworkflow_dir/muscle_per_taxa.nf" \
    addParams( field: 1, outdir: "$params.outdir/up/6-phylum" )
include {align_tax_groups as align_kingdom} from "$subworkflow_dir/muscle_per_taxa.nf" \
    addParams( field: 0, outdir: "$params.outdir/up/7-kingdom" )

include {UPDATE_MSA_WITH_REF as UPDATE_PHYLUM} from "$module_dir/python/biopython/msa/main.nf" \
    addParams( field: 1, outdir: "$params.outdir/down/1-phylum" )
include {UPDATE_MSA_WITH_REF as UPDATE_CLASS} from "$module_dir/python/biopython/msa/main.nf" \
    addParams( field: 2, outdir: "$params.outdir/down/2-class" )
include {UPDATE_MSA_WITH_REF as UPDATE_ORDER} from "$module_dir/python/biopython/msa/main.nf" \
    addParams( field: 3, outdir: "$params.outdir/down/3-order" )
include {UPDATE_MSA_WITH_REF as UPDATE_FAMILY} from "$module_dir/python/biopython/msa/main.nf" \
    addParams( field: 4, outdir: "$params.outdir/down/4-family" )
// include {UPDATE_MSA_WITH_REF as UPDATE_SUBFAMILY} from "$module_dir/python/biopython/msa/main.nf" \
//     addParams( field: 5, outdir: "$params.outdir/down/5-subfamily" )
include {UPDATE_MSA_WITH_REF as UPDATE_GENUS} from "$module_dir/python/biopython/msa/main.nf" \
    addParams( field: 5, outdir: "$params.outdir/down/5-genus" )
include {UPDATE_MSA_WITH_REF as UPDATE_SPECIES} from "$module_dir/python/biopython/msa/main.nf" \
    addParams( field: 6, outdir: "$params.outdir/down/6-species" )

workflow {

    fasta = Channel.fromPath(params.fasta)
    
    species_aln = align_species( fasta )
    genus_aln = align_genus( species_aln.repr )
    // subfamily_aln = align_subfamily( genus_aln.repr )
    // family_aln = align_family( subfamily_aln.repr )
    family_aln = align_family( genus_aln.repr )
    order_aln = align_order( family_aln.repr )
    class_aln = align_class( order_aln.repr )
    phylum_aln = align_phylum( class_aln.repr )
    kingdom_aln = align_kingdom( phylum_aln.repr )    

    phylum_update = UPDATE_PHYLUM(
        phylum_aln.afa.collectFile(name: "phyla.afa"),
        kingdom_aln.afa.collectFile(name: "kingdom.updated.afa", storeDir: "$params.outdir/down")
    )    
    class_update = UPDATE_CLASS(
        class_aln.afa.collectFile(name: "classes.afa"),
        phylum_update.afa.flatten().collectFile(name: "phyla.updated.afa", storeDir: "$params.outdir/down")
    )
    order_update = UPDATE_ORDER(
        order_aln.afa.collectFile(name: "orders.afa"),
        class_update.afa.flatten().collectFile(name: "classes.updated.afa", storeDir: "$params.outdir/down")
    )
    family_update = UPDATE_FAMILY(
        family_aln.afa.collectFile(name: "families.afa"),
        order_update.afa.flatten().collectFile(name: "orders.updated.afa", storeDir: "$params.outdir/down")
    )
    // subfamily_update = UPDATE_SUBFAMILY(
    //     subfamily_aln.afa.collectFile(name: "subfamilies.afa"),
    //     family_update.afa.flatten().collectFile(name: "families.updated.afa", storeDir: "$params.outdir/down")
    // )
    // genus_update = UPDATE_GENUS(
    //     genus_aln.afa.collectFile(name: "genera.afa"),
    //     subfamily_update.afa.flatten().collectFile(name: "subfamilies.updated.afa", storeDir: "$params.outdir/down")
    // )
    genus_update = UPDATE_GENUS(
        genus_aln.afa.collectFile(name: "genera.afa"),
        family_update.afa.flatten().collectFile(name: "families.updated.afa", storeDir: "$params.outdir/down")
    )
    species_update = UPDATE_SPECIES(
        species_aln.afa.collectFile(name: "species.afa"),
        genus_update.afa.flatten().collectFile(name: "genera.updated.afa", storeDir: "$params.outdir/down")
    )
    species_update.afa.flatten().collectFile(name: "species.updated.afa", storeDir: "$params.outdir/down")
}
