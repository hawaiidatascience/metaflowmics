nextflow.enable.dsl = 2

params.options = [:]
module_dir = "../../modules"
subworkflow_dir = "./subworkflows"

include {compute_representatives as get_representatives_species} from "$subworkflow_dir/compute_repr.nf" \
    addParams( level: 7, outdir: "$params.outdir/up/1-species" )
include {compute_representatives as get_representatives_genus} from "$subworkflow_dir/compute_repr.nf" \
    addParams( level: 6, outdir: "$params.outdir/up/2-genus" )
include {compute_representatives as get_representatives_subfamily} from "$subworkflow_dir/compute_repr.nf" \
    addParams( level: 5, outdir: "$params.outdir/up/3-subfamily" )
include {compute_representatives as get_representatives_family} from "$subworkflow_dir/compute_repr.nf" \
    addParams( level: 4, outdir: "$params.outdir/up/4-family" )
include {compute_representatives as get_representatives_order} from "$subworkflow_dir/compute_repr.nf" \
    addParams( level: 3, outdir: "$params.outdir/up/5-order" )
include {compute_representatives as get_representatives_class} from "$subworkflow_dir/compute_repr.nf" \
    addParams( level: 2, outdir: "$params.outdir/up/6-class" )
include {compute_representatives as get_representatives_phylum} from "$subworkflow_dir/compute_repr.nf" \
    addParams( level: 1, outdir: "$params.outdir/up/7-phylum" )
include {compute_representatives as get_representatives_kingdom} from "$subworkflow_dir/compute_repr.nf" \
    addParams( level: 0, outdir: "$params.outdir/up/8-kingdom" )

include {UPDATE_MSA_WITH_REF as UPDATE_PHYLUM} from "$module_dir/python/biopython/msa/main.nf" \
    addParams( level: 1, outdir: "$params.outdir/down/1-phylum" )
include {UPDATE_MSA_WITH_REF as UPDATE_CLASS} from "$module_dir/python/biopython/msa/main.nf" \
    addParams( level: 2, outdir: "$params.outdir/down/2-class" )
include {UPDATE_MSA_WITH_REF as UPDATE_ORDER} from "$module_dir/python/biopython/msa/main.nf" \
    addParams( level: 3, outdir: "$params.outdir/down/3-order" )
include {UPDATE_MSA_WITH_REF as UPDATE_FAMILY} from "$module_dir/python/biopython/msa/main.nf" \
    addParams( level: 4, outdir: "$params.outdir/down/4-family" )
include {UPDATE_MSA_WITH_REF as UPDATE_SUBFAMILY} from "$module_dir/python/biopython/msa/main.nf" \
    addParams( level: 5, outdir: "$params.outdir/down/5-subfamily" )
include {UPDATE_MSA_WITH_REF as UPDATE_GENUS} from "$module_dir/python/biopython/msa/main.nf" \
    addParams( level: 6, outdir: "$params.outdir/down/6-genus" )
include {UPDATE_MSA_WITH_REF as UPDATE_SPECIES} from "$module_dir/python/biopython/msa/main.nf" \
    addParams( level: 7, outdir: "$params.outdir/down/7-species" ) 
include {UPDATE_MSA_WITH_REF as UPDATE_STRAINS} from "$module_dir/python/biopython/msa/main.nf" \
    addParams( level: 8, outdir: "$params.outdir/down/8-strains" ) 


workflow {

    fasta = Channel.fromPath(params.fasta)

    species_aln = get_representatives_species( fasta )
    genus_aln = get_representatives_genus( species_aln.repr )
    subfamily_aln = get_representatives_subfamily( genus_aln.repr )
    family_aln = get_representatives_family( subfamily_aln.repr )
    order_aln = get_representatives_order( family_aln.repr )
    class_aln = get_representatives_class( order_aln.repr )
    phylum_aln = get_representatives_phylum( class_aln.repr )
    kingdom_aln = get_representatives_kingdom( phylum_aln.repr )    

    phylum_update = UPDATE_PHYLUM(
        phylum_aln.afa.collectFile(name: "phyla.afa"),
        kingdom_aln.afa.collectFile(name: "kingdom.updated.afa")
    )    
    class_update = UPDATE_CLASS(
        class_aln.afa.collectFile(name: "classes.afa"),
        phylum_update.afa.flatten().collectFile(name: "phyla.updated.afa")
    )
    order_update = UPDATE_ORDER(
        order_aln.afa.collectFile(name: "orders.afa"),
        class_update.afa.flatten().collectFile(name: "orders.updated.afa")
    )
    family_update = UPDATE_FAMILY(
        family_aln.afa.collectFile(name: "families.afa"),
        order_update.afa.flatten().collectFile(name: "orders.updated.afa")
    )
    subfamily_update = UPDATE_SUBFAMILY(
        subfamily_aln.afa.collectFile(name: "subfamilies.afa"),
        family_update.afa.flatten().collectFile(name: "families.updated.afa")
    )
    genus_update = UPDATE_GENUS(
        genus_aln.afa.collectFile(name: "genera.afa"),
        subfamily_update.afa.flatten().collectFile(name: "subfamilies.updated.afa")
    )
    species_update = UPDATE_SPECIES(
        species_aln.afa.collectFile(name: "species.afa"),
        genus_update.afa.flatten().collectFile(name: "genera.updated.afa")
    )
    strain_updates = UPDATE_STRAINS(
        fasta,
        species_update.afa.flatten().collectFile(name: "species.updated.afa")
    )
}
