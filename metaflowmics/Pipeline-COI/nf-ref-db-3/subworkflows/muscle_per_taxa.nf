#!/usr/bin/env nextflow

nextflow.enable.dsl = 2


process SPLIT_FASTA {
	label "process_low"
	container "nakor/metaflowmics-python:0.0.1"
	conda (params.enable_conda ? "conda-forge::biopython" : null)

	input:
	path fasta

	output:
	path "*.main.faa", optional: true, emit: main
	path "*.others.faa", optional: true, emit: others

	script:
	"""
	split_fasta.py --fasta $fasta \\
	--field $params.field \\
	--sep "$params.sep" \\
	--min-group-size $params.min_group_size
	"""
}

process MUSCLE {
    label "process_high"
    publishDir "${params.outdir}", mode: params.publish_dir_mode
    container "nakor/muscle:5.0.1278"

    input:
    path fasta

    output:
    path "*.afa", emit: afa

    script:
    def arg = fasta.getExtension() == "afa" ? "-profile" : ""
    def prefix = fasta.getBaseName()
    """
    muscle -in $fasta \\
        | awk '/^>/ {printf("\\n%s\\n",\$0);next; } { printf("%s",\$0);}  END {printf("\\n");}' \\
        | tail -n+2 > ${prefix}.afa

    [ -s ${prefix}.afa ] && echo "muscle was successfull" || (echo "Something went wrong. Did muscle segfault?" && exit 1)
    """
}

process COMPUTE_MSA_REPRESENTATIVE {
    label "process_low"
    publishDir "${params.outdir}", mode: params.publish_dir_mode
    container "nakor/metaflowmics-python:0.0.1"

    input:
    path afa

    output:
    path "*repr.afa", emit: repr

    script:
    """	
	compute_consensus.py --fasta $afa
	"""
}

// Main workflow
workflow align_tax_groups {
    take:
    fasta

    main:
    // Split fasta in taxonomic groups
    grouped_taxa = SPLIT_FASTA(
        fasta.collectFile(name: "repr_lvl-${params.field}.faa")
    )

    // Align each taxonomic group
    afa_per_taxa = MUSCLE(grouped_taxa.main.flatten()).afa

    // Get one representative for each aligned fasta to speed up computation
    afa_representatives = COMPUTE_MSA_REPRESENTATIVE(afa_per_taxa).repr

    emit:
    // Re-introduce the sequences for taxa with one sequence
    afa = afa_per_taxa.mix(grouped_taxa.others.flatten())
    repr = afa_representatives.mix(grouped_taxa.others.flatten())
}
