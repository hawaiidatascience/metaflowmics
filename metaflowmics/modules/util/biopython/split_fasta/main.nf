// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

options = initOptions(params.options)


// Split fasta sequences using the header's taxonomy (split in files corresponding to the same taxonomic rank)

process SPLIT_FASTA {
    label "process_low"
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options,
                                        publish_dir:getSoftwareName(task.process)) }

    container "nakor/metaflowmics-python:0.0.1"
    conda (params.enable_conda ? "conda-forge::biopython" : null)

    input:
    path fasta
    val info

    output:
    path "*.faa", emit: faa

    script:
    """
    #!/usr/bin/env python

    import resource
    from Bio.SeqIO.FastaIO import SimpleFastaParser

    resource.setrlimit(
        resource.RLIMIT_NOFILE,
        (int(1e6), int(1e6))
    )

    writers = dict()

    with open("$fasta", "r") as reader:
        for i, (title, seq) in enumerate(SimpleFastaParser(reader)):
            acc = ';'.join(title.split(';')[:-1])
            lineage = title.split(';')[-1].split('$params.sep')

            fname = f'{lineage[$info.field]}.faa'

            if fname not in writers:
                writers[fname] = open(fname, 'w')

            writers[fname].write(f'>{acc};{",".join(lineage)}\\n{seq}\\n')

    for w in writers.values():
        w.close()
    """
}
