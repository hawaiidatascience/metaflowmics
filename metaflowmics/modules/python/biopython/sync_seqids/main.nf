// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

options = initOptions(params.options)


process SYNC_SEQIDS {
	tag "$meta.id"
    label "process_low"
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options,
                                        publish_dir:getSoftwareName(task.process)) }

    container "nakor/metaflowmics-python:0.0.2"
    conda (params.enable_conda ? "conda-forge::biopython conda-forge::pandas" : null)

    input:
    tuple val(meta), path(ref_fa)
    path files

    output:
    tuple val(meta), path("*.sync.{fna,fasta}"), optional: true, emit: fna
    tuple val(meta), path("*.sync.faa"), optional: true, emit: faa
    tuple val(meta), path("*.sync.count_table"), optional: true, emit: count_table

    script:
    """
    #!/usr/bin/env python

    from pathlib import Path
    from Bio import SeqIO
    import pandas as pd

    seq_info = {
      seq.id: seq.description.split("frame=")[-1]
      for seq in SeqIO.parse("$ref_fa", "fasta")
    }

    for f in ["${files.join('","')}"]:
        fpath = Path(f)
        if fpath.suffix == ".count_table":
            table = pd.read_csv(fpath, index_col=0, sep='\\t').loc[seq_info.keys()]
            table.to_csv(f"{fpath.stem}.sync.count_table", sep="\\t")
        else: # we assume it is a fasta formatted file
            seq_data = {seq.id: seq for seq in SeqIO.parse(fpath, "fasta")}

            with open(f"{fpath.stem}.sync{fpath.suffix}", "w") as w:
                for (seq_id, frame) in seq_info.items():
                    seq = seq_data[seq_id]
                    if frame.isdigit() and int(frame) > 3:
                        seq = seq.reverse_complement(id=True, description=True)
                    w.write(seq.format("fasta-2line"))
    """
}
