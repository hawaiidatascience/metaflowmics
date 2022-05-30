// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

options = initOptions(params.options)


process FASTQC {
    publishDir "${params.outdir}/figures", mode: "copy"
    label "high_computation"

    container "quay.io/biocontainers/fastqc:0.11.9--0"

    input:
    path fastqs

    output:
    path "*.html"

    script:
    """
    fastqc -o . --threads $task.cpus $fastqs
    """
}

process GUESS_MATCH_ORDER {
    publishDir "${params.outdir}/interm/guess_matching_order", mode: "copy"
    label "process_low"

    conda (params.enable_conda ? "conda-forge::bash='5.1'" : null)
    container "nakor/bash:5.1.4"
    
    input:
    path index
    path "meta.csv" // rename it to make sure we don't have already the same name

    output:
    path "barcodes_meta.csv"

    script:
    z = index[0].getExtension() == 'gz' ? 'z' : ''
    fwd = index[0]
    rev = index.size() == 2 ? index[1] : ''
    """
    #!/usr/bin/env bash

    # Remove carriage returns at the end of file
    sed -e 's/\\r//g' meta.csv > metadata.csv

    # Reverse complement if set
    if [ "$params.rc_rev_index" == true ]; then
        cat metadata.csv | rev | cut -d, -f1 | tr "ATGC" "TACG" > rev_idx_rc.csv
        cat metadata.csv | rev | cut -d, -f2- | rev > samp_idx1.csv
        paste -d, samp_idx1.csv rev_idx_rc.csv > metadata.csv
    fi

    if [ ${index.size()} -eq 2 ]; then
        reversed=\$([ $params.matching == reversed ] && echo true || echo false)

        # if "auto", we check which matching order looks best
        if [ "$params.matching" == "auto" ]; then
            symbol=\$(${z}cat $fwd | head -c 2)
            
            # extract the most frequent fwd/rev index reads pairs and compare it to the metadata 
            paste -d' ' \\
                <( ${z}cat $fwd | grep "^\$symbol" -A1 | grep -v \$symbol | grep -v "-" ) \\
                <( ${z}cat $rev | grep "^\$symbol" -A1 | grep -v \$symbol | grep -v "-" ) \\
                | grep -v N \\
                | sort | uniq -c | sort -rnk1 | head -n 20 \\
                | sed 's/^[[:space:]]*//g' | sed 's/ /,/g' \\
                | cut -d, -f2,3 > freqs.txt

            awk -F, '{OFS=","} {print \$2,\$1}' freqs.txt > freqs_rc.txt

            # take the order that matches best
            n1=\$(comm -12 <(sort <(cut -d, -f2,3 metadata.csv)) <(sort freqs.txt) | wc -l)
            n2=\$(comm -12 <(sort <(cut -d, -f2,3 metadata.csv)) <(sort freqs_rc.txt) | wc -l)

            reversed=\$([ \$n2 -ge \$n1 ] && echo true || echo false)
        fi

        [ \$reversed==true ] \\
            && awk -F"," '{OFS=","}{print \$1,\$3,\$2}' metadata.csv > barcodes_meta.csv \\
            || mv metadata.csv barcodes_meta.csv
    else
        awk -F"," '{OFS=","}{print \$1,\$2,NaN}' metadata.csv > barcodes_meta.csv
    fi
    """
}

process TO_H5 {
    tag "$split"
    label "process_medium"

    container "nakor/metaflowmics-python:0.0.2"
    conda (params.enable_conda ? "conda-forge::biopython conda-forge::h5py conda-forge::numpy" : null)

    input:
    tuple val(split), file(index)

    output:
    tuple val(split), path("*.h5")

    script:
    """
    #!/usr/bin/env python

    import argparse

    from Bio.SeqIO.QualityIO import FastqGeneralIterator
    import numpy as np
    import h5py

    NUCLS = list('ACGTN')
    MAPPING_BASE5 = str.maketrans(''.join(NUCLS), '01234')
    RC_TRANS = str.maketrans('ACGT', 'TGCA')

    def seq2str(seq):
        if hasattr(seq, 'seq'):
            return str(seq.seq)
        return seq

    def seq2intlist(bc):
        bc_vec = np.fromiter(seq2str(bc).translate(MAPPING_BASE5), dtype='uint8')
        return bc_vec

    def seq2idx(bc):
        bc_idx = int(seq2str(bc).translate(MAPPING_BASE5), 5)
        return bc_idx

    def get_lengths(filename):
        handle = open(filename)
        rid_len = len(next(handle).strip().split()[0])
        bc_len = len(next(handle).strip())

        return (rid_len, bc_len)

    def fastqs_to_h5(fastqs, n_reads, bc_len=None, rid_len=None, split=0, rc=False):
        parsers = [FastqGeneralIterator(open(fastq)) for fastq in fastqs]

        data = [
            {'rid': np.empty(n_reads, dtype='S{}'.format(rid_len+10)),
             'index': np.zeros(n_reads, dtype='uint32'),
             'seq': np.zeros((n_reads, bc_len), dtype='uint8'),
             'qual': np.zeros((n_reads, bc_len), dtype='uint8')}
            for _ in fastqs]

        for k, parser in enumerate(parsers):
            for i, (rid, seq, qual) in enumerate(parser):
                if rc:
                    seq = seq.translate(RC_TRANS)
                data[k]['rid'][i] = rid.split()[0]
                data[k]['index'][i] = seq2idx(seq)
                data[k]['seq'][i] = seq2intlist(seq)
                data[k]['qual'][i] = [ord(score)-33 for score in qual]

                if k == 1:
                    (rid_fwd, rid_rev) = data[0]['rid'][i], data[1]['rid'][i]
                    if rid_fwd != rid_rev:
                        print('Error: {} != {}: forward and reverse read IDs do not match.'
                              .format(rid_fwd, rid_rev))
                        exit(42)

        h5_handle = h5py.File('data_{}.h5'.format(split), 'w')

        options = ['fwd', 'rev'][:len(fastqs)]
        for k, orient in enumerate(options):
            for key, val in data[k].items():
                h5_handle.create_dataset('{}/{}'.format(key, orient), data=val)
        h5_handle.close()
    
    n_reads = sum(1 for _ in open("${index[0]}")) // 4
    (rid_len, bc_len) = get_lengths("${index[0]}")
    
    fastqs_to_h5(
        ["${index.join('", "')}"], 
        n_reads, 
        bc_len=bc_len, rid_len=rid_len, 
        split=$split, rc=${params.rc_rev_index.toString().capitalize()})
    """
}
        
process ERROR_MODEL {
    publishDir "${params.outdir}/interm/error_model", mode: "copy"
    label "process_high"

    container "nakor/metaflowmics-python:0.0.2"
    conda (params.enable_conda ? "conda-forge::pandas conda-forge::h5py conda-forge::scikit-learn conda-forge::bokeh" : null)

    input:
    path h5
    path meta_file

    output:
    path "*.h5", emit: h5
    path "*.html", emit: html

    script:
    """
    #!/usr/bin/env python

    import argparse
    from glob import glob

    import pandas as pd
    import numpy as np
    import h5py
    from sklearn.isotonic import IsotonicRegression

    from bokeh.io import output_file, save
    from bokeh.plotting import figure
    from bokeh.layouts import gridplot

    NUCLS = list('ACGTN')
    MAPPING_BASE5 = str.maketrans(''.join(NUCLS), '01234')


    def seq2str(seq):
        if hasattr(seq, 'seq'):
            return str(seq.seq)
        return seq

    def seq2intlist(bc):
        bc_vec = np.fromiter(seq2str(bc).translate(MAPPING_BASE5), dtype='uint8')
        return bc_vec

    def get_h5_keys():
        h5_file = glob("*.h5")[0]
        handle = h5py.File(h5_file, 'r')

        return list(handle.get('seq').keys())

    def load_h5():
        orientations = get_h5_keys()

        indexes = {orient: [] for orient in orientations}
        sequences = {orient: [] for orient in orientations}
        qualities = {orient: [] for orient in orientations}

        for h5 in sorted(glob("*.h5")):
            handle = h5py.File(h5, 'r')
            for orient in orientations:
                indexes[orient].append(handle.get('index/{}'.format(orient))[:])
                sequences[orient].append(handle.get('seq/{}'.format(orient))[:])
                qualities[orient].append(handle.get('qual/{}'.format(orient))[:])
            handle.close()

        fastq_data = {
            orient: {'index': np.concatenate(indexes[orient]),
                     'seq': np.vstack(sequences[orient]),
                     'qual': np.vstack(qualities[orient])}
            for orient in orientations
        }
        return fastq_data

    def idx2seq(nb):
        repr_b5 = "{:08}".format(int(np.base_repr(nb, 5)))
        nucl = repr_b5.translate(str.maketrans("01234", ''.join(NUCLS)))
        return nucl


    def compute_transition_matrix(fastq_data, barcodes, max_dist=$params.max_mismatches, 
                                  n_bases_max=int($params.n_bases), eps=1e-6):
        '''
        - Barcodes in a numpy array of nucleotide index (base 5)
        -
        '''

        # barcodes_index = np.array([int(''.join(map(str, bc)), 5) for bc in barcodes])
        bc_len = len(barcodes[0])
        transition_matrix = np.zeros((4, 5, 40)) # Transitions (A,C,G,T) to (A,C,G,T,N) and quality_scores
        generator = zip(fastq_data['index'], fastq_data['seq'], fastq_data['qual'])
        n_bases = 0

        mem = {}
        while n_bases < n_bases_max:

            try:
                (code, intlist, qual) = next(generator)
            except StopIteration:
                print('Warning: not enough bases for the error model ({:,}/{:,})'
                      .format(n_bases, n_bases_max))

            if n_bases % 10*bc_len == 0:
                print("{:,}/{:,}".format(n_bases, n_bases_max), end='\\r')

            if code in mem:
                matching_bc = mem[code]
            else:
                matching_bc = barcodes[np.sum(barcodes != intlist[None, :], axis=1) <= max_dist]
                mem[code] = matching_bc

            if len(matching_bc) == 0:
                continue

            for bc in matching_bc:
                tup, counts = np.unique(np.vstack((bc, intlist, qual-1)),
                                        return_counts=True, axis=1)
                transition_matrix[tup[0, :], tup[1, :], tup[2, :]] += counts
                n_bases += bc_len

        sums = transition_matrix.sum(axis=0)
        sums[sums==0] = eps

        transition_matrix /= sums

        return transition_matrix
    
    def regression_model(freqs, deg=2, same=False, method='isotonic'):
        '''
        - qual: all measured quality scores when a transition was observed
        - proportion of transitions for a given quality
        '''

        observed_transitions = (~np.isnan(freqs)) & (freqs>0)

        x = np.arange(1, 41)[observed_transitions]
        y = -np.log10(freqs[observed_transitions])

        if len(x) == 0:
            return np.zeros(40)

        if method == 'polynomial':
            z = np.polyfit(x, y, 3)
            polynom = np.poly1d(z)
            y_interp = 10**-polynom(np.arange(1, 41))
        elif method == 'isotonic':
            ir = IsotonicRegression(y_min=0, out_of_bounds='clip', increasing=not same)
            ir.fit(x, y)
            y_interp = 10**-ir.predict(np.arange(1, 41))
        else:
            print('Unknown method: {}. Aborting.'.format(method))
            exit(1)
        return y_interp

    def compute_error_models(idx, barcodes, orient='fwd'):

        barcodes_int = np.array(
            [list(seq.translate(MAPPING_BASE5)) for seq in barcodes]
        ).astype('uint8')

        transitions = compute_transition_matrix(idx, barcodes_int)

        transitions_interp = np.array(
            [[regression_model(err_prob, same=(i==j))
              for i, err_prob in enumerate(transition_from_orig)]
             for j, transition_from_orig in enumerate(transitions)])

        display_models_bokeh(transitions, transitions_interp, orient=orient)
        print('Error model computed')

        return transitions_interp

    def display_models_bokeh(trans_matrix, trans_matrix_interp, orient='fwd'):
        data = []

        for matrix in [trans_matrix, trans_matrix_interp]:
            # For each matrix, stores non-zero probabilities with transitions and quality
            nz_info = np.nonzero(matrix)
            data.append(np.vstack((matrix[nz_info], *nz_info)).T)

        labels = np.array(['y1']*len(data[0]) + ['y_hat']*len(data[1]))

        data = np.vstack(data)

        transitions = ["{}->{}".format(NUCLS[i], NUCLS[j]) for i, j in data[:, [1, 2]].astype(int)]

        data = pd.DataFrame({'transition': transitions,
                             'quality': 1 + data[:, 3].astype(int),
                             'P_err': data[:, 0],
                             'type': labels})

        data.set_index(['transition', 'type'], inplace=True)

        plots = []
        tools = ['hover', 'box_zoom', 'reset']

        observed_transitions = set(pd.MultiIndex.get_level_values(data.index,'transition'))
        for i, nucl_i in enumerate(NUCLS[:-1]):
            for j, nucl_j in enumerate(NUCLS):
                transition_label = "{}->{}".format(nucl_i, nucl_j)
                if transition_label not in observed_transitions:
                    plots.append(None)
                    continue

                data_s = data.loc[transition_label].sort_values(by='quality')

                p = figure(title="{}->{}".format(nucl_i, nucl_j), tooltips=[], tools=tools,
                           x_range=(0, 40), y_range=(-0.1, 1.1))
                p.circle(x='quality', y='P_err', source=data_s.loc[['y1']], alpha=0.7, size=8)
                p.line(x='quality', y='P_err', source=data_s.loc[['y_hat']], color='red',
                       line_dash="dashed", line_width=2)
                plots.append(p)

        grid = gridplot(plots, ncols=5, plot_width=300, plot_height=300)
        output_file("error_model_{}.html".format(orient))
        save(grid)

    # main script
    barcodes = (pd.read_csv("$meta_file", dtype=str, names=["sample_name", "fwd", "rev"])
                .dropna(axis=1, how='all')
                .drop('sample_name', axis=1))

    fastq_data = load_h5()

    transition_h5 = h5py.File('transition_probs.h5', 'w')
    for (orient, val) in fastq_data.items():
        error_trans = compute_error_models(val, barcodes[orient].unique(), orient=orient)
        transition_h5.create_dataset(orient, data=error_trans, dtype='f4')
    transition_h5.close()
    """
    
}

process MAP_INDEX_TO_SAMPLE {
    publishDir "${params.outdir}/interm/sample_index_mapping", mode: "copy"
    tag "$split"
    label "process_high"
    stageInMode "copy"

    container "nakor/metaflowmics-python:0.0.2"
    conda (params.enable_conda ? "conda-forge::pandas conda-forge::h5py" : null)

    input:
    tuple val(split), path(h5)
    path error_model
    path meta_file

    output:
    tuple val(split), path("demux*.tsv"), emit: tsv
    path "sample_counts*.csv", emit: counts

    script:
    """
    #!/usr/bin/env python

    import pandas as pd
    import numpy as np
    import h5py


    NUCLS = list('ACGTN')
    MAPPING_BASE5 = str.maketrans(''.join(NUCLS), '01234')

    def seq2str(seq):
        if hasattr(seq, 'seq'):
            return str(seq.seq)
        return seq

    def seq2intlist(bc):
        bc_vec = np.fromiter(seq2str(bc).translate(MAPPING_BASE5), dtype='uint8')
        return bc_vec

    def idx2seq(nb):
        repr_b5 = "{:08}".format(int(np.base_repr(nb, 5)))
        nucl = repr_b5.translate(MAPPING_BASE5)
        return nucl

    def from_h5(filename, field=None):
        handle = h5py.File(filename, 'r')
        keys = handle.keys()

        if field is not None:
            keys = ["{}/{}".format(field, key) 
                    for key in handle.get(field).keys()]

        data = {key.split('/')[-1]: np.array(handle.get(key)[:]) for key in keys}

        handle.close()

        return data

    def calc_probs(sequences, qualities, error_model, barcodes, max_mismatches=2):
        '''
        data: index read data (seq index, codes and quals)
        '''

        print('Loading data')
        transition_probs = from_h5(error_model)
        orientations = transition_probs.keys()

        print('Computing likelihood')

        bc_sequences = {orient: np.array([seq2intlist(bc) for bc in barcodes[orient]])
                        for orient in orientations}

        transitions_all = [transition_probs[orient][
            bc_sequences[orient][:, None],
            sequences[orient],
            qualities[orient]-1].prod(axis=2) for orient in orientations]

        if len(transitions_all) == 1:
            assignments = transitions_all[0].argmax(axis=0)
        else:
            assignments = np.multiply(*transitions_all).argmax(axis=0)

        diffs = sum([bc_sequences[orient][assignments, :] != sequences[orient]
                     for orient in orientations]).sum(axis=1)

        assignments = barcodes.iloc[assignments].reset_index()
        assignments.insert(1, 'mismatches', diffs)
        assignments.loc[diffs > max_mismatches, 'sample_name'] = np.nan

        return assignments.values

    def format_result(assignments, rids, seq):

        trans = str.maketrans('01234', 'ACGTN')
        index_string = np.array(
            [[''.join(x) for x in np.core.defchararray.translate(codes.astype('U1'), trans)]
            for codes in seq.values()])

        result = np.vstack([rids[None, :], index_string, assignments[:, [0, 1]].T]).T

        return pd.DataFrame(result)
    
    # main script

    barcodes = pd.read_csv("$meta_file", dtype=str, names=["sample_name", "fwd", "rev"]).dropna(how='all', axis=1)
    barcodes.sample_name = barcodes.sample_name.str.replace("[^a-zA-Z0-9_.]", "_", regex=True)
    barcodes.set_index('sample_name', inplace=True)

    read_ids = from_h5("$h5", field='rid')['fwd'].astype(str)
    sequences = from_h5("$h5", field='seq')
    quals = from_h5("$h5", field='qual')

    assignments = calc_probs(sequences, quals, "$error_model", barcodes, max_mismatches=$params.max_mismatches)
    result = format_result(assignments, read_ids, sequences)

    result.to_csv("demux_info_${split}.tsv", sep='\\t', header=False, index=False)

    sample_sizes = result[result.columns[-2]].astype(str).value_counts()

    (
        sample_sizes[~sample_sizes.index.isnull()]
        .to_csv("sample_counts_${split}.csv", header=False)
    )
    
    
    """
    
}

process SAMPLE_SIZE_DISTRIBUTION {
    publishDir "${params.outdir}/figures", mode: "copy"
    label "process_low"

    container "nakor/metaflowmics-python:0.0.2"
    conda (params.enable_conda ? "conda-forge::pandas conda-forge::bokeh" : null)

    input:
    path counts

    output:
    path "*.html"

    script:
    """
    #!/usr/bin/env python

    from math import log10
    from pathlib import Path

    import numpy as np
    import pandas as pd

    from bokeh.plotting import figure
    from bokeh.io import save, output_file
    from bokeh.models import tickers


    def count_samples(folder='.'):
        '''
        Count the number of read pair/sample in [folder]
        '''

        files = Path(folder).glob("sample_counts*.csv")
        summaries = pd.concat([pd.read_csv(summary, index_col=0, dtype=str, header=None)
                               for summary in files],
                              axis=1, sort=True)

        return summaries.fillna(0).astype(int).sum(axis=1)

    def plot_bokeh(counts):

        hist, edges = np.histogram(np.log10(counts), bins=max(5, len(counts)//10), density=False)
        hist_df = pd.DataFrame({'count': hist,
                                "left": edges[:-1],
                                "right": edges[1:]})
        hist_df["interval"] = ["{:,} - {:,}".format(int(10**left), int(10**right))
                               for left, right in zip(hist_df["left"], hist_df["right"])]

        x_min = int(min(edges))
        x_max = max(4, 1+int(max(edges)))

        p = figure(plot_height=800, plot_width=800,
                   x_range=[x_min, x_max], tools='hover,box_zoom',
                   tooltips=[('Size range', '@interval'),
                             ('#Samples in interval', str("@count"))],
                    title="Sample size distribution",
                    x_axis_label="Sample read count",
                    y_axis_label="Occurrences")

        p.quad(bottom=0, top="count", left="left", 
               right="right", source=hist_df, fill_color="SteelBlue", 
               line_color="black", fill_alpha=0.7,
               hover_fill_alpha=1.0, hover_fill_color="Tan")

        ticks = list(range(x_min, x_max))
        minor_ticks = np.log10([i*10**j for i in range(1, 10) for j in ticks])

        p.xaxis.ticker = tickers.FixedTicker(ticks=ticks, minor_ticks=minor_ticks)
        p.xaxis.major_label_overrides = {tick: "{:,}".format(int(10**tick)) for tick in ticks}
        p.yaxis.minor_tick_line_color = None

        p.axis.major_label_text_font_size = "12pt"
        p.axis.axis_label_text_font_size = "14pt"
        p.title.text_font_size = "18pt"

        output_file('sample_sizes.html')
        save(p)

    # main script
    counts = count_samples()
    counts.to_csv('sample_sizes.csv', header=False)
    plot_bokeh(counts)
    """
    
}

process WRITE_SAMPLES_TO_FASTQ {
    tag "$split"
    label "process_low"

    container "nakor/metaflowmics-python:0.0.2"
    conda (params.enable_conda ? "conda-forge::biopython conda-forge::pandas" : null)

    input:
    tuple val(split), path(fastq), path(mapping)

    output:
    path "*.fastq"

    script:
    """
    #!/usr/bin/env python

    from Bio.SeqIO.QualityIO import FastqGeneralIterator
    import pandas as pd

    demux_info = pd.read_csv("$mapping", header=None, sep="\\t", dtype=str, index_col=0)
    
    if demux_info[2].count() == 0:
        demux_info.drop(columns=[2], inplace=True)
        demux_info.columns = ["fwd", "sample_name", "mismatches"]
    else:
        demux_info.columns = ["fwd", "rev", "sample_name", "mismatches"]

    read_orient = ["fwd", "rev"][:${fastq.size()}]
    
    print("Preparing handles.")
    handles = {}
    for sample in demux_info["sample_name"].unique():
        if not pd.isnull(sample):
            for i, orient in enumerate(read_orient, 1):
                handles[sample+orient] = open(f"{sample}_R{i}.fastq", "w")

    parsers = [FastqGeneralIterator(open(fastq, "r")) for fastq in "$fastq".split()]

    print("Starting demultiplexing")
    for seq_nb, sequences in enumerate(zip(*parsers)):
        ids = [seq[0].split()[0] for seq in sequences]

        if len(ids) > 1:
            if ids[0] != ids[1]:
                print("Sequence #{}: {} (fwd) and {} (rev) do not match. The forward and reverse read files seem to be out of order"
                      .format(seq_nb, *ids))
                exit(42)

        sample_assignment = demux_info.loc[ids[0], "sample_name"]

        if pd.isnull(sample_assignment):
            continue

        for orient, seq in zip(read_orient, sequences):
            handles[sample_assignment + orient].write("@{}\\n{}\\n+\\n{}\\n".format(*seq))

    for sample in demux_info["sample_name"].unique():
        if not pd.isnull(sample):
            for orient in read_orient:
                handles[sample + orient].close()

    print("Demultiplexing finished.")
    """
    
}

process GZIP {
    publishDir "${params.outdir}/fastqs", mode: "copy"
    label "high_computation"
    conda (params.enable_conda ? "conda-forge::bash='5.1'" : null)
    container "nakor/bash:5.1.4"

    input:
	path f

    output:
    path "*.gz"
	
    script:
    """
    ls *.fastq | xargs -P $task.cpus -I % bash -c "gzip -c % > %.gz"
    """	
}

