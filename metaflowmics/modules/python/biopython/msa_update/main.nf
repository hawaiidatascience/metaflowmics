// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

options = initOptions(params.options)


process UPDATE_MSA_WITH_REF {
    label "process_low"
    publishDir params.outdir, mode: params.publish_dir_mode

    container "nakor/metaflowmics-python:0.0.2"
    conda (params.enable_conda ? "conda-forge::biopython conda-forge:pandas" : null)

    input:
    path afa
    path ref
	path updated_ref
	path freqs

    output:
    path "*.updated.afa", emit: afa
	path "backprop.txt", optional true, emit: txt

    script:
    """
    #!/usr/bin/env python

	from functools import lru_cache
	from itertools import groupby

	from Bio.SeqIO.FastaIO import SimpleFastaParser
	import numpy as np
	import pandas as pd


	def make_template(consensus, updated, counts):
	    ref_itv = get_intervals(consensus)
        other_itv = get_intervals(updated)

		bins = pd.cut(counts.index, pd.IntervalIndex(other_itv))
		counts = counts.groupby(bins).agg(lambda x: list(x)[1:])
		counts.index = pd.IntervalIndex(ref_itv)

		(reference_itv, contracted_idx) = map_intervals(ref_itv, other_itv)

	    return (reference_itv, counts, contracted_idx)

	def get_intervals(seq, gap_chars={"-", "."}):
		'''
		Get the location and distance between each non-gap characters
		'''
		indices = [0]+[i for (i, c) in enumerate(seq) if c not in gap_chars] + [len(seq)]
		intervals = pd.arrays.IntervalArray(
			[pd.Interval(i, j) for (i, j) in zip(indices[:-1], indices[1:])],
			closed='left'
		)

		return intervals

	def map_intervals(ref_itv, other_itv, gap_chars={"-", "."}):
		'''
		Map intervals of representative sequences on its aligned version
		Report contracted intervals that need to backpropagate in the other alignments at this level
		'''
		# extend contracted intervals
		contracted = other_itv.length < ref_itv.length
		contracted_idx = np.arange(len(other_itv))[contracted]

		# adjust the interval lengths
		ref_itv_adj = pd.Series(ref_itv.length, index=ref_itv)
		ref_itv_adj.loc[~contracted] = other_itv.length[~contracted]

		return (ref_itv_adj, contracted_idx)

	def map_on_template(template, query, counts, gap_chars={"-", "."}):
		'''
		Map an alignment on the updated template
		The difficulty is that we do not know where new gap were inserted
		among the already existing gaps. There are multiple ways to map.
		We optimize the alignment using the letter frequencies at the 
		undecided positions
		'''
		query_aligned = ""

		for (itv, length) in template.iteritems():

			q_sub = query[itv.left:itv.right]
			length_diff = length-len(q_sub)

			if not q_sub:
				query_aligned += '-'*length
				continue

			query_aligned += q_sub[0]
			q_sub = q_sub[1:]            

			# Case 1: template and query have the same length
			if length_diff == 0:
				query_aligned += q_sub
			# Case 2: template is longer than query
			elif length_diff > 0:
				letters = "".join(c for c in q_sub if c not in gap_chars)

				if not letters: # we only have gaps: we just add enough
					query_aligned += "-" * (length-1)
				if letters: # we need to optimize the placement
					opt_aln = find_optimal_mapping(
						letters, length-1, counts.loc[itv]
					)
					query_aligned += opt_aln
			else: # Should already have been handled when mapping intervals
				raise ValueError("There should not be any contractions at this point")

		return query_aligned

	def find_optimal_mapping(letters, total_len, counts):
		'''
		Recursive function to find the optimal placement of letters
		such that the total string is length <total_len>. The optimal
		location is defined as the one for which the sum of occurences
		for each letter at the chosen column is maximal
		'''
		@lru_cache
		def bfs(letters, ngaps, pos=1):
			'''
			2 options: 1) we put down a letter
					   2) we put down a gap
			We stop when we have no more letters or no more gaps
			'''
			if not letters:
				return (0, [])

			(score1, sol1) = bfs(letters[1:], ngaps, pos+1)
			score1 += counts[letters[0]][pos]

			if ngaps > 0:
				(score2, sol2) = bfs(letters, ngaps-1, pos+1)
				if score2 > score1:
					return (score2, sol2)

			return (score1, [pos]+sol1)

		(score, positions) = bfs(letters, total_len-len(letters)-1, pos=0)
		letter_map = dict(zip(positions, letters))
		opt_aln = "".join(letter_map.get(i, '-') for i in range(total_len))

		return opt_aln


	refs = {}
	with open("$ref") as reader:
		for (title, seq) in SimpleFastaParser(reader):
			lineage = title.split()[1]
			tax = lineage.split("$params.sep")[$params.field]
			refs[tax] = seq

	with open("$afa") as reader:
		for (title, seq) in SimpleFastaParser(reader):
			lineage = title.split()[1]
			tax = lineage.split("$params.sep")[$params.field]

			seq_updated = map_on_ref(refs[tax], seq)

			with open(f"{tax}.updated.afa", "a") as writer:
				writer.write(f">{title}\\n{seq_updated}\\n")    

    """    
}

// TODO
// 1) Loop through the "$ref" and "$updated_ref", make the templates and save the backprop
// 2) Apply the backprop to the templates?
// 3) Construct the frequency table from the updated ref
// 4) Update alignments

