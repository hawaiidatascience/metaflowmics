#!/usr/bin/env bash

pairId=all

for arg in "$@"
do
    case $arg in
    --step=*)
	step="${arg#*=}" ;;
    
    --pairId=*)
	pairId="${arg#*=}" ;;

    --rad=*)
	rad="${arg#*=}" ;;

    --optimize=*)
	optimize="${arg#*=}" ;;

    --criteria=*)
	criteria="${arg#*=}" ;;

    --refAln=*)
	refAln="${arg#*=}" ;;
    
    --refTax=*)
	refTax="${arg#*=}"
	taxRad=`basename $refTax .tax`;;

    --taxaToFilter=*)
	taxaToFilter="${arg#*=}" ;;

    --idThreshold=*)
	idThreshold="${arg#*=}" ;;
    
    --subsamplingNb=*)
	subsamplingNb="${arg#*=}" ;;

    *)
	echo "$arg: Unknown option";;    
    esac
done

# General prefix for output names
out="${pairId}.${step}" 
    
if [ $step == "screening" ]; then
    output_suffix=`echo ${optimize} | sed s/-/_/g`
    out="${out}.${output_suffix}"

    # We need to keep track of the input in case mothur does not create an output in case mothur does not create one.
    inputs_mothur=("all.esv.count_table" "${rad}.fasta") 
    cmd=("screen.seqs(fasta=${rad}.fasta, count=all.esv.count_table, optimize=${optimize}, criteria=${criteria})")
    outputs_mothur=("all.esv.good.count_table" "${rad}.good.fasta")
    outputs_renamed=("${out}.count_table" "${out}.fasta")
	 
elif [ $step == "summary" ]; then
    cmd=("summary.seqs(fasta=${rad}.fasta)")
   
elif [ $step == "MSA" ]; then
    cmd=("align.seqs(fasta=${rad}.fasta, reference=${refAln}) ; "
	 "filter.seqs(fasta=${rad}.align)")
    outputs_mothur=("${rad}.filter.fasta")
    outputs_renamed=("${out}.fasta")

elif [ $step == "chimera" ]; then
    cmd=("chimera.vsearch(fasta=${rad}.fasta, count=${rad}.count_table, dereplicate=t) ; "
	 "remove.seqs(fasta=${rad}.fasta, count=${rad}.count_table, accnos=${rad}.denovo.vsearch.accnos, dups=f)")
    outputs_mothur=("${rad}.pick.fasta" "${rad}.pick.count_table")
    outputs_renamed=("${out}.fasta" "${out}.count_table")

elif [ $step == "taxaFilter" ]; then 
    suffixTax=`echo $taxRad | cut -d. -f2`.wang.taxonomy
    
    cmd=("classify.seqs(fasta=${rad}.fasta, count=${rad}.count_table, template=${refAln}, taxonomy=${refTax})")
    outputs_mothur=("${rad}.${suffixTax}")
    outputs_renamed=("${out}.taxonomy")

    if [ ! -z ${taxaToFilter} ]; then 
        cmd+=("; remove.lineage(taxonomy=${rad}.${suffixTax}, count=${rad}.count_table, fasta=${rad}.fasta, taxon=-${taxaToFilter})")
	outputs_mothur+=("${rad}.pick.count_table" "${rad}.pick.fasta")
	outputs_renamed+=("${out}.count_table" "${out}.fasta")
    fi

elif [ $step == "subsampling" ]; then
    cmd=("sub.sample(persample=true, fasta=${rad}.fasta, count=${rad}.count_table, taxonomy=${rad}.taxonomy, size=${subsamplingNb})")
    outputs_mothur=("${rad}.subsample.fasta" "${rad}.subsample.count_table" "${rad}.subsample.taxonomy")
    outputs_renamed=("${out}.fasta" "${out}.count_table" "${out}.taxonomy")

elif [ $step == "clustering" ]; then
    
    if [ ${idThreshold} = 0 ]; then
	method="unique"
    else
	method="dgc"
    fi
    
    cmd=("cluster(count=${rad}.count_table, fasta=${rad}.fasta, method=${method}, cutoff=${idThreshold}) ; "
	 "make.shared(list=${rad}.${method}.list, count=${rad}.count_table) ; "
	 "get.oturep(count=${rad}.count_table, fasta=${rad}.fasta, list=${rad}.${method}.list, method=abundance, rename=T, label=${idThreshold})")
	
    outputs_mothur=("${rad}.${method}.list" "${rad}.${method}.shared" "${rad}.${method}.${idThreshold}.rep.fasta" "${rad}.${method}.${idThreshold}.rep.count_table")

    outputs_renamed=("${out}.${idThreshold}.list" "${out}.${idThreshold}.shared" "${out}.${idThreshold}.fasta" "${out}.${idThreshold}.count_table" )
    
elif [ $step == "consensusClassification" ]; then
    cmd=("classify.otu(taxonomy=all.subsampling.taxonomy, count=all.subsampling.count_table, list=${rad}.list)")
    outputs_mothur=("all.clustering.${idThreshold}.${idThreshold}.cons.taxonomy" "all.clustering.${idThreshold}.${idThreshold}.cons.tax.summary")
    outputs_renamed=("all.classification.${idThreshold}.cons.taxonomy" "all.classification.${idThreshold}.cons.summary")
fi;

# Process mothur cmd

# Join the instruction array
res=$( IFS=$' '; echo "${cmd[*]}" )

# Execute
mothur "#${res}"

# Rename output files
set -o xtrace

n=$((${#outputs_mothur[@]}-1))
for i in `seq 0 $n`
do
    if [ -e ${outputs_mothur[$i]} ]; then
	mv ${outputs_mothur[$i]} ${outputs_renamed[$i]}
    elif [ $step = "screening" ]; then
	echo "${outputs_mothur[$i]} does not exist. Renaming."
	cp ${inputs_mothur[$i]} ${outputs_renamed[$i]}
    else
	echo "${outputs_mothur[$i]} does not exist. Aborting."
	exit 1
    fi
done
