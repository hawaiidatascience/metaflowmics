#!/usr/bin/env bash

set -o xtrace

getRad() {
    files=$(ls -1 *$1 2>/dev/null)
    if [ ! -z files ]; then
	echo `basename $files $1`
    fi
}

fasta=`getRad .fasta`
shared=`getRad .shared`
count=`getRad .count_table`
tax=`getRad .taxonomy`
list=`getRad .list`

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
    
if [ $step == "MSA" ]; then
    inputs_mothur=("{count}.count_table" "${fasta}.fasta")
    
    output_suffix=`echo ${optimize} | sed s/-/_/g`
    out="${out}.${output_suffix}"
    
    cmd=("align.seqs(fasta=${fasta}.fasta, reference=${refAln}) ; "
	 "filter.seqs(fasta=${fasta}.align) ; "
	 "screen.seqs(fasta=${fasta}.filter.fasta, count=${count}.count_table, optimize=${optimize}, criteria=${criteria}) ; "
	 "summary.seqs(fasta=${fasta}.filter.good.fasta)")
    outputs_mothur=("${fasta}.filter.good.fasta" "${count}.good.count_table")
    outputs_renamed=("${out}.fasta" "${out}.count_table")

elif [ $step == "chimera" ]; then
    method="uchime"
    cmd=("chimera.${method}(fasta=${fasta}.fasta, count=${count}.count_table, dereplicate=t) ; "
	 "remove.seqs(fasta=${fasta}.fasta, count=${count}.count_table, accnos=${fasta}.denovo.${method}.accnos, dups=f)")
    outputs_mothur=("${fasta}.pick.fasta" "${count}.pick.count_table")
    outputs_renamed=("${out}.fasta" "${out}.count_table")

elif [ $step == "taxaFilter" ]; then
    cmd=("remove.lineage(constaxonomy=${tax}.taxonomy, shared=${shared}.shared, taxon='${taxaToFilter}', label=${idThreshold}) ; "
	 "list.seqs(taxonomy=${tax}.pick.taxonomy) ; "
	 "get.seqs(fasta=${fasta}.fasta,accnos=current)")
    outputs_mothur=("${tax}.pick.taxonomy" "${shared}.${idThreshold}.pick.shared" "${fasta}.pick.fasta")
    outputs_renamed=("${out}.taxonomy" "${out}.shared" "${out}.fasta")

elif [ $step == "subsampling" ]; then
    cmd=("sub.sample(persample=true, fasta=${fasta}.fasta, count=${count}.count_table, size=${subsamplingNb})")
    outputs_mothur=("${fasta}.subsample.fasta" "${count}.subsample.count_table")
    outputs_renamed=("${out}.fasta" "${out}.count_table")

elif [ $step == "clustering" ]; then
    
    if [ ${idThreshold} = 0 ]; then
	method="unique"
    else
	method="dgc"
    fi
    
    cmd=("cluster(count=${count}.count_table, fasta=${fasta}.fasta, method=${method}, cutoff=${idThreshold}) ; "
	 "make.shared(list=${fasta}.${method}.list, count=${count}.count_table) ; "
	 "get.oturep(count=${count}.count_table, fasta=${fasta}.fasta, list=${fasta}.${method}.list, method=abundance, rename=T, label=${idThreshold}) ; "
	 "filter.seqs(fasta=${fasta}.${method}.${idThreshold}.rep.fasta, trump=-)")
	
    outputs_mothur=("${fasta}.${method}.list" "${fasta}.${method}.shared" "${fasta}.${method}.${idThreshold}.rep.filter.fasta" "${count}.${method}.${idThreshold}.rep.count_table")
    outputs_renamed=("${out}.${idThreshold}.list" "${out}.${idThreshold}.shared" "${out}.${idThreshold}.fasta" "${out}.${idThreshold}.count_table" )
    
elif [ $step == "consensusClassification" ]; then
    suffixTax=`echo $taxRad | cut -d. -f2`.wang
    cmd=("classify.seqs(fasta=${fasta}.fasta, count=${count}.count_table, template=${refAln}, taxonomy=${refTax}) ; "
	 "classify.otu(taxonomy=${fasta}.${suffixTax}.taxonomy, count=${count}.count_table, list=${list}.list)")
    outputs_mothur=("all.clustering.${idThreshold}.${idThreshold}.cons.taxonomy" "all.clustering.${idThreshold}.${idThreshold}.cons.tax.summary")
    outputs_renamed=("all.classification.${idThreshold}.cons.taxonomy" "all.classification.${idThreshold}.cons.summary")
fi;

# Process mothur cmd

# Join the instruction array
res=$( IFS=$' '; echo "${cmd[*]}" )

# Execute
mothur "#${res}"

# Rename output files

n=$((${#outputs_mothur[@]}-1))
for i in `seq 0 $n`
do
    if [ -e ${outputs_mothur[$i]} ]; then
	mv ${outputs_mothur[$i]} ${outputs_renamed[$i]}
    elif [ $step = "MSA" ]; then
	echo "${outputs_mothur[$i]} does not exist. Renaming."
	cp ${inputs_mothur[$i]} ${outputs_renamed[$i]}
    else
	echo "${outputs_mothur[$i]} does not exist. Aborting."
	exit 1
    fi
done
