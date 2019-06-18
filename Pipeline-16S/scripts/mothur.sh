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
out="${pairId}_${step}"

if [ ! -z $idThreshold ]; then
    if [ $idThreshold -eq 100 ]; then
	mothurThresh=0
    else
	mothurThresh=0.0"$((100-${idThreshold}))"	
    fi
    out="${out}_${idThreshold}"
fi

if [ $step == "MSA" ]; then
    inputs_mothur=("${fasta}.fasta" "${count}.count_table")
    
    # output_suffix=`echo ${optimize} | sed s/-/./g`
    # out="${out}_${output_suffix}"
    
    cmd=("align.seqs(fasta=${fasta}.fasta, reference=${refAln}) ; "
	 "filter.seqs(fasta=${fasta}.align) ; "
	 "screen.seqs(fasta=${fasta}.filter.fasta, count=${count}.count_table, optimize=${optimize}, criteria=${criteria}) ; "
	 "summary.seqs(fasta=${fasta}.filter.good.fasta)")
    
    outputs_mothur=("${fasta}.filter.good.fasta" 
		    "${count}.good.count_table")
    
    outputs_renamed=("${out}.fasta"
		     "${out}.count_table")

elif [ $step == "chimera" ]; then
    method="vsearch"

    cmd=("chimera.${method}(fasta=${fasta}.fasta, count=${count}.count_table, dereplicate=t) ; "
	 "remove.seqs(fasta=${fasta}.fasta, count=${count}.count_table, accnos=${fasta}.denovo.${method}.accnos, dups=f)")
    
    outputs_mothur=("${fasta}.pick.fasta"
		    "${count}.pick.count_table")
    
    outputs_renamed=("${out}.fasta"
		     "${out}.count_table")

elif [ $step == "subsampling" ]; then
    cmd=("sub.sample(persample=true, fasta=${fasta}.fasta, count=${count}.count_table, size=${subsamplingNb})")
    
    outputs_mothur=("${fasta}.subsample.fasta"
		    "${count}.subsample.count_table")
    
    outputs_renamed=("${out}.fasta"
		     "${out}.count_table")

elif [ $step == "clustering" ]; then
    
    if [ ${idThreshold} -eq 100 ]; then
	method="unique"
    else
	method="dgc"
    fi

    cmd=("cluster(count=${count}.count_table, fasta=${fasta}.fasta, method=${method}, cutoff=${mothurThresh}) ; "
	 "make.shared(list=${fasta}.${method}.list, count=${count}.count_table) ; "
	 "get.oturep(count=${count}.count_table, fasta=${fasta}.fasta, list=${fasta}.${method}.list, method=abundance, rename=T, label=${mothurThresh}) ; "
	)
	
    outputs_mothur=("${fasta}.${method}.list"
		    "${fasta}.${method}.shared"
		    "${fasta}.${method}.${mothurThresh}.rep.fasta"
		    "${count}.${method}.${mothurThresh}.rep.count_table")
    
    outputs_renamed=("${out}.list" "${out}.shared" "${out}.fasta" "${out}.count_table" )
    
elif [ $step == "consensusClassification" ]; then
    suffixTax=`echo $taxRad | cut -d. -f2`.wang
    
    cmd=("classify.seqs(fasta=${fasta}.fasta, count=${count}.count_table, template=${refAln}, taxonomy=${refTax}) ; "
	 "classify.otu(taxonomy=${fasta}.${suffixTax}.taxonomy, count=${count}.count_table, list=${list}.list, probs=f)")
    
    outputs_mothur=("${list}.${mothurThresh}.cons.taxonomy"
		    "${list}.${mothurThresh}.cons.tax.summary")
    
    outputs_renamed=("${out}.taxonomy"
		     "${out}.summary")

elif [ $step == "taxaFilter" ]; then
    inputs_mothur=("${tax}.taxonomy" "${shared}.shared" "${fasta}.fasta")
    
    cmd=("remove.lineage(constaxonomy=${tax}.taxonomy, shared=${shared}.shared, taxon='${taxaToFilter}', label=${idThreshold}) ; "
	 "list.seqs(taxonomy=${tax}.pick.cons.taxonomy) ; "
	 "get.seqs(fasta=${fasta}.fasta,accnos=current)")
    
    outputs_mothur=("${tax}.pick.cons.taxonomy"
		    "${shared}.${idThreshold}.pick.shared"
		    "${fasta}.pick.fasta")
    
    outputs_renamed=("annotations_${idThreshold}.taxonomy"
		     "abundance_table_${idThreshold}.shared"
		     "sequences_${idThreshold}.fasta")

elif [ $step == "postprocessing" ]; then
    cmd=("get.relabund(shared=${shared}.shared) ; "
	 "clearcut(fasta=${fasta}.fasta, DNA=T) ; "
	 "count.seqs(shared=${shared}.shared) ; "
	 "unifrac.weighted(tree=current,count=current) ; "
	 "summary.single(shared=${shared}.shared,calc=nseqs-sobs-chao-shannon-shannoneven) ; "
	 "summary.shared(shared=${shared}.shared,calc=braycurtis-thetayc-sharedsobs-sharedchao)")
fi;

# Process mothur cmd

# Join the instruction array
res=$( IFS=$' '; echo "${cmd[*]}" )
echo $res

# Execute
mothur "#${res}"

# Rename output files
n=$((${#outputs_mothur[@]}-1))
for i in `seq 0 $n`
do
    # if the output file exists
    if [ -e ${outputs_mothur[$i]} ]; then
	mv ${outputs_mothur[$i]} ${outputs_renamed[$i]}
    # Special case when screen.seqs (sometime mothur doesnt produce an output file). In this case, just copy the input into the output
    elif [ $step = "MSA" ] || [ $step = "taxaFilter" ]; then
	echo "${outputs_mothur[$i]} does not exist. Renaming."
	cp ${inputs_mothur[$i]} ${outputs_renamed[$i]}
    # Otherwise, raise an error
    else
	echo "${outputs_mothur[$i]} does not exist. Aborting."
	exit 1
    fi
done
