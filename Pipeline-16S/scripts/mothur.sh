#!/usr/bin/env bash

getRad() {
    files=$(ls -1 *$1 2>/dev/null)
    if [ ! -z files ]; then
	echo `basename $files $1`
    fi
}

set -o xtrace

fasta=`getRad .fasta`
shared=`getRad .shared`
count=`getRad .count_table`
tax=`getRad .taxonomy`
list=`getRad .list`

set +o xtrace

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
    inputs_to_copy=("${fasta}.fasta" "${count}.count_table")
    
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

elif [ $step == "preClassification" ]; then
    suffixTax=`echo $taxRad | cut -d. -f2`.wang
    
    cmd=("classify.seqs(fasta=${fasta}.fasta, count=${count}.count_table, template=${refAln}, taxonomy=${refTax})")
    
    outputs_mothur=("${fasta}.${suffixTax}.taxonomy")
    outputs_renamed=("${out}.taxonomy")
    
elif [ $step == "clustering" ]; then
    method=`[ ${idThreshold} -eq 100 ] && echo 'unique' || echo 'dgc'`

    cmd=("cluster(count=${count}.count_table, fasta=${fasta}.fasta, method=${method}, cutoff=${mothurThresh}) ; "
	 "make.shared(list=${fasta}.${method}.list, count=${count}.count_table)")
	
    outputs_mothur=("${fasta}.${method}.shared"
		    "${fasta}.${method}.list")
    
    outputs_renamed=("${out}.shared" "${out}.list")

elif [ $step == "consensusClassification" ]; then
    suffixTax=`echo $taxRad | cut -d. -f2`.wang
    
    cmd=("classify.otu(taxonomy=${tax}.taxonomy, count=${count}.count_table, list=${list}.list, probs=f)")
    
    outputs_mothur=("${list}.${mothurThresh}.cons.taxonomy"
		    "${list}.${mothurThresh}.cons.tax.summary")
    
    outputs_renamed=("${out}.taxonomy"
		     "${out}.summary")

elif [ $step == "otuRepr" ]; then
    cmd=("get.oturep(count=${count}.count_table, fasta=${fasta}.fasta, list=${list}.list, method=abundance, rename=T)")

    outputs_mothur=("${list}.${mothurThresh}.rep.fasta")
    outputs_renamed=("${out}.fasta")

elif [ $step == "subsampling" ]; then
    cmd=("sub.sample(persample=true, fasta=${fasta}.fasta, constaxonomy=${tax}.taxonomy, shared=${shared}.shared, size=${subsamplingNb})")
    inputs_to_copy=("${fasta}.fasta")
    
    outputs_mothur=("${fasta}.subsample.fasta"
		    "${tax}.subsample.taxonomy"
		    "${shared}.${mothurThresh}.subsample.shared")
    
    outputs_renamed=("${out}.fasta"
		     "${out}.taxonomy"
		     "${out}.shared")

elif [ $step == "taxaFilter" ]; then
    # 2 cases depending on whether we do the taxa filtering on the taxonomy or constaxonomy file
    # For the latter case, we need to also process the .shared file 
    if [ -f "${shared}.list" ]; then
	inputs_to_copy=("${tax}.taxonomy" "${list}.list" "${shared}.shared")

        outputs_mothur=("${tax}.pick.cons.taxonomy"
			"${list}.${mothurThresh}.pick.list"
			"${shared}.${mothurThresh}.pick.shared")
			
	outputs_renamed=("${out}.taxonomy"
			 "${out}.list"
			 "${out}.shared")
	
	cmd=("remove.lineage(constaxonomy=${tax}.taxonomy, shared=${shared}.shared, list=${list}.list, taxon='${taxaToFilter}') ; ")
    else
	inputs_to_copy=("${fasta}.fasta" "${shared}.shared" "${tax}.taxonomy")

        outputs_mothur=("${fasta}.pick.fasta"
			"${shared}.pick.shared"
			"${tax}.pick.cons.taxonomy")
	
	outputs_renamed=("${out}.fasta"
			 "${out}.shared"
			 "${out}.taxonomy")
	
	cmd=("remove.lineage(constaxonomy=${tax}.taxonomy, shared=${shared}.shared, taxon='${taxaToFilter}') ; "
	     "list.seqs(taxonomy=${tax}.pick.cons.taxonomy) ; "
	     "get.seqs(fasta=${fasta}.fasta,accnos=current)")

    fi

elif [ $step == "postprocessing" ]; then
    cmd=("get.relabund(shared=${shared}.shared) ; "
	 "clearcut(fasta=${fasta}.fasta, DNA=T) ; "
	 "count.seqs(shared=${shared}.shared) ; "
	 "unifrac.weighted(tree=current,count=current,distance=lt) ; "
	 "summary.single(shared=${shared}.shared,calc=nseqs-sobs-chao-shannon-shannoneven) ; "
	 "summary.shared(shared=${shared}.shared,calc=braycurtis-thetayc-sharedsobs-sharedchao)")
fi;

# Process mothur cmd

# Join the instruction array
res=$( IFS=$' '; echo "${cmd[*]}" )
echo $res

# Execute
set -o xtrace
[ -z $MOTHUR ] && mothur "#${res}" || $MOTHUR/mothur "#${res}"
set +o xtrace

# Rename output files
n=$((${#outputs_mothur[@]}-1))
for i in `seq 0 $n`
do
    # if the output file exists
    if [ -e ${outputs_mothur[$i]} ]; then
	echo "Success: Renaming ${outputs_mothur[$i]} to ${outputs_renamed[$i]}"
	mv ${outputs_mothur[$i]} ${outputs_renamed[$i]}
    # Special case when screen.seqs (sometime mothur doesnt produce an output file). In this case, just copy the input into the output
    elif [ $step = "MSA" ] || [ $step = "taxaFilter" ] || [ $step = "subsampling" ]; then
	echo "WARNING: ${outputs_mothur[$i]} does not exist. Copying input (${inputs_to_copy[$i]})."
	cp ${inputs_to_copy[$i]} ${outputs_renamed[$i]}
    # Otherwise, raise an error
    else
	echo "ERROR: ${outputs_mothur[$i]} does not exist. Aborting."
	exit 1
    fi
done

