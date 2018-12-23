#!/usr/bin/env bash

pairId=all

for arg in "$@"
do
    case $arg in
    --step=*)
	step="${arg#*=}" ;;
    
    --pairId=*)
	pairId="${arg#*=}" ;;

    --fwdFasta=*)
	fwdFasta="${arg#*=}"
	noGzRad=`basename $fwdFasta .gz`
	ext="${noGzRad##*.}"
	faRad=`basename $noGzRad .$ext`;;

    --revFasta=*)
	revFasta="${arg#*=}" ;;

    --inputFasta=*)
	inputFasta="${arg#*=}"
	noGzRad=`basename $inputFasta .gz`
	ext="${noGzRad##*.}"
	faRad=`basename $noGzRad .$ext`;;

    --inputNames=*)
	inputNames="${arg#*=}"
	namesRad=`basename $inputNames .names`;;

    --inputGroups=*)
	inputGroups="${arg#*=}"
	groupsRad=`basename $inputGroups .groups`;;

    --inputList=*)
	inputList="${arg#*=}"
	listRad=`basename $inputList .list`;;

    --inputTax=*)
	inputTax="${arg#*=}"
	taxRad=`basename $inputTax .taxonomy`;;

    --inputCounts=*)
	inputCounts="${arg#*=}"
	countsRad=`basename $inputCounts .count_table`;;
		    
    --optimize=*)
	optimize="${arg#*=}" ;;

    --criteria=*)
	criteria="${arg#*=}" ;;

    --refAln=*)
	refAln="${arg#*=}" ;;
    
    --refTax=*)
	refTax="${arg#*=}"
	refTaxRad=`basename $refTax .tax`;;
    
    --otuMethod=*)
	otuMethod="${arg#*=}" ;;
    
    --clustId=*)
	clustId="${arg#*=}" ;;

    *)
	echo "$arg: Unknown option";;    
    esac
done

out="${pairId}.${step}"

if [ $step == "merging" ]; then
    cmd=("make.contigs(ffasta=${fwdFasta}, rfasta=${revFasta})")
    outputs_mothur=("${faRad}.trim.contigs.fasta")
    outputs_renamed=("${out}.fasta")
    
elif [ $step == "screening" ]; then
    optimize_mod=`echo ${optimize} | sed s/-/_/g`
    if [ -z ${inputNames} ]; then
	cmd=("screen.seqs(fasta=${inputFasta}, optimize=${optimize}, criteria=${criteria})")
	outputs_mothur=("${faRad}.good.fasta")
	outputs_renamed=("${out}.${optimize_mod}.fasta")
    else
	# inputs_mothur=("${inputNames}" "${inputFasta}")
	cmd=("screen.seqs(fasta=${inputFasta}, name=${inputNames}, group=${inputGroups}, optimize=${optimize}, criteria=${criteria})")
	outputs_mothur=("${namesRad}.good.names" "${faRad}.good.fasta" "${groupsRad}.good.groups")
	outputs_renamed=("${out}.${optimize_mod}.names" "${out}.${optimize_mod}.fasta" "${out}.${optimize_mod}.groups")
    fi
	 
elif [ $step == "summary" ]; then
    cmd=("summary.seqs(fasta=${inputFasta})")

elif [ $step == "subsampling" ]; then
    cmd=("sub.sample(persample=true,fasta=${inputFasta},group=${inputGroups},name=${inputNames},taxonomy=${inputTax})")
    outputs_mothur=("${faRad}.subsample.fasta" "${groupsRad}.subsample.groups" "${taxRad}.subsample.names" "${taxRad}.subsample.taxonomy")
    outputs_renamed=("${out}.fasta" "${out}.groups" "${out}.names" "${out}.taxonomy")
    
elif [ $step == "dereplication" ]; then
    cmd=("unique.seqs(fasta=${inputFasta}) ; "
	 "count.seqs(name=${faRad}.names,group=${inputGroups})")
    outputs_mothur=("${faRad}.unique.fasta" "${faRad}.names" "${faRad}.count_table")
    outputs_renamed=("${out}.fasta" "${out}.names" "${out}.count_table")
   
elif [ $step == "MSA" ]; then
    cmd=("align.seqs(fasta=${inputFasta},reference=${refAln}) ; "
	 "filter.seqs(fasta=${faRad}.align)")
    outputs_mothur=("${faRad}.filter.fasta")
    outputs_renamed=("${out}.fasta")

elif [ $step == "chimera" ]; then
    cmd=("chimera.vsearch(fasta=${inputFasta}, name=${inputNames}, group=${inputGroups}, dereplicate=t) ; "
	 "remove.seqs(fasta=${inputFasta}, name=${inputNames}, group=${inputGroups}, accnos=${faRad}.denovo.vsearch.accnos)")
    outputs_mothur=("${faRad}.pick.fasta" "${namesRad}.pick.names" "${groupsRad}.pick.groups")
    outputs_renamed=("${out}.fasta" "${out}.names" "${out}.groups")

elif [ $step == "taxaFilter" ]; then
    suffixTax=`echo $refTaxRad | cut -d. -f2`.wang

    cmd=("classify.seqs(fasta=${inputFasta}, name=${inputNames}, template=${refAln}, taxonomy=${refTax}) ; " 
	 "remove.lineage(taxonomy=${faRad}.${suffixTax}.taxonomy, name=${inputNames}, group=${inputGroups}, fasta=${inputFasta}, taxon=-unknown)")
    outputs_mothur=("${namesRad}.pick.names" "${faRad}.pick.fasta" "${groupsRad}.pick.groups" "${faRad}.${suffixTax}.pick.taxonomy")
    outputs_renamed=("${out}.names" "${out}.fasta" "${out}.groups" "${out}.taxonomy")

elif [ $step == "clustering" ]; then
    cmd=("cluster(name=${inputNames}, fasta=${inputFasta}, method=${otuMethod}, cutoff=${clustId}) ; "
	 "get.oturep(name=${inputNames},fasta=${inputFasta},list=${faRad}.${otuMethod}.list,label=${clustId},method=abundance)")
    outputs_mothur=("${faRad}.${otuMethod}.list" "${faRad}.${otuMethod}.${clustId}.rep.fasta" "${faRad}.${otuMethod}.${clustId}.rep.names")
    outputs_renamed=("${out}.${clustId}.list" "${out}.${clustId}.fasta" "${out}.${clustId}.names")
    
elif [ $step == "classification" ]; then
    cmd=("classify.otu(taxonomy=${inputTax},list=${inputList},name=${inputNames},label=${clustId}) ; "
	  "make.shared(list=${inputList},group=${inputGroups})")
    outputs_mothur=("${listRad}.${clustId}.cons.taxonomy" "${listRad}.${clustId}.cons.tax.summary" "${listRad}.shared")
    outputs_renamed=("${out}.${clustId}.taxonomy" "${out}.${clustId}.summary" "${out}.${clustId}.shared")

fi;

res=$( IFS=$' '; echo "${cmd[*]}" )

mothur "#${res}"

set -o xtrace

n=$((${#outputs_mothur[@]}-1))
for i in `seq 0 $n`
do
    if [ -e ${outputs_mothur[$i]} ]; then
	mv ${outputs_mothur[$i]} ${outputs_renamed[$i]}
    else
	echo "${outputs_mothur[$i]} does not exist. Setting input file as output."
	cp ${inputs_mothur[$i]} ${outputs_renamed[$i]}
    fi
done
