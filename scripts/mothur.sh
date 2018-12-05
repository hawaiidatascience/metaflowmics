#!/usr/bin/env bash

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
		       
    --optimize=*)
	optimize="${arg#*=}" ;;

    --criteria=*)
	criteria="${arg#*=}" ;;

    --refAln=*)
	refAln="${arg#*=}" ;;
    
    --refTax=*)
	refTax="${arg#*=}"
	taxRad=`basename $refTax .tax`;;

    *)
	echo "$arg: Unknown option";;    
    esac
done

out="${pairId}.${step}"

if [ $step == "merging" ]; then
    cmd=("make.contigs(ffasta=${fwdFasta}, rfasta=${revFasta})")
    outputs_mothur=("${faRad}.trim.contigs.fasta")
    outputs_renamed=("${out}.fasta")
    
elif [ $step == "filter" ]; then
    optimize_mod=`echo ${optimize} | sed s/-/_/g`
    if [ -z ${inputNames} ]; then
	cmd=("screen.seqs(fasta=${inputFasta}, optimize=${optimize}, criteria=${criteria})")
	outputs_mothur=("${faRad}.good.fasta")
	outputs_renamed=("${out}.${optimize_mod}.fasta")
    else
	inputs_mothur=("${inputNames}" "${inputFasta}")
	cmd=("screen.seqs(fasta=${inputFasta}, name=${inputNames}, optimize=${optimize}, criteria=${criteria})")
	outputs_mothur=("${namesRad}.good.names" "${faRad}.good.fasta")
	outputs_renamed=("${out}.${optimize_mod}.names" "${out}.${optimize_mod}.fasta")
    fi
	 
elif [ $step == "summary" ]; then
    cmd=("summary.seqs(fasta=${inputFasta})")
    # outputs_mothur=("${faRad}.summary")
    # outputs_renamed=("${faRad}.summary")
    
elif [ $step == "dereplication" ]; then
    cmd=("unique.seqs(fasta=${inputFasta})")
    outputs_mothur=("${faRad}.unique.fasta" "${faRad}.names")
    outputs_renamed=("${out}.fasta" "${out}.names")
   
elif [ $step == "MSA" ]; then
    cmd=("align.seqs(fasta=${inputFasta},reference=${refAln}) ; "
	 "filter.seqs(fasta=${faRad}.align)")
    outputs_mothur=("${faRad}.filter.fasta")
    outputs_renamed=("${out}.fasta")

elif [ $step == "chimera" ]; then
    cmd=("chimera.vsearch(fasta=${inputFasta}, name=${inputNames},dereplicate=t) ; "
	 "remove.seqs(fasta=${inputFasta}, name=${inputNames},accnos=${faRad}.denovo.vsearch.accnos)")
    outputs_mothur=("${faRad}.pick.fasta" "${namesRad}.pick.names")
    outputs_renamed=("${out}.fasta" "${out}.names")

elif [ $step == "taxaFilter" ]; then
    
    suffixTax=`echo $taxRad | cut -d. -f2`.wang.taxonomy
    
    cmd=("classify.seqs(fasta=${inputFasta}, name=${inputNames}, template=${refAln}, taxonomy=${refTax}) ; "
	 "remove.lineage(taxonomy=${faRad}.${suffixTax}, name=${inputNames}, fasta=${inputFasta},taxon=-unknown)")
    outputs_mothur=("${namesRad}.pick.names" "${faRad}.pick.fasta")
    outputs_renamed=("${out}.names" "${out}.fasta")
    
#elif [ $step == "Classification" ]; then
    
#	 "rename.file(input=${faRad},new=${out})"
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
