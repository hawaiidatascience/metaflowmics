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

    --fwdFasta=*)
	fwdFasta="${arg#*=}"
	faRad=`basename $fwdFasta .fasta`;;

    --revFasta=*)
	revFasta="${arg#*=}" ;;
		       		       
    --optimize=*)
	optimize="${arg#*=}" ;;

    --criteria=*)
	criteria="${arg#*=}" ;;

    --refAln=*)
	refAln="${arg#*=}" ;;
    
    --refTax=*)
	refTax="${arg#*=}"
	taxRad=`basename $refTax .tax`;;

    --idThreshold=*)
	idThreshold="${arg#*=}" ;;

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
 
    echo ${rad}.names
    if [ ! -f ${rad}.names ]; then
	cmd=("screen.seqs(fasta=${rad}.fasta, optimize=${optimize}, criteria=${criteria})")
	outputs_mothur=("${rad}.good.fasta")	
	outputs_renamed=("${out}.${optimize_mod}.fasta")
    else
	cmd=("screen.seqs(fasta=${rad}.fasta, name=${rad}.names, group=${rad}.groups, optimize=${optimize}, criteria=${criteria})")
	outputs_mothur=("${rad}.good.names" "${rad}.good.fasta" "${rad}.good.groups")
	outputs_renamed=("${out}.${optimize_mod}.names" "${out}.${optimize_mod}.fasta" "${out}.${optimize_mod}.groups")
    fi
	 
elif [ $step == "summary" ]; then
    cmd=("summary.seqs(fasta=${rad}.fasta)")

elif [ $step == "dereplication" ]; then
    # Merge all fastas
    fasta_files=`ls *.fasta`
    fastas=`tr " " "\n" <<< "$fasta_files" | paste -sd - -`
    ids=`ls -l *.fasta| cut -d\  -f10 | awk -F'.' '{print $1}' | paste -sd - -`
    
    cat *.fasta > all.fasta
    cmd=("make.group(fasta=${fastas}, groups=${ids}) ; "
          "unique.seqs(fasta=all.fasta)")
   
    outputs_mothur=("all.unique.fasta" "all.names" "groups")
    outputs_renamed=("${out}.fasta" "${out}.names" "${out}.groups")
   
elif [ $step == "MSA" ]; then
    cmd=("align.seqs(fasta=${rad}.fasta,reference=${refAln}) ; "
	 "filter.seqs(fasta=${rad}.align)")
    outputs_mothur=("${rad}.filter.fasta")
    outputs_renamed=("${out}.fasta")

elif [ $step == "chimera" ]; then
    cmd=("chimera.vsearch(fasta=${rad}.fasta, name=${rad}.names, dereplicate=t) ; "
	 "remove.seqs(fasta=${rad}.fasta, name=${rad}.names, group=${rad}.groups, accnos=${rad}.denovo.vsearch.accnos)")
    outputs_mothur=("${rad}.pick.fasta" "${rad}.pick.names" "${rad}.pick.groups")
    outputs_renamed=("${out}.fasta" "${out}.names" "${out}.groups")

elif [ $step == "taxaFilter" ]; then    
    suffixTax=`echo $taxRad | cut -d. -f2`.wang.taxonomy
    cmd=("classify.seqs(fasta=${rad}.fasta, name=${rad}.names, group=${rad}.groups, template=${refAln}, taxonomy=${refTax})")
    # "remove.lineage(taxonomy=${rad}.${suffixTax}, name=${rad}.names, fasta=${rad}.fasta, group=${rad}.groups, taxon=-unknown)")
    outputs_mothur=("${rad}.${suffixTax}")
    outputs_renamed=("${out}.taxonomy")
    # outputs_mothur=("${rad}.pick.names" "${rad}.pick.fasta" "${rad}.pick.groups")
    # outputs_renamed=("${out}.names" "${out}.fasta" "${out}.groups")

elif [ $step == "subsampling" ]; then
    cmd=("sub.sample(persample=true,fasta=${rad}.fasta, name=${rad}.names, group=${rad}.groups, taxonomy=${rad}.taxonomy)")
    outputs_mothur=("${rad}.subsample.fasta" "${rad}.subsample.names" "${rad}.subsample.groups" "${rad}.subsample.taxonomy")
    outputs_renamed=("${out}.fasta" "${out}.names" "${out}.groups" "${out}.taxonomy")

elif [ $step == "clustering" ]; then
    if [ ${idThreshold} = 0 ]; then
	cmd=("cluster(name=${rad}.names, fasta=${rad}.fasta, method=unique, cutoff=0) ; "
	     "get.oturep(name=${rad}.names, fasta=${rad}.fasta, list=${rad}.unique.list, group=${rad}.groups, method=abundance) ; "
             "make.shared(list=${rad}.unique.list,group=${rad}.groups)")
	
	outputs_mothur=("${rad}.unique.list" "${rad}.unique.0.rep.fasta" "${rad}.unique.0.rep.names" "${rad}.unique.shared")
	
    else
	cmd=("cluster(name=${rad}.names, fasta=${rad}.fasta, method=dgc, cutoff=${idThreshold}) ; "
             "get.oturep(name=${rad}.names, fasta=${rad}.fasta, list=${rad}.dgc.list, group=${rad}.groups, method=abundance) ; "
	     "make.shared(list=${rad}.dgc.list,group=${rad}.groups)")
	
	outputs_mothur=("${rad}.dgc.list" "${rad}.dgc.${idThreshold}.rep.fasta" "${rad}.dgc.${idThreshold}.rep.names" "${rad}.dgc.shared")
    fi;
    
    outputs_renamed=("${out}.${idThreshold}.list" "${out}.${idThreshold}.fasta" "${out}.${idThreshold}.names" "${out}.${idThreshold}.shared")
    
elif [ $step == "classification" ]; then
    cmd=("classify.otu(taxonomy=all.subsampling.taxonomy, list=${rad}.list)")
    outputs_mothur=("all.clustering.${idThreshold}.${idThreshold}.cons.taxonomy" "all.clustering.${idThreshold}.${idThreshold}.cons.tax.summary")
    outputs_renamed=("all.classification.${idThreshold}.cons.taxonomy" "all.classification.${idThreshold}.cons.summary")
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
