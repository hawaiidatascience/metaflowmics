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

    --minAlnLen=*)
	minAlnLen="${arg#*=}" ;;

    --refAln=*)
	refAln="${arg#*=}" ;;
    
    --refTax=*)
	refTax="${arg#*=}"
	taxRad=`basename $refTax .tax`;;

    --taxaToFilter=*)
	taxaToFilter="${arg#*=}" ;;
    
    --minAbundance=*)
	minAbundance="${arg#*=}" ;;

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
    
    cmd=("align.seqs(fasta=${fasta}.fasta, reference=${refAln})"
		 "filter.seqs(fasta=${fasta}.align)"
		 "screen.seqs(fasta=${fasta}.filter.fasta, count=${count}.count_table, minlength=${minAlnLen})"
		 "screen.seqs(fasta=current, count=current, optimize=${optimize}, criteria=${criteria})"
		 "summary.seqs(fasta=current)")
    
    outputs_mothur=("${fasta}.filter.good.good.fasta" 
					"${count}.good.good.count_table")
    
elif [ $step == "chimera" ]; then
    method="vsearch"

    cmd=("chimera.${method}(fasta=${fasta}.fasta, count=${count}.count_table, dereplicate=t)"
		 "remove.seqs(fasta=${fasta}.fasta, count=${count}.count_table, accnos=${fasta}.denovo.${method}.accnos, dups=f)")
    
    outputs_mothur=("${fasta}.pick.fasta"
					"${count}.pick.count_table")
    
elif [ $step == "preClassification" ]; then
    suffixTax=`echo $taxRad | cut -d. -f2`.wang
    
    cmd=("classify.seqs(fasta=${fasta}.fasta, count=${count}.count_table, template=${refAln}, taxonomy=${refTax})")
    
    outputs_mothur=("${fasta}.${suffixTax}.taxonomy")
    
elif [ $step == "clustering" ]; then
    method=`[ ${idThreshold} -eq 100 ] && echo 'unique' || echo 'dgc'`

    cmd=("cluster(count=${count}.count_table, fasta=${fasta}.fasta, method=${method}, cutoff=${mothurThresh})"
		 "make.shared(list=${fasta}.${method}.list, count=${count}.count_table)")
	
    outputs_mothur=("${fasta}.${method}.shared"
					"${fasta}.${method}.list")
    
elif [ $step == "consensusClassification" ]; then
    suffixTax=`echo $taxRad | cut -d. -f2`.wang
    
    cmd=("classify.otu(taxonomy=${tax}.taxonomy, count=${count}.count_table, list=${list}.list, probs=f)")
    
    outputs_mothur=("${list}.${mothurThresh}.cons.taxonomy"
					"${list}.${mothurThresh}.cons.tax.summary")

elif [ $step == "otuRepr" ]; then
    cmd=("get.oturep(count=${count}.count_table, fasta=${fasta}.fasta, list=${list}.list, method=abundance, rename=T)")

    outputs_mothur=("${list}.${mothurThresh}.rep.fasta")

elif [ $step == 'multipletonsFilter' ]; then
    cmd=("filter.shared(shared=${shared}.shared, mintotal=${minAbundance}, makerare=F)")

    outputs_mothur=("${shared}.${mothurThresh}.filter.shared")

elif [ $step == "subsampling" ]; then
    cmd=("sub.sample(persample=true, fasta=${fasta}.fasta, constaxonomy=${tax}.taxonomy, shared=${shared}.shared, size=${subsamplingNb})")
    
    outputs_mothur=("${fasta}.subsample.fasta"
					"${tax}.subsample.taxonomy"
					"${shared}.${mothurThresh}.subsample.shared")

elif [ $step == "taxaFilter" ]; then
    # 2 cases depending on whether we do the taxa filtering on the taxonomy or constaxonomy file
    # For the latter case, we need to also process the .shared file 
    if [ -f "${shared}.list" ]; then
        outputs_mothur=("${tax}.pick.cons.taxonomy"
						"${list}.${mothurThresh}.pick.list"
						"${shared}.${mothurThresh}.pick.shared")
	
	cmd=("remove.lineage(constaxonomy=${tax}.taxonomy, shared=${shared}.shared, list=${list}.list, taxon='${taxaToFilter}')")
    else
        outputs_mothur=("${fasta}.pick.fasta"
						"${shared}.pick.shared"
						"${tax}.pick.cons.taxonomy")
	
	cmd=("remove.lineage(constaxonomy=${tax}.taxonomy, shared=${shared}.shared, taxon='${taxaToFilter}')"
	     "list.seqs(taxonomy=${tax}.pick.cons.taxonomy)"
	     "get.seqs(fasta=${fasta}.fasta,accnos=current)")
    fi

elif [ $step == "postprocessing" ]; then
	cmd=("count.seqs(shared=${shared}.shared)"
		 "get.relabund(shared=${shared}.shared)")
	outputs_mothur=("${shared}.relabund"
				    "${shared}.${idThreshold}.count_table")
	# "create.database(shared=${shared}.shared,label=${idThreshold},repfasta=${fasta}.fasta,constaxonomy=${tax}.taxonomy)")

elif [ $step == "unifrac" ]; then
    cmd=("clearcut(fasta=${fasta}.fasta, DNA=T)"
		 "unifrac.unweighted(tree=current,count=${count}.count_table,distance=lt)"
		 "unifrac.weighted(tree=current,count=${count}.count_table,distance=lt)")
	
	outputs_mothur=("${fasta}.tre1.unweighted.phylip.dist"
					"${fasta}.tre1.weighted.phylip.dist"
					"${fasta}.uwsummary"
					"${fasta}.tre1.wsummary")
	
	outputs_custom=("${out}.unweighted.dist"
					"${out}.weighted.dist"
					"${out}.unweighted.summary"
					"${out}.weighted.summary")
	
elif [ $step == "alphaDiversity" ]; then
	cmd=("summary.single(shared=${shared}.shared,calc=nseqs-sobs-chao-shannon-shannoneven)")
	outputs_mothur=("${shared}.groups.summary")

elif [ $step == "betaDiversity" ]; then
    cmd=("summary.shared(shared=${shared}.shared,calc=braycurtis-thetayc-sharedsobs-sharedchao)")
	outputs_mothur=("${shared}.summary")
fi;

# Process mothur cmd

# Join the instruction array
res=$( IFS=$';'; echo "${cmd[*]}" )
echo $res

# Execute
set -o xtrace
[ -z $MOTHUR ] && mothur "#${res}" || $MOTHUR/mothur "#${res}"
set +o xtrace

# Rename output files
# for output_mothur in ${outputs_mothur[@]}
for i in ${!outputs_mothur[@]}
do
	output_mothur=${outputs_mothur[i]}
    extension=${output_mothur##*.}
	output_renamed=`[ -z $outputs_custom ] && echo "${out}.${extension}" || echo ${outputs_custom[i]}`
    
    # if the output file exists
    if [ -e $output_mothur ]; then
	echo "Success: Renaming $output_mothur to $output_renamed"
	mv $output_mothur $output_renamed
    # Special case for screen.seqs (sometime mothur doesnt produce an output file).
	# In this case, just copy the input into the output
    elif [ $step = "MSA" ] || [ $step = "taxaFilter" ] || [ $step = "subsampling" ]; then
	filename_to_copy=`ls -t *$extension | head -1`
	echo "WARNING: $output_renamed does not exist. Copying latest file with extension $extension ($filename_to_copy)."
	cp $filename_to_copy $output_renamed
    # Otherwise, raise an error
    else
	echo "ERROR: $output_mothur does not exist. Aborting."
	exit 1
    fi
done

