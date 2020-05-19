#!/usr/bin/env bash

getRad() {
    files=($(ls -t -1 *$1 2>/dev/null))
    if [ ! -z files ]; then
	echo $(basename ${files[0]} $1)
    fi
}

set -o xtrace

fasta=`getRad .fasta`
shared=`getRad .shared`
count=`getRad .count_table`
tax=`getRad .taxonomy`
list=`getRad .list`
relabund=`getRad .relabund`

set +o xtrace

prefix=all

for arg in "$@"
do
    case $arg in
    --step=*)
	step="${arg#*=}" ;;
    
    --prefix=*)
	prefix="${arg#*=}" ;;

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

    --otuId=*)
	otuId="${arg#*=}" ;;
    
    --subsamplingNb=*)
	subsamplingNb="${arg#*=}" ;;

    --rename=*)
	rename="${arg#*=}" ;;

    --uf_mode=*)
	uf_mode="${arg#*=}" ;;	
    *)
	echo "$arg: Unknown option";;    
    esac
done

# General prefix for output names
out="${prefix}_${step}"

if [ ! -z $otuId ]; then
    if [ $otuId -eq 100 ]; then
	mothurThresh=0
    else
	mothurThresh=0.0"$((100-${otuId}))"	
    fi
    out="${out}_${otuId}"
fi

if [ $step == "MSA" ]; then
    
    cmd=("align.seqs(fasta=${fasta}.fasta, reference=${refAln})"
		 "filter.seqs(fasta=${fasta}.align, vertical=T)"
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
    method=`[ ${otuId} -eq 100 ] && echo 'unique' || echo 'dgc'`
	
    cmd=("cluster(count=${count}.count_table, fasta=${fasta}.fasta, method=${method}, cutoff=${mothurThresh})")
	
    outputs_mothur=("${fasta}.${method}.list")
    
elif [ $step == "multipletonsFilter" ]; then
	cmd=("remove.rare(list=${list}.list, count=${count}.count_table, nseqs=${minAbundance})"
		 "make.table(count=${count}.pick.count_table, compress=f)"
		 "list.seqs(count=${count}.pick.full.count_table)"
		 "get.seqs(taxonomy=${tax}.taxonomy, accnos=${count}.pick.full.accnos)"
		 "get.seqs(fasta=${fasta}.fasta, accnos=${count}.pick.full.accnos)")

	outputs_mothur=("${list}.${mothurThresh}.pick.list"
				    "${count}.pick.full.count_table"
				    "${fasta}.pick.fasta"
				    "${tax}.pick.taxonomy")

elif [ $step == "subsampling" ]; then
    cmd=("sub.sample(count=${count}.count_table, taxonomy=${tax}.taxonomy, list=${list}.list, size=${subsamplingNb}, persample=true)"
		 "list.seqs(list=${list}.${mothurThresh}.subsample.list)"
		 "get.seqs(fasta=${fasta}.fasta, accnos=${list}.${mothurThresh}.subsample.accnos)")
    
    outputs_mothur=("${fasta}.pick.fasta"
					"${count}.subsample.count_table"
					"${tax}.subsample.taxonomy"
					"${list}.${mothurThresh}.subsample.list")

elif [ $step == "taxaFilter" ]; then
	cmd=("remove.lineage(taxonomy=${tax}.taxonomy, fasta=${fasta}.fasta, count=${count}.count_table, list=${list}.list, taxon='${taxaToFilter}')")

	outputs_mothur=("${tax}.pick.taxonomy"
				    "${fasta}.pick.fasta"
				    "${count}.pick.count_table"
				    "${list}.${mothurThresh}.pick.list")

elif [ $step == "consensusClassification" ]; then
    suffixTax=`echo $taxRad | cut -d. -f2`.wang
    
    cmd=("classify.otu(taxonomy=${tax}.taxonomy, count=${count}.count_table, list=${list}.list, probs=f)")
    
    outputs_mothur=("${list}.${mothurThresh}.cons.taxonomy"
					"${list}.${mothurThresh}.cons.tax.summary")

elif [ $step == "otuRepr" ]; then
    cmd=("get.oturep(count=${count}.count_table, fasta=${fasta}.fasta, list=${list}.list, method=abundance, rename=${rename})")

    outputs_mothur=("${list}.${mothurThresh}.rep.fasta"
				    "${list}.${mothurThresh}.rep.count_table")

elif [ $step == "abundanceTable" ]; then
	cmd=("make.shared(count=${count}.count_table, list=${list}.list)")
	outputs_mothur=("${list}.shared")
	
elif [ $step == "relabund" ]; then
	cmd=("get.relabund(shared=${shared}.shared)")
	outputs_mothur=("${shared}.relabund")

elif [ $step == "database" ]; then
	cmd=("create.database(relabund=${relabund}.relabund,label=${mothurThresh},repfasta=${fasta}.fasta,count=${count}.count_table,constaxonomy=${tax}.taxonomy)")
	outputs_mothur=("${shared}.${mothurThresh}.database")

elif [ $step == "biom" ]; then
	cmd=("make.biom(shared=${shared}.shared,constaxonomy=${tax}.taxonomy)")
	outputs_mothur=("${shared}.${mothurThresh}.biom")

elif [ $step == "alphaDiversity" ]; then
	cmd=("summary.single(shared=${shared}.shared,calc=nseqs-sobs-chao-shannon-shannoneven)")
	outputs_mothur=("${shared}.groups.summary")

elif [ $step == "betaDiversity" ]; then
    cmd=("summary.shared(shared=${shared}.shared,calc=braycurtis-thetayc-sharedsobs-sharedchao, distance=lt)")
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
# Exclamation mark gets the index of the array element instead of the values
for i in ${!outputs_mothur[@]}
do
	out_file=${outputs_mothur[i]}
    extension=${out_file##*.}
	output_renamed=$([ -z ${outputs_custom} ] && echo "${out}.${extension}" || echo ${outputs_custom[i]})
    
    if [ -f ${out_file} ]; then
		echo "Step successful: Renaming ${out_file} to ${output_renamed}"
		mv ${out_file} ${output_renamed}

    # Special case for screen.seqs (sometime mothur doesnt produce an output file).
	# In this case, just copy the input into the output
    elif [ $step = "MSA" ] || [ $step = "taxaFilter" ] || [ $step = "subsampling" ]; then
		filename_to_copy=$(ls -t *$extension | head -1)
		echo "WARNING: ${output_renamed} does not exist. Copying latest file with extension ${extension} (${filename_to_copy})."
		cp ${filename_to_copy} ${output_renamed}

	else
		echo "ERROR: ${out_file} does not exist. Aborting."
		exit 1
    fi

	if [ "${extension}" == "shared" ]; then
		echo "Calculating summary for ${output_renamed}"
		mothur "#summary.single(shared=${output_renamed}, calc=nseqs-sobs)"
	fi
done

if [ -f "${out}.list" ] && [ -f "${out}.count_table" ] && [ ! -f "${out}.shared" ]; then
	echo "Calculating summary for ${output_renamed} and ${out}.count_table"
	mothur "#make.shared(list=${out}.list, count=${out}.count_table); summary.single(shared=current, calc=nseqs-sobs)"
fi
