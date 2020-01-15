#!/usr/bin/env bash

meta=$1
matching=$2
rc=$3
fwd=$4
rev=$5
fastqs=($4 $5)

sed -e 's/\r//g' ${meta} > metadata.csv

if [ "$rc" == true ]; then
	cat metadata.csv | rev | cut -d, -f1 | tr "ATGC" "TACG" > rev_idx_rc.csv
	cat metadata.csv | rev | cut -d, -f2- | rev > samp_idx1.csv
	paste -d, samp_idx1.csv rev_idx_rc.csv | sed -e 's/\r//g' > metadata.csv
fi

if [ ${#fastqs[@]} -eq 2 ]; then 

	[ ${4##*.} == 'gz' ] && z='z' || z=''

	if [ ${matching} == "auto" ]; then
		symbol=$(${z}cat ${fwd} | head -c 2)

		paste --delimiters=' ' <(\
				${z}grep --no-group-separator "^${symbol}" ${fwd} -A1 | grep -v ${symbol}) <( \
				${z}grep --no-group-separator "^${symbol}" ${rev} -A1 | grep -v ${symbol}) \
			   | grep -v N \
			   | sort | uniq -c | sort -rnk1 | head -n20 \
			   | sed 's/^[[:space:]]*//g' | sed 's/ /,/g' \
			   | cut -d, -f2,3 > pairs_freqs.txt
		
		awk -F, '{OFS=","} {print $2,$1}' pairs_freqs.txt > pairs_freqs_rev.txt

		same1=$(comm -12 <(sort <(cut -d, -f2,3 metadata.csv)) <(sort pairs_freqs.txt) | wc -l)
		same2=$(comm -12 <(sort <(cut -d, -f2,3 metadata.csv)) <(sort pairs_freqs_rev.txt) | wc -l)

		reversed=`[ $same2 -ge $same1 ] && echo true || echo false`
	else
		reversed=`[ ${params.matching} == "reversed" ] && echo true || echo false`
	fi

	if [ $reversed==true ]; then
		awk -F"," '{OFS=","}{print $1,$3,$2}' metadata.csv
	else
		cat metadata.csv
	fi

else
	awk -F"," '{OFS=","}{print $1,$2,NaN}' metadata.csv
fi
