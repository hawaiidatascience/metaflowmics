#!/usr/bin/env bash

meta=$1
matching=$2
fastqs=($3 $4)

if [ ${#fastqs[@]} -eq 2 ]; then 

	[ ${3##*.} == 'gz' ] && z='z' || z=''

	if [ ${matching} == "auto" ]; then
		symbol=$(${z}cat ${fwd} | head -c 2)

		paste --delimiters=' ' <(\
				${z}grep --no-group-separator "^${symbol}" ${fwd} -A1 | grep -v ${symbol}) <( \
				${z}grep --no-group-separator "^${symbol}" ${rev} -A1 | grep -v ${symbol}) \
			   | grep -v N \
			   | sort | uniq -c | sort -rnk1 | head -n20 \
			   | sed 's/^[[:space:]]*//g' | sed 's/ /,/g' | cut -d, -f2,3 \
			   > pairs_freqs.txt

		awk -F, '{OFS=","} {print $2,$1}' pairs_freqs.txt > pairs_freqs_rev.txt

		same1=$(comm -12 <(sort <(cut -d, -f2,3 ${meta})) <(sort pairs_freqs.txt) | wc -l)
		same2=$(comm -12 <(sort <(cut -d, -f2,3 ${meta})) <(sort pairs_freqs_rev.txt) | wc -l)

		reversed=`[ $same2 -ge $same1 ] && echo true || echo false`
	else
		reversed=`[ ${params.matching} == "reversed" ] && echo true || echo false`
	fi

	if [ $reversed==true ]; then
		awk -F"," '{OFS=","}{print $1,$3,$2}' ${meta}
	else
		cat ${meta}
	fi

else
	awk -F"," '{OFS=","}{print $1,$2,NaN}' ${meta}
fi
