#!/usr/bin/env bash

aln=$1
outdir=$(dirname $aln)

ncol=$(grep -v ">" $aln | awk '{L=length($0); if(L>max) max=L} END {print max}')

awk -v L=$ncol '{if($0 !~/^>/ && length($0)<L) $0=$0"-"}1' $aln | sed -r '/^>/s/([^ ]*) (.*)/\1\t100\t\2;/g' > $outdir/db.codons.align

grep '^>' $outdir/db.codons.align | cut -c2- | cut -f1,3 > $outdir/db.codons.tax
