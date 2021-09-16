#!/usr/bin/env bash

n=$(grep -c "^>" Coleoptera.main.faa)

if [ "" == "" ] && [ $n -ge 500 ] ; then
    opt="-maxiters 1 -diags -sv -distance1 kbit20_3"
else
    opt=""
fi

echo $opt

# Run muscle and convert to 2-line fasta
cat Coleoptera.main.faa | muscle  $opt \
    | awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' \
    | tail -n+2 > Coleoptera.main.afa

echo $(muscle -version 2>&1) | cut -d" " -f2  > muscle.version.txt
