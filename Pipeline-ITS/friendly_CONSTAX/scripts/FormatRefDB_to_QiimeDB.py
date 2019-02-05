################################################################################
#	Written by Natalie Vande Pol											   #
#	May 3, 2016																   #
#																			   #
#	*   command line: python generaldb_to_QIIMEdb.py [path/filename(s)]		   #
#	*   accepts unlimited reference database files and splits them to fasta	   #
#		and taxonomy files in the QIIME format								   #
#	*   output files written to original file location, the original filename  #
#		appended "__QIIME.fasta" and "__QIIME_taxonomy.txt"					   #
################################################################################

import sys, os

for arg in sys.argv[1:]:
	gen = open(arg, "r")
	all_lines = gen.readlines()
	gen.close()

	fasta = open(".".join(os.path.splitext(arg)[:-1])+"__QIIME.fasta","w")
	taxon = open(".".join(os.path.splitext(arg)[:-1])+"__QIIME_taxonomy.txt","w")
	for i, line in enumerate(all_lines):
		if line[0]==">":
			temp = line[1:].split("|")
			taxa = temp[4]
			name = [temp[2], temp[1], temp[3]]
			qiime_name = "_".join(name)

			fasta.write(">"+qiime_name+"\n")
			taxon.write(qiime_name+"\t"+taxa)
		else:
			fasta.write(line)
	fasta.close()
	taxon.close()