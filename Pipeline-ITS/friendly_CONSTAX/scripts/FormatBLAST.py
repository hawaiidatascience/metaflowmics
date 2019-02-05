blast_tax_file = open("otu_taxonomy_qiime/combo.ITS1_tax_assignments.txt", "r")
blast_tax = blast_tax_file.readlines()
blast_tax_file.close()

output = open("ITS1_BC_UN_500_otu_BLAST_reformatted.txt", "w")
output.write("OTU_ID\te-value\tKingdom\tPhylum\tClass\tOrder\tFamily\tGenus\tSpecies\n")
for i, line in enumerate(blast_tax):
	temp0 = line.strip().split("\t")
	otu = temp0[0]
	taxonomy = temp0[1].split(";")
	evalue = temp0[2]

	new_tax=[]
	if len(taxonomy)>1 and taxonomy[-1].endswith("_sp"):
			del taxonomy[-1]
	for j, level in enumerate(taxonomy):
		if taxonomy[0]=="No blast hit":
			new_tax=[""]
			evalue="NA"
			break
		elif "unidentified" in level:
			break
		else:
			new_tax.append(level[3:])

	output.write(otu+"\t"+evalue+"\t"+"\t".join(new_tax)+"\n")
