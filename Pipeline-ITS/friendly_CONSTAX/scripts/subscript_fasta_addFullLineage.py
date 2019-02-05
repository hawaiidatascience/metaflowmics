#!/usr/bin/python
import sys, string

def addFullLineage(filebase):
	print("\n\tAdding Full Lineage\n\n")
	
	f1 = open(filebase+"__RDP_taxonomy.txt", 'r').readlines()
	hash = {} #lineage map

	output_file = open(filebase+"__RDP_trained.fasta", 'w')

	for line in f1[1:]:
		cols = line.strip().split('\t')
		lineage = ['Root']
		for node in cols[1:]:
			if not node == '-':
				lineage.append(node)
		ID = cols[0]
		lineage = string.join(lineage, ';')
		hash[ID] = lineage
	f2 = open(filebase+"__RDP.fasta", 'r').readlines()
	for line in f2:
		if line[0] == '>':
			ID = line.strip().replace('>', '')
			try:
				lineage = hash[ID]
			except KeyError:
				print(ID, 'not in taxonomy file')
				sys.exit()
			output_file.write(line.strip() + '\t' + lineage+"\n")
		else:
			output_file.write(line.strip()+"\n")
	output_file.close()
