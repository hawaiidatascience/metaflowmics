# Written by Natalie Vande Pol
# November 15, 2017
#
# Command line: python FormatRefDB.py
# Script will locate, import, and validate the database file specified
# Requires subscript_lineage2taxonomyTrain.py and subscript_fasta_addFullLineage.py
# to be in the same directory
# Creates the following files in the input file location:
#  - a combined taxonomy and fasta database file compatible with the UTAX and
#    SINTAX classifiers
#  - fasta and taxonomy database files compatible with the RDP classifer
#  - training files for the RDP classifier
# Output files are named and located based on the input file



# -*- coding: utf-8 -*-
import sys, unicodedata
from os.path import basename,splitext 
import subscript_lineage2taxonomyTrain, subscript_fasta_addFullLineage

filename = sys.argv[1]
outputdir = sys.argv[2]

try:
        open(filename,"r")
except IOError:
        print("ERROR: file could not be opened.")
        valid = 0
        sys.quit()

input_file = open(filename,"r")
line = input_file.readline()
temp0 = line.split("|")
if line[0]!=">" or len(temp0)!= 5 or "k__" not in temp0[-1]:
        print("Header Format Incompatible. Please Reformat As Below:")
        print(">Peziza_sp|JN102365|SH189857.07FU|reps|k__Fungi;p__Ascomycota;c__Pezizomycetes;o__Pezizales;f__Pezizaceae;g__Peziza;s__Peziza_sp")
        sys.quit()
input_file.close()


print("\n____________________________________________________________________\nReformatting database\n")
filename_base = "{}/{}".format(outputdir, splitext(basename(filename))[0])
# filename_base = ".".join(os.path.splitext(filename)[:-1])
#UTAX output file
fastatax = open(filename_base+"__UTAX.fasta","w")

#RDP output files
fasta = open(filename_base+"__RDP.fasta","w")
taxon = open(filename_base+"__RDP_taxonomy.txt","w")
taxon.write("\t".join(["Seq_ID","Kingdom","Phylum","Class","Order","Family","Genus","Species\n"]))

num = 0
with open(filename) as database:
        for line in database:
                if line[0]==">":
                        #correct umlauts or special letters
                        unico_line = unicode(line, '1252')
                        ascii_line = unicodedata.normalize('NFKD', unico_line).encode('ASCII', 'ignore')
                        temp = ascii_line[1:].split("|")

                        #UTAX file
                        utax_name = temp[1]+"|"+temp[2]
                        utax_taxa = temp[4][1:].strip().replace("__",":").replace(";",",")
                        temp_utax = utax_taxa.strip().split(",")
                        if temp_utax[0].endswith("Fungi"):
                                temp_utax2 = [x for x in temp_utax if "unidentified" not in x]
                                temp_utax3 = [y for y in temp_utax2 if "Incertae_sedis" not in y]
                                temp_utax4 = [z for z in temp_utax3 if "unknown" not in z]
                                new_utax_taxa = ",".join(temp_utax4)
                                fastatax.write(">"+utax_name+";tax=d"+new_utax_taxa+";\n")

                        #RDP files
                        name = str(temp[1])
                        temp2 = temp[4].strip().split("__")
                        to_genus = [ item[:-2] for item in temp2[1:-1] ]

                        if "Incertae_sedis" in to_genus:
                                indices = [i for i,x in enumerate(to_genus) if x == "Incertae_sedis"]
                                for j in indices:
                                        if "Incertae_sedis" not in to_genus[j-1]:
                                                to_genus[j] = str(to_genus[j-1])+"_Incertae_sedis"
                                        else:
                                                to_genus[j] = str(to_genus[j-1])
                        if "unidentified" in to_genus:
                                indices = [i for i,x in enumerate(to_genus) if x == "unidentified"]
                                for j in indices:
                                        to_genus[j] = "-"

                        if to_genus[0] != "-":
                                species = str(temp2[-1])
                                if "Incertae" in species:
                                        species = "unidentified_sp"
                                elif to_genus[-1] not in species:
                                        temp=species.split("_")
                                        species = temp[0]+"_unidentified_"+temp[1]
                                if species.endswith("sp"):
                                        species+= "_"+str(num)
                                        num += 1

                                taxonomy = name+"\t"+"\t".join(to_genus)+"\t"+species+"\n"
                                fasta.write(">"+name+"\n")
                                taxon.write(taxonomy)
                                seq = next(database)
                                fastatax.write(seq)
                                fasta.write(seq)

        fasta.close()
        taxon.close()
        fastatax.close()

        subscript_lineage2taxonomyTrain.lin2tax(filename_base)
        subscript_fasta_addFullLineage.addFullLineage(filename_base)

print("Database formatting complete\n____________________________________________________________________\n\n")
