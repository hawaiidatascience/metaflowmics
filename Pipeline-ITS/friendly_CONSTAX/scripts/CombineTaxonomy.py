# Written by Natalie Vande Pol
# November 15, 2017
#
# Command line: python CombineTaxonomy.py
# Script will locate the taxonomy file from each classifier and validate file
# Creates a taxonomy file in standardized format for each input file and
# generates a consensus taxonomy file based on the three input files
# Reformatted taxonomy files are named based on the input file
# All output files are placed in dedicated outputs folder

import sys, itertools

################################################################################
def reformat_RDP(rdp_file,outputdir="."):
        input = open(rdp_file, "U")
        all_lines = input.readlines()
        input.close()

        output_file = outputdir + "/outputs/otu_taxonomy_rdp_final.txt"
        output = open(output_file,"w")
        output.write("\t".join(["OTU_ID","OTU_Score","Kingdom","K_score","Phylum","P_score","Class","C_score"]))
        output.write("\t".join(["Order","O_score","Family","F_score","Genus","G_score","Species","S_score\n"]))

        for i, line in enumerate(all_lines):
                # capture confidence level at genus before altering line
                temp = line.strip().split("\t")
                confi = temp[7:-2][::3]
                confi.append(temp[-1])
                taxon = temp[5:][::3]

                # remove any taxonomic levels after first "unidentified"
                j=0
                new_taxon = []
                while j<len(taxon):
                        if taxon[j].endswith("Incertae_sedis"):
                                taxon[j] = "Incertae_sedis"
                        if "unidentified" in taxon[j] or float(confi[j])<confidence:
                                del confi[j:]
                                break
                        else: new_taxon.append(taxon[j].capitalize())
                        j+=1

                # remove "_sp" species classificaitons
                if len(new_taxon)>0 and " sp" in new_taxon[-1]:
                        del new_taxon[-1]
                        del confi[-1]
                # remove terminal Incertae_sedis
                while len(new_taxon)>0 and "Incertae_sedis" in new_taxon[-1]:
                        del new_taxon[-1]
                        del confi[-1]


                if confi == []:
                        score = "NA"
                else:
                        score = confi[-1]

                iters = [iter(new_taxon), iter(confi)]
                tax_confi = list(str(it.next()) for it in itertools.cycle(iters))

                output.write(temp[0]+"\t"+score+"\t"+"\t".join(tax_confi)+"\n")

        output.close()
        return output_file

################################################################################
def reformat_UTAX(utax_file,outputdir="./"):
        input = open(utax_file, "U")
        all_lines = input.readlines()
        input.close()

        output_file = outputdir + "/outputs/otu_taxonomy_utax_final.txt"
        output = open(output_file,"w")
        output.write("\t".join(["OTU_ID","OTU_Score","Kingdom","K_score","Phylum","P_score","Class","C_score"]))
        output.write("\t".join(["Order","O_score","Family","F_score","Genus","G_score","Species","S_score\n"]))
        
        for i, line in enumerate(all_lines):
                #remove unwanted third column and convert "(" and ")" to "*"
                temp = line.replace("(", "*").replace(")","*").split()
                line = temp[0]+"\t"+temp[1]

                temp0 = line.split("*")
                confid = temp0[1:][::2]
                line2 = "".join(temp0[0:-2][::2])

                temp1 = line2.split(",")
                j=1
                new_line = [temp1[0]]
                while j<len(temp1):
                        if  "unidentified" in temp1[j] or float(confid[j-1]) < confidence:
                                del confid[j-1:]
                                break
                        else:
                                new_line.append(temp1[j])
                        j+=1

                line2 = ",".join(new_line)

                temp2 = line2.split(",")
                # OTU_1328        d:Fungi, p:Ascomycota, c:Archaeorhizomycetes, o:Archaeorhizomycetales, f:Archaeorhizomycetaceae, g:Archaeorhizomyces, s:Archaeorhizomyces_sp
                temp3 = temp2[0].split()
                # OTU_1328, d:Fungi
                temp4 = []
                for item in temp2[1:]:
                        temp4.append(item[2:].capitalize())
                # Ascomycota, Archaeorhizomycetes, Archaeorhizomycetales, Archaeorhizomycetaceae, Archaeorhizomyces, Archaeorhizomyces_sp

                # remove "_sp" species classificaitons
                if len(temp4)>0 and temp4[-1].endswith("_sp"):
                        del temp4[-1]
                        del confid[-1]
                        
                # remove terminal Incertae_sedis
                while len(temp4)>0 and temp4[-1]=="Incertae_sedis":
                        del temp4[-1]
                        del confid[-1]

                confid = [str(x) for x in confid]
                confid.insert(0, "NA")
                score = confid[-1]
                if len(temp3) == 0:
                        continue
                final_line = temp3[0]+"\t"+score+"\t"+temp3[1][2:].capitalize()+"\t"

                iters = [iter(confid), iter(temp4)]
                tax_confi = list(str(it.next()) for it in itertools.cycle(iters))

                output.write(final_line+"\t".join(tax_confi)+"\n")

        output.close()
        return output_file

################################################################################
def reformat_SINTAX(sintax_file,outputdir="./"):
        input = open (sintax_file, "U")
        all_lines = input.readlines()
        input.close()

        output_file = outputdir + "/outputs/otu_taxonomy_sintax_final.txt"
        output = open(output_file, "w")
        output.write("\t".join(["OTU_ID","OTU_Score","Kingdom","K_score","Phylum","P_score","Class","C_score"]))
        output.write("\t".join(["Order","O_score","Family","F_score","Genus","G_score","Species","S_score\n"]))
        
        for i, line in enumerate(all_lines):
                # remove unwanted third column and convert "(" and ")" to "*"
                temp = line.replace("(", "*").replace(")","*").split("\t")
                del temp[2:]
                # temp =>  'OTU_999'        'd:Fungi*1.0000*,p:Ascomycota*1.0000*,c:Leotiomycetes*1.0000*,o:Helotiales*1.0000*,s:Helotiales_sp*1.0000*'

                temp0 = temp[1].split("*")
                # temp0 =>  'd:Fungi,' '1.0000' ',p:Ascomycota' '1.0000' ',c:Leotiomycetes' '1.0000' ',o:Helotiales' '1.0000' ',s:Helotiales_sp' '1.0000'
                confid = temp0[1:][::2]
                # confid =>  '1.0000' '1.0000' '1.0000' '1.0000' '1.0000'
                temp_line = "".join(temp0[0:-2][::2])
                # temp_line =>  "d:Fungi,p:Ascomycota,c:Leotiomycetes,o:Helotiales,s:Helotiales_sp"

                # fix missing taxonomic levels
                temp1 = temp_line.split(",")
                # temp1 =>  'd:Fungi' 'p:Ascomycota' 'c:Leotiomycetes' 'o:Helotiales' 's:Helotiales_sp'
                levels = ["d:", "p:", "c:", "o:", "f:", "g:", "s:"]
                if len(temp1)<len(levels):
                        if "g:" in temp1[-2]:
                                for k, level in enumerate(temp1):
                                        if levels[k] not in temp1[k]:
                                                temp1.insert(k, levels[k]+"Incertae_sedis")
                                                confid.insert(k, 9)
                        else:
                                for k, level in enumerate(temp1):
                                        if levels[k] not in temp1[k]:
                                                temp1.insert(k, levels[k]+"unidentified")
                                                confid.insert(k, 0)

                j=0
                temp2 = []
                while j<len(temp1):
                        if  "unidentified" in temp1[j]:
                                del confid[j:]
                                break
                        elif confid[j]==9 and float(confid[j+1])<confidence:
                                del confid[j:]
                                break
                        elif float(confid[j])<confidence:
                                del confid[j:]
                                break
                        else:
                                temp2.append(temp1[j].capitalize())
                        j+=1

                # remove "_sp" species classificaitons
                if len(temp2)>0 and temp2[-1].endswith("_sp"):
                        del temp2[-1]
                        del confid[-1]
                # remove terminal Incertae_sedis
                while len(temp2)>0 and "Incertae_sedis" in temp2[-1]:
                        del temp2[-1]
                        del confid[-1]

                confid = [str(x) if x!=9 else "NA" for x in confid]

                new_taxonomy = []
                for item in temp2:
                        new_taxonomy.append(item[2:].capitalize())

                if confid == []:
                        score = "NA"
                else:
                        score = confid[-1]

                iters = [iter(new_taxonomy), iter(confid)]
                tax_confi = list(str(it.next()) for it in itertools.cycle(iters))

                output.write(temp[0]+"\t"+str(score)+"\t"+"\t".join(tax_confi)+"\n")


        output.close()
        return output_file

################################################################################
def build_dict(filename):
        file = open(filename, "r")
        all_lines = file.readlines()
        file.close()

        dict = {}
        for i, line in enumerate(all_lines[1:]):
                temp = line.replace(" ", "_").strip().split()
                dict[temp[0]] = []
                if len(temp)>2:
                        dict[temp[0]]=temp[2:]
                if len(dict[temp[0]])<14:
                        while len(dict[temp[0]])<14:
                                dict[temp[0]].append("")
                # strip numbers from species identifications
                species = filter(lambda c: not c.isdigit(), dict[temp[0]][-2])
                dict[temp[0]][-2] = species.replace("_"," ")
        return dict

################################################################################
def vote(rdp, sin, uta):
        winner = ""
        taxa = [rdp[0], sin[0], uta[0]]
        scores = [rdp[1], sin[1], uta[1]]
        # tally = ["0","0","0"]
        duplicates_notempty = [i for i, x in enumerate(taxa) if x!= "" and taxa.count(x) > 1]
        unique_2empty = [i for i, x in enumerate(taxa) if x!="" and taxa.count("") > 1]
        unique = [i for i, x in enumerate(taxa) if x!="" and taxa.count("") == 1]
        for j in range(0,3):
                if taxa[j]!="":
                        if j in duplicates_notempty:
                                winner = taxa[j]
                                break
                        elif j in unique_2empty:
                                winner = taxa[j]
                                break
                        elif j in unique:
                                scores = [float(x) if x!="NA" and x!="" else 0 for x in scores]
                                winner = taxa[scores.index(max([scores[x] for x in unique]))]
                else:
                        winner = taxa[j]
        return winner

################################################################################
def count_classifications(filenames,outputdir="./"):
        # rdp, utax, sintax, consensus
        file_num = 0
        unique_dict = {0:{}, 1:{}, 2:{}, 3:{}, 4:{}, 5:{}, 6:{}}
        for i, file in enumerate(filenames):
                input = open(file, "r")
                all_lines = input.readlines()
                input.close()
                count_y = [0,0,0,0,0,0,0]
                count_n = [0,0,0,0,0,0,0]
                output1 = open(outputdir + "/outputs/otu_taxonomy_CountClassified.txt", "w")
                if i<3:        #first 3 files have scores
                        start= 2
                        freq = 2
                else:         #consensus file does not have scores
                        start= 1
                        freq = 1
                for j, line in enumerate(all_lines[1:]):
                        temp = line.strip().split("\t")
                        taxonomy = temp[start::freq]
                        if len(taxonomy)==7:
                                # strip numbers from species identifications
                                species = filter(lambda c: not c.isdigit(), taxonomy[-1])
                                taxonomy[-1] = species.replace("_"," ")

                        for k in range(0,7):
                                if k<len(taxonomy):
                                        count_y[k]+=1
                                        if taxonomy[k] not in unique_dict[k]:
                                                unique_dict[k][taxonomy[k]] = [0,0,0,0]
                                                unique_dict[k][taxonomy[k]][file_num]+= 1
                                        else:
                                                unique_dict[k][taxonomy[k]][file_num]+= 1
                                else:
                                        taxonomy.append("Unidentified")
                                        count_n[k]+=1
                                        if taxonomy[k] not in unique_dict[k]:
                                                unique_dict[k][taxonomy[k]] = [0,0,0,0]
                                                unique_dict[k][taxonomy[k]][file_num]+= 1
                                        else:
                                                unique_dict[k][taxonomy[k]][file_num]+= 1


                for l, level in enumerate(count_y):
                        count_y[l] = str(level)
                for m, level in enumerate(count_n):
                        count_n[m] = str(level)
                output1.write("\t".join(["Kingdom","Phylum","Class","Order","Family","Genus","Species\n"]))
                output1.write("Classified\t"+"\t".join(count_y)+"\n")
                output1.write("Unclassified\t"+"\t".join(count_n))
                output1.close()
                file_num+=1

        output2 = open(outputdir + "/outputs/Classification_Summary.txt", "w")
        output2.write("\t".join(["Classification","RDP","UTAX","SINTAX","Consensus\n"]))
        for l in range(0, 7):
                key_list = unique_dict[l].keys()
                key_list.sort()
                for key in key_list:
                        output2.write(key+"\t"+"\t".join(str(x) for x in unique_dict[l][key])+"\n")
        output2.close()

###############################################################################
# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
###############################################################################
#Example RDP input:
#OTU_2836        _        Root        rootrank        1.0        Fungi        Kingdom        0.98        Zygomycota        Phylum        0.05        Zygomycota_Incertae_sedis        Class        0.05        Mucorales        Order        0.04        Syncephalastraceae        Family        0.01        Fennellomyces        Genus        0.01        Fennellomyces linderi        Species        0.01
#Criteria: temp = line.split("\t"), temp[3]=="rootrank" len(temp)> 21

#Example UTAX input
#OTU_2545        d:Fungi,p:Basidiomycota(0.9941),c:Agaricomycetes(0.8548),o:Agaricales(0.8523),f:Mycenaceae(0.5664),g:Mycena(0.2776),s:Mycena_epipterygia(0.0882)        p:Basidiomycota,c:Agaricomycetes,o:Agaricales        -
#Criteria: "),c:" in line, temp = line.split("\t"), temp[1].startswith("d:")

#Example SINTAX input
#OTU_5443        d:Fungi(1.0000),p:Ascomycota(0.9700),c:Pezizomycetes(0.8000),o:Pezizales(0.7900),f:Sarcosomataceae(0.7700),g:Pseudoplectania(0.3700),s:Pseudoplectania_nigrella(0.3700)        +        d:Fungi,p:Ascomycota,c:Pezizomycetes
#Criteria: "),p:" in line, temp = line.split("\t"), temp[1].startswith("d:")

confidence = float(sys.argv[1])
outputdir = sys.argv[2]

three_classifiers = ["rdp", "utax", "sintax"]

for classifier in three_classifiers:
        file_name = outputdir+"/taxonomy_assignments/otu_taxonomy."+classifier
        try:
                open(file_name,"r")
        except IOError:
                print("ERROR: "+classifier.upper()+" file could not be opened.")
                quit()
        input_file = open(file_name,"r")
        line = input_file.readline()
        temp0 = line.split("\t")
        if classifier == "rdp":
                if len(temp0)<10 or temp0[3]!="rootrank":
                        print("Input file not in RDP format. Please Reformat As Below:")
                        print("OTU_###        _        Root        rootrank        1.0        Fungi        Kingdom        0.98        Zygomycota        Phylum        0.05        Zygomycota_Incertae_sedis        Class        0.05        Mucorales        Order        0.04        Syncephalastraceae        Family        0.01        Fennellomyces        Genus        0.01        Fennellomyces linderi        Species        0.01")
                        quit()
        elif classifier == "utax":
                if ",s:" not in temp0[1] or not temp0[1].startswith("d:"):
                        print("Input file not in UTAX format. Please Reformat As Below:")
                        print("OTU_###        d:Fungi,p:Ascomycota(0.9700),c:Pezizomycetes(0.8000),o:Pezizales(0.7900),f:Sarcosomataceae(0.7700),g:Pseudoplectania(0.3700),s:Pseudoplectania_nigrella(0.3700)        +        d:Fungi,p:Ascomycota,c:Pezizomycetes")
                        quit()
        else:
                if "),s:" not in temp0[1] or not temp0[1].startswith("d:"):
                        print("Input file not in SINTAX format. Please Reformat As Below:")
                        print("OTU_###        d:Fungi(1.0000),p:Ascomycota(0.9700),c:Pezizomycetes(0.8000),o:Pezizales(0.7900),f:Sarcosomataceae(0.7700),g:Pseudoplectania(0.3700),s:Pseudoplectania_nigrella(0.3700)        +        d:Fungi,p:Ascomycota,c:Pezizomycetes")
                        quit()

        input_file.close()

        if classifier == "rdp":
                print("\n____________________________________________________________________\nReformatting "+classifier.upper()+" file\n")
                rdp_file = reformat_RDP(file_name,outputdir=outputdir)
                rdp_dict = build_dict(rdp_file)
        elif classifier == "utax":
                print("\nReformatting "+classifier.upper()+" file\n")
                uta_file = reformat_UTAX(file_name,outputdir=outputdir)
                uta_dict = build_dict(uta_file)
        else:
                print("\nReformatting "+classifier.upper()+" file\n")
                sin_file = reformat_SINTAX(file_name,outputdir=outputdir)
                sin_dict = build_dict(sin_file)
        print("\tDone\n")

print("\nGenerating consensus taxonomy & combined taxonomy table\n")
consensus_file = outputdir + "/outputs/consensus_taxonomy.txt"
consensus = open(consensus_file, "w")
consensus.write("\t".join(["OTU_ID","Kingdom","Phylum","Class","Order","Family","Genus","Species\n"]))

combined = open(outputdir+"/outputs/combined_taxonomy.txt", "w")
combined.write("\t".join(["OTU_ID","Kingdom_RDP","Kingdom_SINTAX","Kingdom_UTAX","Kingdom_Consensus","Phylum_RDP","Phylum_SINTAX","Phylum_UTAX",
                          "Phylum_Consensus","Class_RDP","Class_SINTAX","Class_UTAX","Class_Consensus","Order_RDP","Order_SINTAX","Order_UTAX",
                          "Order_Consensus","Family_RDP","Family_SINTAX","Family_UTAX","Family_Consensus","Genus_RDP","Genus_SINTAX","Genus_UTAX",
                          "Genus_Consensus","Species_RDP","Species_SINTAX","Species_UTAX","Species_Consensus\n"]))

for otu in rdp_dict.keys():
        consensus.write(otu+"\t")
        combined.write(otu)
        levels = []

        for m in range(0,14,2):
                rdp_val = ["",""]
                sin_val = ["",""]
                uta_val = ["",""]
                
                if otu in rdp_dict:
                        rdp_val = rdp_dict[otu][m:m+2]
                if otu in sin_dict:
                        sin_val = sin_dict[otu][m:m+2]
                if otu in uta_dict:
                        uta_val = uta_dict[otu][m:m+2]
                        
                level = vote(rdp_val,sin_val,uta_val)
                # level = vote(rdp_dict[otu][m:m+2], sin_dict[otu][m:m+2], uta_dict[otu][m:m+2])
                combined.write("\t"+rdp_val[0]+"\t"+sin_val[0]+"\t"+uta_val[0]+"\t")                
                # combined.write("\t"+rdp_dict[otu][m]+"\t"+sin_dict[otu][m]+"\t"+uta_dict[otu][m]+"\t")
                if level != "":
                        levels.append(level)
                combined.write(level)
        consensus.write("\t".join(levels)+"\n")
        combined.write("\n")
print("\tDone\n")

consensus.close()
combined.close()

print("\nGenerating classification counts & summary table\n")
count_classifications([rdp_file, uta_file, sin_file, consensus_file], outputdir=outputdir)
print("\tDone\n\n")
print("____________________________________________________________________\n")
