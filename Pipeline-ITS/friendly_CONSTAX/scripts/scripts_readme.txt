### Interested in modifying CONSTAX?
Please see descriptions of all custom python scripts implemented in CONSTAX. We've included these descriptions in anticipation that users would like to use CONSTAX with other reference databases. All of the scripts described below can be found within the *scripts* directory. 

#### *FormatRefDB.py*
The reference formatting step takes the UNITE general release fasta, and changes the information in the header line, and generates additional files (for the RDP Classifier) so all three classifier programs begin from a reference that contains an identical number of species. The series of scripts *FormatRefDB.py*, *subscript_fasta_addFullLineage.py*, and *subscript_lineage2taxonomyTrain.py* are all called by the *config* file.

#### ConsensusTaxonomy.py
This is the workhorse script. It takes the classifier-specific taxonomy outputs, combines them using the rules described in the methods, and generates a new consensus taxonomy table. The output is filtered at the cutoff specified in the *config* file. This script also creates a combined taxonomy table which shows each classifier's output side-by-side with the consensus output. 

#### FormatBLAST.py
We've also included a python script that formats BLAST output to match the format of CONSTAX output (when used with Qiime's *parallel_assign_taxonomy.py*). Although we do not recommend use of BLAST for taxonomy assignment of ribosomal sequences. This script is NOT automated through *constax.sh*. 
