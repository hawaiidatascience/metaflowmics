#!/bin/bash

echo "###############*** CONSTAX v1. modified ***######################"
echo "#See https://github.com/natalie-vandepol/compare_taxonomy for the original repository#\n"

echo "________________Parameters________________"
set -o xtrace

RDPTools_path=$HOME/.local/bin/RDPTools

SCRIPT_DIR="./" # $PWD/friendly_CONSTAX
ROOT_DIR=`dirname $1`

otuPath=$PWD/$1
constaxOutput=$ROOT_DIR/CONSTAX_outputs
taxOutput=$constaxOutput/taxonomy_assignments
trainOutput=$constaxOutput/training_files
dbOutput=$constaxOutput/DB
globalOutput=$constaxOutput/outputs

ref_database="uniteDB_01-12-2017_nodup.fasta"
base=${ref_database%.fasta}

threads=10
conf_threshold="0.8"

set +o xtrace

mkdir -p $taxOutput
mkdir -p $trainOutput
mkdir -p $globalOutput
mkdir -p $dbOutput

# Execute the python script, passing as the first argument the value of the variable ref_database declared in the config file
python2 $SCRIPT_DIR/scripts/FormatRefDB.py $SCRIPT_DIR/DB/$ref_database $dbOutput

echo "__________________________________________________________________________"
echo "Training UTAX Classifier"

usearch8 -utax_train $dbOutput/${base}__UTAX.fasta \
	 -report $trainOutput/utax_db_report.txt \
	 -taxconfsout $trainOutput/utax.tc \
	 -utax_splitlevels NVpcofgs \
	 -utax_trainlevels kpcofgs \
	 -log $trainOutput/utax_train.log \
	 -report $trainOutput/utax_report.txt

usearch8 -makeudb_utax $dbOutput/${base}__UTAX.fasta \
	 -taxconfsin $trainOutput/utax.tc \
	 -output $trainOutput/utax.db \
	 -log $trainOutput/make_udb.log \
	 -report $trainOutput/utax_report.txt

usearch8 -utax $otuPath \
	 -db $trainOutput/utax.db \
	 -strand both \
	 -utaxout $taxOutput/otu_taxonomy.utax \
	 -utax_cutoff $conf_threshold \
	 -threads $threads

echo "__________________________________________________________________________"
echo "Training RDP Classifier"

java -Xmx10g -jar $RDPTools_path/classifier.jar train \
     -o $trainOutput \
     -s $dbOutput/${base}__RDP_trained.fasta \
     -t $dbOutput/${base}__RDP_taxonomy_trained.txt

cp $RDPTools_path/classifier/samplefiles/rRNAClassifier.properties $trainOutput/. 

java -Xmx10g -jar $RDPTools_path/classifier.jar classify $otuPath \
               --conf $conf_threshold \
	       --format allrank \
	       -o $taxOutput/otu_taxonomy.rdp \
	       --train_propfile $trainOutput/rRNAClassifier.properties 

echo "__________________________________________________________________________"
echo "Training SINTAX Classifier"

usearch10 -makeudb_sintax $dbOutput/${base}__UTAX.fasta \
	  -output $trainOutput/sintax.db

usearch10 -sintax $otuPath \
	  -db $trainOutput/sintax.db \
	  -tabbedout $taxOutput/otu_taxonomy.sintax \
	  -strand both \
	  -sintax_cutoff $conf_threshold

echo "__________________________________________________________________________"
echo "Creating Consensus Taxonomy"

python2 $SCRIPT_DIR/scripts/CombineTaxonomy.py $conf_threshold $constaxOutput

# plot graphs in R
Rscript $SCRIPT_DIR/R/ComparisonBars.R $constaxOutput








