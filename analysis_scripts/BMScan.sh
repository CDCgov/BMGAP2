#! /bin/bash
#$ -cwd

###################################################
# Defining variables
FASTA_DIR=${1} # Directory of each sample
NAME=${2}      # Sample name (M#)
# Characterization file with location of fasta and results
INPUT="$FASTA_DIR/characterization"
###################################################

###################################################
# Python script for BMScan
python analysis_scripts/SpeciesDB/bin/identify_species.py -d $INPUT -o $INPUT/BMScan -t 10

echo "BMScan complete for $NAME."
###################################################

