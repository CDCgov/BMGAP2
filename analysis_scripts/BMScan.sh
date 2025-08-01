#! /bin/bash
#$ -cwd

###################################################
# Defining variables
FASTA_DIR=${1} # Directory of each sample
NAME=${2}      # Sample name (M#)
# Characterization file with location of fasta and results
INPUT="$FASTA_DIR/characterization"
NSLOTS=${NSLOTS:1} # Number of slots, default to 1 if not set
###################################################

###################################################
# Python script for BMScan
python analysis_scripts/SpeciesDB/bin/identify_species.py -d $INPUT -o $INPUT/BMScan -t $NSLOTS

echo "BMScan complete for $NAME."
###################################################

