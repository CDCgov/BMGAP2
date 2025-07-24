#! /bin/bash
#$ -cwd

###################################################
# Defining variables
###################################################
SAMPLE_DIR=${1} # Directory of each sample
NAME=${2}       # Sample name (M#)
# Characterization file with location of fasta and results
INPUT="$SAMPLE_DIR/characterization/PMGA/json"
OUTPUT="$SAMPLE_DIR/characterization/AMR_$NAME"
mkdir $OUTPUT
ORIGINAL_DIR=$(pwd)
###################################################

###################################################
##SPECIES EXTRACTION##
# Script to extract species from BMScan result for input to AMR 
# Move to BMScan folder to extract species
cd "$SAMPLE_DIR/characterization/BMScan" || { echo "BMScan directory not found"; exit 1; }

bmscan_output=$(find . -maxdepth 1 -name 'species_analysis*.csv' -print -quit)

# Use awk to get the result from the second row (first row of results)
species_name=$(awk -F ',' 'NR==2 {print $2}' "$bmscan_output")

# Based on species, code for runAST.py
case $species_name in
    "Haemophilus influenzae")
        amr_code="Hi"
        ;;
    "Neisseria meningitidis")
        amr_code="Nm"
        ;;
    *)
        echo "Species is other or not found; AMR analysis will not be run."
        exit 0 # Exit successfully since it is not an error
        ;;
esac
###################################################

###################################################
##AMR TESTING##
# Run AMR for each isolate, using the amr_code for specific species
# Move back to original directory
cd "$ORIGINAL_DIR"
python3 analysis_scripts/amr_variants_component/runAST.py -i $INPUT -s $amr_code -x -o $OUTPUT

echo "AMR completed for $NAME"


