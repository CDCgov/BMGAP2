#! /bin/bash
#$ -cwd

###################################################
# Defining variables
###################################################
FASTA_DIR=${1} # Directory of each sample
NAME=${2}      # Sample name (M#)

# Characterization file with location of fasta and results
INPUT="$FASTA_DIR/characterization"
###################################################
##SPECIES IDENTIFICATION##

BMSCAN_file=$(find "$FASTA_DIR/characterization/BMScan" -maxdepth 1 -name 'species_analysis*.csv' -print -quit)

# Extract the value from the second column
species_name=$(awk -F ',' 'NR==2 {print $2}' "$BMSCAN_file")

# Check if species is either H.influenzae or N.meningitidis
if [ -z "$species_name" ]; then
    echo "Error: species is empty, unable to determine species."
    exit 1
fi

if [[ "$species_name" == "Haemophilus influenzae" || "$species_name" == "Neisseria meningitidis" ]]; then
    echo "Species is $species_name, continue with testing" 

    # PMGA run on cleaned fasta file 
    python3 analysis_scripts/PMGA/blast_pubmlst.py -sc -a -o $INPUT/PMGA -sg -p -d $INPUT -t 16 -fa

else
    echo "Species is not H.influenzae or N.meningitidis, species is \"$species_name\", done." 
    exit 0 # Not an error

fi
###################################################
echo "PMGA completed for $NAME"

