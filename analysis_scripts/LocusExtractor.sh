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
# get species
BMSCAN_file=$(find "$FASTA_DIR/characterization/BMScan" -maxdepth 1 -name 'species_analysis*.csv' -print -quit)

# Extract the value from the second column
species_name=$(awk -F ',' 'NR==2 {print $2}' "$BMSCAN_file")

# Check if species is either H.influenzae or N.meningitidis
if [ -z "$species_name" ]; then
    echo "Error: species is empty, unable to continue testing."
    exit 1
fi

if [[ "$species_name" == "Haemophilus influenzae" || "$species_name" == "Neisseria meningitidis" ]]; then
    echo "Species is $species_name, continue with testing" 

    # Python script to extract MLST and Vaccine antigens for samples
    python3 analysis_scripts/locusextractor/LocusExtractor_RHEL8.py --no_update -p $NAME $INPUT

    # Move result folder to sample directory
    mv LE*${NAME}* $INPUT

    echo "Locus Extractor completed for $NAME"

else
    echo "Species is not H.influenzae or N.meningitidis, species is \"$species_name\", done." 
    exit 0 # Not an error

fi
###################################################

