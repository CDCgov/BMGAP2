#!/bin/bash
#$ -cwd 

set -e

# Check if the correct number of arguments is provided
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <result_directory> <LAB_id>"
    exit 1
fi
#Analysis output directoy 
RESULT_DIR="${1}"
echo "Analysis Directory is $RESULT_DIR"
Lab_ID="${2}"

# Check if the parent directory exists
if [ ! -d "$RESULT_DIR" ]; then
    echo "Result directory does not exist."
    exit 1
fi

#Directory to store results to share
DEST_DIR="$RESULT_DIR/shareFiles"
echo "destination dir is $DEST_DIR"
mkdir -pv "$DEST_DIR"

# Find the species_analysis_datestamp_.json file that has the latest date.
# I believe that the numbers in the filename are either timestamps or job numbers but that still gives us an ordering.
species_analysis_json=$(find "$RESULT_DIR" -type f -name "species_analysis_*.json" | sort | tail -n 1)
# same sort of thing for csv
species_analysis_csv=$(find "$RESULT_DIR" -type f -name "species_analysis_*.csv" | sort | tail -n 1)
# e.g., 
# $RESULT_DIR/characterization/BMScan/species_analysis_1754783735.3227649.json
# $RESULT_DIR/characterization/BMScan/species_analysis_1754809463.2148802.json
# $RESULT_DIR/characterization/BMScan/species_analysis_1754867544.7031925.json <-- latest

# find latest locus extractor directory. Timestamp is in the LE directory
LE_dir=$(\ls -d "$RESULT_DIR"/characterization/LE_* | sort | tail -n 1)
# find latest serogroup_predictions_ .tab file in the PMGA serogroup directory
PMGA_serogroup_tab=$(find "$RESULT_DIR"/characterization/PMGA/serogroup -type f -name "serogroup_predictions_*.tab" | sort | tail -n 1)

cp -vf "$RESULT_DIR"/characterization/*_cleaned.fasta "$DEST_DIR/assembly_cleaned.fasta"
cp -vf "$species_analysis_json" "$DEST_DIR/bmscan_species_analysis.json"
cp -vf "$species_analysis_csv" "$DEST_DIR/bmscan_species_analysis.csv"
cp -vf "$RESULT_DIR"/characterization/AMR_*/*_amr_data.json "$DEST_DIR/amr_data.json"
cp -vf "$LE_dir"/molecular_data_*.json "$DEST_DIR/le_molecular_data.json"
cp -vf "$LE_dir"/Results_text/molecular_data_*.csv "$DEST_DIR/le_molecular_data.csv"
cp -vf "$RESULT_DIR"/characterization/PMGA/scheme_counts.json "$DEST_DIR/pmga_scheme_counts.json"
cp -vf "$RESULT_DIR"/characterization/PMGA/json/*_cleaned_final_results.json "$DEST_DIR/pmga_cleaned_final_results.json"
cp -vf "$RESULT_DIR"/characterization/PMGA/serogroup/serogroup_results.json "$DEST_DIR/pmga_serogroup_results.json"
cp -vf "$PMGA_serogroup_tab" "$DEST_DIR/pmga_serogroup_predictions.tab"

# Make sure the permissions are sane
find "$DEST_DIR" -type f -exec chmod 644 {} \;

# Create the tar archive
TAR_FILE="$Lab_ID"-"share.tgz"
echo "zip file name $TAR_FILE"
tar -cvzf "$TAR_FILE" "$DEST_DIR"

# Check if the tar command was successful
if [ $? -eq 0 ]; then
    echo "Successfully created tar archive: $TAR_FILE"
else
    echo "Error: Failed to create tar archive."
    exit 1
fi
