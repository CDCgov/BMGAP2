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

# Loop through each subdirectory in the parent directory and copying assembly fasta file to new folder
find "$RESULT_DIR" "$RESULT_DIR/characterization" -maxdepth 1 -type d | \
  grep -v shareFiles | \
  xargs -n 1 bash -c '
    # If there are json files, then cp them
    files=$(\ls $0/*json 2>/dev/null || true)
    for file in $files; do
        if [ ! -f "$file" ]; then
            continue;
        fi
        cp -vf $file "'$DEST_DIR'"
    done
  '

cp -v "$RESULT_DIR"/characterization/*_cleaned.fasta "$DEST_DIR/cleaned.fasta"
cp -v "$RESULT_DIR"/characterization/BMScan/species_analysis_*.csv "$DEST_DIR/species_analysis.csv"

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

# delete the staging directory and just be verbose on the rmdir
#rm -rf "$DEST_DIR/"*
#rmdir -v "$DEST_DIR"
