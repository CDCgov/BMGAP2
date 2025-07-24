#!/bin/bash
#$ -cwd 

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
mkdir  "$DEST_DIR"

# Loop through each subdirectory in the parent directory and copying assembly fasta file to new folder
for SUBDIR in "$RESULT_DIR"/*; do
    # Check if it's a directory
    if [ -d "$SUBDIR" ]; then
        base="${SUBDIR##*/}"
        if [ "$base" == "shareFiles" ]; then
           echo "skipping $base"
        else
            SRC_DIR="$RESULT_DIR/$base/characterization"
#            echo "destination sub dir $SRC_DIR"
            cp -r "$SRC_DIR"/*fasta "$DEST_DIR/"
        fi
    fi
done
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
