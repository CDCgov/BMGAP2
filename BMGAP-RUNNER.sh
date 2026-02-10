#! /bin/bash
#$ -cwd

# Written by ujz6 on 02/25/2025

###########################################
# Function to display the usage message
usage() {
    echo "Usage: $0 <FASTQ_DIR> <ANALYSIS_DIRECTORY>"
    echo
    echo "Required Arguments:"
    echo "[FASTQ_DIR] 		Directory with FASTQ locations."
    echo "[ANALYSIS_DIRECTORY] 	Directory for results."
    echo
}

# Script versioning (BMGAP Version, found in GitHub)
SCRIPT_VERSION="2"
SCRIPT_SUBVERSION="2.0"
export BMGAP_VERSION="$SCRIPT_VERSION.$SCRIPT_SUBVERSION"

# Check if no arguments or insufficient arguments are provided
if [[ $# -lt 2 ]]; then
    usage
    exit 1
fi
###########################################

###########################################
#Defining directories containing FASTQ files and analysis
FASTQ_DIR=${1}
ANALYSIS_DIRECTORY=${2}

PATH2="$(pwd)"
ANALYSIS_SCRIPTS="$PATH2/analysis_scripts"

# Check if FASTQ_DIR exists and is a directory
if [[ ! -d "$FASTQ_DIR" ]]; then
    echo "Error: The directory '$FASTQ_DIR' does not exist or is not a valid directory."
    exit 1
fi

# Check if FASTQ_DIR contains .fastq.gz files
if ! ls "$FASTQ_DIR"/*.fastq.gz &>/dev/null; then
    echo "Error: No fastq.gz files found in the directory '$FASTQ_DIR'."
    exit 1
fi

mkdir -p $ANALYSIS_DIRECTORY
###########################################

# Prepare run report
RUN_REPORT="$ANALYSIS_DIRECTORY/run_report.txt"
START_TIME=$(date +"%Y-%m-%d %H:%M:%S")

{
echo "Current path: $PATH2"
echo "Run report generated: $START_TIME"
echo "BMGAP Version: $SCRIPT_VERSION.$SCRIPT_SUBVERSION"
echo "FASTQ input directory: $FASTQ_DIR"
echo "Analysis output directory: $ANALYSIS_DIRECTORY"
echo "Analysis scripts directory: $ANALYSIS_SCRIPTS"
echo
} > "$RUN_REPORT"

LOG_DIR="$ANALYSIS_DIRECTORY/log"
mkdir $LOG_DIR
CTRL_FILE="$ANALYSIS_DIRECTORY/fastq.fofn"
touch $CTRL_FILE
###########################################

echo "Beginning analysis of isolates"
echo ""
#Loop through all R1 files in the Directory
for r1_file in "$FASTQ_DIR"/*R1*.fastq.gz; do
        #Check if file exists
	if [[ ! -f $r1_file ]]; then
		echo "No R1 files found (expected: $r1_file) ."
		        echo "No R1 files found in $FASTQ_DIR, please use a different directory" >> "$RUN_REPORT"
		exit 1
	fi
	
	#Derive the corresponding R2 filename from R1 filename
	r2_file="${r1_file/R1/R2}"  # Replace 'R1' with 'R2'
       	
	#Check if the corresponding R2 file exists
	if [[ -f $r2_file ]]; then
		echo -e "$r1_file\t$r2_file"
	else
		echo "Warning: Corresponding R2 file not found for $r1_file. Skipping this pair." >&2
	fi
done 

END_TIME=$(date +"%Y-%m-%d %H:%M:%S")

{
echo
echo "Run start time: $START_TIME"
echo "Run end time: $END_TIME"
echo "Full results directory: $ANALYSIS_DIRECTORY"
echo "Note: Individual sample output directories are under the analysis directory."
} >> "$RUN_REPORT"
