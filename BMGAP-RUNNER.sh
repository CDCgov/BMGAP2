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

echo "Current path: $PATH2"
echo "FASTQ directory: $FASTQ_DIR"
echo "Analysis directory: $ANALYSIS_DIRECTORY"
echo "Analysis scripts directory: $ANALYSIS_SCRIPTS"

mkdir -p $ANALYSIS_DIRECTORY
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
		echo "No R1 files found." >&2
		exit 1
	fi
	
	#Derive the corresponding R2 filename from R1 filename
	r2_file="${r1_file/R1/R2}"  # Replace 'R1' with 'R2'
       	
	#Check if the corresponding R2 file exists
	if [[ -f $r2_file ]]; then
		echo -e "$r1_file\t$r2_file"
	fi
done > $CTRL_FILE

# Use -tc 2 to limit the number of concurrent jobs to 2 especially for pubmlst interactions
qsub -N bmgap -o $ANALYSIS_DIRECTORY/log -j y -V -cwd -tc 2 -t 1-$(cat $CTRL_FILE | wc -l) \
  -v "ANALYSIS_SCRIPTS=$ANALYSIS_SCRIPTS" -v "ANALYSIS_DIRECTORY=$ANALYSIS_DIRECTORY" -v "CTRL_FILE=$CTRL_FILE" <<- "END_OF_SCRIPT"

	r1_file=$(sed -n ${SGE_TASK_ID}p $CTRL_FILE | awk '{print $1}')
	r2_file=$(sed -n ${SGE_TASK_ID}p $CTRL_FILE | awk '{print $2}')

	#Extract base name without path
	base_name=$(basename "$r1_file" .fastq.gz)
	base_name=${base_name%%_*} #Gets only M#
	OUTPUT_DIR="$ANALYSIS_DIRECTORY/$base_name"
	mkdir -p "$OUTPUT_DIR"
	
	# Define error and output file paths for isolate
	error_file="$OUTPUT_DIR/${base_name}_error.e"
	output_file="$OUTPUT_DIR/${base_name}_output.o"
	
	# Align PE-FASTQ files and remove human DNA
	echo "Alignment of $base_name." | tee --append "$output_file" "$error_file"
	bash $ANALYSIS_SCRIPTS/Alignment.sh "$r1_file" "$r2_file" "$OUTPUT_DIR" "$base_name" \
		>> "$output_file" 2>> "$error_file"
	
	# Assembly cleanup, check QC 
	echo "AssemblyCleanup for $base_name." | tee --append "$output_file" "$error_file"
	$ANALYSIS_SCRIPTS/cleanupSingle.sh "$OUTPUT_DIR" "$base_name" \
		>> "$output_file" 2>> "$error_file"
	
	# Run BMScan on each new fasta file
	echo "BMSCAN for $base_name." | tee --append "$output_file" "$error_file"
	bash $ANALYSIS_SCRIPTS/BMScan.sh "$OUTPUT_DIR" "$base_name" \
		>> "$output_file" 2>> "$error_file"
	
	# Run PMGA on each new fasta file 
	echo "PMGA for $base_name." | tee --append "$output_file" "$error_file"
	bash $ANALYSIS_SCRIPTS/PMGA.sh "$OUTPUT_DIR" "$base_name" \
		>> "$output_file" 2>> "$error_file"
	
	# Run LocusExtractor on each new fasta file
	echo "LocusExtractor for $base_name." | tee --append "$output_file" "$error_file"
	bash $ANALYSIS_SCRIPTS/LocusExtractor.sh "$OUTPUT_DIR" "$base_name" \
		>> "$output_file" 2>> "$error_file"
	
	# Run AMR with species code for each sample
	echo "AMR for $base_name." | tee --append "$output_file" "$error_file"
	bash $ANALYSIS_SCRIPTS/AMR.sh "$OUTPUT_DIR" "$base_name" \
		>> "$output_file" 2>> "$error_file"
	
END_OF_SCRIPT