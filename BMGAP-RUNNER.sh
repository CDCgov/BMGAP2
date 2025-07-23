#Written by ujz6 on 02/25/2025
#! /bin/bash
#$ -cwd

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
###########################################

echo "Beginning analysis of isolates"
echo ""
#Loop through all R1 files in the Directory
for r1_file in "$FASTQ_DIR"/*R1*.fastq.gz; do
        #Check if file exists
	if [[ ! -f $r1_file ]]; then
		echo "No R1 files found."
		exit 1
	fi
	
	#Derive the corresponding R2 filename from R1 filename
	r2_file="${r1_file/R1/R2}"  # Replace 'R1' with 'R2'
       	
	#Check if the corresponding R2 file exists
	if [[ -f $r2_file ]]; then
		#Extract base name without path
		base_name=$(basename "$r1_file" .fastq.gz)
		base_name=${base_name%%_*} #Gets only M#
		OUTPUT_DIR="$ANALYSIS_DIRECTORY/$base_name"
		
		# Define error and output file paths for isolate
		error_file="$OUTPUT_DIR/${base_name}_error.e"
		output_file="$OUTPUT_DIR/${base_name}_output.o"
		
		# Align PE-FASTQ files and remove human DNA
		echo "Beginning alignment of $base_name."
		align_id=$(qsub -V -e "$error_file" -o "$output_file" $ANALYSIS_SCRIPTS/Alignment.sh "$r1_file" "$r2_file" "$OUTPUT_DIR" "$base_name" | awk '{print $3}')
		
		# Assembly cleanup, check QC 
		cleanup_id=$(qsub -V -hold_jid "$align_id" -e "$error_file" -o "$output_file" $ANALYSIS_SCRIPTS/cleanupSingle.sh "$OUTPUT_DIR" "$base_name" | awk '{print $3}')
		echo "Submitted AssemblyCleanup, waiting for Alignment job $align_id to complete before beginning AssemblyCleanup."
		
		# Run BMScan on each new fasta file
		bmscan_id=$(qsub -V -hold_jid "$cleanup_id" -e "$error_file" -o "$output_file" $ANALYSIS_SCRIPTS/BMScan.sh "$OUTPUT_DIR" "$base_name" | awk '{print $3}')
		
		# Run PMGA on each new fasta file 
		pmga_id=$(qsub -V -hold_jid "$bmscan_id" -e "$error_file" -o "$output_file" $ANALYSIS_SCRIPTS/PMGA.sh "$OUTPUT_DIR" "$base_name" | awk '{print $3}')
		sleep 200
		
		# Run LocusExtractor on each new fasta file
		le_id=$(qsub -V -hold_jid "$bmscan_id" -e "$error_file" -o "$output_file" $ANALYSIS_SCRIPTS/LocusExtractor.sh "$OUTPUT_DIR" "$base_name" | awk '{print $3}')
		echo "Submitted LocusExtractor, BMScan and PMGA, waiting for AssemblyCeanup job $cleanup_id to complete before running."
		
		# Run AMR with species code for each sample
		amr_id=$(qsub -V -hold_jid "$pmga_id" -e "$error_file" -o "$output_file" $ANALYSIS_SCRIPTS/AMR.sh "$OUTPUT_DIR" "$base_name" | awk '{print $3}')
		echo "Submitted AMR, waiting for PMGA job $pmga_id to complete before running."	
		sleep 200
	else
		echo "Warning: corresponding R2 file not found for $r1_file"
	fi
done 
