#! /bin/bash
#$ -cwd

###################################################
# Defining variables
OUT_DIR=${1}     # Output directory of each sample
NAME=${2}        # Sample name (M#)
###################################################

###################################################
# Create new folders for characterization with fasta
mkdir $OUT_DIR/characterization
mkdir $OUT_DIR/AssemblyCleanup

# Python script for AssemblyCleanup, makes new contig.fasta file
python analysis_scripts/AssemblyCleanup/cleanupSingle.py -c 0 -p 0.1 -b $NAME $OUT_DIR/${NAME}_SPAdes/contigs.fasta $OUT_DIR/AssemblyCleanup/${NAME}_DIS.fasta 

# Move the new, cleaned fasta file to new directory for analysis
cp "$OUT_DIR"/AssemblyCleanup/${NAME}_DIS.fasta $OUT_DIR/characterization/${NAME}_cleaned.fasta

echo "AssemblyCleanup complete and newly cleaned file ready for characterization"
###################################################

