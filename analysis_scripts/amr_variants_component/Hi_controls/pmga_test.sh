#!/bin/bash -l

### This script runs PMGA (blast_pubmlst) on the control sequences. Make sure you are using the same installation of PMGA that you use in your pipeline.

#$ -cwd
#$ -q dbd.q,all.q 
#$ -pe smp 8-16

ASSEMBLY_DIR=sequences/    
module load Python/3.7
module load ncbi-blast+/LATEST
module load Mash/1.1
python3 /scicomp/groups/OID/NCIRD/DBD/MVPDB/ML/tools/PMGA/blast_pubmlst.py -p -o pmga -x test_files.xlsx -t $NSLOTS #Note: this will run BMSCAN ... then default to Nm.

