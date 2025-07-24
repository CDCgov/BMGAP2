#!/bin/bash -l

### This script runs the AMR genotyping script on the PMGA output folder created by the pmga script

#$ -pe smp 8-16
#$ -M ujz6@cdc.gov
#$ -m abe
#$ -l h_vmem=10G

module load miniconda3
conda activate base_plus2
#module load Python/3.4

##First move each JSON into a directory
shopt -s nullglob
for json in pmga/json/*final_results.json; 
do 
    mydir=${json%_final_results.json}/;
    mkdir $mydir;
    mv $json $mydir;
done;

##Make an output directory
mkdir AMR_test4

##Then run AMR on each of those directories
for json_dir in pmga/json/*/;
do
    result_dir='AMR_test4/'${json_dir##pmga/json/}
    python3 ../runAST.py -i $json_dir -s Nm -o $result_dir;
done;

