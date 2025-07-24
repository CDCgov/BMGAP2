#!/bin/bash -l

### This script runs the AMR genotyping script on the PMGA output folder created by the pmga script

#$ -cwd
#$ -q dbd.q,all.q 
#$ -pe smp 1

  
module load Python/3.4

##First move each JSON into a directory
shopt -s nullglob
for json in pmga/json/*final_results.json; 
do 
    mydir=${json%_final_results.json}/;
    mkdir $mydir;
    mv $json $mydir;
done;

##Make an output directory
mkdir AMR

##Then run AMR on each of those directories
for json_dir in pmga/json/*/;
do
    result_dir='AMR/'${json_dir##pmga/json/}
    python3 ../runAST.py -i $json_dir -s Hi -o $result_dir;
done;

