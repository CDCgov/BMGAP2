#! /bin/bash
#$ -cwd

NSLOTS=${NSLOTS:-1}  # Number of slots available for parallel processing

###################################################
# Definining variables
R1=${1}        #R1 fastq
R2=${2}        #R2 fastq
OUT_DIR=${3}   #Sample name directory
NAME=${4}      #Sample name (M#)
###################################################

###################################################
# Remove human DNA
bowtie2 -p $NSLOTS --local -t -x analysis_scripts/hg38/genome --un-conc-gz $OUT_DIR/${NAME}_R%_dedup_NoHuman.fastq.gz -1 $R1 -2 $R2 > \
  $OUT_DIR/${NAME}.human.sam
samtools view -bhS $OUT_DIR/${NAME}.human.sam > $OUT_DIR/${NAME}.human.bam
rm -v $OUT_DIR/${NAME}.human.sam

###################################################
# adapter, and quality trimming, filter for Phred quality 20
cutadapt -o $OUT_DIR/${NAME}_R1_dedup_NoHuman_cutTruSeq_trim.fastq.gz -p $OUT_DIR/${NAME}_R2_dedup_NoHuman_cutTruSeq_trim.fastq.gz -n 5 --trim-n -m 50 -q 15 -a Prefix=AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -A Universal_rc=AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT $OUT_DIR/${NAME}_R1_dedup_NoHuman.fastq.gz $OUT_DIR/${NAME}_R2_dedup_NoHuman.fastq.gz

###################################################
# assemble the cleaned PE fastq reads
# spades.py -t $NSLOTS -m 32 --pe1-1 $OUT_DIR/${NAME}_R1_dedup_NoHuman_cutTruSeq_trim.fastq.gz --pe1-2 $OUT_DIR/${NAME}_R2_dedup_NoHuman_cutTruSeq_trim.fastq.gz -o $OUT_DIR/${NAME}_SPAdes
# assemble the downsampled PE fastq reads
spades.py -t $NSLOTS -m 32 --pe1-1 $OUT_DIR/${NAME}_R1_dedup_NoHuman_cutTruSeq_trim_downsampled.fastq.gz --pe1-2 $OUT_DIR/${NAME}_R2_dedup_NoHuman_cutTruSeq_trim_downsampled.fastq.gz -o $OUT_DIR/${NAME}_SPAdes

echo "Completed alignment of $NAME!"
###################################################
