#!/bin/bash


declare -a names=();


for i in ${names[@]};do

	PREFIX=$i
    FILE1=$PREFIX"_R1_001.fastq.gz"
    FILE2=$PREFIX"_R2_001.fastq.gz"  
    FILE1_P=$PREFIX"_paired_R1_001.fastq.gz"
    FILE1_NP=$PREFIX"_unpaired_R1_001.fastq.gz"
    FILE2_P=$PREFIX"_paired_R2_001.fastq.gz"
    FILE2_NP=$PREFIX"_unpaired_R2_001.fastq.gz"
    java -jar trimmomatic-0.38.jar PE -threads 10 $FILE1 $FILE2  $FILE1_P $FILE1_NP $FILE2_P $FILE2_NP ILLUMINACLIP:readthrough_primers_Firm5.fa:2:3:10 LEADING:28 TRAILING:28 SLIDINGWINDOW:4:15 MINLEN:90
done

