#!/bin/bash


declare -a names=();


for i in ${names[@]};do

	PREFIX=$i
    FILE1=$PREFIX"paired_L001_R1_001.fastq.gz"
    FILE2=$PREFIX"paired_L001_R2_001.fastq.gz"  
	/Software/BioTools/PEAR/bin/pear-0.9.6-bin-64 -f $FILE1 -r $FILE2 -o $i -m 290 -n 284 -j 4 -q 26 -v 10 -b 33
done

