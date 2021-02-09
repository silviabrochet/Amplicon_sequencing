#!/bin/bash


for i in $( ls *_paired_L001_R*_001.fastq.gz); do

    ~/..Software/FastQC/fastqc $i >fastqc_post_trimmo.log

done
