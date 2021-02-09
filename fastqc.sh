#!/bin/bash


for i in $( ls *.fastq.gz); do

    ~/../Software/FastQC/fastqc $i >fastqc.log

done
