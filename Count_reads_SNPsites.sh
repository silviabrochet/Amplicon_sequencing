#!/bin/bash

declare -a names=();


for i in ${names[@]};do

	PREFIX=$i
    FILE1=$PREFIX".assembled.fastq"
    perl count_reads_SNPs_ampliseq.pl -r $FILE1 -a 183_184_185_186_amplicon1.snps -i index_phasing_list.txt -pf 44 -pr 41 -l 199 -o $i"_strain_assignment.results"

done
