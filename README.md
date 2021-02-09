# Amplicon_sequencing

Analysis Amplicon sequencing data for Lactobacillus Firm5 (see Brochet et al 2021). 


#Quality check 
	
	bash fastqc.sh


#Trim and filter 
	
	bash trimmomatic.sh


#Quality check on trimmed and filtered fastq
	
	bash fastqc_postTrimmo.sh


#Assemble fw and rev read

Min/max amplicon size and overlap should be changed if the script is run on an amplicon of different size
	
	bash pear.sh

#Determine SNP sites in your reference sequences, 
	
	perl do.extract.SNPfreq_v_1.0.pl -i 183_184_185_186_amplicon.txt.fasta -o 183_184_185_186_amplicon1.snps -m 3

#Assign reads to strains according to SNP variants 

Here, a user-defined cutoff score is used, which was set to 0.75 within the script in our case (75% of the SNPs need to be present in a read)

This script needs as input:
	
SNP profiles generated in the previous steps
index_phasing_list.txt: file containing the phasing introduced for each primer (phasing is needed on the MiniSeq platform for sequencing low-complexity samples, due to the Illumina 2-channel technology). The position of SNPs will be shifted depending in the phasing (0 to +3)
Fw primer length (here 41); RV primer length (here 44) 
Length of amplicon 199 
Extension for output files (here: strain_assignment.results)
	
	bash Count_reads_SNPsites.sh

#Summarize number of reads per strain over all files in the current directory

A file name needs to be given. Two files are being produced: 
1. number of reads assigned to each strain in each sample
2. % of reads with correct length and % binned reads

		perl summarize_strain_counts.pl FILENAME

Once the number of reads are obtained, the proportions of each strain per sample are calculated using R:

	data[, -1] <- lapply( data[ , -1], function(x) x/sum(x, na.rm=TRUE) )
	library(reshape)
	prop <- melt(data, id.vars = "Strain")

and then normalised by the total number of CFUs. The limit of detection (LOD) is calculated by dividing the total number of reads by the total number of CFUs for each sample. 

