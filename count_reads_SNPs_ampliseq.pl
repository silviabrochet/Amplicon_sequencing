#!usr/bin/perl

use strict;

##read input from command line, print program requirements if no or missing input
our %opts;
&getopt(\%opts);

if(exists($opts{'h'})||(scalar(keys(%opts))==0)|| ! $opts{'r'} || ! -e $opts{'r'}) 
    #HELP OUTPUT
{ 
	print "usage: $0 [options]\n";
	print " -r Input FASTQ file (REQUIRED)\n";
	print " -a Input SNP file (REQUIRED)\n";
	print " -i Input index and phasing list\n";
	print " -pf total lenght of fw primer annealing sites (REQUIRED)\n";
	print " -pr total lenght of rv primer annealing sites (REQUIRED)\n";
	print " -l amplicon length without primer annealing sites (REQUIRED)\n";
    print " -o Output File (Default: STDOUT)\n";
	exit;
}

my $count_snps;
my %snp_hash;
my @snp_site;
my @strains;

##store SNP patterns for each strain in a hash by reading SNP file
open(READFILE,'<'.$opts{a});
while(<READFILE>){
    chomp $_;
    $count_snps += 1;
    #separate position and SNP variants
    @snp_site = split (/[\t,]+/, $_);
    #go through each SNP variant
    for (my $i = 1; $i<scalar(@snp_site); $i++){
            #separate base/freq/strain info about each SNP variant
            my @variant = split(':', $snp_site[$i]);
            #store base in variable
            my $nt = $variant[0];
            # separate different strains that show this SNP variant, and go through each of them
            @strains = split('\|', $variant[2]);
            for (my $s = 0; $s<scalar(@strains); $s++){
            	#store SNP pos, strain, and base at this position for this strain in a hash of hashes
            	$snp_hash{$snp_site[0]}{$strains[$s]} = $nt;
            	}
        }
    }

close READFILE;

#store info about which index is linked to which phasing
my %phasing_hash;
open(PHASELIST,'<'.$opts{i});
while(<PHASELIST>){
	chomp $_;
	my @phase_array = split('\t', $_);
	$phasing_hash{$phase_array[0]} = $phase_array[1];
	}
 close PHASELIST;	


#open OUTPUT
if($opts{o}){
    my $var = $opts{o};
    if(-d $opts{o}){
        print STDERR "WARNING: You gave us a folder, instead of a file. Printing on STDOUT instead\n";
    }else{
        my @path = split /\//, $opts{o};
        my $final = pop @path;
        my $path = join '/', @path;
        $path .='/';
        if( -d $path){
            my $tmp = $opts{o};
            qx/touch $tmp/;
            open STDOUT, ">", $opts{o};
        }else{
            print STDERR "WARNING: The folder where you are trying to print your output does not exist. Printing on STDOUT instead\n";
        }
    }
    
    
}

print STDOUT "Total SNP sites: $count_snps\n";



#get sequence lengths of different parts of the assembled read pair
my $fw_primer_len = $opts{pf};
my $rev_primer_len = $opts{pr};
my $amplicon_len = $opts{l};
my %score;
my @each_nt_read;
my $base;
my $count_fastq_reads = 0;
my $count_correct_fastq_reads = 0;
my %best_scores;
my %secondbest_scores;
my $count_total_goodreads;
my %count_strain_goodreads;

my @file_name_array = split('_', $opts{r});
my $phasing_fw = $phasing_hash{$file_name_array[7]};
my $phasing_rev = $phasing_hash{$file_name_array[6]};
my $phasing_tot = $phasing_fw + $phasing_rev;


open(FASTQ,'<'.$opts{r});
while(<FASTQ>){

    chomp $_;
    #identify the line with the actual sequence in the FASTQ file
	if($_ =~ /\A [ACGTN]+ \z/ix){
		$count_fastq_reads += 1;
    	#we will only consider assembled contigs that have the expected size
    	if (length($_) == $fw_primer_len + $rev_primer_len + $amplicon_len + $phasing_tot){
    		$count_correct_fastq_reads += 1;
    		@each_nt_read = split('', $_);
    		foreach my $snp_pos (keys %snp_hash) {
    			my $read_pos = $phasing_fw + $fw_primer_len + $snp_pos - 1;
    			$base = $each_nt_read[$read_pos];
  				foreach my $strain (keys %{$snp_hash{$snp_pos}}) {
    			
    			if ($snp_hash{$snp_pos}{$strain} eq $base){
    				$score{$strain} += 1;
    				
    					}
    				}
  				}
    		
		}
	}



my %rel_scores;



foreach my $strains (keys %score) { 
	$rel_scores{$strains} = $score{$strains}/$count_snps;
	#print STDOUT "$strains:\n  Total score: $score{$strains}\n  Relative score: $rel_scores{$strains}\n\n";
	}

my $count = 0;
foreach my $strains (sort { $rel_scores{$b} <=> $rel_scores{$a} } keys %rel_scores) {
#foreach my $strains	(sort {$b <=> $a} values %rel_scores){
	if ($count <2){

		if ($count == 0){
			$best_scores{$count_fastq_reads} = $rel_scores{$strains};
			#print "hey\t$rel_scores{$strains}\t$strains\t$count_fastq_reads\t$best_scores[$count_fastq_reads]\t";
			
			if ($rel_scores{$strains} >= 0.75){
				$count_total_goodreads += 1;
				$count_strain_goodreads{$strains} += 1;
				}
			}
		elsif($count == 1){
			$secondbest_scores{$count_fastq_reads} = $rel_scores{$strains}
			}
		$count += 1;
		}
	}
	
%score = ();
%rel_scores = ();	
	
}

print STDOUT "read no\tbest strain score\t 2nd best strain score\n";
foreach (sort { $a <=> $b } keys(%best_scores) ){
	print STDOUT "$_\t$best_scores{$_}\t$secondbest_scores{$_}\n";
	}

print STDOUT "Total reads\t$count_fastq_reads\nTotal reads of correct length\t$count_correct_fastq_reads\nTotal binned reads analyzed\t$count_total_goodreads\n\n";
print STDOUT "Binned reads per strain:\n";

foreach (sort { $a <=> $b } keys(%count_strain_goodreads)){
	print STDOUT "$_\t$count_strain_goodreads{$_}\n";
	}


close FASTQ;
close STDOUT;

#Now report the overall results, number of reads, good reads, reads per strain, relative per strain, best vs secondbest scores.

exit;


sub getopt (;$$) {
    my ($hash) = @_;
    my ($first,$rest);
    local $_;
    my @EXPORT;

    while (@ARGV && ($_ = $ARGV[0]) =~ /^-(.*)/) {
	($first) = ($1);
         $rest='';
	if (/^--$/) {	# early exit if --
	    shift @ARGV;
	    last;
	}
	if ($first ne '') {
	    if ($rest ne '') {
		shift(@ARGV);
	    }
	    else {
		shift(@ARGV);
		$rest = shift(@ARGV);
	    }
	    if (ref $hash) {
	        $$hash{$first} = $rest;
	    }
	    else {
	        ${"opt_$first"} = $rest;
	        push( @EXPORT, "\$opt_$first" );
	    }
	}
	else {
	    if (ref $hash) {
	        $$hash{$first} = 1;
	    }
	    else {
	        ${"opt_$first"} = 1;
	        push( @EXPORT, "\$opt_$first" );
	    }
	    if ($rest ne '') {
		$ARGV[0] = "-$rest";
	    }
	    else {
		shift(@ARGV);
	    }
	}
    }
    unless (ref $hash) { 
	local $Exporter::ExportLevel = 1;
	import Getopt::Std;
    }
}
