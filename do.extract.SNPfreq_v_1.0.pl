#!usr/bin/perl

use strict;

our %opts;
&getopt(\%opts);

if(exists($opts{'h'})||(scalar(keys(%opts))==0)|| ! $opts{'i'} || ! -e $opts{'i'}) 
    #HELP OUTPUT
{ 
	print "usage: $0 [options]\n";
	print " -i Input FASTA file (REQUIRED)\n";
    print " -o Output File (Default: STDOUT)\n";
	print ' -m {1|2|3} output method (Optional). Default == 3' . "\n";
	print "\t -m 1 Only the variants\n";
	print "\t -m 2 Variants and Frequencies\n";
	print "\t -m 3 Variants, Frequencies and Sequencies\n";
    exit;
}

my %aln;
#load alignment

open FILE, $opts{i};


my $name;
while(<FILE>){
    chomp;
    if(/>(\S+)/){
        $name = $1;
    }else{
        $aln{$name} .= $_;
    }
}

close FILE;

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
#checking output format
if(!$opts{m}){
    print STDERR "WARNING: No output format was defined. Using format 1, instead\n";
    $opts{m} = 1;
}elsif($opts{m}>3 || $opts{m}<1){
    print STDERR "WARNING: wrong output format defined. Using format 1, instead\n";
    $opts{m} = 1;
}
#keep the number of seqs in mind
my $dash_seqs = scalar keys %aln;
#check size;
my @keys = sort {$a cmp $b} keys %aln;
my %sizes;
foreach(@keys){
    my $name = $_;
    my $size = length $aln{$name};
    $sizes{sizes}{$size}++;
    $sizes{names}{$size} .= " " . $name;
}

if (scalar keys %{$sizes{sizes}}>1 ){
    my $scalar = scalar keys %{$sizes{sizes}};
    print STDERR "ERROR: Samples do not seem to be aligned. There seems to be $scalar groups of samples:\n";
    my @sizes = sort {$sizes{sizes}{$b} <=> $sizes{sizes}{$a} } keys %{$sizes{sizes}};
    while (@sizes){
        my $size = shift @sizes;
        print STDERR "SIZE: $size. Samples: " . $sizes{names}{$size} . "\n";
    }
    exit;
}else{
    #start checking for SNPs;
    my %snps;
    my @seqs = sort {$a cmp $b} keys %aln;
    my $nseqs = scalar @seqs;
    while(@seqs){
        my $id = shift @seqs;
        my @seq = split '', $aln{$id};
        my $scalar = scalar @seq;
        for (my $i = 0; $i< $scalar; $i++){
            my $nt = shift @seq;
            $snps{$i}{variants}{$nt}++;
            $snps{$i}{seqs}{$nt} .= " ". $id;
        }
    }
    #print OUTPUT in STDOUT 
    my @positions = sort { $a <=> $b } keys %snps;
    my $scalar = scalar @positions;
    print "$nseqs SEQS $scalar POSITIONS\n";
    while(@positions){
        my $position = shift @positions;
        my $line = ($position + 1) . "\t";
        if (scalar keys $snps{$position}{variants} > 1){
            my @nts = sort { $a cmp $b } keys %{$snps{$position}{variants}};
            while(@nts){
                my $snp = shift @nts;
                if ($opts{m} == 1){
                    $line .= $snp;
                }elsif($opts{m} == 2){
                    $line .= $snp;
                    my $num = $snps{$position}{variants}{$snp};
                    my $freq = sprintf "%.3f", $num / $dash_seqs;
                    $line .=  ':' . $freq;
                }elsif($opts{m} == 3){
                    $line .= $snp;
                    my $num = $snps{$position}{variants}{$snp};
                    my $freq = sprintf "%.3f", $num / $dash_seqs;
                    $line .=  ':' . $freq;
                    my @names = split " ", $snps{$position}{seqs}{$snp};
                    shift @names if $names[0] eq '';
                    my $names;
                    while (@names){
                        $names .= shift @names;
                        $names .= '|' if scalar @names != 0;
                    }
                    $line .= ':' . $names;
                }
                $line .= ',' if scalar @nts != 0;
            }
            print STDOUT $line . "\n";
        }else{
            my $wow;
        }
        
    }
    
}

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